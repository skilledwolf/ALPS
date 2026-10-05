// SPDX-License-Identifier: MIT
#include <alps/mcbase.hpp>
#include <alps/check_schedule.hpp>
#include <alps/alea/hdf5.hpp>
#include <alps/alea/checkpoint.hpp>

#include <boost/filesystem.hpp>

#include <iostream>
#include <stdexcept>
#include <string>

namespace {
void require(bool condition, char const* message) {
    if (!condition) throw std::runtime_error(message);
}

class simulation : public alps::mcbase {
public:
    simulation() : mcbase(alps::params{}) {}
    using mcbase::save;
    using mcbase::load;
    void update() override {}
    void measure() override {}
    double fraction_completed() const override { return 1; }
    void save(alps::hdf5::archive& ar) const override {
        mcbase::save(ar);
        ar["state"] << state;
        if (fail) throw std::runtime_error("checkpoint hook failure");
    }
    void load(alps::hdf5::archive& ar) override {
        mcbase::load(ar);
        ar["state"] >> state;
    }
    int state = 0;
    bool fail = false;
};

bool has_temporary(boost::filesystem::path const& directory, std::string const& target) {
    for (auto const& entry : boost::filesystem::directory_iterator(directory))
        if (entry.path().filename().string().find(target + ".tmp.") == 0) return true;
    return false;
}

void save_failure(simulation const& sim, boost::filesystem::path const& target) {
    bool rejected = false;
    try { sim.save(target); }
    catch (std::runtime_error const& error) {
        require(std::string(error.what()) == "checkpoint hook failure", "save lost the original exception");
        rejected = true;
    }
    require(rejected, "failed checkpoint was accepted");
    require(!has_temporary(target.parent_path(), target.filename().string()), "failed checkpoint left temporary data");
}

alps::alea::batch_result<double> read_result(alps::hdf5::archive& ar, std::string const& path) {
    alps::alea::hdf5_serializer serializer(ar, path);
    alps::alea::batch_result<double> result;
    alps::alea::deserialize(serializer, "", result);
    return result;
}

void failed_base_load(simulation& sim, boost::filesystem::path const& filename) {
    auto const parameters = sim.get_parameters();
    auto const handle = sim.measurement("value");
    auto const result = handle->result();
    auto random = sim.get_random();
    bool rejected = false;
    try { sim.load(filename); } catch (std::exception const&) { rejected = true; }
    require(rejected, "invalid base checkpoint was accepted");
    require(sim.get_parameters()["retained"].as<int>() == parameters["retained"].as<int>()
            && sim.measurement("value") == handle && handle->result() == result
            && sim.get_random()() == random(), "failed base load changed prior state");
}

}

int main() {
    alps::check_schedule schedule(60., 60.);
    require(schedule.check_interval()==0. && schedule.pending(), "Initial progress check must be immediate");
    schedule.update(0.25);
    require(schedule.check_interval()==60. && !schedule.pending(), "Updated progress check must respect interval");
    auto const directory = boost::filesystem::unique_path("test_mcbase_checkpoint.%%%%-%%%%");
    boost::filesystem::create_directory(directory);
    try {
        boost::variate_generator<boost::mt19937,boost::uniform_01<double>> mt(boost::mt19937(37),{});
        boost::random::lagged_fibonacci<uint32_t,48,607,273> fib(37);
        alps::random01 native_mt(37), native_fib(37,"lagged_fibonacci607");
        for (int i=0;i<2000;++i) {
            require(native_mt()==mt(), "MT19937 stream changed");
            require(native_fib()==boost::uniform_01<double>()(fib), "wrong released Fibonacci engine");
        }
        auto const checkpoint = directory / "state.h5";
        simulation sim;
        sim.state = 1;
        sim.save(checkpoint);
        require(!has_temporary(directory, "state.h5"), "successful save left temporary data");
        simulation restored;
        restored.load(checkpoint);
        require(restored.state == 1, "checkpoint did not round trip");

        sim.state = 2;
        sim.fail = true;
        save_failure(sim, checkpoint);
        restored.load(checkpoint);
        require(restored.state == 1, "failed save changed the prior checkpoint");
        auto const missing = directory / "missing.h5";
        save_failure(sim, missing);
        require(!boost::filesystem::exists(missing), "failed save published a new checkpoint");

        sim.fail = false;
        sim.save(checkpoint);
        restored.load(checkpoint);
        require(restored.state == 2, "successful save did not replace the prior checkpoint");

        // Force publication to fail after serialization and close. Existing
        // destination contents must survive and no temporary data may remain.
        auto const blocked = directory / "blocked.h5";
        boost::filesystem::create_directory(blocked);
        auto const sentinel = blocked / "untouched.h5";
        {
            alps::hdf5::archive ar(sentinel, "w");
            ar["/sentinel"] << 7;
        }
        bool rejected = false;
        try { sim.save(blocked); }
        catch (boost::filesystem::filesystem_error const&) { rejected = true; }
        require(rejected, "publication failure was not reported");
        require(!has_temporary(directory, "blocked.h5"), "publication failure left temporary data");
        {
            alps::hdf5::archive ar(sentinel, "r");
            int value = 0;
            ar["/sentinel"] >> value;
            require(value == 7, "publication failure changed unrelated data");
        }

        auto const results_file = directory / "results.h5";
        {
            alps::hdf5::archive ar(results_file, "w");
            ar["/unrelated"] << 8;
        }
        sim.get_measurements().emplace("value", std::make_shared<alps::alea::batch_acc<double>>(1, 4));
        *sim.measurement("value") << alps::alea::make_adapter(1.);
        auto const results = sim.collect_results();
        alps::save_results(results, sim.get_parameters(), results_file, "/first/results");
        alps::save_results(results, sim.get_parameters(), results_file, "/second/results");
        {
            alps::hdf5::archive ar(results_file);
            int value = 0;
            ar["/unrelated"] >> value;
            require(value == 8 && ar.is_group("/first/results") && ar.is_group("/second/results"),
                    "saving results at a path erased unrelated groups");
        }
        sim.get_parameters()["empty_save"] = "updated";
        alps::save_results(simulation::results_type{}, sim.get_parameters(), results_file, "/first/results");
        {
            alps::hdf5::archive ar(results_file);
            require(ar.is_group("/first/results") &&
                    ar.list_children("/first/results").empty(),
                    "saving empty results retained old entries");
            require(read_result(ar, "/second/results/value").count() == 1,
                    "saving empty results changed an unrelated result collection");
            alps::params parameters;
            ar["/parameters"] >> parameters;
            require(parameters["empty_save"].as<std::string>() == "updated",
                    "saving empty results skipped the parameters");
            int value = 0;
            ar["/unrelated"] >> value;
            require(value == 8, "saving empty results changed unrelated data");
        }

        // Base loads stage parameters, every native accumulator and RNG.
        sim.get_parameters()["retained"] = 17;
        for (int i=0; i<41; ++i)
            *sim.measurement("value") << alps::alea::make_adapter(double(i));
        sim.save(checkpoint);
        auto retained = sim.measurement("value");
        simulation resumed;
        resumed.load(checkpoint);
        require(resumed.measurement("value")->result() == retained->result(),
                "native MC accumulator checkpoint lost state");
        auto const invalid = directory / "invalid.h5";
        for (int damage=0; damage<4; ++damage) {
            sim.save(invalid);
            {
                alps::hdf5::archive ar(invalid, "a");
                auto changed = sim.get_parameters();
                changed["retained"] = 99;
                ar["/parameters"] << changed;
                auto const clone = "/simulation/realizations/0/clones/0/";
                if (damage == 0)
                    ar[std::string(clone)+"checkpoint/engine/engine"] << std::string("invalid engine");
                else if (damage == 1)
                    ar.delete_group(std::string(clone)+"measurements/value");
                else if (damage == 2) {
                    alps::alea::batch_acc<double> wrong_shape(2,4);
                    alps::alea::hdf5_serializer serializer(ar, std::string(clone)+"measurements/value");
                    alps::alea::serialize(serializer, "", wrong_shape);
                } else
                    ar[std::string(clone)+"measurements/value/@kind"] << std::uint64_t(5);
            }
            failed_base_load(sim, invalid);
        }
        auto null_handle = sim.measurement("value");
        sim.measurement("value").reset();
        bool invalid_measurement = false;
        try { sim.save(checkpoint); }
        catch (std::invalid_argument const&) { invalid_measurement = true; }
        require(invalid_measurement, "null MC measurement was accepted");
        sim.measurement("value") = null_handle;
        resumed.load(checkpoint);
        require(resumed.measurement("value")->result() == null_handle->result(),
                "invalid measurement save changed prior checkpoint");
        boost::filesystem::remove_all(directory);
        return 0;
    } catch (std::exception const& error) {
        std::cerr << error.what() << '\n';
        boost::system::error_code cleanup_error;
        boost::filesystem::remove_all(directory, cleanup_error);
        return 1;
    }
}
