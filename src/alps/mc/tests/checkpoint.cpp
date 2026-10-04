// SPDX-License-Identifier: MIT
#include <alps/mcbase.hpp>

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

template <typename Collection>
void collection_replacement(boost::filesystem::path const& filename, Collection collection) {
    auto const original = collection;
    alps::hdf5::archive ar(filename, "w");
    ar["/unrelated"] << 8;
    ar["/collection"] << collection;
    ar["/collection/@tag"] << 9;

    collection.erase("gone");
    ar["/collection/orphan"] << 10;
    ar["/collection"] << collection;
    Collection restored;
    ar["/collection"] >> restored;
    require(restored.size() == 1 && restored.has("kept") && !restored.has("gone"),
            "collection save retained a removed entry");
    require(!ar.is_data("/collection/orphan"), "collection save retained an orphan dataset");
    int value = 0;
    ar["/unrelated"] >> value;
    require(value == 8, "collection replacement changed an unrelated group");
    ar["/collection/@tag"] >> value;
    require(value == 9, "collection replacement removed the group attribute");

    collection.clear();
    ar["/collection"] << collection;
    Collection empty;
    ar["/collection"] >> empty;
    require(ar.is_group("/collection") && ar.list_children("/collection").empty() && empty.empty(),
            "empty collection save retained old entries");
    ar["/collection/@tag"] >> value;
    require(value == 9, "empty collection save removed the group attribute");

    // Root is also a valid collection context. Replace its children without
    // trying to unlink root, then persist an explicitly empty collection.
    ar["/"] << original;
    ar["/@tag"] << 11;
    Collection root;
    ar["/"] >> root;
    require(root.size() == 2 && root.has("kept") && root.has("gone"),
            "root collection save did not replace its children");
    ar["/"] << collection;
    Collection empty_root;
    ar["/"] >> empty_root;
    require(ar.list_children("/").empty() && empty_root.empty(),
            "empty root collection save retained old entries");
    ar["/@tag"] >> value;
    require(value == 11, "root collection replacement removed the root attribute");
    ar.close();
}
}

int main() {
    auto const directory = boost::filesystem::unique_path("test_mcbase_checkpoint.%%%%-%%%%");
    boost::filesystem::create_directory(directory);
    try {
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
        sim.get_measurements() << alps::accumulator::RealObservable("value");
        sim.get_measurements()["value"] << 1.;
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
            simulation::results_type empty, preserved;
            ar["/first/results"] >> empty;
            ar["/second/results"] >> preserved;
            require(empty.empty() && ar.is_group("/first/results") &&
                    ar.list_children("/first/results").empty(),
                    "saving empty results retained old entries");
            require(preserved.size() == 1 && preserved.has("value") && preserved["value"].count() == 1,
                    "saving empty results changed an unrelated result collection");
            alps::params parameters;
            ar["/parameters"] >> parameters;
            require(parameters["empty_save"].as<std::string>() == "updated",
                    "saving empty results skipped the parameters");
            int value = 0;
            ar["/unrelated"] >> value;
            require(value == 8, "saving empty results changed unrelated data");
        }

        alps::mcobservables observables;
        observables << alps::accumulator::RealObservable("kept") << alps::accumulator::RealObservable("gone");
        observables["kept"] << 1.;
        observables["gone"] << 2.;
        alps::mcresults collection_results;
        collection_results.insert("kept", alps::mcresult(observables["kept"]));
        collection_results.insert("gone", alps::mcresult(observables["gone"]));
        collection_replacement(directory / "observables.h5", observables);
        collection_replacement(directory / "collection_results.h5", collection_results);
        boost::filesystem::remove_all(directory);
        return 0;
    } catch (std::exception const& error) {
        std::cerr << error.what() << '\n';
        boost::system::error_code cleanup_error;
        boost::filesystem::remove_all(directory, cleanup_error);
        return 1;
    }
}
