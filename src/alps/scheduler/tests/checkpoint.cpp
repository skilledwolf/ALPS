// SPDX-License-Identifier: MIT
#include <alps/scheduler/task.h>
#include <alps/osiris/xdrdump.h>

#include <boost/filesystem.hpp>

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {
void require(bool condition, char const* message) {
    if (!condition) throw std::runtime_error(message);
}

alps::Parameters parameters() {
    alps::Parameters result;
    result["SEED"] = 17;
    return result;
}

class worker : public alps::scheduler::Worker {
public:
    worker() : Worker(parameters()) {}
    void save(alps::hdf5::archive& ar) const override {
        ar["state"] << state;
        if (fail) throw std::runtime_error("HDF5 hook failure");
    }
    void save_worker(alps::ODump& dump) const override { dump << state; }
    int state = 1;
    bool fail = false;
};

class task : public alps::scheduler::Task {
public:
    task() : Task(alps::ProcessList{}, parameters()) {}
    void dostep() override {}
    void save(alps::hdf5::archive& ar) const override {
        ar["state"] << state;
        if (fail) throw std::runtime_error("HDF5 hook failure");
    }
    int state = 1;
    bool fail = false;
    bool fail_xml = false;
protected:
    void write_xml_header(alps::oxstream&) const override {}
    void write_xml_trailer(alps::oxstream&) const override {}
    void write_xml_body(alps::oxstream& out, boost::filesystem::path const&, bool) const override {
        out << alps::start_tag("STATE") << alps::attribute("value", state) << alps::end_tag("STATE");
        if (fail_xml) throw std::runtime_error("XML hook failure");
    }
};

int hdf5_state(boost::filesystem::path const& filename) {
    alps::hdf5::archive ar(filename);
    int value = 0;
    ar["/state"] >> value;
    return value;
}

int xdr_state(boost::filesystem::path const& filename) {
    alps::IXDRFileDump dump(filename);
    int value = 0;
    dump >> value;
    return value;
}

std::string contents(boost::filesystem::path const& filename) {
    std::ifstream file(filename.string());
    require(bool(file), "missing XML snapshot");
    std::ostringstream text;
    text << file.rdbuf();
    return text.str();
}

template <typename Save>
void expect_failure(Save save, char const* expected) {
    try { save(); }
    catch (std::runtime_error const& error) {
        require(error.what() == std::string(expected), "checkpoint changed the hook exception");
        return;
    }
    throw std::runtime_error("failed checkpoint was accepted");
}
}

int main() {
    auto const directory = boost::filesystem::unique_path("test_scheduler_checkpoint.%%%%-%%%%");
    boost::filesystem::create_directory(directory);
    try {
        worker run;
        auto const worker_xdr = directory / "worker.xdr";
        auto const worker_h5 = directory / "worker.h5";
        run.save_to_file(worker_xdr, worker_h5);
        run.state = 2;
        run.fail = true;
        expect_failure([&] { run.save_to_file(worker_xdr, worker_h5); }, "HDF5 hook failure");
        require(hdf5_state(worker_h5) == 1 && xdr_state(worker_xdr) == 1,
                "failed HDF5 save changed a prior worker snapshot");
        expect_failure([&] { run.save_to_file(directory / "missing.xdr", directory / "missing.h5"); },
                       "HDF5 hook failure");
        require(!boost::filesystem::exists(directory / "missing.h5") &&
                !boost::filesystem::exists(directory / "missing.xdr"),
                "failed initial HDF5 save published worker files");
        run.fail = false;
        run.save_to_file(worker_xdr, worker_h5);
        require(hdf5_state(worker_h5) == 2 && xdr_state(worker_xdr) == 2,
                "successful worker save did not replace both snapshots");

        task simulation;
        auto const task_xml = directory / "task.xml";
        auto const task_h5 = directory / "task.h5";
        simulation.checkpoint(task_xml);
        auto const original_xml = contents(task_xml);
        simulation.state = 2;
        simulation.fail = true;
        expect_failure([&] { simulation.checkpoint(task_xml); }, "HDF5 hook failure");
        require(hdf5_state(task_h5) == 1 && contents(task_xml) == original_xml,
                "failed HDF5 stage changed a prior task snapshot");
        simulation.fail = false;
        simulation.fail_xml = true;
        expect_failure([&] { simulation.checkpoint(task_xml); }, "XML hook failure");
        require(hdf5_state(task_h5) == 1 && contents(task_xml) == original_xml &&
                hdf5_state(directory / "task.h5.bak") == 2,
                "XML failure published the closed HDF5 stage early");
        simulation.fail_xml = false;
        simulation.checkpoint(task_xml);
        require(hdf5_state(task_h5) == 2 && contents(task_xml) != original_xml &&
                !boost::filesystem::exists(directory / "task.h5.bak"),
                "successful task save did not publish the staged pair");
        simulation.state = 3;
        simulation.checkpoint_hdf5(task_xml);
        require(hdf5_state(task_h5) == 3, "HDF5-only task save did not replace its snapshot");
        for (auto const& entry : boost::filesystem::directory_iterator(directory))
            require(entry.path().filename().string().find(".tmp.") == std::string::npos,
                    "scheduler save left temporary HDF5 data");

        boost::filesystem::remove_all(directory);
        return 0;
    } catch (std::exception const& error) {
        std::cerr << error.what() << '\n';
        boost::system::error_code cleanup_error;
        boost::filesystem::remove_all(directory, cleanup_error);
        return 1;
    }
}
