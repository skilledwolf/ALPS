// SPDX-License-Identifier: MIT
// Characterize the legacy task layer (WorkerTask and MCSimulation) that a
// configuration front end has to drive: which runs a task creates or loads,
// how their progress adds up, how their measurements are pooled, which
// checkpoint files a task writes, and what a resumed task does with changed
// parameters. Only the single-process paths run here. Remote runs and
// add_process need MPI messaging and are not exercised.
//
// The probe run counts its steps and records values derived from its SEED,
// without drawing random numbers, so the results do not depend on the RNG.
// With ALPS_TASK_DUMP_DIR set, every test copies the HDF5 files it wrote
// there, for comparing builds.
#include <alps/alea.h>
#include <alps/hdf5.hpp>
#include <alps/osiris/comm.h>
#include <alps/osiris/dump.h>
#include <alps/parameter.h>
#include <alps/parser/xmlstream.h>
#include <alps/parser/xslt_path.h>
#include <alps/scheduler/montecarlo.h>
#include <alps/scheduler/options.h>
#include <alps/scheduler/scheduler.h>
#include <boost/filesystem.hpp>
#include <gtest/gtest.h>

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <memory>
#include <regex>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

namespace fs = boost::filesystem;

alps::Parameters parameters(const std::string& text) {
    alps::Parameters p;
    std::istringstream input(text);
    p.parse(input);
    return p;
}

class counting_run : public alps::scheduler::MCRun {
public:
    counting_run(const alps::ProcessList& where, const alps::Parameters& p, int node)
        : alps::scheduler::MCRun(where, p, node), steps_done_(0) {
        // BINNING=0 records Value without binning, for the one merge that
        // cannot take binned data (GetAndRemoveObservableCountsSavedResultsAgain).
        if (static_cast<int>(parms.value_or_default("BINNING", 1)))
            measurements << alps::RealObservable("Value");
        else
            measurements << alps::SimpleRealObservable("Value");
    }
    void dostep() {
        ++steps_done_;
        measurements["Value"] << value(steps_done_);
    }
    bool is_thermalized() const { return true; }
    double work_done() const {
        return double(steps_done_) / static_cast<double>(parms["STEPS"]);
    }
    bool change_parameter(const std::string& name, const alps::StringValue&) {
        return name == "STEPS";
    }
    void save(alps::ODump& dump) const { dump << steps_done_; }
    void load(alps::IDump& dump) { dump >> steps_done_; }

    uint32_t steps_done() const { return steps_done_; }
    const alps::Parameters& parameters() const { return parms; }
    const alps::ObservableSet& own_measurements() const { return measurements; }

private:
    // Distinct values per run, derived from the run's SEED and OFFSET.
    double value(uint32_t step) const {
        return static_cast<double>(parms.value_or_default("OFFSET", 0.)) +
               0.37 * static_cast<uint32_t>(parms["SEED"]) + 0.1 * step;
    }
    uint32_t steps_done_;
};

class probe_simulation : public alps::scheduler::MCSimulation {
public:
    probe_simulation(const alps::ProcessList& where, const alps::Parameters& p)
        : alps::scheduler::MCSimulation(where, p) {}
    probe_simulation(const alps::ProcessList& where, const fs::path& file)
        : alps::scheduler::MCSimulation(where, file) {}

    std::vector<RunStatus> status() const { return workerstatus; }
    const counting_run& run(std::size_t i) const {
        return dynamic_cast<const counting_run&>(*runs.at(i));
    }
    void steps(int n) {
        start();
        for (int i = 0; i < n; ++i)
            run();
    }
    using alps::scheduler::MCSimulation::run;
};

typedef alps::scheduler::WorkerTask::RunStatus RunStatus;
const RunStatus local_run = alps::scheduler::WorkerTask::LocalRun;
const RunStatus on_dump = alps::scheduler::WorkerTask::RunOnDump;

const alps::ProcessList here(1, alps::local_process());

double mean(const alps::ObservableSet& set, const std::string& name) {
    return alps::RealObsevaluator(set[name]).mean();
}

double count(const alps::ObservableSet& set, const std::string& name) {
    return alps::RealObsevaluator(set[name]).count();
}

std::string read_file(const fs::path& file) {
    std::ifstream in(file.string().c_str());
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

// The (format, file) pairs of the CHECKPOINT elements in a task file, in order.
std::vector<std::pair<std::string, std::string>> checkpoints(const fs::path& task_file) {
    const std::string text = read_file(task_file);
    const std::regex element(R"re(<CHECKPOINT\s+format="([^"]*)"\s+file="([^"]*)")re");
    std::vector<std::pair<std::string, std::string>> out;
    for (std::sregex_iterator it(text.begin(), text.end(), element), end; it != end; ++it)
        out.emplace_back((*it)[1], (*it)[2]);
    return out;
}

std::vector<std::string> files_in(const fs::path& dir) {
    std::vector<std::string> out;
    for (fs::directory_iterator it(dir), end; it != end; ++it)
        out.push_back(it->path().filename().string());
    std::sort(out.begin(), out.end());
    return out;
}

class TaskLayer : public ::testing::Test {
protected:
    TaskLayer()
        : scheduler(options, factory),
          dir(fs::temp_directory_path() / fs::unique_path("alps-task-%%%%-%%%%-%%%%")) {
        fs::create_directories(dir);
    }

    ~TaskLayer() {
        if (const char* dump = std::getenv("ALPS_TASK_DUMP_DIR")) {
            const ::testing::TestInfo* test =
                ::testing::UnitTest::GetInstance()->current_test_info();
            const fs::path target = fs::path(dump) / test->name();
            fs::create_directories(target);
            for (fs::directory_iterator it(dir), end; it != end; ++it)
                if (it->path().extension() == ".h5")
                    fs::copy_file(it->path(), target / it->path().filename(),
                                  fs::copy_options::overwrite_existing);
        }
        fs::remove_all(dir);
    }

    // Runs a fresh one-run task for a number of steps and checkpoints it as
    // <name>.xml. Returns the run's XDR and HDF5 checkpoint files.
    std::pair<fs::path, fs::path> source_run(const std::string& name, const std::string& text,
                                             int steps) {
        probe_simulation task(here, parameters(text));
        task.steps(steps);
        task.checkpoint(dir / (name + ".xml"));
        return {dir / (name + ".run1"), dir / (name + ".run1.h5")};
    }

    // Writes <name>.xml with the given parameters and one MCRUN element per
    // run, as a task with several runs writes it: the run checkpoints are
    // copied to <name>.run<N> and <name>.run<N>.h5. The header is Task's; the
    // parser needs the stylesheet line after the XML declaration.
    fs::path compose(const std::string& name, const std::string& text,
                     const std::vector<std::pair<fs::path, fs::path>>& runs) {
        const fs::path file = dir / (name + ".xml");
        alps::oxstream out(file);
        out << alps::header("UTF-8") << alps::stylesheet(alps::xslt_path("ALPS.xsl"))
            << alps::start_tag("SIMULATION") << parameters(text);
        for (std::size_t i = 0; i < runs.size(); ++i) {
            const std::string run = name + ".run" + std::to_string(i + 1);
            fs::copy_file(runs[i].first, dir / run);
            fs::copy_file(runs[i].second, dir / (run + ".h5"));
            out << alps::start_tag("MCRUN") << alps::start_tag("CHECKPOINT")
                << alps::attribute("format", "osiris") << alps::attribute("file", run)
                << alps::end_tag("CHECKPOINT") << alps::start_tag("CHECKPOINT")
                << alps::attribute("format", "hdf5") << alps::attribute("file", run + ".h5")
                << alps::end_tag("CHECKPOINT") << alps::end_tag("MCRUN");
        }
        out << alps::end_tag("SIMULATION");
        return file;
    }

    // Source runs with SEED 11, 22, ... and 2, 3, ... steps, out of STEPS=16,
    // and any further parameters.
    // Every run has at least two values, since the mean of one is NaN.
    std::vector<std::pair<fs::path, fs::path>> source_runs(int n, const std::string& more = "") {
        std::vector<std::pair<fs::path, fs::path>> out;
        for (int i = 1; i <= n; ++i)
            out.push_back(source_run("source" + std::to_string(i),
                                     "STEPS=16;SEED=" + std::to_string(11 * i) + more, i + 1));
        return out;
    }

    alps::scheduler::SimpleMCFactory<counting_run> factory;
    alps::scheduler::NoJobfileOptions options;
    alps::scheduler::Scheduler scheduler;
    fs::path dir;
};

// One process gives a task one local run. The task's work is
// WORK_FACTOR * (1 - work done), and the task is finished when that is no
// longer positive, also after the run overshoots.
TEST_F(TaskLayer, CreatesOneLocalRunAndCountsItsWork) {
    probe_simulation task(here, parameters("STEPS=4;SEED=1"));
    ASSERT_EQ(task.runs.size(), 1u);
    EXPECT_EQ(task.status(), std::vector<RunStatus>{local_run});
    EXPECT_DOUBLE_EQ(task.work(), 1.);
    EXPECT_FALSE(task.finished_notime());

    task.steps(1);
    EXPECT_DOUBLE_EQ(task.work_done(), 0.25);
    EXPECT_DOUBLE_EQ(task.work(), 0.75);
    task.steps(3);
    EXPECT_DOUBLE_EQ(task.work(), 0.);
    EXPECT_TRUE(task.finished_notime());
    task.steps(2);
    EXPECT_DOUBLE_EQ(task.work_done(), 1.5);
    EXPECT_DOUBLE_EQ(task.work(), -0.5);
    EXPECT_TRUE(task.finished_notime());

    probe_simulation scaled(here, parameters("STEPS=4;SEED=1;X=1.5;WORK_FACTOR=\"2*X\""));
    scaled.steps(1);
    EXPECT_DOUBLE_EQ(scaled.work(), 2.25);
}

// The percentage the scheduler reports is 1 - work, clipped to [0, 1], so it
// includes WORK_FACTOR and is not the fraction of work done.
TEST_F(TaskLayer, ReportsPercentageFromScaledWork) {
    probe_simulation task(here, parameters("STEPS=4;SEED=1;WORK_FACTOR=0.5"));
    task.steps(1);
    double more_time = -1., percentage = -1.;
    EXPECT_FALSE(task.finished(more_time, percentage));
    EXPECT_DOUBLE_EQ(percentage, 0.625);

    probe_simulation heavy(here, parameters("STEPS=4;SEED=1;WORK_FACTOR=3"));
    heavy.steps(1);
    more_time = -1.;
    percentage = -1.;
    EXPECT_FALSE(heavy.finished(more_time, percentage));
    EXPECT_DOUBLE_EQ(percentage, 0.);
}

// A checkpoint <base>.xml writes <base>.h5 and, per run, <base>.run<N> (XDR)
// and <base>.run<N>.h5, named relative to the task file. Later checkpoints
// keep the run names, also when the task file is renamed.
TEST_F(TaskLayer, NamesRunCheckpointsAfterTheFirstTaskFile) {
    probe_simulation task(here, parameters("STEPS=4;SEED=1"));
    task.steps(1);
    task.checkpoint(dir / "job.task1.xml");
    EXPECT_EQ(files_in(dir), (std::vector<std::string>{"job.task1.h5", "job.task1.run1",
                                                       "job.task1.run1.h5", "job.task1.xml"}));
    using checkpoint = std::pair<std::string, std::string>;
    EXPECT_EQ(checkpoints(dir / "job.task1.xml"),
              (std::vector<checkpoint>{{"osiris", "job.task1.run1"},
                                       {"hdf5", "job.task1.run1.h5"}}));

    task.steps(1);
    task.checkpoint(dir / "job.task1.xml");
    EXPECT_EQ(files_in(dir).size(), 4u);

    task.checkpoint(dir / "renamed.xml");
    EXPECT_EQ(checkpoints(dir / "renamed.xml"),
              (std::vector<checkpoint>{{"osiris", "job.task1.run1"},
                                       {"hdf5", "job.task1.run1.h5"}}));
    EXPECT_EQ(files_in(dir), (std::vector<std::string>{"job.task1.h5", "job.task1.run1",
                                                       "job.task1.run1.h5", "job.task1.xml",
                                                       "renamed.h5", "renamed.xml"}));
}

// A resumed task on one process runs its first dumped run and keeps the
// others on dump, in the order of the MCRUN elements. It creates no run, so
// its SEED is not advanced, and every run keeps the SEED it was created with.
TEST_F(TaskLayer, ResumesTheFirstDumpedRunAndKeepsTheOthersOnDump) {
    probe_simulation task(here, compose("resume", "STEPS=16;SEED=5", source_runs(3)));
    ASSERT_EQ(task.runs.size(), 3u);
    EXPECT_EQ(task.status(), (std::vector<RunStatus>{local_run, on_dump, on_dump}));
    EXPECT_EQ(std::string(task.get_parameters()["SEED"]), "5");
    for (std::size_t i = 0; i < 3; ++i) {
        EXPECT_EQ(std::string(task.run(i).parameters()["SEED"]), std::to_string(11 * (i + 1)));
        EXPECT_EQ(task.run(i).steps_done(), i + 2);
    }
}

// Work done is summed over all runs, including runs on dump, so STEPS is a
// budget for the whole task rather than for each run.
TEST_F(TaskLayer, SumsWorkOverLocalAndDumpedRuns) {
    probe_simulation task(here, compose("resume", "STEPS=16;SEED=5", source_runs(3)));
    EXPECT_DOUBLE_EQ(task.work_done(), 9. / 16.);
    task.steps(1);
    EXPECT_DOUBLE_EQ(task.work_done(), 10. / 16.);
    EXPECT_FALSE(task.finished_notime());
    task.steps(6);
    EXPECT_EQ(task.run(0).steps_done(), 9u);
    EXPECT_EQ(task.run(1).steps_done(), 3u);
    EXPECT_TRUE(task.finished_notime());
}

// The runs' measurements are pooled into one observable, not kept as separate
// runs. Runs are added in order and merged like a binary counter: two groups
// of equal size merge, and the remaining groups merge from the right.
TEST_F(TaskLayer, PoolsRunMeasurementsInRunOrder) {
    const auto sources = source_runs(7);
    const auto pooled_runs = [&](std::size_t n) {
        probe_simulation task(here, compose("pool" + std::to_string(n), "STEPS=16;SEED=5",
                                            {sources.begin(), sources.begin() + n}));
        std::vector<alps::ObservableSet> runs;
        for (std::size_t i = 0; i < n; ++i)
            runs.push_back(task.run(i).own_measurements());
        task.checkpoint(dir / ("pool" + std::to_string(n) + ".xml")); // for comparing builds
        return std::make_pair(task.get_measurements(), runs);
    };
    const auto merged = [](alps::ObservableSet a, const alps::ObservableSet& b) {
        a << b;
        return a;
    };
    for (std::size_t n : {1u, 2u, 3u, 4u, 5u, 7u}) {
        SCOPED_TRACE(n);
        const auto result = pooled_runs(n);
        const alps::ObservableSet& pooled = result.first;
        const std::vector<alps::ObservableSet>& r = result.second;
        alps::ObservableSet expected;
        switch (n) {
        case 1: expected = r[0]; break;
        case 2: expected = merged(r[0], r[1]); break;
        case 3: expected = merged(merged(r[0], r[1]), r[2]); break;
        case 4: expected = merged(merged(r[0], r[1]), merged(r[2], r[3])); break;
        case 5:
            expected = merged(merged(merged(r[0], r[1]), merged(r[2], r[3])), r[4]);
            break;
        case 7:
            expected = merged(merged(merged(r[0], r[1]), merged(r[2], r[3])),
                              merged(merged(r[4], r[5]), r[6]));
            break;
        }
        EXPECT_EQ(pooled["Value"].number_of_runs(), 1u);
        EXPECT_EQ(count(pooled, "Value"), double(n * (n + 3) / 2));
        EXPECT_EQ(mean(pooled, "Value"), mean(expected, "Value"));
        EXPECT_EQ(alps::RealObsevaluator(pooled["Value"]).error(),
                  alps::RealObsevaluator(expected["Value"]).error());
    }

    // Merging keeps each run's data in run order and computes the mean from
    // it, so the grouping does not change the result but the run order does:
    // with offsets of opposite sign, the reversed order rounds differently.
    std::vector<std::pair<fs::path, fs::path>> offset;
    const char* offsets[] = {"0", "4e15", "-3e15", "0"};
    for (int i = 0; i < 4; ++i)
        offset.push_back(source_run("offset" + std::to_string(i + 1),
                                    "STEPS=16;SEED=" + std::to_string(11 * (i + 1)) +
                                        ";OFFSET=" + offsets[i],
                                    i + 2));
    probe_simulation task(here, compose("offset", "STEPS=16;SEED=5", offset));
    std::vector<alps::ObservableSet> r;
    for (std::size_t i = 0; i < 4; ++i)
        r.push_back(task.run(i).own_measurements());
    const double binary = mean(merged(merged(r[0], r[1]), merged(r[2], r[3])), "Value");
    const double from_left = mean(merged(merged(merged(r[0], r[1]), r[2]), r[3]), "Value");
    const double reversed = mean(merged(merged(r[3], r[2]), merged(r[1], r[0])), "Value");
    EXPECT_EQ(binary, from_left);
    EXPECT_NE(binary, reversed);
    EXPECT_EQ(mean(task.get_measurements(), "Value"), binary);
    task.checkpoint(dir / "offset.xml"); // for comparing builds
}

// Observables added to the simulation itself only fill in names that no run
// has, and a task without runs returns no measurements at all.
TEST_F(TaskLayer, AddsSimulationObservablesOnlyForNamesTheRunsLack) {
    probe_simulation task(here, parameters("STEPS=4;SEED=1"));
    task.steps(2);
    alps::RealObservable extra("Extra");
    extra << 4.5;
    extra << 5.5;
    alps::RealObservable shadowed("Value");
    shadowed << 1000.;
    task << extra << shadowed;
    const alps::ObservableSet pooled = task.get_measurements();
    EXPECT_DOUBLE_EQ(mean(pooled, "Extra"), 5.);
    EXPECT_EQ(count(pooled, "Value"), 2.);
    EXPECT_DOUBLE_EQ(mean(pooled, "Value"), 0.37 + 0.15);

    probe_simulation empty(alps::ProcessList(), parameters("STEPS=4;SEED=1"));
    EXPECT_TRUE(empty.runs.empty());
    empty << extra;
    EXPECT_TRUE(empty.get_measurements().empty());
}

// A resumed task saves the pooled measurements as /simulation/results of its
// own HDF5 file. Checkpointed under its own name, it saves the local run and
// leaves the runs on dump as they are, and it resumes from there.
TEST_F(TaskLayer, SavesPooledResultsAndResumesFromItsCheckpoint) {
    probe_simulation task(here, compose("resume", "STEPS=16;SEED=5", source_runs(3)));
    const std::string dumped = read_file(dir / "resume.run2");
    task.steps(1);
    task.checkpoint(dir / "resume.xml");
    using checkpoint = std::pair<std::string, std::string>;
    EXPECT_EQ(checkpoints(dir / "resume.xml"),
              (std::vector<checkpoint>{{"osiris", "resume.run1"},
                                       {"hdf5", "resume.run1.h5"},
                                       {"osiris", "resume.run2"},
                                       {"hdf5", "resume.run2.h5"},
                                       {"osiris", "resume.run3"},
                                       {"hdf5", "resume.run3.h5"}}));
    EXPECT_EQ(read_file(dir / "resume.run2"), dumped);

    alps::ObservableSet saved;
    {
        alps::hdf5::archive ar((dir / "resume.h5").string());
        ar["/simulation/results"] >> saved;
    }
    const alps::ObservableSet pooled = task.get_measurements();
    EXPECT_EQ(count(saved, "Value"), 10.);
    EXPECT_EQ(mean(saved, "Value"), mean(pooled, "Value"));

    probe_simulation resumed(here, dir / "resume.xml");
    EXPECT_EQ(resumed.status(), (std::vector<RunStatus>{local_run, on_dump, on_dump}));
    EXPECT_EQ(resumed.run(0).steps_done(), 3u);
    EXPECT_DOUBLE_EQ(resumed.work_done(), 10. / 16.);
}

// get_and_remove_observable pools differently from get_measurements: it merges
// the runs from the left and then merges the simulation's own observable of
// that name too. A task resumed from its checkpoint holds the previously
// pooled results there, so they are counted again. Nothing in ALPS calls it.
//
// Value is unbinned here. Results loaded from an HDF5 checkpoint have bin
// size 0, since the bin size is not saved, and merging them into binned data
// divides by that bin size: SIGFPE on x86, a silent double count on arm64.
// That undefined case is not pinned.
TEST_F(TaskLayer, GetAndRemoveObservableCountsSavedResultsAgain) {
    {
        probe_simulation task(here, compose("resume", "STEPS=16;SEED=5;BINNING=0",
                                            source_runs(3, ";BINNING=0")));
        task.steps(1);
        task.checkpoint(dir / "resume.xml");
    }
    probe_simulation resumed(here, dir / "resume.xml");
    EXPECT_EQ(count(resumed.get_measurements(), "Value"), 10.);
    const alps::ObservableSet removed = resumed.get_and_remove_observable("Value");
    EXPECT_EQ(count(removed, "Value"), 20.);
    EXPECT_FALSE(resumed.get_measurements().has("Value"));
}

// Checkpointed under a new name, the runs are renamed after the new task file:
// the local run is saved and the runs on dump are copied, so the new task
// resumes. Later checkpoints into another directory copy them again.
TEST_F(TaskLayer, CopiesDumpedRunsUnderANewName) {
    probe_simulation task(here, compose("resume", "STEPS=16;SEED=5", source_runs(3)));
    task.steps(1);
    task.checkpoint(dir / "moved.xml");
    using checkpoint = std::pair<std::string, std::string>;
    EXPECT_EQ(checkpoints(dir / "moved.xml"),
              (std::vector<checkpoint>{{"osiris", "moved.run1"},
                                       {"hdf5", "moved.run1.h5"},
                                       {"osiris", "moved.run2"},
                                       {"hdf5", "moved.run2.h5"},
                                       {"osiris", "moved.run3"},
                                       {"hdf5", "moved.run3.h5"}}));
    for (const char* run : {"run2", "run3"}) {
        EXPECT_EQ(read_file(dir / ("moved." + std::string(run))),
                  read_file(dir / ("resume." + std::string(run))));
        EXPECT_EQ(read_file(dir / ("moved." + std::string(run) + ".h5")),
                  read_file(dir / ("resume." + std::string(run) + ".h5")));
    }
    probe_simulation resumed(here, dir / "moved.xml");
    EXPECT_EQ(resumed.status(), (std::vector<RunStatus>{local_run, on_dump, on_dump}));
    EXPECT_DOUBLE_EQ(resumed.work_done(), 10. / 16.);

    fs::create_directories(dir / "elsewhere");
    task.checkpoint(dir / "elsewhere" / "moved.xml");
    EXPECT_EQ(read_file(dir / "elsewhere" / "moved.run2"), read_file(dir / "resume.run2"));
    EXPECT_EQ(read_file(dir / "elsewhere" / "moved.run3.h5"),
              read_file(dir / "resume.run3.h5"));
    probe_simulation elsewhere(here, dir / "elsewhere" / "moved.xml");
    EXPECT_EQ(elsewhere.run(2).steps_done(), 4u);
    EXPECT_DOUBLE_EQ(elsewhere.work_done(), 10. / 16.);
}

// A changed parameter that the runs accept reaches every loaded run, also
// runs on dump. Their checkpoints are not rewritten and keep the old value,
// and the task file carries the new one, so the change is applied again on
// each resume.
// A parameter the runs did not have is added to them. A changed parameter
// the runs do not accept fails the task, and a changed SEED is ignored.
TEST_F(TaskLayer, AppliesChangedParametersToEveryLoadedRun) {
    const auto sources = source_runs(2);
    probe_simulation task(here, compose("longer", "STEPS=32;SEED=99", sources));
    EXPECT_EQ(std::string(task.run(0).parameters()["STEPS"]), "32");
    EXPECT_EQ(std::string(task.run(1).parameters()["STEPS"]), "32");
    EXPECT_DOUBLE_EQ(task.work_done(), 5. / 32.);
    EXPECT_EQ(std::string(task.run(1).parameters()["SEED"]), "22");

    task.checkpoint(dir / "longer.xml");
    alps::Parameters dumped;
    {
        alps::hdf5::archive ar((dir / "longer.run2.h5").string());
        ar["/parameters"] >> dumped;
    }
    EXPECT_EQ(std::string(dumped["STEPS"]), "16");
    probe_simulation again(here, dir / "longer.xml");
    EXPECT_EQ(std::string(again.run(1).parameters()["STEPS"]), "32");

    probe_simulation added(here, compose("added", "STEPS=16;SEED=5;EXTRA=1", sources));
    EXPECT_EQ(std::string(added.run(1).parameters()["EXTRA"]), "1");

    const auto fixed = source_run("fixed", "STEPS=16;SEED=11;X=1", 2);
    try {
        probe_simulation changed(here, compose("changed", "STEPS=16;SEED=11;X=2", {fixed}));
        FAIL() << "expected an exception";
    } catch (const std::runtime_error& error) {
        EXPECT_STREQ(error.what(), "Cannot change parameter X");
    }
}

// A checkpointed task has its parameters both in <task>.h5 and in the task
// file. The task file is read last, so editing it changes the parameters,
// and the change reaches every run.
TEST_F(TaskLayer, ReadsParametersEditedInTheTaskFile) {
    {
        probe_simulation task(here, compose("edit", "STEPS=16;SEED=5", source_runs(2)));
        task.steps(1);
        task.checkpoint(dir / "edit.xml");
    }
    ASSERT_TRUE(fs::exists(dir / "edit.h5"));
    std::string text = read_file(dir / "edit.xml");
    const std::string from = "<PARAMETER name=\"STEPS\">16</PARAMETER>";
    ASSERT_NE(text.find(from), std::string::npos);
    text.replace(text.find(from), from.size(), "<PARAMETER name=\"STEPS\">32</PARAMETER>");
    {
        std::ofstream out((dir / "edit.xml").string().c_str());
        out << text;
    }
    probe_simulation resumed(here, dir / "edit.xml");
    EXPECT_EQ(std::string(resumed.get_parameters()["STEPS"]), "32");
    EXPECT_EQ(std::string(resumed.run(0).parameters()["STEPS"]), "32");
    EXPECT_EQ(std::string(resumed.run(1).parameters()["STEPS"]), "32");
    EXPECT_DOUBLE_EQ(resumed.work_done(), 6. / 32.);
}

} // namespace
