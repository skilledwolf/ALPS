// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Checkpoints of the classic DMRG solver (Options::checkpoint_directory,
// resume, resume_point). A callback that throws interrupts a solve after a
// chosen iteration of the regular sweeps; the resumed solve must return
// what the uninterrupted one returns, bit for bit: energies, every
// iteration row, the measurements and the truncation error. The states are
// kept below exact, so truncation shapes the result. Bad resume input is
// rejected and leaves the checkpoint usable. The dmrg application's
// dmrg_resume test kills a run outright and resumes it.
#include "support.hpp"
#include <gtest/gtest.h>

#include <alps/dmrg/molecular.hpp>

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace classic_test;
namespace fs = std::filesystem;

namespace {

struct Directory {
    fs::path path;
    Directory() {
        std::random_device random;
        path = fs::temp_directory_path() / ("alps-classic-dmrg-checkpoint-" + std::to_string(random()));
        fs::create_directories(path);
    }
    ~Directory() {
        std::error_code ignored;
        fs::remove_all(path, ignored);
    }
};

std::string text(fs::path const& file) {
    std::ifstream in(file);
    return std::string(std::istreambuf_iterator<char>(in), {});
}

struct Stop {};

int const length = 10;
std::size_t const half = length - 4; // iterations per half sweep

mdl::Hamiltonian chain(double jz = 1.0) {
    return mdl::xxz(lat::chain(length, 1.0, lat::Boundary::open), 1.0, jz);
}

// Three sweeps below exact (16 states of 32 at the centre).
cd::Options truncated(std::size_t states = 12, std::size_t count = 1) {
    cd::Options options;
    options.sweeps = {{6, 8}, {10, 12}, {states, states}};
    options.warmup_states = 6;
    options.num_eigenvalues = count;
    options.sector = {{"Sz", 0.}};
    options.lanczos_tolerance = 1e-12;
    return options;
}

cd::Measurements requests(mdl::Hamiltonian const& model) {
    cd::Measurements m;
    m.terms["Sz"] = cd::site_terms(model, "Sz");
    m.terms["SzSz"] = cd::pair_terms(model, "Sz", "Sz");
    m.expectation["H"] = mdl::compile_hamiltonian(model).op;
    return m;
}

// Bitwise equality of everything a solve reports.
void expect_same(cd::Result const& a, cd::Result const& b) {
    EXPECT_EQ(a.energies, b.energies);
    EXPECT_EQ(a.final_truncation_error, b.final_truncation_error);
    ASSERT_EQ(a.iterations.size(), b.iterations.size());
    for (std::size_t i = 0; i < a.iterations.size(); ++i) {
        auto const& x = a.iterations[i];
        auto const& y = b.iterations[i];
        EXPECT_TRUE(x.direction == y.direction && x.index == y.index && x.energy == y.energy
                    && x.truncation_error == y.truncation_error && x.entropy == y.entropy)
            << "iteration row " << i << " differs";
    }
    ASSERT_EQ(a.measurements.size(), b.measurements.size());
    for (auto const& [name, measurement] : a.measurements) {
        auto const& other = b.measurements.at(name);
        EXPECT_EQ(measurement.values, other.values) << name;
        EXPECT_EQ(measurement.sites, other.sites) << name;
    }
}

// The rows of the warm-up (three passes): all rows less those of the regular
// sweeps and of the final sweep, which reports its right-to-left half and
// its closing configuration. With stored products (several targets and
// two-site measurements, as requests() has) that half runs one site further.
std::size_t warmup_rows(cd::Result const& reference, cd::Options const& options) {
    std::size_t const sweeps = options.sweeps.size();
    std::size_t const final_rows = half + (options.num_eigenvalues > 1 ? 1 : 0) + 1;
    std::size_t const regular = sweeps * 2 * half;
    EXPECT_GT(reference.iterations.size(), regular + final_rows);
    return reference.iterations.size() - regular - final_rows;
}

// The point after regular-sweep iteration k (from one) completes.
cd::SweepPoint after(std::size_t k, std::size_t sweeps) {
    std::size_t const sweep = (k - 1) / (2 * half) + 1;
    bool const right_to_left = ((k - 1) % (2 * half)) < half;
    std::size_t const iteration = (k - 1) % half + 1;
    if (iteration < half)
        return {sweep, right_to_left ? cd::Direction::right_to_left : cd::Direction::left_to_right, iteration + 1};
    if (right_to_left) return {sweep, cd::Direction::left_to_right, 1};
    return {sweep + 1, cd::Direction::right_to_left, 1}; // sweep sweeps + 1 is the final sweep
}

// Interrupts a checkpointed solve when the callback sees row `stop` (from one).
void interrupt(mdl::Hamiltonian const& model, cd::Options options, cd::Measurements const& measurements,
               fs::path const& directory, std::size_t stop) {
    options.checkpoint_directory = directory;
    std::size_t seen = 0;
    EXPECT_THROW(cd::solve(model, options, measurements, [&](cd::Iteration const&) {
        if (++seen == stop) throw Stop{};
    }), Stop);
}

std::string next_line(fs::path const& directory) {
    std::istringstream manifest(text(directory / "checkpoint.txt"));
    std::string line;
    while (std::getline(manifest, line))
        if (line.rfind("next ", 0) == 0) return line;
    return {};
}

std::string encoded(cd::SweepPoint const& p) {
    return "next " + std::to_string(p.sweep) + " " + (p.direction == cd::Direction::right_to_left ? "1" : "0") + " "
           + std::to_string(p.iteration);
}

template<class F>
void expect_rejection(F&& f, std::string const& message) {
    try {
        f();
        ADD_FAILURE() << "accepted; expected: " << message;
    } catch (std::invalid_argument const& error) {
        EXPECT_NE(std::string(error.what()).find(message), std::string::npos) << error.what();
    }
}

void check_resume(cd::Options const& options, std::vector<std::size_t> const& stops) {
    auto const model = chain();
    auto const measurements = requests(model);
    auto const reference = cd::solve(model, options, measurements);
    std::size_t const sweeps = options.sweeps.size();
    std::size_t const warmup = warmup_rows(reference, options);
    // The rows of the regular sweeps follow the sweep order that after() assumes.
    for (std::size_t k = 1; k <= sweeps * 2 * half; ++k) {
        auto const& row = reference.iterations.at(warmup + k - 1);
        bool const right_to_left = ((k - 1) % (2 * half)) < half;
        ASSERT_EQ(row.direction == cd::Direction::right_to_left, right_to_left) << "row " << warmup + k - 1;
        ASSERT_EQ(row.index, (k - 1) % half + 1) << "row " << warmup + k - 1;
    }
    for (auto k : stops) {
        SCOPED_TRACE("interrupted after regular iteration " + std::to_string(k));
        Directory directory;
        interrupt(model, options, measurements, directory.path, warmup + k);
        auto const point = after(k, sweeps);
        ASSERT_EQ(next_line(directory.path), encoded(point));
        auto resumed_options = options;
        resumed_options.checkpoint_directory = directory.path;
        resumed_options.resume = true;
        resumed_options.resume_point = point;
        std::size_t observed = 0;
        auto const resumed = cd::solve(model, resumed_options, measurements,
                                       [&](cd::Iteration const&) { ++observed; });
        expect_same(resumed, reference);
        EXPECT_EQ(observed, reference.iterations.size() - (warmup + k)) << "the callback saw earlier rows";
        EXPECT_TRUE(fs::is_empty(directory.path)) << "the spent checkpoint was not removed";
    }
}

} // namespace

TEST(ClassicDmrgCheckpoint, CheckpointedSolveMatchesAndCleansUp) {
    auto const model = chain();
    auto const measurements = requests(model);
    auto const options = truncated();
    auto const reference = cd::solve(model, options, measurements);
    EXPECT_GT(reference.final_truncation_error, 1e-12) << "the test needs truncation";
    Directory directory;
    auto checkpointed = options;
    checkpointed.checkpoint_directory = directory.path;
    bool files = false;
    auto const result = cd::solve(model, checkpointed, measurements, [&](cd::Iteration const&) {
        files = files || fs::exists(directory.path / "checkpoint.txt");
    });
    EXPECT_TRUE(files) << "no checkpoint was written";
    expect_same(result, reference);
    EXPECT_TRUE(fs::is_empty(directory.path));
}

// Inside a right-to-left half sweep, at its last iteration (the turn reads
// the left-hand wave function), inside left to right, at the end of a sweep,
// and at the end of the last regular sweep (only the final sweep remains).
TEST(ClassicDmrgCheckpoint, InterruptedSolveResumesBitwise) {
    check_resume(truncated(), {3, half, half + 2, 2 * half, 2 * half + 1, 2 * 3 * half});
}

// Two targets with two-site measurements store products in the blocks.
TEST(ClassicDmrgCheckpoint, StoredProductsResumeBitwise) {
    check_resume(truncated(12, 2), {half + 3, 2 * half + half});
}

TEST(ClassicDmrgCheckpoint, ComplexArithmeticAndTruncationErrorResumeBitwise) {
    auto options = truncated();
    options.arithmetic = cd::Arithmetic::complex;
    options.truncation_error = 1e-6;
    options.truncation_max_states = 12;
    check_resume(options, {4, 3 * half + 1});
}

TEST(ClassicDmrgCheckpoint, ResumeWithoutPointAndTwiceInterrupted) {
    auto const model = chain();
    auto const measurements = requests(model);
    auto const options = truncated();
    auto const reference = cd::solve(model, options, measurements);
    std::size_t const warmup = warmup_rows(reference, options);
    Directory directory;
    interrupt(model, options, measurements, directory.path, warmup + 2);
    auto again = options;
    again.resume = true;
    // Interrupted once more after resuming; the second resume takes the
    // checkpoint's own point.
    interrupt(model, again, measurements, directory.path, 5);
    EXPECT_EQ(next_line(directory.path), encoded(after(2 + 5, options.sweeps.size())));
    again.checkpoint_directory = directory.path;
    expect_same(cd::solve(model, again, measurements), reference);
    EXPECT_TRUE(fs::is_empty(directory.path));
}

TEST(ClassicDmrgCheckpoint, BadResumeInputIsRejected) {
    auto const model = chain();
    auto const measurements = requests(model);
    auto const options = truncated();
    auto const reference = cd::solve(model, options, measurements);
    std::size_t const warmup = warmup_rows(reference, options);
    std::size_t const sweeps = options.sweeps.size();

    auto resume = options;
    resume.resume = true;
    expect_rejection([&] { cd::solve(model, resume, measurements); }, "resume needs the checkpoint_directory");
    auto point_only = options;
    point_only.resume_point = cd::SweepPoint{};
    expect_rejection([&] { cd::solve(model, point_only, measurements); }, "resume_point requires resume");
    auto missing = options;
    missing.checkpoint_directory = fs::temp_directory_path() / "alps-classic-dmrg-no-such-directory";
    expect_rejection([&] { cd::solve(model, missing, measurements); }, "is not a directory");

    Directory empty;
    resume.checkpoint_directory = empty.path;
    expect_rejection([&] { cd::solve(model, resume, measurements); }, "holds no checkpoint");
    std::ofstream(empty.path / "unrelated") << "keep me";
    auto fresh = options;
    fresh.checkpoint_directory = empty.path;
    expect_rejection([&] { cd::solve(model, fresh, measurements); }, "must be empty");
    EXPECT_EQ(text(empty.path / "unrelated"), "keep me");

    Directory directory;
    std::size_t const k = half + 2;
    interrupt(model, options, measurements, directory.path, warmup + k);
    auto const point = after(k, sweeps);
    fresh.checkpoint_directory = directory.path;
    expect_rejection([&] { cd::solve(model, fresh, measurements); }, "holds the checkpoint of an interrupted solve");

    resume.checkpoint_directory = directory.path;
    auto wrong = resume;
    wrong.resume_point = cd::SweepPoint{point.sweep, point.direction, point.iteration + 1};
    expect_rejection([&] { cd::solve(model, wrong, measurements); },
                     "checkpoint continues at sweep 1, left-to-right iteration 3, not at sweep 1, "
                     "left-to-right iteration 4");
    wrong.resume_point = cd::SweepPoint{0, cd::Direction::right_to_left, 1};
    expect_rejection([&] { cd::solve(model, wrong, measurements); }, "sweep must lie between 1 and 4");
    wrong.resume_point = cd::SweepPoint{sweeps + 2, cd::Direction::right_to_left, 1};
    expect_rejection([&] { cd::solve(model, wrong, measurements); }, "sweep must lie between 1 and 4");
    wrong.resume_point = cd::SweepPoint{1, cd::Direction::right_to_left, 0};
    expect_rejection([&] { cd::solve(model, wrong, measurements); }, "iteration must lie between 1 and 6");
    wrong.resume_point = cd::SweepPoint{1, cd::Direction::right_to_left, half + 1};
    expect_rejection([&] { cd::solve(model, wrong, measurements); }, "iteration must lie between 1 and 6");
    wrong.resume_point = cd::SweepPoint{sweeps + 1, cd::Direction::left_to_right, 1};
    expect_rejection([&] { cd::solve(model, wrong, measurements); }, "must be the start of the final sweep");
    wrong.resume_point = cd::SweepPoint{1, static_cast<cd::Direction>(7), 1};
    expect_rejection([&] { cd::solve(model, wrong, measurements); }, "unknown classic DMRG resume_point direction");

    auto other_model = chain(0.5);
    expect_rejection([&] { cd::solve(other_model, resume, requests(other_model)); }, "belongs to another solve");
    auto other_states = resume;
    other_states.sweeps.back().left_to_right = 11;
    expect_rejection([&] { cd::solve(model, other_states, measurements); }, "belongs to another solve");
    auto other_measurements = measurements;
    other_measurements.terms.erase("Sz");
    expect_rejection([&] { cd::solve(model, resume, other_measurements); }, "belongs to another solve");

    cd::Evolution evolution;
    evolution.steps = 1;
    expect_rejection([&] { cd::evolve(model, model, evolution, resume, measurements); },
                     "evolve takes no checkpoint_directory");
    cd::MolecularInput molecule;
    molecule.orbitals = 2;
    molecule.alpha_electrons = molecule.beta_electrons = 1;
    molecule.one_body.assign(4, 0.);
    molecule.two_body.assign(16, 0.);
    cd::Options molecular_options;
    molecular_options.checkpoint_directory = directory.path;
    expect_rejection([&] { cd::solve_molecular(molecule, molecular_options); },
                     "solve_molecular takes no checkpoint_directory");

    // None of these touched the checkpoint, which still resumes exactly.
    resume.resume_point = point;
    expect_same(cd::solve(model, resume, measurements), reference);
    EXPECT_TRUE(fs::is_empty(directory.path));
}

// A solve that throws where it cannot be resumed removes its files. A
// killed process leaves them, with the stage it reached in the manifest,
// which a resume refuses; so does a damaged manifest.
TEST(ClassicDmrgCheckpoint, WarmupAndFinalSweepCannotBeResumed) {
    auto const model = chain();
    auto const measurements = requests(model);
    auto const options = truncated();
    auto const reference = cd::solve(model, options, measurements);
    std::size_t const warmup = warmup_rows(reference, options);
    auto resume = options;
    resume.resume = true;

    Directory during_warmup;
    interrupt(model, options, measurements, during_warmup.path, warmup);
    EXPECT_TRUE(fs::is_empty(during_warmup.path));
    Directory during_final;
    interrupt(model, options, measurements, during_final.path, reference.iterations.size() - 2);
    EXPECT_TRUE(fs::is_empty(during_final.path));

    auto killed_at = [&](fs::path const& directory, std::string const& stage) {
        interrupt(model, options, measurements, directory, warmup + 3);
        std::string manifest = text(directory / "checkpoint.txt");
        manifest.replace(manifest.find("stage sweeps"), 12, "stage " + stage);
        std::ofstream(directory / "checkpoint.txt") << manifest;
    };
    Directory killed_in_warmup;
    killed_at(killed_in_warmup.path, "warmup");
    resume.checkpoint_directory = killed_in_warmup.path;
    expect_rejection([&] { cd::solve(model, resume, measurements); }, "taken during the warm-up");
    Directory killed_in_final;
    killed_at(killed_in_final.path, "final");
    resume.checkpoint_directory = killed_in_final.path;
    expect_rejection([&] { cd::solve(model, resume, measurements); }, "taken during the final sweep");

    Directory damaged;
    interrupt(model, options, measurements, damaged.path, warmup + 3);
    std::ofstream(damaged.path / "checkpoint.txt") << "alps classic_dmrg checkpoint 1\nfingerprint zz\n";
    resume.checkpoint_directory = damaged.path;
    expect_rejection([&] { cd::solve(model, resume, measurements); }, "not readable");
}

// The first iteration of a half sweep reads the density matrix that the
// previous half sweep in its direction left on its side and then writes its
// own there. A process killed after that write but before the manifest's
// update leaves the manifest of the previous iteration with the new matrix;
// the copy taken before the iteration started restores the old one. The
// manifest saved at the turn and put back after the iteration emulates this.
TEST(ClassicDmrgCheckpoint, KilledInsideTheFirstIterationOfAHalfSweep) {
    auto const model = chain();
    auto const measurements = requests(model);
    auto const options = truncated();
    auto const reference = cd::solve(model, options, measurements);
    std::size_t const warmup = warmup_rows(reference, options);
    for (std::size_t k : {half, 2 * half}) { // before a left-to-right and a right-to-left first iteration
        SCOPED_TRACE("turn after regular iteration " + std::to_string(k));
        Directory directory;
        interrupt(model, options, measurements, directory.path, warmup + k);
        auto const point = after(k, options.sweeps.size());
        ASSERT_EQ(point.iteration, 1u);
        std::string const turn = text(directory.path / "checkpoint.txt");
        auto resume = options;
        resume.resume = true;
        interrupt(model, resume, measurements, directory.path, 1); // the first iteration completes
        ASSERT_EQ(next_line(directory.path), encoded(after(k + 1, options.sweeps.size())));
        std::ofstream(directory.path / "checkpoint.txt") << turn;
        resume.checkpoint_directory = directory.path;
        resume.resume_point = point;
        expect_same(cd::solve(model, resume, measurements), reference);
        EXPECT_TRUE(fs::is_empty(directory.path));
    }
}

// A uniform warm-up that misplaces a domain wall of fields (48 sites, the
// wall at site 12, as classic_dmrg_warmup) fails its check, and the solve
// is repeated from the product state in a private directory: checkpointed,
// and interrupted and resumed, the solve returns what the plain one does.
TEST(ClassicDmrgCheckpoint, RepeatAfterAFailedWarmupCheck) {
    int const sites = 48, wall = 12;
    auto const xx = mdl::xxz(lat::chain(sites, 1.0, lat::Boundary::open), 1.0, 0.0);
    auto op = cd::site_terms(xx, "Sz");
    for (auto& term : op.terms) term.coefficient *= term.factors[0].site < wall ? -1. : 1.;
    auto const basis = mdl::physical_basis(xx);
    cd::Measurements measurements;
    measurements.terms["Sz"] = cd::site_terms(xx, "Sz");
    cd::Options options;
    options.sweeps = cd::fixed_sweeps(4, 32);
    options.sector = {{"Sz", wall - sites / 2.}};
    auto const reference = cd::solve(op, basis, options, measurements);
    ASSERT_TRUE(reference.warmup_check && reference.warmup_check->misplaced);

    Directory directory;
    auto checkpointed = options;
    checkpointed.checkpoint_directory = directory.path;
    expect_same(cd::solve(op, basis, checkpointed, measurements), reference);
    EXPECT_TRUE(fs::is_empty(directory.path));

    std::size_t seen = 0;
    EXPECT_THROW(cd::solve(op, basis, checkpointed, measurements, [&](cd::Iteration const&) {
        if (++seen == 3 * (sites - 4) + 5) throw Stop{}; // past the warm-up's three passes
    }), Stop);
    ASSERT_TRUE(fs::exists(directory.path / "checkpoint.txt"));
    auto resume = checkpointed;
    resume.resume = true;
    expect_same(cd::solve(op, basis, resume, measurements), reference);
    EXPECT_TRUE(fs::is_empty(directory.path));
}
