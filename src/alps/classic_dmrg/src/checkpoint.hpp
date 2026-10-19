// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// The checkpoint of a solve (Options::checkpoint_directory and resume),
// written here, not taken from alpsx or from the dmrg application's former
// copy of dmtk. That copy read START_SWEEP, START_DIR and START_ITER and
// called System::resume on the scratch files of an earlier run, which could
// not work from a new process: every scratch file had an mkstemp suffix
// known only to the process that made it, so the new process read fresh
// empty files, and System::resume handed START_DIR and START_ITER to every
// remaining sweep. Here the files keep fixed names (dmtk::FileList with a
// registry), this manifest records where the solve continues, and the
// driver in classic.cpp starts its sweeps there.
//
// The directory holds dmtk's files, the registry of their names
// (checkpoint-files.txt), the iterations completed so far
// (checkpoint-iterations.txt, one per line, appended), a copy of the file
// that the first iteration of a half sweep reads and then overwrites
// (checkpoint-turn-<sweep>-<direction>.dat, below) and the manifest
// (checkpoint.txt, replaced by a rename after every iteration):
//
//   alps classic_dmrg checkpoint 1
//   fingerprint <16 hex digits of the problem and the options>
//   stage warmup | sweeps | final
//   next <sweep> <direction 0 left-to-right, 1 right-to-left> <iteration>
//   iterations <the leading lines of checkpoint-iterations.txt that count>
//   random <real> <complex>   (dmtk's Lanczos start-vector streams)
//
// Only stage `sweeps` can be resumed: its blocks, density matrices and
// wave functions on disk are those after the last completed iteration,
// which the next iteration reads. An iteration overwrites none of the files
// it reads, except the first of a half sweep: it transforms its seed with
// the density matrix of the two-site block on its side from the previous
// half sweep in its direction, and then writes its own there. That file is
// copied before the iteration starts and put back when a solve resumes
// there, so that a process ended between that write and the manifest's
// update resumes exactly. The warm-up has not written every file yet, and
// the final sweep's own start is such an iteration, so neither is resumed.
#pragma once

#include <alps/dmrg/classic.hpp>

#include <algorithm>
#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

namespace alps::dmrg::classic::detail {

// FNV-1a over a canonical byte stream: stable across processes and runs of
// one build, unlike std::hash.
class Fingerprint {
public:
    void add(std::uint64_t v) {
        for (int i = 0; i < 8; ++i) byte(static_cast<unsigned char>(v >> (8 * i)));
    }
    void add(double v) {
        std::uint64_t bits = 0;
        static_assert(sizeof bits == sizeof v);
        std::memcpy(&bits, &v, sizeof v);
        add(bits);
    }
    void add(std::string const& s) {
        add(static_cast<std::uint64_t>(s.size()));
        for (char c : s) byte(static_cast<unsigned char>(c));
    }
    std::uint64_t value() const noexcept { return value_; }

private:
    void byte(unsigned char c) {
        value_ ^= c;
        value_ *= 1099511628211ull;
    }
    std::uint64_t value_ = 1469598103934665603ull;
};

class Checkpoint {
public:
    // Opens `options.checkpoint_directory` for a solve of `num_sites` sites
    // and the fingerprint of its problem and options: an empty directory to
    // start in, or, with options.resume, the checkpoint to continue. Writes
    // nothing; begin() does, once the invocation owns the engine.
    Checkpoint(Options const& options, std::uint64_t fingerprint, int num_sites)
        : sweeps_(options.sweeps.size()), half_(static_cast<std::size_t>(num_sites - 4)),
          fingerprint_(fingerprint), resuming_(options.resume) {
        namespace fs = std::filesystem;
        std::error_code error;
        directory_ = fs::absolute(options.checkpoint_directory, error);
        if (error || !fs::is_directory(directory_))
            throw std::invalid_argument("classic DMRG checkpoint_directory is not a directory: "
                                        + options.checkpoint_directory.string());
        bool const held = fs::exists(manifest());
        if (!resuming_) {
            if (held)
                throw std::invalid_argument("classic DMRG checkpoint_directory holds the checkpoint of an "
                                            "interrupted solve; resume it, or remove its files to start afresh: "
                                            + directory_.string());
            if (!fs::is_empty(directory_))
                throw std::invalid_argument("classic DMRG checkpoint_directory must be empty to start a solve: "
                                            + directory_.string());
            return;
        }
        if (!held)
            throw std::invalid_argument("classic DMRG checkpoint_directory holds no checkpoint to resume: "
                                        + directory_.string());
        read_manifest();
        if (stage_ == "warmup")
            throw std::invalid_argument("classic DMRG checkpoint was taken during the warm-up, which cannot be "
                                        "resumed; remove its files and start afresh");
        if (stage_ == "final")
            throw std::invalid_argument("classic DMRG checkpoint was taken during the final sweep, which rewrites "
                                        "the files a resume needs; remove its files and start afresh");
        if (options.resume_point) {
            auto const& point = *options.resume_point;
            if (point.sweep < 1 || point.sweep > sweeps_ + 1)
                throw std::invalid_argument("classic DMRG resume_point sweep must lie between 1 and "
                                            + std::to_string(sweeps_ + 1) + " (the final sweep)");
            if (point.direction != Direction::left_to_right && point.direction != Direction::right_to_left)
                throw std::invalid_argument("unknown classic DMRG resume_point direction");
            if (point.iteration < 1 || point.iteration > std::max<std::size_t>(half_, 1))
                throw std::invalid_argument("classic DMRG resume_point iteration must lie between 1 and "
                                            + std::to_string(std::max<std::size_t>(half_, 1)));
            if (point.sweep == sweeps_ + 1 && (point.direction != Direction::right_to_left || point.iteration != 1))
                throw std::invalid_argument("classic DMRG resume_point in sweep " + std::to_string(sweeps_ + 1)
                                            + " must be the start of the final sweep (right to left, iteration 1)");
            if (point != next_)
                throw std::invalid_argument("classic DMRG checkpoint continues at " + describe(next_) + ", not at "
                                            + describe(point));
        }
        read_rows();
    }

    // A fresh solve marks the directory as its checkpoint (stage warmup); a
    // resumed one drops rows past the manifest's count and puts back the
    // copy of the file that its first iteration reads.
    void begin() {
        namespace fs = std::filesystem;
        if (!resuming_) return write_manifest();
        std::ofstream out(iterations(), std::ios::trunc);
        for (auto const& row : rows_) out << line(row);
        if (!out.flush()) fail("rewrite the iterations of");
        if (next_.iteration == 1 && next_.sweep <= sweeps_) {
            auto const copy = turn_copy(next_.sweep, next_.direction);
            std::error_code error;
            if (fs::exists(copy)) {
                fs::copy_file(copy, directory_ / turn_file(next_.direction), fs::copy_options::overwrite_existing,
                              error);
                if (error) fail("restore a density matrix of");
                turn_ = copy;
            }
        }
    }

    bool resuming() const noexcept { return resuming_; }
    // A solve that ends here can be resumed.
    bool resumable() const noexcept { return stage_ == "sweeps"; }
    SweepPoint const& next() const noexcept { return next_; }
    std::vector<Iteration> const& rows() const noexcept { return rows_; }
    long real_random() const noexcept { return random_[0]; }
    long complex_random() const noexcept { return random_[1]; }
    std::filesystem::path const& directory() const noexcept { return directory_; }
    std::filesystem::path registry() const { return directory_ / "checkpoint-files.txt"; }

    // The warm-up is complete: the sweeps start at sweep 1.
    void after_warmup(std::vector<Iteration> const& rows, long real_random, long complex_random) {
        record(half_ ? SweepPoint{1, Direction::right_to_left, 1} : final_start(), rows, real_random,
               complex_random);
    }

    // The first iteration of a half sweep of the regular sweeps starts and
    // has read nothing yet: copy the density matrix that it reads and then
    // overwrites, by a rename, and drop the previous copy.
    void before_first_iteration(std::size_t sweep, Direction direction) {
        namespace fs = std::filesystem;
        auto const source = directory_ / turn_file(direction);
        if (!fs::exists(source)) return;
        auto const copy = turn_copy(sweep, direction), staged_copy = directory_ / "checkpoint-turn.new";
        std::error_code error;
        fs::copy_file(source, staged_copy, fs::copy_options::overwrite_existing, error);
        if (!error) fs::rename(staged_copy, copy, error);
        if (error) fail("copy a density matrix into");
        if (!turn_.empty() && turn_ != copy) fs::remove(turn_, error);
        turn_ = copy;
    }

    // The iteration (sweep, direction, iteration) of the regular sweeps is
    // complete and its files are written.
    void after(std::size_t sweep, Direction direction, std::size_t iteration, std::vector<Iteration> const& rows,
               long real_random, long complex_random) {
        SweepPoint next{sweep, direction, iteration + 1};
        if (iteration >= half_) {
            if (direction == Direction::right_to_left) next = {sweep, Direction::left_to_right, 1};
            else next = sweep < sweeps_ ? SweepPoint{sweep + 1, Direction::right_to_left, 1} : final_start();
        }
        record(next, rows, real_random, complex_random);
    }

    // The final sweep starts; from here on nothing can be resumed.
    void before_final_sweep() {
        stage_ = "final";
        write_manifest();
    }

    // Removes the manifest, the iteration list and the copies (dmtk's files
    // go with tmp_files.remove_registered()).
    void remove() noexcept {
        namespace fs = std::filesystem;
        std::error_code ignored;
        fs::remove(manifest(), ignored);
        fs::remove(staged(), ignored);
        fs::remove(iterations(), ignored);
        fs::remove(directory_ / "checkpoint-turn.new", ignored);
        std::vector<fs::path> copies;
        for (fs::directory_iterator it(directory_, ignored), end; !ignored && it != end; it.increment(ignored))
            if (it->path().filename().string().rfind("checkpoint-turn-", 0) == 0) copies.push_back(it->path());
        for (auto const& copy : copies) fs::remove(copy, ignored);
    }

    static std::string describe(SweepPoint const& point) {
        return "sweep " + std::to_string(point.sweep) + ", "
               + (point.direction == Direction::right_to_left ? "right-to-left" : "left-to-right") + " iteration "
               + std::to_string(point.iteration);
    }

private:
    std::filesystem::path manifest() const { return directory_ / "checkpoint.txt"; }
    std::filesystem::path staged() const { return directory_ / "checkpoint.txt.new"; }
    std::filesystem::path iterations() const { return directory_ / "checkpoint-iterations.txt"; }

    SweepPoint final_start() const { return {sweeps_ + 1, Direction::right_to_left, 1}; }

    // dmtk::System::write_rho's name for the density matrix of the two-site
    // block on the side a half sweep grows, of classic.cpp's System "ALPS".
    static std::string turn_file(Direction direction) {
        return direction == Direction::right_to_left ? "rho_ALPS_2_r.dat" : "rho_ALPS_2_l.dat";
    }
    std::filesystem::path turn_copy(std::size_t sweep, Direction direction) const {
        return directory_ / ("checkpoint-turn-" + std::to_string(sweep) + "-"
                             + (direction == Direction::right_to_left ? "1" : "0") + ".dat");
    }

    // Hexadecimal floating point, which strtod reads back exactly.
    static std::string line(Iteration const& row) {
        char text[160];
        std::snprintf(text, sizeof text, "%d %zu %a %a %a\n", row.direction == Direction::right_to_left ? 1 : 0,
                      row.index, row.energy, row.truncation_error, row.entropy);
        return text;
    }

    void record(SweepPoint const& next, std::vector<Iteration> const& rows, long real_random, long complex_random) {
        if (rows.size() < written_) throw std::logic_error("classic DMRG checkpoint lost iterations");
        if (rows.size() > written_) {
            std::ofstream out(iterations(), std::ios::app);
            for (std::size_t i = written_; i < rows.size(); ++i) out << line(rows[i]);
            if (!out.flush()) fail("write the iterations of");
            written_ = rows.size();
        }
        stage_ = "sweeps";
        next_ = next;
        random_[0] = real_random;
        random_[1] = complex_random;
        write_manifest();
    }

    void write_manifest() {
        {
            std::ofstream out(staged(), std::ios::trunc);
            char hex[17];
            std::snprintf(hex, sizeof hex, "%016llx", static_cast<unsigned long long>(fingerprint_));
            out << "alps classic_dmrg checkpoint 1\n"
                << "fingerprint " << hex << '\n'
                << "stage " << stage_ << '\n'
                << "next " << next_.sweep << ' ' << (next_.direction == Direction::right_to_left ? 1 : 0) << ' '
                << next_.iteration << '\n'
                << "iterations " << written_ << '\n'
                << "random " << random_[0] << ' ' << random_[1] << '\n';
            if (!out.flush()) fail("write");
        }
        std::error_code error;
        std::filesystem::rename(staged(), manifest(), error);
        if (error) fail("replace the manifest of");
    }

    [[noreturn]] void fail(char const* what) const {
        throw std::runtime_error(std::string("classic DMRG could not ") + what + " the checkpoint in "
                                 + directory_.string());
    }

    [[noreturn]] void unreadable() const {
        throw std::invalid_argument("classic DMRG checkpoint is not readable or was written by another version: "
                                    + manifest().string());
    }

    void read_manifest() {
        std::ifstream in(manifest());
        std::string line, key, hex;
        if (!std::getline(in, line) || line != "alps classic_dmrg checkpoint 1") unreadable();
        int direction = -1;
        unsigned long long stored = 0;
        bool complete = (in >> key >> hex) && key == "fingerprint" && (in >> key >> stage_) && key == "stage"
                        && (in >> key >> next_.sweep >> direction >> next_.iteration) && key == "next"
                        && (in >> key >> written_) && key == "iterations"
                        && (in >> key >> random_[0] >> random_[1]) && key == "random";
        if (!complete || hex.size() != 16 || (direction != 0 && direction != 1)
            || (stage_ != "warmup" && stage_ != "sweeps" && stage_ != "final"))
            unreadable();
        char* end = nullptr;
        errno = 0;
        stored = std::strtoull(hex.c_str(), &end, 16);
        if (errno || end != hex.c_str() + hex.size()) unreadable();
        next_.direction = direction ? Direction::right_to_left : Direction::left_to_right;
        if (stored != fingerprint_)
            throw std::invalid_argument("classic DMRG checkpoint belongs to another solve: the Hamiltonian, "
                                        "physical basis, measurements or options (sweeps, states, sector, "
                                        "tolerances, arithmetic, num_eigenvalues, warm-up, keep_mps) differ");
        bool const valid_point = stage_ != "sweeps"
                                 || (next_.sweep >= 1 && next_.iteration >= 1
                                     && ((next_.sweep <= sweeps_ && next_.iteration <= half_)
                                         || next_ == final_start()));
        if (!valid_point) unreadable();
    }

    // The leading `written_` lines; a process ended between appending and
    // renaming the manifest leaves more, which are dropped.
    void read_rows() {
        std::ifstream in(iterations());
        std::vector<Iteration> rows;
        std::string text;
        while (rows.size() < written_ && std::getline(in, text)) {
            std::istringstream fields(text);
            int direction = -1;
            std::string index, energy, truncation, entropy;
            if (!(fields >> direction >> index >> energy >> truncation >> entropy) || (direction != 0 && direction != 1))
                unreadable();
            auto number = [&](std::string const& text) {
                char* end = nullptr;
                double const value = std::strtod(text.c_str(), &end);
                if (end != text.c_str() + text.size()) unreadable();
                return value;
            };
            char* end = nullptr;
            unsigned long long const at = std::strtoull(index.c_str(), &end, 10);
            if (end != index.c_str() + index.size()) unreadable();
            rows.push_back(Iteration{direction ? Direction::right_to_left : Direction::left_to_right,
                                     static_cast<std::size_t>(at), number(energy), number(truncation),
                                     number(entropy)});
        }
        if (rows.size() != written_) unreadable();
        rows_ = std::move(rows);
    }

    std::filesystem::path directory_;
    std::size_t sweeps_;
    std::size_t half_; // iterations per half sweep, L - 4
    std::uint64_t fingerprint_;
    std::string stage_ = "warmup";
    SweepPoint next_{};
    std::size_t written_ = 0;
    long random_[2] = {42, 42};
    std::vector<Iteration> rows_;
    bool resuming_;
    std::filesystem::path turn_; // the current copy, if any
};

} // namespace alps::dmrg::classic::detail
