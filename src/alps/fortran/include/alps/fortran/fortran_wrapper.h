// Copyright (C) 2011 Synge Todo; 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include <alps/mcbase.hpp>
#include <exception>

namespace alps {
// Each callback receives this instance as an opaque ISO_C_BINDING C_PTR.
// Application state belongs in context, never in shared Fortran module globals.
class fortran_wrapper : public mcbase {
public:
    fortran_wrapper(params const&, std::size_t bins, std::size_t chain);
    ~fortran_wrapper() override;
    fortran_wrapper(fortran_wrapper const&) = delete;
    fortran_wrapper& operator=(fortran_wrapper const&) = delete;
    void update() override;
    void measure() override {} // Fortran records observations during alps_run.
    double fraction_completed() const override;
    std::uint64_t completed_sweeps() const { return steps_; }
    static params checkpoint_parameters(params p) { p.erase("SWEEPS"); return p; }
    void save(hdf5::archive&) const override;
    void load(hdf5::archive&) override;
    void* context = nullptr;
    bool collecting() const { return collecting_; }
    std::size_t bins() const { return bins_; }
    std::string field(int type, std::size_t count, bool writing);
    hdf5::archive& archive() { return *archive_; }
    template<class F> void guard(F f) noexcept {
        if (!error_) try { f(); } catch (...) { error_ = std::current_exception(); }
    }
    bool failed() const { return bool(error_); }
    void check() const { if (error_) std::rethrow_exception(error_); }
    static int main(int argc, char** argv, char const* application, char const* schema);
private:
    std::size_t bins_, chain_;
    std::uint64_t steps_ = 0;
    bool collecting_ = false;
    mutable hdf5::archive* archive_ = nullptr;
    mutable std::size_t field_ = 0;
    mutable bool writing_ = false;
    mutable std::exception_ptr error_;
};
}
