// Copyright (C) 2010-2013, 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#pragma once
#include <alps/ngs/config.hpp>
#ifdef ALPS_HAVE_MPI
#include <alps/alea/mpi.hpp>
#include <alps/ngs/boost_mpi.hpp>
#include <alps/check_schedule.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <functional>
#include <limits>

namespace alps {

template<class Base, class ScheduleChecker = alps::check_schedule>
class mcmpiadapter : public Base {
public:
    mcmpiadapter(typename Base::parameters_type const& parameters,
                 boost::mpi::communicator const& comm,
                 ScheduleChecker const& check = ScheduleChecker())
        : Base(parameters, comm.rank()), communicator(comm), schedule_checker(check) {}

    double fraction_completed() const { return fraction; }

    bool run(std::function<bool()> const& stop_callback) {
        bool initial = true, stopped = false;
        std::string error;
        for (;;) {
            if (initial || !error.empty() || schedule_checker.pending()) {
                agree_failure(error);
                double local_fraction = 0;
                try {
                    stopped = stop_callback();
                    local_fraction = stopped ? 1. : Base::fraction_completed();
                } catch (std::exception const& failure) { error = failure.what(); }
                catch (...) { error = "MC progress callback failed"; }
                agree_failure(error);
                fraction = boost::mpi::all_reduce(communicator, local_fraction, std::plus<double>());
                stopped = boost::mpi::all_reduce(communicator, stopped, std::logical_or<bool>());
                schedule_checker.update(fraction);
                if (fraction >= 1. || stopped) return !stopped;
                initial = false;
            }
            try { this->update(); this->measure(); }
            catch (std::exception const& failure) { error = failure.what(); }
            catch (...) { error = "MC update or measurement failed"; }
        }
    }

    typename Base::results_type collect_results() const {
        return collect_results(this->result_names());
    }

    typename Base::results_type collect_results(typename Base::result_names_type const& names) const {
        auto root_names = names;
        boost::mpi::broadcast(communicator, root_names, 0);
        agree_failure(names == root_names ? "" : "MC result names must agree on every rank");
        typename Base::results_type results;
        std::string error;
        try { results = Base::collect_results(names); }
        catch (std::exception const& failure) { error = failure.what(); }
        catch (...) { error = "MC result snapshot failed"; }
        agree_failure(error);
        alps::alea::mpi_reducer reduction(communicator);
        for (auto& entry : results) entry.second.reduce(reduction);
        if (communicator.rank()) results.clear();
        return results;
    }

protected:
    boost::mpi::communicator communicator;
    ScheduleChecker schedule_checker;
    double fraction = 0.;

private:
    void agree_failure(std::string const& error) const {
        auto origin = boost::mpi::all_reduce(communicator,
            error.empty() ? std::numeric_limits<int>::max() : communicator.rank(),
            boost::mpi::minimum<int>());
        if (origin == std::numeric_limits<int>::max()) return;
        auto message = error;
        boost::mpi::broadcast(communicator, message, origin);
        throw std::runtime_error(message);
    }
};
}
#endif
