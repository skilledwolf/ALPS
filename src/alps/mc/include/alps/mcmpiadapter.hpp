/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2013 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <alps/ngs/config.hpp>

#if !defined(ALPS_MCMPIADAPTER_HPP) && defined(ALPS_HAVE_MPI)
#define ALPS_MCMPIADAPTER_HPP

#include <alps/ngs.hpp>
#include <alps/ngs/boost_mpi.hpp>
#include <alps/check_schedule.hpp>
#include <alps/ngs/detail/mpi_agreement.hpp>
#include <optional>
#include <limits>

namespace alps {

    template<typename Base, typename ScheduleChecker = alps::check_schedule> class mcmpiadapter : public Base {

        public:

            mcmpiadapter(
                  typename Base::parameters_type const & parameters
                , boost::mpi::communicator const & comm
                , ScheduleChecker const & check = ScheduleChecker()
            )
                : Base(parameters, comm.rank())
                , communicator(comm)
                , schedule_checker(check)
                , clone(comm.rank())
#ifndef ALPS_NGS_USE_NEW_ALEA
                , binnumber(parameters["BINNUMBER"] | 128)
#endif
            {}

            double fraction_completed() const {
                return fraction;
            }

            bool run(boost::function<bool ()> const & stop_callback) {
                bool done = false, stopped = false;
                do {
                    this->update();
                    this->measure();
                    if (stopped || schedule_checker.pending()) {
                        stopped = stop_callback(); 
                        double local_fraction = stopped ? 1. : Base::fraction_completed();
                        schedule_checker.update(fraction = boost::mpi::all_reduce(communicator, local_fraction, std::plus<double>()));
                        done = fraction >= 1.;
                    }
                } while(!done);
                return !stopped;
            }

            typename Base::results_type collect_results() const {
                return collect_results(this->result_names());
            }

            typename Base::results_type collect_results(typename Base::result_names_type const & names) const {
#ifdef ALPS_NGS_USE_NEW_ALEA
                typename Base::results_type partial_results;
                for(typename Base::result_names_type::const_iterator it = names.begin(); it != names.end(); ++it) {
                        if (communicator.rank() == 0) {
                            if (this->measurements[*it].count()) {
                                typename Base::observable_collection_type::value_type merged = this->measurements[*it];
                                merged.collective_merge(communicator, 0);
                                partial_results.insert(*it, merged.result());
                            } else
                                partial_results.insert(*it, this->measurements[*it].result());
                        } else if (this->measurements[*it].count())
                            this->measurements[*it].collective_merge(communicator, 0);
                }
                return partial_results;
#else
                // The retained reducer assumes identical observable layouts and
                // cannot represent an empty contribution. Agree participation
                // before entering it: a local count()/bin branch can deadlock.
                alps::detail::mpi_group group(communicator);
                auto agreed_names = names;
                boost::mpi::broadcast(communicator, agreed_names, 0);
                group.checked_local([&] {
                    if (names != agreed_names)
                        throw std::invalid_argument("MPI result names/order must agree on every rank");
                });
                typename Base::results_type partial_results;
                for (const auto &name : names) {
                    std::optional<alps::mcresult> result;
                    bool populated = false;
                    group.checked_local([&] {
                        result.emplace(this->measurements[name]);
                        populated = result->count() != 0;
                    });
                    const bool any_populated = group.any(populated);
                    const bool all_populated = group.all(populated);
                    group.checked_local([&] {
                        if (any_populated != all_populated)
                            throw std::invalid_argument("legacy MPI reduction requires measurements on every rank: " + name);
                    });
                    if (!any_populated) {
                        group.checked_local([&] { partial_results.insert(name, *result); });
                        continue;
                    }

                    // Type/extent and optional statistics control the reducer's
                    // collective sequence and message lengths. Counts may differ.
                    std::vector<std::uint64_t> layout;
                    group.checked_local([&] {
                        if (binnumber == 0 || binnumber > std::uint64_t(std::numeric_limits<int>::max()))
                            throw std::invalid_argument("legacy MPI BINNUMBER must be positive and fit an int");
                        std::uint64_t kind = 0, width = 0;
                        if (result->template is_type<double>()) {
                            kind = 1;
                            width = 1;
                            (void)result->template mean<double>();
                            (void)result->template error<double>();
                        } else if (result->template is_type<std::vector<double>>()) {
                            kind = 2;
                            width = result->template mean<std::vector<double>>().size();
                            if (result->template error<std::vector<double>>().size() != width)
                                throw std::invalid_argument("legacy MPI result error extent differs: " + name);
                            if (result->has_variance() && result->template variance<std::vector<double>>().size() != width)
                                throw std::invalid_argument("legacy MPI result variance extent differs: " + name);
                            if (result->has_tau() && result->template tau<std::vector<double>>().size() != width)
                                throw std::invalid_argument("legacy MPI result tau extent differs: " + name);
                            for (const auto &bin : result->template bins<std::vector<double>>())
                                if (bin.size() != width)
                                    throw std::invalid_argument("legacy MPI result bin extent differs: " + name);
                        } else {
                            throw std::invalid_argument("unsupported legacy MPI result type: " + name);
                        }
                        layout = {kind, width, result->has_variance(), result->has_tau(),
                                  result->bin_number() != 0, binnumber};
                    });
                    auto agreed_layout = layout;
                    boost::mpi::broadcast(communicator, agreed_layout, 0);
                    group.checked_local([&] {
                        if (layout != agreed_layout)
                            throw std::invalid_argument("legacy MPI result layouts must agree: " + name);
                    });
                    // This call itself is collective and must stay outside
                    // checked_local. MPI/process failures remain unrecoverable.
                    auto reduced = result->reduce(communicator, binnumber);
                    group.checked_local([&] { partial_results.insert(name, reduced); });
                }
                return partial_results;
#endif
            }

        protected:

            boost::mpi::communicator communicator;

            ScheduleChecker schedule_checker;
            double fraction;
            int clone;
#ifndef ALPS_NGS_USE_NEW_ALEA
            std::size_t binnumber;
#endif
    };
}

#endif
