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

#ifndef ALPS_NGS_MCBASE_HPP
#define ALPS_NGS_MCBASE_HPP

#include <alps/ngs.hpp>

#include <boost/filesystem/path.hpp>

#include <alps/alea/result.hpp>
#include <functional>
#include <map>
#include <memory>
#include <vector>
#include <string>

namespace alps {

    class ALPS_DECL mcbase {

        public:

            using observable_type = std::variant<
                std::shared_ptr<alea::mean_acc<double>>,
                std::shared_ptr<alea::mean_acc<std::complex<double>>>,
                std::shared_ptr<alea::var_acc<double>>,
                std::shared_ptr<alea::var_acc<std::complex<double>>>,
                std::shared_ptr<alea::var_acc<std::complex<double>,alea::elliptic_var>>,
                std::shared_ptr<alea::cov_acc<double>>,
                std::shared_ptr<alea::cov_acc<std::complex<double>>>,
                std::shared_ptr<alea::cov_acc<std::complex<double>,alea::elliptic_var>>,
                std::shared_ptr<alea::autocorr_acc<double>>,
                std::shared_ptr<alea::autocorr_acc<std::complex<double>>>,
                std::shared_ptr<alea::batch_acc<double>>,
                std::shared_ptr<alea::batch_acc<std::complex<double>>>>;
            using observable_collection_type = std::map<std::string,observable_type>;

            typedef alps::params parameters_type;
            typedef std::vector<std::string> result_names_type;

            using results_type = std::map<std::string, alea::result::variant_type>;

            mcbase(parameters_type const & parms, std::size_t seed_offset = 0);
            virtual ~mcbase();

            virtual void update() = 0;
            virtual void measure() = 0;
            virtual double fraction_completed() const = 0;
            bool run(std::function<bool ()> const & stop_callback);

            result_names_type result_names() const;
            results_type collect_results() const;
            results_type collect_results(result_names_type const & names) const;

            template<class A=alea::batch_acc<double>> auto& measurement(std::string const& name) {
                return std::get<std::shared_ptr<A>>(measurements.at(name));
            }
            template<class A=alea::batch_acc<double>> auto const& measurement(std::string const& name) const {
                return std::get<std::shared_ptr<A>>(measurements.at(name));
            }
            template<class R> std::map<std::string,R> collect_results_as(result_names_type const& names = {}) const {
                auto all = names.empty() ? collect_results() : collect_results(names);
                std::map<std::string,R> typed;
                for (auto& [name,value]:all) typed.emplace(name,std::get<R>(std::move(value)));
                return typed;
            }

            void save(boost::filesystem::path const & filename) const;
            void load(boost::filesystem::path const & filename);
            virtual void save(alps::hdf5::archive & ar) const;
            virtual void load(alps::hdf5::archive & ar);

            // Non-virtual accessors for language bindings and downstream
            // exporters. Keeping these on the actual base class avoids
            // assuming that every derived simulation is a Python trampoline.
            alps::random01 & get_random() { return random; }
            parameters_type & get_parameters() { return parameters; }
            observable_collection_type & get_measurements() { return measurements; }

        protected:

            parameters_type parameters;
            alps::random01 mutable random;
            observable_collection_type measurements;
    };

}

#endif
