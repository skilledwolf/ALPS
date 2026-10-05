/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2013 by Lukas Gamper <gamperl@gmail.com>                   *
 *                              Matthias Troyer <troyer@comp-phys.org>             *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_NGS_RANDOM01_HPP
#define ALPS_NGS_RANDOM01_HPP

#include <alps/hdf5/archive.hpp>

#include <boost/random.hpp>

#include <string>
#include <sstream>
#include <stdexcept>
#include <variant>
#include <cstdint>

namespace alps {

    class random01 {
        // Use the exact engine registered by the released ALPS RNG factory.
        using fibonacci = boost::random::lagged_fibonacci<uint32_t,48,607,273>;
        std::variant<boost::mt19937,fibonacci> engine_;
    public:
        using result_type = double;
        random01(int seed = 42, std::string const& name = "mt19937") {
            if (name == "lagged_fibonacci607") engine_.emplace<fibonacci>();
            else if (name != "mt19937") throw std::invalid_argument("Unknown RNG: " + name);
            this->seed(seed);
        }
        std::string name() const { return engine_.index() ? "lagged_fibonacci607" : "mt19937"; }
        void seed(uint32_t value) { std::visit([&](auto& engine) { engine.seed(value); }, engine_); }
        template<class F> decltype(auto) with_engine(F&& f) { return std::visit(std::forward<F>(f),engine_); }
        double operator()() {
            return with_engine([](auto& engine) { return boost::uniform_01<double>()(engine); });
        }
        friend bool operator==(random01 const& a, random01 const& b) { return a.engine_ == b.engine_; }
        void save(alps::hdf5::archive & ar) const {
            std::ostringstream os;
            std::visit([&](auto const& engine) { os << engine; },engine_);
            ar["name"] << name();
            ar["engine"] << os.str();
        }
        void load(alps::hdf5::archive & ar) {
            std::string state, type;
            ar["name"] >> type;
            ar["engine"] >> state;
            random01 restored(0,type);
            std::istringstream is(state + " ");
            restored.with_engine([&](auto& engine) {
                if (!(is >> engine) || !is.eof())
                    throw std::runtime_error("invalid random01 checkpoint");
            });
            *this = std::move(restored);
        }
    };

}

#endif 
