/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2012 by Lukas Gamper <gamperl@gmail.com>,                  *
 *                              Synge Todo <wistaria@comp-phys.org>                *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_NGS_CALLBACK_HPP
#define ALPS_NGS_CALLBACK_HPP

#include <alps/ngs/config.hpp>
#include <alps/ngs/signal.hpp>

#include <boost/chrono.hpp>

namespace alps {

	class ALPS_DECL stop_callback {
		public:
		    stop_callback(std::size_t timelimit);
		    bool operator()();
		private:
		    boost::chrono::duration<std::size_t> limit;
		    alps::ngs::signal signals;
		    boost::chrono::high_resolution_clock::time_point start;
	};

}

#endif
