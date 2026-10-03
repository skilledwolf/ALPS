/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_NGS_SIGNAL_HPP
#define ALPS_NGS_SIGNAL_HPP

#include <alps/ngs/config.hpp>
#include <alps/utilities_export.h>

#include <boost/array.hpp>

namespace alps {
  namespace ngs {

    class ALPS_UTILITIES_DECL signal{

    public:

      /*!
Listen to the following posix signals SIGINT, SIGTERM, SIGXCPU, SIGQUIT,
SIGUSR1 and SIGUSR2. These signals can be checked by empty, top and pop.
      */
      signal();

      /*!
Returns if a signal has been captured.
      */
      bool empty();

      /*!
Returns the last signal that has been captured .
      */
      int top();


      /*!
Pops a signal form the stack.
       */
      void pop();

      static void slot(int signal);

    private:

      static std::size_t begin_;
      static std::size_t end_;
      static boost::array<int, 0x20> signals_;
    };
  }
}

#endif
