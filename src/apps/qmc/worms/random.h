/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2001-2004 by Matthias Troyer <troyer@comp-phys.org>,
*                            Simon Trebst <trebst@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <limits>
#include <cmath>

#ifndef ALPS_APPLICATIONS_WORM_RANDOM_H
#define ALPS_APPLICATIONS_WORM_RANDOM_H

class new_finite_exponential_random {
public:
  new_finite_exponential_random(double x, double l, double t)
   : random_x(x), lambda(l), time(t) {}

  double operator()() {
    if(fabs(lambda*time)<1e-8)
      return time*random_x;
    else if(lambda<0)
    {   
      if(lambda*time > std::log(std::numeric_limits<double>::min()))
      { 
        double factor=std::exp(lambda*time);
        return 1./lambda*std::log(factor+(1.-factor)*random_x);
      }   
      else
      {
        double t= 1./lambda*std::log(random_x);
        return (t < time ? t : time-std::numeric_limits<double>::epsilon());
       }
    }
    else 
    {
      if(-lambda*time > std::log(std::numeric_limits<double>::min()))
      {
        double factor=std::exp(-lambda*time);
        return time - 1./(-lambda)*std::log(factor+(1.-factor)*random_x);
        }
      else
      {
        double t = time -1./(-lambda)*std::log(random_x);
        return (t>0 ? t : std::numeric_limits<double>::epsilon());
      } 
    }     
  }

private:
  double random_x;
  double lambda;
  double time; 
};   // new_finite_exponential_random

#endif
