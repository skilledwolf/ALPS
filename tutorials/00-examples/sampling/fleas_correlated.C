/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2007 - 2010 by Matthias Troyer <troyer@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include "fleas.h"
#include <alps/alea/autocorr.hpp>
#include <alps/alea/computed.hpp>
#include <boost/random.hpp>
#include <iostream>
#include <vector>


int main()
{
  const int N=50; // total number of fleas
  
  int M; // number of hops
  std::cout << "How many hops? ";
  std::cin >> M;
  unsigned int seed;
  std::cout << "Random number seed? ";
  std::cin >> seed;
  int n=N; // all fleas on left dog
    
  typedef boost::mt19937 engine_type;
  typedef boost::uniform_int<> dist_type;
  
  engine_type engine(seed);
  dist_type dist(1,N);
  boost::variate_generator<engine_type,dist_type> rng(engine,dist);
  
  // Successive hops are correlated; the autocorrelation accumulators
  // estimate errors and integrated autocorrelation times by binning.
  std::vector<double> current_distribution(N+1,0.);
  alps::alea::autocorr_acc<double> histogram(N+1), number;
  
  // equilibration
  for (int i=0;i<M/5;++i) {
    if (rng() <= n )
     --n;
    else
      ++n;
  }
  
  for (int i=0;i<M;++i) {
    if (rng() <= n )
     --n;
    else
      ++n;

    current_distribution[n]=1;
    histogram << alps::alea::make_adapter(current_distribution);
    number << alps::alea::make_adapter(double(n));
    current_distribution[n]=0;  
  }
  
  // Autocorrelation times need enough batches; short runs report none.
  const auto count=number.finalize();
  std::cout << "Mean number on Anik: " << count.mean()[0] << " +/- " << count.stderror()[0];
  if (count.tau_available())
    std::cout << " (tau = " << count.tau()[0] << ")";
  std::cout << "\n";

  const auto distribution=histogram.finalize();
  const auto mean=distribution.mean(), error=distribution.stderror(), tau=distribution.tau();
  for (int i=0;i<=N;++i) {
    std::cout << i << "\t" << probability(N,i) << "\t" << mean[i] << "\t" << error[i];
    if (distribution.tau_available())
      std::cout << "\t" << tau[i];
    std::cout << "\n";
  }
  
  return 0;
}
