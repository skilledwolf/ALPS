/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2002-2003 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/fixed_capacity_vector.h>
#include <iostream>
#include <ctime>
#include <vector>

const std::size_t n = (2<<20);
const std::size_t m = 4;

int main()
{
  typedef std::vector<int>::iterator iterator;
  typedef alps::fixed_capacity_vector<int, m> sv_type0;
  typedef std::vector<int> sv_type1;

  std::cout << "allocating " << n << " short vectors of length " << m
            << std::endl;

  const std::clock_t t0 = std::clock();
  std::vector<sv_type0> vec0;
  for (std::size_t i = 0; i < n; ++i) {
    vec0.push_back(sv_type0());
    for (std::size_t j = 0; j < m; ++j) {
      vec0.back().push_back(i);
    }
  }
  std::cout << "fixed_capacity_vector  "
            << static_cast<double>(std::clock() - t0) / CLOCKS_PER_SEC << " sec\n";

  const std::clock_t t1 = std::clock();
  std::vector<sv_type1> vec1;
  for (std::size_t i = 0; i < n; ++i) {
    vec1.push_back(sv_type1());
    for (std::size_t j = 0; j < m; ++j) {
      vec1.back().push_back(i);
    }
  }
  std::cout << "std::vector            "
            << static_cast<double>(std::clock() - t1) / CLOCKS_PER_SEC << " sec\n";
  
  return 0;
}
