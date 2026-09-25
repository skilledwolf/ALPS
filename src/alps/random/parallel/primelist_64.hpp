/* 
 * Copyright Matthias Troyer 2006
 * SPDX-License-Identifier: MIT
*
 */

#define PRIMELISTSIZE1 1000
#define STEP 10000
#define PRIMELISTSIZE2 15613

unsigned int prime_list_64[PRIMELISTSIZE2] = 
{
#include <alps/random/parallel/detail/primelist_64.ipp>
};

