/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include "../classical.hpp"

int main(int argc,char** argv) {return classical_main<classical_walker<1>>(argc,argv,"ising",{"ising","ising; exchange"});}
