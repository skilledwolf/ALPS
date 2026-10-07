/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1994-2004 by Matthias Troyer <troyer@itp.phys.ethz.ch>
* Modifications (C) 2026 ALPS Collaboration
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

// ALPS sparse diagonalization application; see A.F. Albuquerque et al.,
// J. of Magn. and Magn. Materials 310, 1187 (2007).
#include "sparsediag.h"
#include "../application.hpp"
#include "schema.hpp"

int main(int argc, char** argv)
{
  return diag::main<SparseDiagMatrix>(argc, argv, "sparsediag", diag_schema);
}
