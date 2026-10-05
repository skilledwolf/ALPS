// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include <alps/fortran/fortran_wrapper.h>
#include "schema.hpp"
int main(int argc, char** argv) { return alps::fortran_wrapper::main(argc, argv, application, schema); }
