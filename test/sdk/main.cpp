// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#include <alps/ngs/params.hpp>
#include <alps/version.h>
int main() {
    alps::params parameters;
    parameters["value"] = 42;
    return parameters["value"].cast<int>() == 42 ? 0 : 1;
}
