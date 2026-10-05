// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include "WRun.h"
#include "schema.hpp"
#include "../native_driver.hpp"
int main(int argc,char** argv) {
    return native_qmc::main<WRun>(argc,argv,"worm",qmc_common_schema,qmc_application_schema);
}
