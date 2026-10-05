// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
#include "sse.h"
#include "schema.hpp"
#include "../native_driver.hpp"
int main(int argc,char** argv) {
    return native_qmc::main<SSE_run>(argc,argv,"dirloop_sse",qmc_common_schema,qmc_application_schema);
}
