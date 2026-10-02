# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
add_executable(run_config_contract contract.cpp)
target_link_libraries(run_config_contract PRIVATE ALPS::run_config)
alps_add_test(run_config_contract)
set_property(TEST run_config_contract PROPERTY LABELS "params;configuration")
