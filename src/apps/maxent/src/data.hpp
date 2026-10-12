// SPDX-License-Identifier: MIT
#pragma once
#include <alps/maxent.hpp>
#include <alps/ngs/params.hpp>
namespace maxent_detail {
alps::maxent::data legacy_data(const alps::params&);
int data_count(const alps::maxent::data&);
}
