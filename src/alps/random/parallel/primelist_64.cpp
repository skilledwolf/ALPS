/* 
 * Copyright Matthias Troyer 2006
 * SPDX-License-Identifier: MIT
*
 */

#include <alps/random/parallel/detail/primelist_64.hpp>

namespace alps { namespace random { namespace detail {

unsigned int primelist_64[primelist_size2] = 
{
#include <alps/random/parallel/detail/primelist_64.ipp>
};

} } } // namespace alps::random::detail
