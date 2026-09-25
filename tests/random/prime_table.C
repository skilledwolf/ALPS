// SPDX-License-Identifier: MIT
#include <alps/random/parallel/primelist_64.hpp>
#include <alps/random/parallel/detail/get_prime.hpp>
#include <cstdint>
#include <initializer_list>
#include <stdexcept>
#include <type_traits>

int main() {
  static_assert(std::is_same_v<decltype(prime_list_64), unsigned int[15613]>);
  std::uint64_t hash = 14695981039346656037ull;
  for (auto prime : prime_list_64) hash = (hash ^ prime) * 1099511628211ull;
  if (hash != 1281122858450542406ull)
    throw std::runtime_error("parallel RNG prime table changed");
  using alps::random::detail::get_prime_64;
  for (unsigned i = 0; i < PRIMELISTSIZE1; ++i)
    if (get_prime_64(i) != prime_list_64[i])
      throw std::runtime_error("parallel RNG direct prime lookup changed");
  for (unsigned i : {1000u, 1001u, 2000u, 10000u, 15612u})
    if (get_prime_64(999 + (i - 999) * STEP) != prime_list_64[i])
      throw std::runtime_error("parallel RNG sparse prime lookup changed");
  auto first = prime_list_64[0];
  ++prime_list_64[0];
  if (get_prime_64(0) != first)
    throw std::runtime_error("public and internal prime arrays must remain independent");
  prime_list_64[0] = first;
}
