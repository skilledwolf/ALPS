// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
#include <alps/solvers.hpp>
#include <alps/cthyb.hpp>
#include <alps/ctint.hpp>

// Keep references to all three entry points, including in optimized builds:
// the installed archives must be complete and link together without collisions.
int main() {
  using solver = void (*)(alps::run_configuration const&);
  solver volatile entrypoints[] = {
      &alps::solvers::cthyb, &alps::solvers::ctint};
  auto volatile maxent = static_cast<bool (*)(alps::params const&, const alps::maxent::data&,
                                              std::string const&, int, bool)>(&alps::solvers::maxent);
  auto volatile configured_maxent = static_cast<bool (*)(alps::run_configuration const&)>(
      &alps::solvers::maxent);
  if (!maxent || !configured_maxent || alps::cthyb::schema().empty() ||
      alps::ctint::schema().empty()) return 1;
  for (auto const& entrypoint : entrypoints)
    if (!entrypoint) return 1;
}
