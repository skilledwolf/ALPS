// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
#include <alps/solvers.hpp>

// Keep references to all three entry points, including in optimized builds:
// the installed archives must be complete and link together without collisions.
int main() {
  using solver = void (*)(alps::params const&, std::string const&);
  solver volatile entrypoints[] = {
      &alps::solvers::cthyb, &alps::solvers::ctint};
  auto volatile maxent = &alps::solvers::maxent;
  if (!maxent) return 1;
  for (auto const& entrypoint : entrypoints)
    if (!entrypoint) return 1;
}
