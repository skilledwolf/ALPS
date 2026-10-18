# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
#
# Adapted from skilledwolf/alpsx (private repository), branch
# agent/v3-execution-and-api-audit at d37e51227,
# python/qdk-chemistry-alps/examples/rdm_observables.py. Changes: the RDMs
# are those pyalps.classic_dmrg.solve_molecular reads from the final state,
# on the integrals hamiltonian_to_integrals extracts from the lattice model,
# instead of contractions of an alpsx MPS; a local-moment identity checks
# the spin blocks against the density.
"""Fourteen-site Hubbard observables without constructing a CI vector.

An ALPS-native route: pyalps.classic_dmrg.solve_molecular returns the final
state's reduced density matrices, so this runs with released qdk-chemistry,
which has no RdmContainer. Returning them as a QDK RdmContainer is
large_active_space.py, with a qdk-chemistry build that provides it.
"""

import time
from math import comb

import numpy as np
from demo_support import finish, parser
from pyalps import classic_dmrg as dmrg
from pyalps import models

from qdk_chemistry_alps.interop import hamiltonian_to_integrals

SITES, T, U = 14, 1., 4.
SWEEPS, STATES = 8, 64


def plot(plt, data):
    figure, axes = plt.subplots(1, 2, figsize=(9, 3.5), layout="constrained")
    axes[0].plot(data["natural_occupations"], "o-", color="#087e8b")
    axes[0].set(xlabel="Natural orbital (descending occupation)", ylabel="Occupation", ylim=(0, 2),
                title="14-site Hubbard: one-particle RDM")
    axes[1].plot(range(SITES), data["double_occupancy"], "o-", color="#d68c22")
    axes[1].set(xlabel="Site", ylabel="⟨n↑ n↓⟩", title="Local two-particle information")
    return figure


def main():
    args = parser(__doc__).parse_args()
    started = time.perf_counter()
    h = models.hubbard(models.lattice.chain(SITES, bc="open"), t=T, U=U)
    h1, eri, core = hamiltonian_to_integrals(h)
    half = SITES // 2
    result = dmrg.solve_molecular(h1, eri, half, half, core_energy=core,
                                  sweeps=dmrg.fixed_sweeps(SWEEPS, STATES))
    gamma = np.asarray(result.one_rdm)  # gamma[i, j] = sum_s <c+_js c_is>
    doubles = dmrg.local_expectations(result, "n_updn")
    # For Hubbard only these two-particle diagonal elements are needed.
    from_rdms = -T * (np.trace(gamma, offset=1) + np.trace(gamma, offset=-1)) + U * doubles.sum()
    np.testing.assert_allclose(from_rdms, result.energy, atol=1e-9, rtol=0)
    np.testing.assert_allclose(np.trace(gamma), SITES, atol=1e-9, rtol=0)
    # <S_i^2> = 3/4 (<n_i> - 2 <n_i,up n_i,dn>) ties the spin blocks to the density.
    moments = np.diag(dmrg.spin_spin_correlations(result))
    np.testing.assert_allclose(moments, .75 * (np.diag(gamma) - 2 * doubles), atol=1e-9, rtol=0)
    natural = np.linalg.eigvalsh(gamma)[::-1]
    assert np.min(natural) > -1e-9 and np.max(natural) < 2 + 1e-9
    determinants = comb(SITES, half)**2
    error = abs(from_rdms - result.energy)
    print("ALPS-native classic DMRG; no QDK RdmContainer and no CI vector")
    print(f"14-site Hubbard: E/t={result.energy:.10f}; RDM energy error={error:.2e}")
    print(f"{determinants:,} determinants would need {8*determinants/1e6:.1f} MB for real amplitudes alone.")
    print("This checks energy consistency and particle count, not DMRG convergence.")
    finish(args, {"model": {"sites": SITES, "t": T, "U": U, "sector": {"N": SITES, "Sz": 0}, "boundary": "open"},
                  "controls": {"states": STATES, "sweeps": SWEEPS},
                  "workflow": "ALPS-native observables; not a QDK RdmContainer",
                  "energy": result.energy, "energy_from_rdms": float(from_rdms),
                  "final_truncation_error": result.final_truncation_error,
                  "natural_occupations": natural, "double_occupancy": doubles, "local_moments": moments,
                  "one_rdm": gamma, "determinants": determinants}, started, plot)


if __name__ == "__main__":
    main()
