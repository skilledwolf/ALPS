# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT

"""Exact diagonalization of compiled models (the C++ ``ALPS::ed`` component).

>>> from pyalps import models, ed
>>> model = models.heisenberg(models.lattice.chain(12))
>>> problem = ed.prepare_product(model, sector={"Sz": 0.0})     # dimension 924
>>> result = problem.eigenstates(count=3)
>>> result.eigenvalues[0]                                        # doctest: +SKIP
-5.387390917445...
>>> spectrum = ed.model_spectrum(model, all_momenta=True, count=2)  # Sz and k resolved
>>> spin = models.Spin()
>>> ground = ed.ground_state(model, sector={"Sz": 0.0},
...                          observables={"SzSz": spin.Sz(0) * spin.Sz(1), "Sz": "Sz"})
>>> ground.energy, ground.observables["SzSz"]                    # doctest: +SKIP

``prepare_product``, ``prepare_projected``, ``prepare_finite``,
``model_spectrum`` and ``finite_spectrum`` also take an XML catalog model
(``models.CatalogModel``) in place of a ``Hamiltonian``: its compiled
Hamiltonian runs on its physical basis, in the sector of its ``*_total``
parameters unless ``sector=`` says otherwise (``{}`` is the whole space).
A momentum request (``momentum=``, ``all_momenta=True``) uses the catalog
lattice's own translations (``CatalogModel.translations()``, the generators
sparsediag resolves ``TOTAL_MOMENTUM`` by) unless ``translations=`` gives
others; an open lattice or a plain ``GRAPH`` has none, and such a request is
refused.

>>> ring = models.CatalogModel({"MODEL": "spin", "LATTICE": "chain lattice", "L": 12,
...                             "CONSERVED_QUANTUMNUMBERS": "Sz", "Sz_total": 0})
>>> by_k = ed.model_spectrum(ring, all_momenta=True, count=2)       # Sz=0, every k

Eigenvectors are returned as ``(dimension, count)`` arrays whose column ``j``
belongs to eigenvalue ``j``; sparse matrices carry ``to_scipy()``. Every solve
releases the GIL. ``lanczos(matrix, ...)`` runs the deflated Lanczos driver on
a scipy CSR matrix.
"""

from . import models as _models  # registers the model types first  # noqa: F401

try:
    from ._ext import ed_c as _ed_c
except ImportError as error:  # pragma: no cover - build without PYALPS_BUILD_MODELS
    raise ImportError(
        "pyalps was built without the exact diagonalization bindings (PYALPS_BUILD_MODELS=OFF "
        "or an SDK without ALPS::ed)"
    ) from error

from ._ext.ed_c import *  # noqa: F401,F403

_CatalogModel = getattr(_models, "CatalogModel", None)  # absent without the legacy model library
_CatalogTranslation = getattr(_models, "CatalogTranslation", None)


def _is_catalog(model):
    return _CatalogModel is not None and isinstance(model, _CatalogModel)


def _catalog_inputs(catalog, sector):
    """A catalog model's compiled Hamiltonian, physical basis and sector (its own unless given)."""
    physical = catalog.physical_basis()
    return catalog.compile_hamiltonian(), physical, catalog.sector() if sector is None else dict(sector)


def _catalog_translations(catalog, translations, momentum, all_momenta=False):
    """The SpectrumTranslation values of a catalog request and their lattice axes (None: positional).

    ``translations=None`` takes the catalog lattice's own when a momentum is
    requested and none otherwise; explicit translations pass through, a
    ``CatalogTranslation`` converted.
    """
    requested = bool(len(momentum)) or all_momenta
    if translations is None:
        if not requested:
            return [], None
        translations = catalog.translations()
        if not translations:
            raise ValueError(
                f"momenta need lattice translations, and the catalog lattice of {catalog!r} has none: it has no "
                "periodic direction, is a plain GRAPH, or is disordered or depleted (translations= supplies them)")
    elif requested and not len(translations):
        raise ValueError("momenta need lattice translations, and translations= is empty")
    converted, directions = [], []
    for translation in translations:
        if _CatalogTranslation is not None and isinstance(translation, _CatalogTranslation):
            converted.append(SpectrumTranslation(list(translation.sites), translation.period))
            if directions is not None:
                directions.append(translation.direction)
        else:  # a SpectrumTranslation, which carries no axis; the binding checks the type
            converted.append(translation)
            directions = None
    return converted, directions


def _catalog_product(catalog, sector=None, max_full_dimension=10_000_000):
    op, physical, sector = _catalog_inputs(catalog, sector)
    return _ed_c.prepare_product(op, physical, sector, max_full_dimension)


def _catalog_projected(catalog, sector=None, momentum=(), translations=None):
    translations, _ = _catalog_translations(catalog, translations, momentum)
    op, physical, sector = _catalog_inputs(catalog, sector)
    return _ed_c.prepare_projected(op, projected_spaces(physical), translations, sector, list(momentum))


def _catalog_finite(catalog, irrep, sector=None):
    op, physical, sector = _catalog_inputs(catalog, sector)
    return _ed_c.prepare_finite(op, projected_spaces(physical), irrep, sector)


def _catalog_model_spectrum(catalog, sector=None, momentum=(), all_momenta=False, count=None, which="SA",
                            tolerance=None, max_iterations=None, translations=None):
    translations, directions = _catalog_translations(catalog, translations, momentum, all_momenta)
    op, physical, sector = _catalog_inputs(catalog, sector)
    spectrum = _ed_c.model_spectrum(op, projected_spaces(physical), translations, sector, list(momentum),
                                    all_momenta, count, which, tolerance, max_iterations)
    if directions is not None and translations:
        spectrum.directions = directions  # the lattice axes, as for a Hamiltonian
    return spectrum


def _catalog_finite_spectrum(catalog, irreps, sector=None, count=None, which="SA", tolerance=None,
                             max_iterations=None):
    op, physical, sector = _catalog_inputs(catalog, sector)
    return _ed_c.finite_spectrum(op, projected_spaces(physical), irreps, sector, count, which, tolerance,
                                 max_iterations)


def _accepting_catalogs(name, catalog_form, signature):
    """``name`` from ed_c, with a CatalogModel first argument routed to ``catalog_form``."""
    native = getattr(_ed_c, name)

    def dispatch(model, *args, **kwargs):
        _models.require_bound(model)
        if _is_catalog(model):
            return catalog_form(model, *args, **kwargs)
        return native(model, *args, **kwargs)

    dispatch.__name__ = dispatch.__qualname__ = name
    dispatch.__doc__ = (f"{name}{signature}\n\nA CatalogModel runs its compiled Hamiltonian on its physical "
                        f"basis; ``sector`` defaults to the catalog's own and a momentum request to the catalog "
                        f"lattice's translations.\n\nThe other forms:\n\n{native.__doc__}")
    return dispatch


prepare_product = _accepting_catalogs(
    "prepare_product", _catalog_product, "(catalog, sector=None, max_full_dimension=10000000)")
prepare_projected = _accepting_catalogs(
    "prepare_projected", _catalog_projected, "(catalog, sector=None, momentum=(), translations=None)")
prepare_finite = _accepting_catalogs("prepare_finite", _catalog_finite, "(catalog, irrep, sector=None)")
model_spectrum = _accepting_catalogs(
    "model_spectrum", _catalog_model_spectrum,
    "(catalog, sector=None, momentum=(), all_momenta=False, count=None, which='SA', tolerance=None, "
    "max_iterations=None, translations=None)")
finite_spectrum = _accepting_catalogs(
    "finite_spectrum", _catalog_finite_spectrum,
    "(catalog, irreps, sector=None, count=None, which='SA', tolerance=None, max_iterations=None)")


class GroundState:
    """The lowest eigenpair of a model and the observables measured on it.

    ``energy``, ``vector`` (in the basis of ``problem``), ``observables``
    (name -> expectation value: a float when the imaginary part vanishes to
    rounding, else complex), and the full ``result`` and its ``problem``.
    """

    __slots__ = ("energy", "vector", "observables", "result")

    def __init__(self, energy, vector, observables, result):
        self.energy, self.vector, self.observables, self.result = energy, vector, observables, result

    @property
    def problem(self):
        return self.result.problem

    def __repr__(self):
        return f"GroundState(energy={self.energy!r}, dimension={len(self.vector)}, observables={sorted(self.observables)})"


def _value(z):
    z = complex(z)
    return z.real if abs(z.imag) <= 1e-12 * max(1.0, abs(z.real)) else z


def ground_state(model, *, sector=None, momentum=None, translations=None, observables=None, method="auto",
                 tolerance=None, max_iterations=None):
    """The ground state of a ``Hamiltonian`` or ``CatalogModel`` in one call.

    ``sector`` selects charges (default: the whole space, or a catalog's
    ``*_total`` sector); ``momentum`` selects a translation-projected basis
    (on a catalog, by its lattice's translations unless ``translations=``
    gives others). ``observables`` maps names to anything
    ``models.compile_observable`` accepts (located terms, operator names,
    compiled operators), measured on the ground state; under a
    momentum they must preserve the symmetry sector. ``method`` and the
    Lanczos controls are those of ``eigenstates``.
    """
    _models.require_bound(model)
    sector_args = {} if sector is None else {"sector": sector}
    if momentum is None:
        problem = prepare_product(model, **sector_args)
    elif _is_catalog(model):
        problem = prepare_projected(model, momentum=list(momentum), translations=translations, **sector_args)
    else:
        problem = prepare_projected(model, momentum=list(momentum), **sector_args)
    result = problem.eigenstates(count=1, method=method, tolerance=tolerance, max_iterations=max_iterations)
    measured = {name: _value(result.expectation_values(_models.compile_observable(model, observable))[0])
                for name, observable in (observables or {}).items()}
    return GroundState(float(result.eigenvalues[0]), result.eigenvectors[:, 0], measured, result)


def lanczos(matrix, count=1, which="SA", tolerance=0.0, max_iterations=0, want_vectors=True):
    """Deflated Lanczos on a scipy sparse matrix (real symmetric or complex Hermitian).

    The matrix is converted to CSR; ``which`` is ``SA``, ``LA`` or ``LM`` and
    ``tolerance``/``max_iterations`` at or below zero select the driver defaults.
    Degenerate eigenvalues come back with multiplicity.
    """
    import scipy.sparse

    csr = scipy.sparse.csr_matrix(matrix)
    if csr.shape[0] != csr.shape[1]:
        raise ValueError("lanczos needs a square matrix")
    return _ed_c.lanczos_csr(csr.indptr, csr.indices, csr.data, csr.shape[0], count, which,
                             tolerance, max_iterations, want_vectors)


__all__ = [name for name in dir(_ed_c) if not name.startswith("_")] + ["GroundState", "ground_state", "lanczos"]
