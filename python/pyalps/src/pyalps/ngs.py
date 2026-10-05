 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 #                                                                                 #
 # ALPS Project: Algorithms and Libraries for Physics Simulations                  #
 #                                                                                 #
 # ALPS Libraries                                                                  #
 #                                                                                 #
 # Copyright (C) 2010 - 2013 by Lukas Gamper <gamperl@gmail.com>                   #
 #                      2012 by Troels F. Roennow <tfr@nanophysics.dk>             #
 #                                                                                 #
# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
 #                                                                                 #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from collections.abc import MutableMapping
from .cxx.pyngsparams_c import params
from .cxx.pyngsbase_c import mcbase
from .cxx.pyngsapi_c import collectResults, saveResults
from .cxx.pyngsrandom01_c import random01

# Extension types cannot inherit the Python mapping mixin directly.
MutableMapping.register(params)
for _method in ("keys", "values", "items", "get", "popitem", "clear", "update",
                "setdefault", "__eq__", "__ne__", "__hash__"):
    if getattr(params, _method, None) is getattr(object, _method, None):
        setattr(params, _method, getattr(MutableMapping, _method))

_MAPPING_POP_MARKER = object()


def _mapping_pop(self, key, default=_MAPPING_POP_MARKER):
    if key in self:
        value = self[key]
        del self[key]
        return value
    if default is _MAPPING_POP_MARKER:
        raise KeyError(key)
    return default


params.pop = _mapping_pop
