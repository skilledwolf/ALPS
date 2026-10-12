# ****************************************************************************
# 
# ALPS Project: Algorithms and Libraries for Physics Simulations
# 
# ALPS Libraries
# 
# Copyright (C) 1994-2010 by Bela Bauer <bauerb@phys.ethz.ch>
#                            Ping Nang Ma
#
# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
# 
# ****************************************************************************

import numpy as np
from pyalps.alea import *
import math as pm

def _dispatcher(operation):
    numpy_name = {
        "asin": "arcsin", "acos": "arccos", "atan": "arctan",
        "asinh": "arcsinh", "acosh": "arccosh", "atanh": "arctanh",
    }.get(operation, operation)

    def apply(obj):
        if isinstance(obj, (float, int, np.ndarray)):
            if operation == "sq":
                return obj ** 2
            if operation == "cb":
                return obj ** 3
            library = np if isinstance(obj, np.ndarray) else pm
            name = numpy_name if library is np else operation
            return getattr(library, name)(obj)
        return getattr(obj, operation)()

    apply.__name__ = operation
    return apply


for _operation in ["sq", "sqrt", "cb", "cbrt", "exp", "log", "sin", "cos",
                   "tan", "asin", "acos", "atan", "sinh", "cosh", "tanh",
                   "asinh", "acosh", "atanh"]:
    globals()[_operation] = _dispatcher(_operation)
