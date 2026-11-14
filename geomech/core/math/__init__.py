"""Coefficient extraction, ``collect`` and integration by parts.

``collect`` is defined in ``geomech.core.math.collecting`` (a module name distinct from
the function, so neither hides the other): ``from geomech.core.math import collect`` gives
the function, ``import geomech.core.math.collecting as m`` the module.
"""

from geomech.core.math.collecting import collect
from geomech.core.math.extract import extract_linear_coeff
from geomech.core.math.ibp import integrate_by_parts

__all__ = [
    "extract_linear_coeff",
    "collect",
    "integrate_by_parts",
]
