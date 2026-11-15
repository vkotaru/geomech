"""Finite sums over a symbolic number of bodies (``Σ_{i=1}^{n}``).

See ``geomech.indexed.core`` for the representation and ``geomech.indexed.eom`` for the
derivation of the equation of the i-th member of a family.
"""

from geomech.indexed.core import (
    Family,
    Index,
    IndexedExpr,
    Size,
    Sum,
    Term,
    as_indexed,
    indices,
    n,
    rename_indices,
)

__all__ = [
    "Family",
    "Index",
    "IndexedExpr",
    "Size",
    "Sum",
    "Term",
    "as_indexed",
    "indices",
    "n",
    "rename_indices",
]

from geomech.indexed.eom import (  # noqa: E402
    FREE,
    IndexedEOM,
    IndexedEquation,
    IndexedVariables,
    compute_indexed_eom,
    unroll,
)

__all__ += [
    "FREE",
    "IndexedEOM",
    "IndexedEquation",
    "IndexedVariables",
    "compute_indexed_eom",
    "unroll",
]
