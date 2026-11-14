"""Tree rewrites: ``expand``, ``pull``, ``vector_rules``, ``simplify``/``full_simplify`` and
``apply_manifold_rules``.

Each function lives in a module with a different name, so neither hides the other:
``expand`` in ``expansion``, ``pull`` in ``pulling``, ``vector_rules`` in
``vector_identities``, ``simplify`` and ``full_simplify`` in ``simplification``, and
``apply_manifold_rules`` in ``manifold_rules`` (all under
``geomech.core.transformations``).

``from geomech.core.transformations import simplify`` gives the function (the same object
as ``geomech.simplify``); ``import geomech.core.transformations.simplification as m``
gives the module, whose docstring lists its named rules (set breakpoints or patch helpers
there).
"""

from geomech.core.transformations.expansion import expand
from geomech.core.transformations.manifold_rules import apply_manifold_rules
from geomech.core.transformations.pulling import pull
from geomech.core.transformations.simplification import full_simplify, simplify
from geomech.core.transformations.utils import has_nested_add, is_leaf
from geomech.core.transformations.vector_identities import vector_rules

__all__ = [
    "apply_manifold_rules",
    "expand",
    "full_simplify",
    "has_nested_add",
    "is_leaf",
    "pull",
    "simplify",
    "vector_rules",
]
