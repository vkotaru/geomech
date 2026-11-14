r"""Readable LaTeX names for the leaves of an expression.

The engine names symbols the way they are typed (``J1``, ``e3``, ``\Omega_{R1}``), and
``to_latex`` prints a leaf by its name, so a formula shows ``J1 \dot{\Omega_{R1}}``.
``latex_names`` builds a map from each leaf name to a subscripted form, which
``to_latex(expr, names=...)`` applies to the leaves while it prints them. Nothing is ever
rewritten in an assembled formula: a pattern applied to a whole formula could match
across symbol boundaries (``\Omega_{R1}`` inside ``\dot{\Omega_{R1}}``, a digit of a
number) and corrupt it.

Rules, applied to one leaf name at a time (``pretty_name``):

====================================  ===========================  =======================
Leaf                                  Example                      Printed
====================================  ===========================  =======================
letters (or one ``\command``) + digits ``J1``, ``e3``, ``m12``     ``J_1``, ``e_3``, ``m_{12}``
tangent/variation vector of a point    ``\Omega_{R1}``, ``\eta_{R2}``  ``\Omega_1``, ``\eta_2``
named letters + digits
anything else                          ``m_L``, ``\rho_1``, ``q``, ``\omega_{q}``  unchanged
====================================  ===========================  =======================

The map must stay one-to-one: when two different leaves would print alike (``J1`` and a
leaf already named ``J_1``; ``\Omega_{R1}`` and ``\Omega_{Q1}``), neither is renamed.
"""

from __future__ import annotations

import re
from typing import Iterable

__all__ = ["pretty_name", "latex_names"]

# one leaf name, never a formula: letters or a single LaTeX command, then digits
_INDEXED = re.compile(r"^([A-Za-z]+|\\[A-Za-z]+)(\d+)$")
# the tangent/variation vector names of an S2/SO3 point: \omega_{q}, \xi_{q}, \Omega_{R}, \eta_{R}
_TANGENT = re.compile(r"^(\\(?:omega|xi|Omega|eta))_\{(.+)\}$")


def _subscript(base: str, index: str) -> str:
    return f"{base}_{index}" if len(index) == 1 else f"{base}_{{{index}}}"


def pretty_name(name: str, point: str | None = None) -> str:
    r"""The LaTeX of one leaf name (see the module docstring for the rules).

    *point* is the name of the S² or SO(3) point when the leaf is its tangent or variation
    vector (``\Omega_{R1}`` of ``R1``); the vector then takes the point's index
    (``\Omega_1``). Without *point* such a name is left as it is.
    """
    if point is not None:
        m = _TANGENT.match(name)
        p = _INDEXED.match(point)
        if m and m.group(2) == point and p:
            return _subscript(m.group(1), p.group(2))
        return name
    m = _INDEXED.match(name)
    if m:
        return _subscript(m.group(1), m.group(2))
    return name


def _leaves(expr, out: dict):
    nodes = getattr(expr, "nodes", None)
    if nodes:
        for n in nodes:
            _leaves(n, out)
        return
    name = getattr(expr, "name", None)
    if isinstance(name, str):
        out.setdefault(name, expr)


def _point_of(leaf) -> str | None:
    from geomech.core.base.expressions import TS2, TSO3

    point = None
    if isinstance(leaf, TS2):
        point = leaf.S2
    elif isinstance(leaf, TSO3):
        point = leaf.SO3
    return getattr(point, "name", None) if point is not None else None


def _gather(obj, out: dict, depth: int = 0):
    """Leaves of an expression, a ``compute_eom``/``to_standard_form`` dict, a
    ``StandardFormEquation``, a ``(variation, equation)`` pair or a list of these."""
    from geomech.core.base.expressions import Expr
    from geomech.dynamics.standard_form import StandardFormEquation

    if isinstance(obj, Expr):
        _leaves(obj, out)
    elif isinstance(obj, StandardFormEquation):
        for block in [*obj.M.values(), obj.f, *obj.G.values()]:
            if block is not None:
                _gather(block, out, depth + 1)
    elif isinstance(obj, dict):
        for value in obj.values():
            _gather(value, out, depth + 1)
    elif isinstance(obj, (list, tuple)) and depth < 6:
        for item in obj:
            _gather(item, out, depth + 1)


def latex_names(*objs: Iterable) -> dict[str, str]:
    r"""``{leaf name: LaTeX}`` for every leaf of *objs* whose printed form changes.

    *objs* are expressions, the dicts of ``compute_eom`` or ``to_standard_form``, the
    ``SystemVariables`` lists or any list of these. The tangent and variation vectors of
    each S² and SO(3) point in *objs* are included. A name is left out when its pretty
    form would collide with another leaf's name or pretty form, so the map is one-to-one
    on the leaves it sees. Pass the result to ``to_latex(expr, names=...)``.
    """
    from geomech.core.base.expressions import S2, SO3

    leaves: dict = {}
    for obj in objs:
        _gather(obj, leaves)
    # the velocity and variation vectors of every manifold point, whether or not they occur
    for leaf in list(leaves.values()):
        if isinstance(leaf, (S2, SO3)):
            for vec in (leaf.get_tangent_vector(), leaf.get_variation_vector()):
                leaves.setdefault(str(vec), vec)
    proposed = {name: pretty_name(name, _point_of(leaf)) for name, leaf in leaves.items()}
    # one-to-one: a target shared by two leaves, or equal to another leaf's own name, is
    # dropped for every leaf involved
    count: dict[str, int] = {}
    for target in proposed.values():
        count[target] = count.get(target, 0) + 1
    out = {}
    for name, target in proposed.items():
        if target == name:
            continue
        if count[target] > 1 or (target in leaves and target != name):
            continue
        out[name] = target
    return out
