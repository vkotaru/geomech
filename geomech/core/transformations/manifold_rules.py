"""Manifold kinematic substitutions applied during EOM derivation.

For S2 (unit sphere) with manifold point q, tangent vector ω, variation vector ξ:
  - Variation(ω) → ξ̇ - 2 ω × ξ        (Lee et al. eq. 5.13 on the constraint set ξ ⊥ q)

For SO3 (rotation group) with tangent vector Ω, variation vector η:
  - Variation(Ω) → η̇ + Ω × η          (body-frame angular velocity convention)

Derivation for S2.  The kinematics are q̇ = ω × q with ω ⊥ q (so ω = q × q̇) and the
variation is δq = ξ × q with ξ ⊥ q.  Lee et al. eq. (5.13) give the variation of the
angular velocity as

    δω = -S(ω) ξ + (I - q qᵀ) ξ̇ = ξ̇ - ω × ξ - q (q · ξ̇).

The last term is not zero.  Differentiating the constraint ξ · q = 0 along the curve,

    q · ξ̇ = -ξ · q̇ = -ξ · (ω × q) = q · (ω × ξ),

and because ω and ξ both lie in the tangent plane, ω × ξ is parallel to q, so
q (q · ξ̇) = q (q · (ω × ξ)) = ω × ξ.  Hence, on the constraint set,

    δω = ξ̇ - 2 ω × ξ.

δω is *not* tangent to the sphere: its normal component is q · δω = -q · ξ̇ = -q · (ω × ξ),
which is what keeps δ(ω · q) = 0.  The earlier rule ξ̇ - ω × ξ is only the tangential
projection (I - q qᵀ) δω; it is correct whenever ∂L/∂ω ⊥ q (every Lagrangian that depends
on ω only through q̇ = ω × q) and wrong otherwise: with ∂L/∂ω = m l² ω + c the gyroscopic
term 2 ω × c came out as ω × c, i.e. at half strength.

For SO3, with Ṙ = R·Hat(Ω) and δR = R·Hat(η), equating the mixed derivatives
d/dt(δR) = δ(Ṙ) gives δΩ = η̇ + Ω × η exactly (Lee et al. eq. 6.9).
"""

from __future__ import annotations

from geomech.core.base.expressions import S2, SO3, TS2, TSO3, Scalar
from geomech.core.operations.addition import VAdd
from geomech.core.operations.calculus import TimeDerivative, Variation
from geomech.core.operations.geometry import Cross
from geomech.core.operations.multiplication import SVMul
from geomech.dynamics.variables import SystemVariables

_NEG2 = Scalar("(-2)", value=-2, attr=["Constant"])


def apply_manifold_rules(expr, variables: SystemVariables):
    """Apply manifold kinematic substitutions to an expression.

    Walks the expression tree and replaces Variation(TS2) and
    Variation(TSO3) nodes with their kinematic expansions.
    """
    tangent_map = {}

    # S2 manifolds (vectors)
    for vec in variables.vectors:
        if isinstance(vec, S2):
            omega = vec.get_tangent_vector()
            xi = vec.get_variation_vector()
            tangent_map[str(omega)] = ("S2", vec, omega, xi)

    # SO3 manifolds (matrices)
    for mat in variables.matrices:
        if isinstance(mat, SO3):
            Omega = mat.get_tangent_vector()
            eta = mat.get_variation_vector()
            tangent_map[str(Omega)] = ("SO3", mat, Omega, eta)

    if not tangent_map:
        return expr

    return _substitute(expr, tangent_map)


def _substitute(expr, tangent_map):
    """Recursively substitute manifold variations."""

    # Variation(TS2) or Variation(TSO3) → kinematic expansion
    if isinstance(expr, Variation):
        inner = expr.expr
        if isinstance(inner, (TS2, TSO3)) and str(inner) in tangent_map:
            kind, parent, omega, var_vec = tangent_map[str(inner)]
            if kind == "S2":
                # δ(ω) = ξ̇ - 2 ω × ξ   (see the module docstring)
                return VAdd(
                    TimeDerivative(var_vec),
                    SVMul(Cross(omega, var_vec), _NEG2),
                )
            elif kind == "SO3":
                # δ(Ω) = η̇ + Ω × η
                return VAdd(
                    TimeDerivative(var_vec),
                    Cross(omega, var_vec),
                )
            else:
                raise ValueError(f"Unknown manifold kind: {kind}")
        # Variation of non-tangent or unknown — recurse into inner
        inner_sub = _substitute(inner, tangent_map)
        if inner_sub is not inner:
            return Variation(inner_sub)
        return expr

    # Recurse into children
    nodes = getattr(expr, "nodes", None)
    if nodes is None or len(nodes) == 0:
        return expr

    new_nodes = [_substitute(n, tangent_map) for n in nodes]
    if all(n is o for n, o in zip(new_nodes, nodes)):
        return expr  # nothing changed

    # Rebuild the node with new children
    return _rebuild(expr, new_nodes)


def _rebuild(expr, new_nodes):
    """Rebuild an expression node with new children."""
    from geomech.core.operations.addition import Add, MAdd, VAdd
    from geomech.core.operations.calculus import TimeDerivative, TimeIntegral, Variation
    from geomech.core.operations.geometry import Cross, Dot, Hat, Transpose, Vee
    from geomech.core.operations.multiplication import (
        MMMul,
        Mul,
        MVMul,
        SMMul,
        SVMul,
        VVMul,
    )

    cls = type(expr)

    # N-ary (Add, VAdd, MAdd)
    if isinstance(expr, (Add, VAdd, MAdd)):
        return cls(*new_nodes)

    # Binary ops
    if isinstance(expr, (Mul, SVMul, SMMul, MVMul, MMMul, VVMul, Dot, Cross)):
        return cls(new_nodes[0], new_nodes[1])

    # Unary ops
    if isinstance(expr, (Hat, Vee, Transpose, Variation, TimeDerivative, TimeIntegral)):
        return cls(new_nodes[0])

    return expr
