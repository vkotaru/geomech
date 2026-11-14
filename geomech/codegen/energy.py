"""Total energy ``E = Σ (∂L/∂v)·v − L`` of a Lagrangian, derived with the engine.

For every configuration variable the velocity ``v`` is its velocity slot: ``ṡ``, ``ẋ``,
``ω`` (S²) or ``Ω`` (SO(3)). ``Σ (∂L/∂v)·v`` is the derivative of ``L(q, λv)`` at ``λ = 1``
(for ``L = L₂ + L₁ + L₀``, homogeneous of degree 2, 1, 0 in the velocities, it is
``2L₂ + L₁``, so ``E = L₂ − L₀``). It does not depend on the velocity coordinates as long as
they are linear in the rates: on S², ``q̇ = ω × q`` with ``ω ⊥ q``, so scaling ``ω`` scales
``q̇``; on SO(3), ``Ṙ = R Ω̂``. ``E`` is therefore the energy of Lee et al. (the Legendre
transform on ``TS²`` and ``TSO(3)``), conserved when ``δW = 0``.

``∂L/∂v`` is the coefficient of the velocity perturbation in the first variation of ``L``:
``δL = Σ (∂L/∂v)·δv + (configuration terms)``, with the velocity perturbations of
``linearize`` (``δṡ = d/dt(δs)``, ``δẋ = d/dt(δx)``, ``δω``, ``δΩ``). ``δL`` is expanded,
its terms are partitioned by perturbation symbol exactly as ``linearize`` does, and
``extract_linear_coeff`` gives the coefficient of each velocity perturbation.
"""

from __future__ import annotations

from typing import Callable

from geomech.codegen.layout import StateLayout
from geomech.codegen.source import function_source, get_target, namespace
from geomech.core.base.expressions import S2, SO3, Expr
from geomech.core.base.types import ExprType
from geomech.core.math.extract import extract_linear_coeff
from geomech.core.operations.geometry import Dot
from geomech.core.operations.multiplication import Mul
from geomech.core.transformations.expansion import expand
from geomech.core.transformations.simplification import full_simplify
from geomech.core.transformations.utils import sum_of
from geomech.dynamics.eom import check_variables
from geomech.dynamics.variables import SystemVariables
from geomech.dynamics.vbl import _partition_terms, _perturbation_variables

__all__ = ["momenta", "energy_expr", "energy_function", "energy_source"]


def _velocity(var):
    if isinstance(var, (S2, SO3)):
        return var.get_tangent_vector()
    return var.t_diff()


def momenta(L: Expr, variables: SystemVariables) -> dict[str, tuple[Expr, Expr]]:
    """``{str(v): (v, ∂L/∂v)}`` for the velocity ``v`` of every configuration variable.

    ``∂L/∂v`` is a scalar for a scalar variable and a column vector otherwise (for S², the
    full gradient in ``ω``, which may have a component along ``q``).
    """
    check_variables(variables, caller="momenta")
    cfg_vars, cfg_perts, vel_perts = _perturbation_variables(variables)
    dL = expand(full_simplify(L.delta()))
    buckets = _partition_terms(dL, vel_perts + cfg_perts, "L")
    out = {}
    for var, pert, terms in zip(cfg_vars, vel_perts, buckets):
        coeff = extract_linear_coeff(sum_of(terms, ExprType.SCALAR), pert)
        v = _velocity(var)
        out[str(v)] = (v, full_simplify(coeff))
    return out


def energy_expr(L: Expr, variables: SystemVariables) -> Expr:
    """The total energy ``Σ (∂L/∂v)·v − L`` as an expression tree (simplified)."""
    terms: list[Expr] = []
    for v, p in momenta(L, variables).values():
        if p.is_zero:
            continue
        terms.append(Mul(p, v) if v.type == ExprType.SCALAR else Dot(p, v))
    return full_simplify(sum_of(terms, ExprType.SCALAR) - L)


def energy_source(
    L: Expr, layout: StateLayout, target: str = "numpy", name: str = "energy", *, expr=None
) -> str:
    """Source of ``def energy(state, p)`` returning the total energy (see ``energy_expr``);
    *expr* is the energy if already derived."""
    E = energy_expr(L, layout.variables) if expr is None else expr
    doc = "Total energy E = sum((dL/dv) . v) - L at the state."
    kinds = {"configuration", "velocity", "param"}
    return function_source(
        name,
        ("state", "p"),
        doc,
        layout,
        get_target(target),
        [E],
        lambda em, codes: [f"return {codes[0].text}"],
        kinds=kinds,
        what="the energy",
    )


def energy_function(L: Expr, layout: StateLayout, target: str = "numpy") -> Callable:
    """The total energy of *L* as a callable ``energy(state, p)`` (generated code).

    ``E = Σ (∂L/∂v)·v − L`` over the velocities of ``layout.variables``, derived
    symbolically (``energy_expr``) and generated through the same code path as every
    other block. The callable carries its expression in ``.expr`` and its source in
    ``.source``. Every leaf of ``L`` other than the configuration and velocity symbols
    must be a parameter of *layout*.
    """
    E = energy_expr(L, layout.variables)
    src = "_PARAMS = " + repr(tuple(layout.param_names)) + "\n\n\n"
    src += energy_source(L, layout, target, expr=E)
    ns = namespace(target)
    exec(compile(src, "<geomech.codegen.energy>", "exec"), ns)
    fn = ns["energy"]
    fn.source = src
    fn.expr = E
    return fn
