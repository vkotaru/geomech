"""Standard form extraction: M(q)*a + f(q,dq) + G(q,dq)*u = 0.

to_standard_form(eom_dict, variables, inputs)
    Decompose EOM into mass matrix, nonlinear terms, and input influence.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from geomech.core.base.expressions import S2, SO3, Expr, Zero, ZeroVector
from geomech.core.base.types import ExprType
from geomech.core.math.extract import _has_algebraic, extract_linear_coeff
from geomech.core.operations.addition import Add, VAdd
from geomech.core.operations.calculus import TimeDerivative
from geomech.core.transformations.expansion import expand
from geomech.core.transformations.simplification import full_simplify
from geomech.dynamics.variables import SystemVariables
from geomech.utils.errors import AlgebraicError, NonLinearError
from geomech.utils.trace import debug_enabled, stage_recorder


@dataclass
class StandardFormEquation:
    """Standard form for one EOM equation.

    ``Σ M[a_j] * a_j  +  f  +  Σ G[u_k] * u_k  =  0``

    Attributes
    ----------
    M : dict[str, Expr]
        Mass/inertia coefficient for each acceleration variable.
    f : Expr
        Nonlinear terms (Coriolis, gravity, etc.) — scalar or vector.
    G : dict[str, Expr]
        Input influence coefficient for each input variable.
    """

    M: dict = field(default_factory=dict)
    f: Expr = None
    G: dict = field(default_factory=dict)


def to_standard_form(eom_dict, variables: SystemVariables, inputs=None):
    """Decompose EOM into standard manipulator form.

    Parameters
    ----------
    eom_dict : dict
        Output of ``compute_eom()``:
        ``{str(variation_vec): (variation_vec, eom_expr)}``.
    variables : SystemVariables
        Configuration variables (used to derive acceleration variables).
    inputs : list[Expr], optional
        Input variables (forces, torques). Default: empty.

    Returns
    -------
    dict[str, StandardFormEquation]
        Keyed by the same variation vector names as *eom_dict*.

    Notes
    -----
    On S² the equation ``E = 0`` of a variable ``q`` is only meaningful on the tangent
    plane, ``(I − qqᵀ) E = 0``: ``E`` is the coefficient of ``ξ ⊥ q`` and in general has
    a component along ``q`` (from the forcing, from ``ω × c``, from ``A ω̇``, …). The
    blocks returned here are the terms of ``E`` as derived, not projected. To solve for
    ``ω̇`` solve ``(I − qqᵀ)(M ω̇ + f + Σ G u) = 0`` with ``ω̇ ⊥ q``; the unprojected
    3×3 solve gives ``q · ω̇ ≠ 0`` and a trajectory that leaves the sphere. The same
    holds for the rows of ``linearize``.

    Errors name the equation and the term: a term that holds an acceleration only under
    a further time derivative (a jerk ``d/dt(ẍ)``) raises ``AlgebraicError`` saying the
    equation is not second order; a term bilinear in two accelerations raises
    ``NonLinearError``. With ``GEOMECH_DEBUG=1`` the blocks of each equation are logged
    on the ``geomech`` logger and ``check_tree`` runs on every equation and block.
    """
    if inputs is None:
        inputs = []

    accel_vars = _get_accel_vars(variables)
    debug = debug_enabled()
    if debug:
        from geomech.utils.invariants import check_tree
    record = stage_recorder(None, "to_standard_form")

    result = {}
    for key, (var_vec, eom_expr) in eom_dict.items():
        if debug:
            check_tree(eom_expr)
        sf = _extract_single(eom_expr, accel_vars, inputs, key)
        if record is not None:
            for name, block in sf.M.items():
                record(f"M[{key}][{name}]", block)
            record(f"f[{key}]", sf.f)
            for name, block in sf.G.items():
                record(f"G[{key}][{name}]", block)
        if debug:
            for block in (*sf.M.values(), sf.f, *sf.G.values()):
                check_tree(block)
        result[key] = sf

    return result


def _get_accel_vars(variables: SystemVariables):
    """Derive acceleration variables from config variables.

    - Scalar x  →  d²x/dt²  =  TimeDerivative(TimeDerivative(x))
    - Vector x  →  d²x/dt²  =  TimeDerivative(TimeDerivative(x))
    - S2 q      →  dω/dt    =  TimeDerivative(ω_q)
    - SO3 R     →  dΩ/dt    =  TimeDerivative(Ω_R)
    """
    accels = []
    for s in variables.scalars:
        accels.append(s.t_diff().t_diff())
    for v in variables.vectors:
        if isinstance(v, S2):
            accels.append(v.get_tangent_vector().t_diff())
        else:
            accels.append(v.t_diff().t_diff())
    for m in variables.matrices:
        if isinstance(m, SO3):
            accels.append(m.get_tangent_vector().t_diff())
        else:
            accels.append(m.t_diff())
    return accels


def _flatten_addends(expr):
    """Return a flat list of addends from a scalar or vector expression."""
    if isinstance(expr, (Add, VAdd)):
        return list(expr.nodes)
    return [expr]


def _zero_for(eom_expr):
    """Return the appropriate zero element for the EOM's type."""
    if getattr(eom_expr, "type", None) == ExprType.SCALAR:
        return Zero
    return ZeroVector


def _combine(terms, eom_expr):
    """Combine a list of terms into a single expression matching *eom_expr*'s type."""
    if not terms:
        return _zero_for(eom_expr)
    if len(terms) == 1:
        return terms[0]
    if getattr(eom_expr, "type", None) == ExprType.SCALAR:
        return Add(*terms)
    return VAdd(*terms)


def _calculus_holding(term, query):
    """The outermost calculus node of *term* that contains *query* (for error messages)."""
    stack = [term]
    while stack:
        node = stack.pop()
        if isinstance(node, TimeDerivative) and node.has(query):
            return node
        stack.extend(getattr(node, "nodes", None) or ())
    return None


def _has_in_term(term, query, key, kind):
    """``_has_algebraic(term, query)`` with the equation and term named on failure."""
    try:
        return _has_algebraic(term, query)
    except AlgebraicError as err:
        node = _calculus_holding(term, query)
        where = f" inside '{node}'" if node is not None else ""
        if kind == "acceleration":
            what = (
                f"the acceleration '{query}' occurs only{where}, a third or higher time"
                " derivative of the configuration (a jerk): the equation is not second"
                " order, so it has no standard form M a + f + G u = 0"
            )
        else:
            what = (
                f"the input '{query}' occurs only{where}, a time derivative of the input:"
                " an input enters the standard form linearly and undifferentiated"
            )
        raise AlgebraicError(
            f"to_standard_form: in equation '{key}', term '{term}': {what} (engine: {err})"
        ) from err


def _extract(expr, query, key):
    """``extract_linear_coeff`` with the equation named in a re-raised engine error."""
    try:
        return extract_linear_coeff(expr, query)
    except (AlgebraicError, NonLinearError) as err:
        raise type(err)(
            f"to_standard_form: in equation '{key}', extracting the coefficient of"
            f" '{query}': {err}"
        ) from err


def _extract_single(eom_expr, accel_vars, inputs, key="?"):
    """Extract M, f, G from a single EOM expression (scalar or vector)."""
    eom_expr = expand(eom_expr)
    terms = _flatten_addends(eom_expr)

    accel_terms = []
    input_terms = []
    residual_terms = []

    # Partition with the same algebraic test extraction uses, so no term can be
    # classified as acceleration/input and then extracted as zero; a variable
    # hidden under a higher derivative (e.g. jerk) raises AlgebraicError here.
    for term in terms:
        has_accel = any(_has_in_term(term, a, key, "acceleration") for a in accel_vars)
        has_input = any(_has_in_term(term, u, key, "input") for u in inputs)

        if has_accel:
            accel_terms.append(term)
        elif has_input:
            input_terms.append(term)
        else:
            residual_terms.append(term)

    M = {}
    if accel_terms:
        accel_expr = _combine(accel_terms, eom_expr)
        for a in accel_vars:
            if _has_algebraic(accel_expr, a):
                coeff = _extract(accel_expr, a, key)
                coeff = full_simplify(coeff)
                # M must be free of accelerations: a term bilinear in two
                # accelerations (c ẍ₁ ẍ₂) has no standard-form coefficient.
                for b in accel_vars:
                    if _has_algebraic(coeff, b):
                        raise NonLinearError(
                            f"to_standard_form: in equation '{key}', the coefficient of"
                            f" '{a}' still contains the"
                            f" acceleration '{b}' in '{coeff}'; the equation is not linear"
                            " in the accelerations"
                        )
                M[str(a)] = coeff

    G = {}
    if input_terms:
        input_expr = _combine(input_terms, eom_expr)
        for u in inputs:
            if _has_algebraic(input_expr, u):
                coeff = _extract(input_expr, u, key)
                coeff = full_simplify(coeff)
                G[str(u)] = coeff

    if not residual_terms:
        f = _zero_for(eom_expr)
    else:
        f = full_simplify(_combine(residual_terms, eom_expr))

    return StandardFormEquation(M=M, f=f, G=G)
