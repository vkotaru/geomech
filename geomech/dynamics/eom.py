"""Equations of motion via the principle of least action.

compute_eom(lagrangian, inf_work, variables)
    variation → simplify → IBP → expand → extract coefficients

separate_variations(expr, variation_vectors)
    extract coefficient of each variation vector from the action integral
"""

from __future__ import annotations

from typing import Sequence

from geomech.core.base.expressions import S2, SO3, TS2, TSO3, Expr, Matrix, Scalar, Vector, Zero
from geomech.core.base.types import ExprType
from geomech.core.math.extract import extract_linear_coeff
from geomech.core.math.ibp import integrate_by_parts
from geomech.core.operations.calculus import TimeDerivative, Variation
from geomech.core.transformations.expand import expand
from geomech.core.transformations.manifold_rules import apply_manifold_rules
from geomech.core.transformations.simplify import full_simplify
from geomech.dynamics.variables import SystemVariables, check_group_is_sequence


def _tangent_names(variables: SystemVariables):
    """``{printed tangent vector (ω_q, Ω_R): printed point}`` of the listed S2/SO3 entries."""
    names = {}
    for group in (variables.vectors, variables.matrices):
        for var in group:
            if isinstance(var, (S2, SO3)):
                names[str(var.get_tangent_vector())] = str(var)
    return names


def _check_variables(variables: SystemVariables, caller: str = "compute_eom"):
    """Every configuration variable must be a distinct symbol with a non-vanishing δq.

    Rejected, naming the entry: a group that is not a list (``vectors=x``,
    ``TypeError``); a variable of the wrong type (a Vector in ``scalars``), a tangent or
    variation vector (``ω``, ``ξ``, ``Ω``, ``η``) or a ``TimeDerivative``/``Variation``
    node listed as a configuration variable (``TypeError``); a Constant/numeric leaf,
    whose δq = 0 would yield an empty equation (``ValueError``); a plain ``Matrix`` in
    ``matrices``, since only SO3 is supported (``TypeError``); and the same printed name
    listed twice, across groups or as a Scalar and a Vector, or an entry printed like
    another entry's variation or tangent vector, which would silently merge two symbols
    into one (``ValueError``). Messages start with *caller* (``compute_eom``/``linearize``).
    """
    groups = (
        ("scalars", variables.scalars, ExprType.SCALAR),
        ("vectors", variables.vectors, ExprType.VECTOR),
        ("matrices", variables.matrices, ExprType.MATRIX),
    )
    for group_name, group, _ in groups:
        try:
            check_group_is_sequence(group_name, group)
        except TypeError as e:
            raise TypeError(f"{caller}: {e}") from None
    tangent_names = _tangent_names(variables)  # printed ω_q/Ω_R → printed point
    seen_names = {}  # printed name → group it was listed in
    seen_variations = {}  # printed variation vector → printed variable
    for group_name, group, expected_type in groups:
        for var in group:
            where = f"{caller}: SystemVariables.{group_name} entry"
            if getattr(var, "type", None) != expected_type:
                raise TypeError(
                    f"{where} '{var}' has type {getattr(var, 'type', None)}, expected"
                    f" {expected_type}"
                )
            if isinstance(var, (TS2, TSO3)):
                parent = var.S2 if isinstance(var, TS2) else var.SO3
                hint = f"; list the manifold point '{parent}' instead" if parent else ""
                raise TypeError(
                    f"{where} '{var}' is a tangent or variation vector, not a configuration"
                    f" variable{hint}"
                )
            if isinstance(var, (TimeDerivative, Variation)):
                raise TypeError(
                    f"{where} '{var}' is a {type(var).__name__} node, not a configuration"
                    f" variable; list the underlying variable '{var.expr}' instead"
                )
            if var.is_constant or var.is_numeric or var.is_zero or var.delta().is_zero:
                raise ValueError(
                    f"{where} '{var}' is constant: its variation vanishes, so it cannot be a"
                    " configuration variable"
                )
            if group_name == "matrices" and not isinstance(var, SO3):
                raise TypeError(
                    f"{where} '{var}' is a general {type(var).__name__}; only SO3 rotation"
                    " matrices are supported as matrix configuration variables"
                )
            name = str(var)
            variation = str(
                var.delta() if expected_type == ExprType.SCALAR else var.get_variation_vector()
            )
            if name in seen_names:
                groups_text = (
                    f"twice in SystemVariables.{group_name}"
                    if seen_names[name] == group_name
                    else f"in SystemVariables.{seen_names[name]} and {group_name}"
                )
                raise ValueError(
                    f"{caller}: SystemVariables: the name '{name}' is listed more than once"
                    f" ({groups_text}); configuration variables must be distinct symbols"
                )
            if name in tangent_names:
                raise ValueError(
                    f"{where} '{name}' has the printed name of the tangent vector of"
                    f" '{tangent_names[name]}'; configuration variables must be distinct"
                    " symbols (rename it)"
                )
            if name in seen_variations:
                raise ValueError(
                    f"{where} '{name}' has the printed name of the variation vector of"
                    f" '{seen_variations[name]}'; configuration variables must be distinct"
                    " symbols"
                )
            if variation in seen_names:
                raise ValueError(
                    f"{where} '{name}' has the variation vector '{variation}', which is also"
                    f" the name of the entry listed in {seen_names[variation]};"
                    " configuration variables must be distinct symbols"
                )
            if variation in seen_variations:
                raise ValueError(
                    f"{caller}: SystemVariables: entries '{seen_variations[variation]}' and"
                    f" '{name}' share the variation vector '{variation}'; configuration"
                    " variables must be distinct symbols"
                )
            seen_names[name] = group_name
            seen_variations[variation] = name


def _check_scalar_input(label, obj):
    """``lagrangian``/``inf_work`` must be scalar expression trees."""
    if not isinstance(obj, Expr):
        raise TypeError(
            f"compute_eom: {label} must be a scalar expression, got {type(obj).__name__} {obj!r}"
        )
    if obj.type != ExprType.SCALAR:
        kind = getattr(obj.type, "name", obj.type)
        raise TypeError(
            f"compute_eom: {label} must be a scalar expression, got a {kind} expression '{obj}'"
        )


def _check_params(params, variables):
    """Prescribed parameters: distinct, non-constant plain leaves that are not variables.

    *params* must be an iterable of leaves (list, tuple, set, generator), not one leaf. A
    param printed like a configuration variable or like the tangent or variation vector of
    a listed ``S2``/``SO3`` point (``ω_q``, ``Ω_R``, ``ξ_q``, ``η_R``) raises: the tree
    would hold two different leaves with one printed name, which every str-keyed consumer
    merges.
    """
    if isinstance(params, (Expr, str)) or not hasattr(params, "__iter__"):
        hint = f"; write params=[{params}]" if isinstance(params, Expr) else ""
        raise TypeError(
            f"compute_eom: params must be a sequence of leaves, got {type(params).__name__}"
            f" {str(params)!r}{hint}"
        )
    reserved = {}  # printed tangent/variation vector of a listed point → description
    for group in (variables.vectors, variables.matrices):
        for var in group:
            if isinstance(var, (S2, SO3)):
                reserved[str(var.get_tangent_vector())] = f"tangent vector of '{var}'"
                reserved[str(var.get_variation_vector())] = f"variation vector of '{var}'"
    config_names = {}
    for group_name, group in (
        ("scalars", variables.scalars),
        ("vectors", variables.vectors),
        ("matrices", variables.matrices),
    ):
        for var in group:
            config_names[str(var)] = group_name
    seen = set()
    for p in params:
        if isinstance(p, (TS2, TSO3)):
            raise TypeError(
                f"compute_eom: params entry '{p}' is a tangent or variation vector, not a"
                " prescribed parameter"
            )
        if isinstance(p, (S2, SO3)):
            raise TypeError(
                f"compute_eom: params entry '{p}' is a manifold point, not a prescribed"
                " parameter; list it in SystemVariables"
            )
        if not isinstance(p, (Scalar, Vector, Matrix)):
            raise TypeError(
                f"compute_eom: params entry {p!r} must be a plain Scalar, Vector or Matrix"
                f" leaf, got {type(p).__name__}"
            )
        if p.is_constant or p.is_numeric:
            raise ValueError(
                f"compute_eom: params entry '{p}' is constant; a prescribed parameter varies"
                " in time (drop it from params or remove the 'Constant' attribute)"
            )
        name = str(p)
        if name in config_names:
            raise ValueError(
                f"compute_eom: params entry '{name}' is also a configuration variable"
                f" (SystemVariables.{config_names[name]}); a symbol is either prescribed or"
                " a degree of freedom"
            )
        if name in reserved:
            raise ValueError(
                f"compute_eom: params entry '{name}' has the printed name of the"
                f" {reserved[name]}; rename the parameter (two different symbols with one"
                " printed name are merged when the equation is printed or evaluated)"
            )
        if name in seen:
            raise ValueError(f"compute_eom: params entry '{name}' is listed more than once")
        seen.add(name)
    return seen


def _stray_variation_symbols(expr, allowed, param_names, wrapped=None):
    """Variation symbols in *expr* that belong to no listed variable or declared param.

    Returns ``{printed symbol: printed underlying symbol}`` for every ``Variation`` node
    whose printed form is not in *allowed* (the variation vectors of ``SystemVariables``)
    and whose operand is not a declared param, and for every ``TS2``/``TSO3`` leaf that is
    the variation vector (ξ, η) of a manifold point other than the listed ones. Such a
    symbol left in δS after extraction belongs to a degree of freedom that was not listed
    (or to a parameter that is neither ``Constant`` nor declared in ``params``): its terms
    were dropped from every equation without a message.

    A ``Variation`` of a ``TimeDerivative`` (a hand-built ``Variation(x.t_diff())``) is
    reported under the differentiated symbol ``x``; if *wrapped* is a dict, it also
    receives ``{printed symbol: x}`` for those, so the message can point to the variation
    of ``x`` itself.
    """
    stray = {}

    def walk(node):
        if isinstance(node, Variation):
            inner = node.expr
            if str(node) in allowed or str(inner) in param_names:
                return
            if isinstance(inner, TimeDerivative):
                # a hand-built δ(ẋ): name the differentiated symbol, not 'd/dt(x)'
                while isinstance(inner, TimeDerivative):
                    inner = inner.expr
                if wrapped is not None:
                    wrapped[str(node)] = str(inner)
            if isinstance(inner, (TS2, TSO3)):
                parent = inner.S2 if isinstance(inner, TS2) else inner.SO3
                stray[str(node)] = str(parent) if parent is not None else str(inner)
            else:
                stray[str(node)] = str(inner)
            return
        if isinstance(node, (TS2, TSO3)):
            parent = node.S2 if isinstance(node, TS2) else node.SO3
            if (
                parent is not None
                and node.name == parent.manifold_info.variation_vector_name
                and str(node) not in allowed
            ):
                stray[str(node)] = str(parent)
            return
        for child in getattr(node, "nodes", ()):
            walk(child)

    walk(expr)
    return stray


def _derivatives_of_variations(expr, allowed):
    """Outermost ``TimeDerivative`` nodes of *expr* whose subtree holds a listed variation.

    Returns ``{printed derivative: printed variation vector}``; *allowed* is the set of
    printed variation vectors of ``SystemVariables`` (``δx``, ``ξ_q``, ``η_R``).
    """
    found = {}

    def contained(node):
        if str(node) in allowed and isinstance(node, (Variation, TS2, TSO3)):
            return str(node)
        for child in getattr(node, "nodes", ()):
            hit = contained(child)
            if hit is not None:
                return hit
        return None

    def walk(node):
        if isinstance(node, TimeDerivative):
            hit = contained(node.expr)
            if hit is not None:
                found[str(node)] = hit
                return
        for child in getattr(node, "nodes", ()):
            walk(child)

    walk(expr)
    return found


def _check_virtual_work(inf_work, variables, allowed):
    """δW is linear in the variation vectors themselves, never in their time derivatives.

    ``Dot(x.t_diff().delta(), f)`` builds ``d/dt(δx)·f``; integration by parts would turn
    it into ``−δx·ḟ`` and the equation would carry ``−ḟ`` instead of ``f``. The manifold
    rules are applied to δW alone first, so ``δ(q̇)`` on S² (``δω → ξ̇ − …``) is caught as
    well. A velocity-dependent force (``Dot(x.delta(), ẋ c)``) has no variation under the
    derivative and passes. δL is not checked: ``ξ̇``/``η̇`` are legitimate there.
    """
    work = full_simplify(apply_manifold_rules(inf_work, variables))
    bad = _derivatives_of_variations(work, allowed)
    if bad:
        term, var = sorted(bad.items())[0]
        others = f" (and {len(bad) - 1} more)" if len(bad) > 1 else ""
        if term not in str(inf_work):
            others += " after the manifold kinematic substitutions"
        raise ValueError(
            f"compute_eom: inf_work contains '{term}'{others}, the time derivative of the"
            f" variation vector '{var}'; virtual work must be linear in the variations"
            " themselves, not their time derivatives (write Dot(x.delta(), f), not"
            " Dot(x.t_diff().delta(), f))"
        )


def compute_eom(
    lagrangian,
    inf_work,
    variables: SystemVariables,
    params: Sequence[Expr] = (),
):
    """Compute equations of motion using the principle of least action.

    1. Take the variation of the Lagrangian: δL
    2. Form the infinitesimal action integral: δS = δL + δW
    3. Apply manifold kinematic substitutions (e.g. δ(ω) for S2)
    4. Simplify
    5. Gather variation vectors and their time derivatives
    6. Integration by parts to move time derivatives off variation vectors
    7. Expand
    8. Extract coefficients of each independent variation vector

    ``inf_work`` may be ``None`` or ``0`` for a conservative system. It must be linear in
    the variation vectors themselves: a time derivative of a listed variation vector in
    ``inf_work`` (``Dot(x.t_diff().delta(), f)``, ``Dot(ξ.t_diff(), f)``) raises
    ``ValueError``. ``variables`` must list at least one configuration variable
    (``ValueError`` otherwise); ``params`` must be a sequence of leaves (``TypeError``).

    Every symbol in ``lagrangian`` and ``inf_work`` is one of three things: a
    configuration variable (listed in ``variables``, varied by Hamilton's principle), a
    constant (``attr=["Constant"]``), or a prescribed time-varying parameter declared in
    ``params`` (a plain non-constant ``Scalar``/``Vector``/``Matrix`` leaf such as a cable
    length ``l(t)``: its rate ``l̇`` appears in the equations, its variation ``δl`` is not
    free and is dropped). A variation symbol that survives extraction and belongs to none
    of these (``δx`` of an unlisted ``x``, ``ξ``/``δω`` of an unlisted ``S2`` point, …)
    raises ``ValueError`` naming it: before, such terms were dropped silently and the
    system lost a degree of freedom. ``numeric.check_eom`` takes the same ``params``.
    ``to_standard_form`` accepts the result (a param's rate stays in ``f``, the param
    itself in ``M``); ``linearize`` does not support prescribed parameters yet and raises
    ``AlgebraicError`` on the ``δl`` term.
    """
    _check_variables(variables)
    param_names = _check_params(params, variables)
    if isinstance(inf_work, bool):
        raise TypeError(
            f"compute_eom: inf_work must be a scalar expression, None or 0, got bool {inf_work}"
        )
    if inf_work is None or (isinstance(inf_work, (int, float)) and inf_work == 0):
        inf_work = Zero
    _check_scalar_input("lagrangian", lagrangian)
    _check_scalar_input("inf_work", inf_work)

    # gather variation vectors and their time derivatives
    variation_vectors = []
    variation_vector_dots = []

    for s in variables.scalars:
        ds = s.delta()  # Variation(scalar) — scalar type
        variation_vectors.append(ds)
        variation_vector_dots.append(ds.t_diff())

    for vec in variables.vectors:
        x = vec.get_variation_vector()
        variation_vectors.append(x)
        variation_vector_dots.append(x.t_diff())

    for mat in variables.matrices:
        x = mat.get_variation_vector()
        variation_vectors.append(x)
        variation_vector_dots.append(x.t_diff())

    allowed = {str(v) for v in variation_vectors}
    _check_virtual_work(inf_work, variables, allowed)

    # variation of the Lagrangian
    dL = lagrangian.delta()

    # infinitesimal action integral
    dS = dL + inf_work

    # apply manifold kinematic substitutions
    dS = apply_manifold_rules(dS, variables)

    dS = full_simplify(dS)

    # integration by parts
    dS = integrate_by_parts(dS, variation_vector_dots)

    # expand
    dS = expand(dS)

    # extract equations of motion
    eom = separate_variations(dS, variation_vectors)

    # a variation symbol that survives extraction belongs to no listed variable or param
    wrapped = {}  # δ(d/dt(x)) built by hand → x
    stray = _stray_variation_symbols(dS, allowed, param_names, wrapped)
    for _, eq in eom.values():
        stray.update(_stray_variation_symbols(eq, allowed, param_names, wrapped))
    if stray:
        symbols = ", ".join(f"'{s}'" for s in sorted(stray))
        names = ", ".join(f"'{n}'" for n in sorted(set(stray.values())))
        message = (
            f"compute_eom: δS still contains the variation symbols {symbols} of the symbols"
            f" {names}, which are neither entries of SystemVariables nor declared in params"
        )
        if any(sym not in wrapped for sym in stray):
            message += (
                "; for each symbol, list it in SystemVariables, declare it in params=..., or"
                " mark it 'Constant'"
            )
        for sym, leaf in sorted(wrapped.items()):
            message += (
                f"; '{sym}' is the variation of a time derivative, which the engine never"
                f" builds: write the variation of the configuration variable itself,"
                f" {leaf}.delta() (the variation of the velocity is {leaf}.t_diff().delta(),"
                " built as d/dt of the variation, and is not allowed in inf_work)"
            )
        raise ValueError(message)
    if not variation_vectors:
        raise ValueError(
            "compute_eom: no configuration variables listed in SystemVariables; a system"
            " with no degree of freedom has no equation of motion"
        )
    return eom


def separate_variations(inf_action_integral, variation_vectors):
    """Extract the equation of motion for each variation vector.

    Returns a dict mapping str(variation_vector) → (variation_vector, equation).
    """
    eom = {}
    for vec in variation_vectors:
        dyn_eqn = extract_linear_coeff(inf_action_integral, vec)
        dyn_eqn = full_simplify(dyn_eqn)
        eom[str(vec)] = (vec, dyn_eqn)
    return eom
