"""Equations of motion via the principle of least action.

compute_eom(lagrangian, inf_work, variables, params=(), *, trace=None)
    δL + δW → apply_manifold_rules → full_simplify → integrate_by_parts (expand,
    full_simplify, then per target: collect → IBP → full_simplify) → expand →
    extract_linear_coeff per variation vector → full_simplify.
    ``docs/pipeline.md`` describes each stage; ``trace=EOMTrace()`` keeps them all.

separate_variations(expr, variation_vectors)
    extract coefficient of each variation vector from the action integral
"""

from __future__ import annotations

from typing import Sequence

from geomech.core.base.expressions import S2, SO3, TS2, TSO3, Expr, Matrix, Scalar, Vector, Zero
from geomech.core.base.types import ExprType
from geomech.core.math.extract import extract_linear_coeff
from geomech.core.math.ibp import _integrate_by_parts
from geomech.core.operations.calculus import TimeDerivative, Variation
from geomech.core.operations.geometry import Dot
from geomech.core.operations.multiplication import VVMul
from geomech.core.transformations.expansion import expand
from geomech.core.transformations.manifold_rules import apply_manifold_rules
from geomech.core.transformations.simplification import full_simplify
from geomech.core.transformations.utils import addends, flatten_mul
from geomech.dynamics.variables import SystemVariables, check_group_is_sequence
from geomech.utils.errors import AlgebraicError, NonLinearError
from geomech.utils.invariants import check_tree
from geomech.utils.trace import EOMTrace, debug_enabled, stage_recorder


def _tangent_names(variables: SystemVariables):
    """``{printed tangent vector (ω_q, Ω_R): printed point}`` of the listed S2/SO3 entries."""
    names = {}
    for group in (variables.vectors, variables.matrices):
        for var in group:
            if isinstance(var, (S2, SO3)):
                names[str(var.get_tangent_vector())] = str(var)
    return names


def check_variables(variables: SystemVariables, caller: str = "compute_eom") -> None:
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
    seen_names: dict[str, str] = {}  # printed name → group it was listed in
    seen_variations: dict[str, str] = {}  # printed variation vector → printed variable
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
            if getattr(var, "nodes", None):
                raise TypeError(
                    f"{where} '{var}' is a {type(var).__name__} expression, not a single"
                    " symbol; a configuration variable is one Scalar, Vector, S2 or SO3 leaf"
                    " (list the symbols it is built from instead)"
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


def _symbol_key(leaf):
    """What makes two leaves of one printed name the same symbol: class, flags, size."""
    return type(leaf), getattr(leaf, "flags", None), getattr(leaf, "size", None)


def _describe_leaf(leaf):
    """``constant Scalar`` / ``Vector`` / ``Vector with attr=['UnitNorm']`` for a message."""
    kind = type(leaf).__name__
    if leaf.is_constant or leaf.is_numeric:
        return f"constant {kind}"
    attr = getattr(leaf, "attr", None)
    return f"{kind} with attr={list(attr)}" if attr else kind


def check_symbol_clash(
    symbols: Sequence[Expr],
    exprs: Sequence[object],
    caller: str,
    role: str = "input",
    what: str = "an input",
    where: str = "the equations",
) -> None:
    """Reject a declared symbol printed like a *different* leaf of *exprs*.

    *symbols* are the inputs of ``to_standard_form``/``linearize`` or the ``params`` of
    ``compute_eom``; *exprs* the equations (or the Lagrangian and the virtual work). The
    engine identifies a leaf by its printed name, so a declared ``Scalar("g")`` would be
    matched against, or silently miss, a constant ``Scalar("g", attr=['Constant'])`` in
    the equations. A declared symbol whose printed name equals that of a leaf of another
    class, other flags (``Constant``, ``UnitNorm``, ...) or another size raises
    ``ValueError`` naming both; the same symbol (the normal case) passes. *role*, *what*
    and *where* word the message (``params entry``, ``a prescribed parameter``, ``the
    Lagrangian or the virtual work``).
    """
    leaves: dict[str, list] = {}
    stack = list(exprs)
    while stack:
        node = stack.pop()
        children = getattr(node, "nodes", None)
        if children:
            stack.extend(children)
            continue
        same_name = leaves.setdefault(str(node), [])
        if all(_symbol_key(node) != _symbol_key(other) for other in same_name):
            same_name.append(node)
    for sym in symbols:
        name = str(sym)
        for leaf in leaves.get(name, ()):
            if _symbol_key(leaf) == _symbol_key(sym):
                continue
            if (leaf.is_constant or leaf.is_numeric) and type(leaf) is type(sym):
                raise ValueError(
                    f"{caller}: {role} '{name}' has the name of the constant '{name}' in"
                    f" {where}; a symbol is either a constant parameter or {what} (rename"
                    f" one, or drop attr=['Constant'] if '{name}' is {what})"
                )
            raise ValueError(
                f"{caller}: {role} '{name}' ({_describe_leaf(sym)}) has the name of the"
                f" {_describe_leaf(leaf)} '{name}' in {where}, a different symbol;"
                " leaves are matched by printed name, so the two would be confused (rename"
                f" one, or build the {role} with the class and attr it has in {where})"
            )


def _input_state_symbols(variables: SystemVariables) -> set[str]:
    """Printed perturbation states of *variables* and their time derivatives.

    The configuration perturbations (``δx``, ``δs``, ``ξ_q``, ``η_R``), the velocity
    perturbations (``δẋ``, ``δṡ``, ``δω_q``, ``δΩ_R``) and the d/dt of each: the columns
    of ``linearize`` other than the input variations. A general ``Matrix`` (rejected by
    ``check_variables``) contributes nothing.
    """
    states: list[Expr] = []
    for s in variables.scalars:
        states += [s.delta(), s.t_diff().delta()]
    for v in variables.vectors:
        if isinstance(v, S2):
            states += [v.get_variation_vector(), Variation(v.get_tangent_vector())]
        else:
            states += [v.delta(), v.t_diff().delta()]
    for m in variables.matrices:
        if isinstance(m, SO3):
            states += [m.get_variation_vector(), Variation(m.get_tangent_vector())]
    return {str(x) for x in states} | {str(x.t_diff()) for x in states}


def check_inputs(
    inputs: Sequence[Expr] | None,
    variables: SystemVariables,
    exprs: Sequence[object],
    caller: str,
) -> list[Expr]:
    """Validate the ``inputs`` of ``to_standard_form``/``linearize``; return them as a list.

    ``None`` means no inputs. *inputs* must be a sequence of leaves, not one leaf
    (``TypeError``). Each input must be a plain non-manifold ``Scalar`` or ``Vector`` leaf
    (``TypeError`` for an ``S2``/``SO3`` point, a ``Matrix`` or an expression such as
    ``d/dt(x)`` or ``F + u``), and is rejected with ``ValueError`` if it is a tangent or
    variation vector (``ω``, ``ξ``, ``Ω``, ``η``), constant (``attr=['Constant']``, a
    number, zero), a configuration variable, a symbol whose variation prints like a
    perturbation state (a plain ``Vector`` named like ``Ω``), listed twice, or printed
    like a different leaf of *exprs* (``check_symbol_clash``). Messages start with
    *caller*. Both functions call this, so they accept and reject the same inputs.
    """
    if inputs is None:
        return []
    if isinstance(inputs, (Expr, str)) or not hasattr(inputs, "__iter__"):
        hint = f"; write inputs=[{inputs}]" if isinstance(inputs, Expr) else ""
        raise TypeError(
            f"{caller}: inputs must be a sequence of Scalar/Vector leaves, got"
            f" {type(inputs).__name__} {str(inputs)!r}{hint}"
        )
    inputs = list(inputs)
    config = {str(v) for v in variables.scalars + variables.vectors + variables.matrices}
    states = _input_state_symbols(variables)
    seen = set()
    for u in inputs:
        if not isinstance(u, (Scalar, Vector)) or isinstance(u, S2):
            raise TypeError(
                f"{caller}: input '{u}' must be a Scalar or Vector leaf, got {type(u).__name__}"
            )
        if isinstance(u, (TS2, TSO3)):
            raise ValueError(
                f"{caller}: input '{u}' is a tangent or variation vector of a configuration"
                " variable; its variation is part of the perturbation state, not an input"
            )
        if u.is_constant or u.is_numeric or u.is_zero:
            raise ValueError(
                f"{caller}: input '{u}' is constant: its variation vanishes,"
                " so it cannot be an input"
            )
        if str(u) in config:
            raise ValueError(f"{caller}: input '{u}' is also a configuration variable")
        if str(u.delta()) in states:
            raise ValueError(
                f"{caller}: the variation '{u.delta()}' of input '{u}' is already a"
                f" perturbation state (is '{u}' a tangent vector of a configuration"
                " variable?)"
            )
        if str(u) in seen:
            raise ValueError(f"{caller}: input '{u}' is listed more than once")
        seen.add(str(u))
    check_symbol_clash(inputs, exprs, caller)
    return inputs


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
    lagrangian: Expr,
    inf_work: Expr | int | float | None,
    variables: SystemVariables,
    params: Sequence[Expr] = (),
    *,
    trace: EOMTrace | None = None,
) -> dict[str, tuple[Expr, Expr]]:
    """Compute equations of motion using the principle of least action.

    Stages, in the order the code runs them (``docs/pipeline.md`` has one paragraph per
    stage; the bracketed names are the ``EOMTrace`` stage names):

    1. ``[dL]`` the variation of the Lagrangian, ``[dW]`` the virtual work, and
       ``[dS]`` the infinitesimal action integral δS = δL + δW
    2. ``[manifold_rules]`` ``apply_manifold_rules``: δω and δΩ replaced by their
       kinematic expansions (S², SO(3))
    3. ``[simplify]`` ``full_simplify`` (expand once, then simplify to a fixed point)
    4. integration by parts against the time derivative of each variation vector:
       ``[ibp.prepare]`` expand + full_simplify, then per target ``[ibp.collect[t]]``
       (rotate each term so the target is the left operand of a ``Dot``),
       ``[ibp.integrate[t]]`` (one IBP pass) and ``[ibp.simplify[t]]``; ``[ibp]`` is
       the result
    5. ``[expand]``
    6. per variation vector: ``extract_linear_coeff`` ``[raw[ξ]]``, then
       ``full_simplify`` ``[equation[ξ]]``

    ``trace`` (keyword only) is an ``EOMTrace``: it is filled with every stage as the
    pipeline runs, so after an exception it holds the stages before the failure and
    ``trace.failed_at`` names the stage that raised (``trace.summary()``,
    ``trace.show("simplify")``). The trace is cleared first, so a reused trace holds
    this run only. Without a trace, and with the
    ``geomech`` logger not at ``DEBUG``, nothing is recorded. With ``GEOMECH_DEBUG=1``
    each stage is logged and ``check_tree`` runs on the inputs and the equations.

    An extraction failure is re-raised with the same exception type and the pipeline
    context: the variation vector, and when a time derivative of a variation vector
    survived integration by parts, the term it survived in.

    ``inf_work`` may be ``None`` or ``0`` for a conservative system. It must be linear in
    the variation vectors themselves: a time derivative of a listed variation vector in
    ``inf_work`` (``Dot(x.t_diff().delta(), f)``, ``Dot(ξ.t_diff(), f)``) raises
    ``ValueError``. ``variables`` must list at least one configuration variable
    (``ValueError`` "no configuration variables" otherwise; this check runs last, so a
    Lagrangian that holds a symbol first fails the stray-variation check below);
    ``params`` must be a sequence of leaves (``TypeError``).

    Every symbol in ``lagrangian`` and ``inf_work`` is one of three things: a
    configuration variable (listed in ``variables``, varied by Hamilton's principle), a
    constant (``attr=["Constant"]``), or a prescribed time-varying parameter declared in
    ``params`` (a plain non-constant ``Scalar``/``Vector``/``Matrix`` leaf such as a cable
    length ``l(t)``: its rate ``l̇`` appears in the equations, its variation ``δl`` is not
    free and is dropped). A variation symbol that survives extraction and belongs to none
    of these (``δx`` of an unlisted ``x``, ``ξ``/``δω`` of an unlisted ``S2`` point, …)
    raises ``ValueError`` naming it: before, such terms were dropped silently and the
    system lost a degree of freedom. ``numeric.check_eom`` takes the same ``params``. A
    param printed like a different leaf of ``lagrangian`` or ``inf_work`` (a constant
    ``l`` of ``L`` declared as ``Scalar("l")``) raises ``ValueError``: it would have no
    effect.
    ``to_standard_form`` accepts the result (a param's rate stays in ``f``, the param
    itself in ``M``); ``linearize`` does not support prescribed parameters yet and raises
    ``AlgebraicError`` on the ``δl`` term.
    """
    if trace is not None:
        trace._reset()
    stage = "validation"
    try:
        check_variables(variables)
        if not isinstance(params, (Expr, str)) and hasattr(params, "__iter__"):
            params = list(params)  # a generator is read twice (here and the clash check)
        param_names = _check_params(params, variables)
        if isinstance(inf_work, bool):
            raise TypeError(
                "compute_eom: inf_work must be a scalar expression, None or 0, got bool"
                f" {inf_work}"
            )
        if inf_work is None or (isinstance(inf_work, (int, float)) and inf_work == 0):
            inf_work = Zero
        _check_scalar_input("lagrangian", lagrangian)
        _check_scalar_input("inf_work", inf_work)
        check_symbol_clash(
            params,
            [lagrangian, inf_work],
            "compute_eom",
            role="params entry",
            what="a prescribed parameter",
            where="the Lagrangian or the virtual work",
        )

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

        debug = debug_enabled()
        if debug:
            check_tree(lagrangian)
            check_tree(inf_work)
        record = stage_recorder(trace, "compute_eom")
        if trace is not None:
            trace.variation_vectors = list(variation_vectors)
            trace.ibp_targets = list(variation_vector_dots)

        # variation of the Lagrangian
        stage = "dL"
        dL = lagrangian.delta()

        # infinitesimal action integral
        dS = dL + inf_work
        if record is not None:
            record("dL", dL)
            record("dW", inf_work)
            record("dS", dS)

        # apply manifold kinematic substitutions
        stage = "manifold_rules"
        dS = apply_manifold_rules(dS, variables)
        if record is not None:
            record("manifold_rules", dS)

        stage = "simplify"
        dS = full_simplify(dS)
        if record is not None:
            record("simplify", dS)

        # integration by parts
        stage = "ibp"
        try:
            dS = _integrate_by_parts(dS, variation_vector_dots, record)
        except (AlgebraicError, NonLinearError) as err:
            raise type(err)(
                f"compute_eom: while integrating by parts (stage ibp, targets"
                f" {', '.join(f'{t}' for t in variation_vector_dots)}): {err}"
            ) from err
        if record is not None:
            record("ibp", dS)

        # expand
        stage = "expand"
        dS = expand(dS)
        if record is not None:
            record("expand", dS)

        # extract equations of motion
        stage = "extract"
        try:
            eom = _separate_variations(dS, variation_vectors, record)
        except _ExtractionFailure as failure:
            # retold with L when it depends on accelerations; otherwise the δS message
            raise failure.retold(dS, variation_vectors, lagrangian, variables) from (
                failure.engine
            )
        if debug:
            for _, eq in eom.values():
                check_tree(eq)

        stage = "stray_check"
        # a variation symbol that survives extraction belongs to no listed variable or param
        wrapped: dict[str, str] = {}  # δ(d/dt(x)) built by hand → x
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
    except Exception as err:
        if trace is not None:
            trace._fail(_failed_stage(trace, stage, err), err)
        raise


def separate_variations(
    inf_action_integral: Expr, variation_vectors: Sequence[Expr]
) -> dict[str, tuple[Expr, Expr]]:
    """Extract the equation of motion for each variation vector.

    Returns a dict mapping str(variation_vector) → (variation_vector, equation). An
    ``AlgebraicError``/``NonLinearError`` from the extraction is re-raised (same type)
    naming the variation vector and, if a time derivative of a variation vector is left
    in δS, the term that holds it.
    """
    try:
        return _separate_variations(inf_action_integral, variation_vectors)
    except _ExtractionFailure as failure:
        raise failure.retold(inf_action_integral, variation_vectors) from failure.engine


class _ExtractionFailure(Exception):
    """An extraction failure of ``_separate_variations``, with what is needed to explain it.

    *engine* is the ``AlgebraicError``/``NonLinearError`` of ``extract_linear_coeff`` and
    *vec* the variation vector it was extracting. The callers turn it into an error of the
    engine's type with ``retold``: ``separate_variations`` from δS alone, ``compute_eom``
    with the Lagrangian too (so the message can say that L depends on accelerations).
    It never leaves this module.
    """

    def __init__(self, engine, vec):
        super().__init__(str(engine))
        self.engine = engine
        self.vec = vec

    def retold(self, dS, variation_vectors, lagrangian=None, variables=None):
        return _extraction_error(
            self.engine, dS, self.vec, variation_vectors, lagrangian, variables
        )


def _separate_variations(dS, variation_vectors, record=None):
    """``separate_variations`` with the stage recorder; raises ``_ExtractionFailure``."""
    eom = {}
    for vec in variation_vectors:
        try:
            raw = extract_linear_coeff(dS, vec)
        except (AlgebraicError, NonLinearError) as err:
            raise _ExtractionFailure(err, vec) from err
        if record is not None:
            record(f"raw[{vec}]", raw)
        dyn_eqn = full_simplify(raw)
        if record is not None:
            record(f"equation[{vec}]", dyn_eqn)
        eom[str(vec)] = (vec, dyn_eqn)
    return eom


def _failed_stage(trace, stage, err):
    """The stage name a failed ``compute_eom`` run reports in ``trace.failed_at``.

    *stage* is the coarse step that was running; for integration by parts and the
    extraction it is refined from the stages already recorded to the one that raised
    (``ibp.collect[t]``, ``raw[ξ]`` ...), which is the first one missing.
    """
    names = set(trace.stage_names())
    if stage == "ibp":
        if "ibp.prepare" not in names:
            return "ibp.prepare"
        for target in trace.ibp_targets or ():
            for step in ("collect", "integrate", "simplify"):
                if f"ibp.{step}[{target}]" not in names:
                    return f"ibp.{step}[{target}]"
    if stage == "extract":
        for vec in trace.variation_vectors or ():
            for step in ("raw", "equation"):
                if f"{step}[{vec}]" not in names:
                    return f"{step}[{vec}]"
    return stage


def _nodes(expr):
    stack = [expr]
    while stack:
        node = stack.pop()
        yield node
        stack.extend(getattr(node, "nodes", None) or ())


def _derivative_order(derivative, var):
    """Most ``d/dt`` nodes on a path from *derivative* down to the variation *var*."""
    best = 0
    stack = [(derivative, 0)]
    while stack:
        node, depth = stack.pop()
        if isinstance(node, TimeDerivative):
            depth += 1
        if str(node) == var and isinstance(node, (Variation, TS2, TSO3)):
            best = max(best, depth)
            continue
        stack.extend((child, depth) for child in getattr(node, "nodes", None) or ())
    return best


def _is_exposed(term, derivative, scalar):
    """True if *derivative* already has the shape IBP integrates in *term*: the left
    operand of a ``Dot`` (vector target) or a direct factor of a product (scalar)."""
    if scalar:
        return any(str(f) == derivative for f in flatten_mul(term))
    return any(isinstance(n, Dot) and str(n.left) == derivative for n in _nodes(term))


def _accelerations_in(lagrangian, variables):
    """Printed accelerations of configuration variables that occur in *lagrangian*:
    ``ẍ`` of a scalar or vector variable, ``ω̇``/``Ω̇`` of an S² or SO(3) point."""
    if lagrangian is None or variables is None:
        return []
    wanted = set()
    for var in (*variables.scalars, *variables.vectors, *variables.matrices):
        if isinstance(var, (S2, SO3)):
            wanted.add(str(var.get_tangent_vector().t_diff()))
        else:
            wanted.add(str(var.t_diff().t_diff()))
    return sorted({str(n) for n in _nodes(lagrangian) if str(n) in wanted})


def _extraction_error(err, dS, vec, variation_vectors, lagrangian=None, variables=None):
    """*err* (from ``extract_linear_coeff(dS, vec)``) rebuilt with the pipeline context.

    When a time derivative of a variation vector is left in δS, the message first tells
    the cause apart:

    - the derivative is of order 2 or more (``d/dt(d/dt(δx))``), or the Lagrangian holds
      an acceleration of a configuration variable: the Lagrangian depends on
      accelerations, which is not supported (integration by parts runs once per variation
      vector, against its first time derivative); ``_higher_order_error``;
    - a first derivative already in the shape IBP integrates (the left operand of its
      ``Dot``; for a scalar variation, a direct factor of its product): the IBP pass left
      it, and the message points to ``ibp.integrate[t]``; ``_ibp_left_error``;
    - a first derivative that is not yet in that shape: ``collect`` could not expose it,
      and the message points to the ``ibp.collect[t]`` stage; ``_collect_error``.

    Any other failure keeps the engine message after the context.
    """
    allowed = {str(v) for v in variation_vectors}
    by_name = {str(v): v for v in variation_vectors}
    found = []  # (term, derivative, variation vector), terms holding vec's derivative first
    for term in addends(dS):
        for derivative, var in sorted(_derivatives_of_variations(term, allowed).items()):
            found.append((var != str(vec), term, derivative, var))
    accelerations = _accelerations_in(lagrangian, variables)
    for _, term, derivative, var in sorted(found, key=lambda item: item[0]):
        node = next(
            n for n in _nodes(term) if isinstance(n, TimeDerivative) and str(n) == derivative
        )
        order = _derivative_order(node, var)
        target = str(by_name[var].t_diff())  # the key of the ibp.* trace stages
        if order >= 2 or accelerations:
            return _higher_order_error(
                err, vec, term, derivative, var, order, target, accelerations
            )
        scalar = by_name[var].type == ExprType.SCALAR
        if _is_exposed(term, derivative, scalar):
            return _ibp_left_error(err, vec, term, derivative, var, target, scalar)
        return _collect_error(err, vec, term, derivative, var, target, scalar)
    return type(err)(
        f"compute_eom: while extracting the equation for '{vec}' from δS (extraction, after"
        f" the stage 'expand'): {err}"
    )


def _higher_order_error(err, vec, term, derivative, var, order, target, accelerations):
    """The Lagrangian depends on accelerations (or the derivative has order 2 or more)."""
    if accelerations:
        why = "L contains " + ", ".join(f"'{a}'" for a in accelerations)
    else:
        why = f"'{derivative}' is a time derivative of order {order} of '{var}'"
    return type(err)(
        f"compute_eom: while extracting the equation for '{vec}': the Lagrangian"
        f" depends on accelerations ({why}); higher-order Lagrangians are not"
        " supported. Integration by parts runs once per variation vector, against"
        f" its first time derivative '{target}', so the term '{_cut(term)}' of δS"
        f" keeps '{derivative}' and has no coefficient. Write L in the"
        " configuration variables and their velocities only; run"
        " compute_eom(..., trace=EOMTrace()) and look at the stage 'dL' to see"
        f" where the acceleration enters. (engine: {err})"
    )


def _ibp_left_error(err, vec, term, derivative, var, target, scalar):
    """A first derivative already in the integrable shape survived the IBP pass."""
    shape = f"{derivative} g" if scalar else f"Dot({derivative}, f)"
    return type(err)(
        f"compute_eom: while extracting the equation for '{vec}': the time"
        f" derivative '{derivative}' of the variation vector '{var}' survived"
        f" integration by parts in the term '{_cut(term)}' of δS, although it is"
        f" already in the shape {shape} that the integration-by-parts pass"
        " integrates. Run compute_eom(..., trace=EOMTrace()) and compare the stages"
        f" 'ibp.collect[{target}]' and 'ibp.integrate[{target}]'. (engine: {err})"
    )


def _collect_error(err, vec, term, derivative, var, target, scalar):
    """A first derivative that ``collect`` could not bring into the integrable shape."""
    hint = ""
    if _has_inner_product_row(term):
        hint = (
            " The term holds a row-vector product (v.T() * A * w); collect turns an inner"
            " product that holds the derivative into a Dot, but this one could not be"
            " exposed: write it as Dot(v, A * w)."
        )
    if scalar:
        where = f"a direct factor of a product ({derivative} g)"
    else:
        where = f"the left operand of a Dot (Dot({derivative}, f))"
    return type(err)(
        f"compute_eom: while extracting the equation for '{vec}': the time derivative"
        f" '{derivative}' of the variation vector '{var}' survived integration by parts"
        f" in the term '{_cut(term)}' of δS; collect could not expose it as {where},"
        f" so the term was not integrated by parts and has no coefficient.{hint}"
        " Inspect the term with print_tree(term, style='indent'), or run"
        " compute_eom(..., trace=EOMTrace()) and look at the stage"
        f" 'ibp.collect[{target}]'. (engine: {err})"
    )


def _has_inner_product_row(term):
    stack = [term]
    while stack:
        node = stack.pop()
        if isinstance(node, VVMul) and node.type == ExprType.SCALAR:
            return True
        stack.extend(getattr(node, "nodes", None) or ())
    return False


def _cut(expr, limit=300):
    text = str(expr)
    return text if len(text) <= limit else text[: limit - 3] + "..."
