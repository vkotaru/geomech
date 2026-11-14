"""Variation-based linearization (VBL) of the equations of motion.

``linearize(eom, variables, inputs)`` takes the equations returned by
``compute_eom`` and produces their linearization along an arbitrary trajectory
in the first-order form

    M δẋ = A δx + B δu

where the perturbation state δx is built from manifold-respecting perturbation
variables rather than local coordinates (Wu & Sreenath 2015; Lee et al.
Appendix B). The blocks ``M``, ``A``, ``B`` are symbolic functions of the
trajectory ``(q, v, v̇, u)``.

Perturbation variables
----------------------

Each configuration variable contributes one configuration perturbation and one
velocity perturbation. Velocity perturbations are the ``Variation`` of the
velocity symbol and are **not** substituted by the EOM manifold rules
(``apply_manifold_rules`` is never called here):

    configuration      config. perturbation          velocity perturbation
    s : Scalar (R)     δs  = Variation(s)            δṡ = TimeDerivative(Variation(s))
    x : Vector (R³)    δx  = Variation(x)            δẋ = TimeDerivative(Variation(x))
    q : S2             ξ   = q.get_variation_vector  δω = Variation(ω),  ω = q.get_tangent_vector()
    R : SO3            η   = R.get_variation_vector  δΩ = Variation(Ω),  Ω = R.get_tangent_vector()

δ(constant) = 0, δq = ξ × q and δR = R·Hat(η) are produced eagerly by the leaf
types, and δ(v̇) is canonicalised to ``TimeDerivative(Variation(v))``.

Kinematic rows (one per configuration variable)
------------------------------------------------

- R and R³: d/dt(δs) = δṡ and d/dt(δx) = δẋ (identities between the same
  symbols, included so that the block structure is uniform).
- S² with q̇ = ω × q, ω ⊥ q, δq = ξ × q, ξ ⊥ q: equating δ(q̇) with d/dt(δq)
  and using the tangency constraint d/dt(ξ·q) = 0 gives

      ξ̇ = q qᵀ Hat(ω) ξ + (I − q qᵀ) δω

  so ``M[ξ̇] = I``, ``A[ξ] = q qᵀ Hat(ω)``, ``A[δω] = I − q qᵀ``.
- SO(3) with Ṙ = R·Hat(Ω), δR = R·Hat(η): δΩ = η̇ + Ω × η rearranged,

      η̇ = −Hat(Ω) η + δΩ

  so ``M[η̇] = I``, ``A[η] = −Hat(Ω)``, ``A[δΩ] = I``.

Dynamic rows (one per equation of ``compute_eom``) and sign convention
---------------------------------------------------------------------

``compute_eom`` returns, for each variation vector, an expression
``E(q, v, v̇, u)`` with ``E = 0``. ``E`` is the coefficient of the variation in
the post-IBP integrand, i.e. the *negative* of the Euler–Lagrange expression,
which is why ``to_standard_form`` reports ``M = −(physical inertia)``. VBL takes
the first variation δE, which is linear in the perturbation variables, and
defines the blocks so that the physical inertia comes out positive:

    δE = −M δv̇ + A δx + B δu

    M[row][δv̇_j] = −coeff(δv̇_j, δE)
    A[row][δx_j] =  coeff(δx_j, δE)      (δx_j over configuration and velocity perturbations)
    B[row][δu_j] =  coeff(δu_j, δE)

Free rigid body (Lee et al. Eq. 6.16): E = −JΩ̇ − Ω × JΩ + M_in gives
δE = −J δΩ̇ + (Hat(JΩ) − Hat(Ω) J) δΩ + δM, so ``M = J``,
``A[δΩ] = Hat(JΩ) − Hat(Ω) J`` and ``B = I``.

Completeness rule
-----------------

After extraction every additive term of ``expand(δE)`` must contain exactly one
perturbation symbol (configuration, velocity, d/dt of velocity, or input
variation). A term with none, or with a perturbation symbol only inside a
non-canonical calculus node, means the linearization is incomplete and
``linearize`` raises ``AlgebraicError`` naming the term. Nothing is silently
dropped or set to zero.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from geomech.core.base.expressions import (
    S2,
    SO3,
    TS2,
    TSO3,
    Expr,
    IdentityMatrix,
    One,
    Scalar,
    Vector,
    Zero,
    ZeroMatrix,
    ZeroVector,
)
from geomech.core.base.types import ExprType
from geomech.core.math.extract import _has_algebraic, _same, extract_linear_coeff
from geomech.core.operations.addition import Add, MAdd, VAdd
from geomech.core.operations.calculus import TimeDerivative, TimeIntegral, Variation
from geomech.core.operations.geometry import Hat, Transpose
from geomech.core.operations.multiplication import MMMul, Mul, SMMul, VVMul
from geomech.core.transformations.expansion import expand
from geomech.core.transformations.simplification import full_simplify
from geomech.dynamics.eom import _check_variables
from geomech.dynamics.variables import SystemVariables
from geomech.utils.errors import AlgebraicError
from geomech.utils.trace import debug_enabled, stage_recorder

# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class LinearizedSystem:
    """Linearized dynamics ``M δẋ = A δx + B δu`` as symbolic blocks.

    Attributes
    ----------
    rows : list[str]
        Ordered row labels: kinematic rows first (``"kin:" + str(var)``, one per
        configuration variable in ``SystemVariables`` order), then the dynamic rows
        labelled by the ``compute_eom`` keys.
    states : list[Expr]
        Perturbation symbols of δx: configuration perturbations in
        ``SystemVariables`` order, then velocity perturbations in the same order.
    dstates : list[Expr]
        ``[s.t_diff() for s in states]``.
    inputs : list[Expr]
        ``[u.delta() for u in inputs]``.
    M, A, B : dict[str, dict[str, Expr]]
        ``M[row][str(dstate)]``, ``A[row][str(state)]``, ``B[row][str(input_variation)]``.
        Dense: every (row, column) pair is present. Block shapes follow
        row type × column type (vector × vector → 3×3 ``Matrix``, scalar × vector →
        1×3 ``Transpose(Vector)``, vector × scalar → 3×1 ``Vector``, scalar × scalar
        → ``Scalar``); zero blocks are the library zero constants of that shape.
    row_types, col_types : dict[str, ExprType]
        Expression type of each row and of each column symbol (used by ``block_shape``).
    """

    rows: list = field(default_factory=list)
    states: list = field(default_factory=list)
    dstates: list = field(default_factory=list)
    inputs: list = field(default_factory=list)
    M: dict = field(default_factory=dict)
    A: dict = field(default_factory=dict)
    B: dict = field(default_factory=dict)
    row_types: dict = field(default_factory=dict, repr=False)
    col_types: dict = field(default_factory=dict, repr=False)

    def block_shape(self, row, col) -> tuple[int, int]:
        """Shape of the block at (*row*, *col*); *col* may be a symbol or its ``str``."""
        row_key = row if isinstance(row, str) else str(row)
        col_key = col if isinstance(col, str) else str(col)
        if row_key not in self.row_types:
            raise KeyError(f"LinearizedSystem: unknown row '{row_key}'")
        if col_key not in self.col_types:
            raise KeyError(f"LinearizedSystem: unknown column '{col_key}'")
        n_rows = 1 if self.row_types[row_key] == ExprType.SCALAR else 3
        n_cols = 1 if self.col_types[col_key] == ExprType.SCALAR else 3
        return (n_rows, n_cols)

    def __str__(self):
        lines = ["LinearizedSystem: M d/dt(dx) = A dx + B du"]
        lines.append("  states : [" + ", ".join(str(s) for s in self.states) + "]")
        lines.append("  inputs : [" + ", ".join(str(u) for u in self.inputs) + "]")
        for row in self.rows:
            lines.append(f"  row {row}")
            for name, blocks in (("M", self.M[row]), ("A", self.A[row]), ("B", self.B[row])):
                for col, block in blocks.items():
                    if not _is_zero_block(block):
                        lines.append(f"    {name}[{col}] = {block}")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def _failed_row_stage(trace, stage):
    """``trace.failed_at`` for ``linearize``: the first stage of the failing row that is
    not recorded (``dE[row]``, then ``bucket``/``raw``/``block`` of the last symbol), or
    *stage* (``validation``, ``row[label]``) when that cannot be told."""
    if not stage.startswith("row["):
        return stage
    label = stage[4:-1]
    names = [n for n in trace.stage_names() if n.split("[", 1)[-1].startswith(label + "]")]
    if f"dE[{label}]" not in names:
        return f"dE[{label}]"
    last = names[-1]
    kind, _, rest = last.partition("[")
    following = {"bucket": "raw", "raw": "block"}.get(kind)
    if following is not None:
        return f"{following}[{rest}"
    return f"bucket[{label}]"


def linearize(
    eom, variables: SystemVariables, inputs=None, params=(), *, trace=None
) -> LinearizedSystem:
    """Linearize the equations of motion along an arbitrary trajectory.

    Parameters
    ----------
    eom : dict
        Output of ``compute_eom``: ``{str(variation_vec): (variation_vec, eom_expr)}``.
    variables : SystemVariables
        The configuration variables used for ``compute_eom``.
    inputs : list[Scalar | Vector], optional
        Input symbols (forces, torques, thrust). Their variations form δu.
    params : list, optional
        The prescribed parameters passed to ``compute_eom(..., params=...)``. Prescribed
        parameters are not supported by ``linearize`` yet: a term carrying the variation
        of a declared param raises ``AlgebraicError`` saying so (without ``params`` the
        message cannot tell a param from a missing input or a missing ``Constant``).
    trace : LinearizationTrace, optional
        Keyword only. Filled per dynamic row with δE, the partition buckets, the raw
        coefficients (before ``full_simplify`` and before the mass blocks are negated)
        and the final blocks; see ``LinearizationTrace``. Cleared when the call
        starts; after an exception ``trace.failed_at`` names the stage that raised.

    Returns
    -------
    LinearizedSystem
        Blocks of ``M δẋ = A δx + B δu`` (see the module docstring for the sign
        convention and the kinematic rows).

    Raises
    ------
    TypeError
        Invalid ``variables`` or ``inputs`` (not plain ``Scalar``/``Vector`` leaves),
        a non-``SO3`` matrix configuration, or an equation whose type differs from
        that of its variation vector (a scalar equation under a vector key).
    ValueError
        Constant or duplicate inputs, an input that is also a configuration variable
        or a tangent vector (``Ω``, ``ω``) whose variation is already a state, or
        ``eom`` keys that do not match ``variables``. A stray variation symbol in an
        equation also surfaces as the engine's ``ValueError`` (δ of a variation is
        not representable).
    AlgebraicError
        A term of δE carries no perturbation symbol, more than one (a nonlinear or
        incompletely distributed term), or one only inside a non-canonical calculus
        node, or a coefficient could not be factored.
    """
    if trace is not None:
        trace._reset()
    stage = "validation"
    try:
        _check_variables(variables, caller="linearize")
        cfg_vars, cfg_perts, vel_perts = _perturbation_variables(variables)
        _check_eom(eom, variables)
        param_names = _param_names(params)
        debug = debug_enabled()
        if debug:
            from geomech.utils.invariants import check_tree

            for _, eq in eom.values():
                check_tree(eq)
        record = stage_recorder(trace, "linearize")
        states = cfg_perts + vel_perts
        dstates = [s.t_diff() for s in states]
        inputs = _check_inputs(inputs, variables, states + dstates)
        input_perts = [u.delta() for u in inputs]
        dvel_perts = dstates[len(cfg_perts) :]

        col_types = {}
        for sym in states + dstates + input_perts:
            col_types[str(sym)] = sym.type

        rows, row_types = [], {}
        M, A, B = {}, {}, {}

        # ---- kinematic rows -------------------------------------------------
        for i, var in enumerate(cfg_vars):
            label = "kin:" + str(var)
            cfg, vel = cfg_perts[i], vel_perts[i]
            row_type = cfg.type
            rows.append(label)
            row_types[label] = row_type
            M[label] = {str(d): _zero_block(row_type, d.type) for d in dstates}
            A[label] = {str(s): _zero_block(row_type, s.type) for s in states}
            B[label] = {str(u): _zero_block(row_type, u.type) for u in input_perts}

            M[label][str(dstates[i])] = One if row_type == ExprType.SCALAR else IdentityMatrix
            if isinstance(var, S2):
                omega = var.get_tangent_vector()
                qqT = VVMul(var, Transpose(var))
                A[label][str(cfg)] = MMMul(qqT, Hat(omega))
                A[label][str(vel)] = MAdd(IdentityMatrix, SMMul(qqT, -1))
            elif isinstance(var, SO3):
                A[label][str(cfg)] = SMMul(Hat(var.get_tangent_vector()), -1)
                A[label][str(vel)] = IdentityMatrix
            else:
                A[label][str(vel)] = One if row_type == ExprType.SCALAR else IdentityMatrix

        # ---- dynamic rows ---------------------------------------------------
        # Priority order for assigning a term to its perturbation symbol: the most
        # deeply wrapped symbols first, so that a term in d/dt(δv) is never tested
        # for δv (which would raise inside _has_algebraic).
        priority = dvel_perts + vel_perts + cfg_perts + input_perts
        all_perts = priority + [d for d in dstates[: len(cfg_perts)]]
        n_mass = len(dvel_perts)
        n_state = n_mass + len(vel_perts) + len(cfg_perts)

        for label, (var_vec, eq) in eom.items():
            stage = f"row[{label}]"
            row_type = var_vec.type
            rows.append(label)
            row_types[label] = row_type
            M[label] = {str(d): _zero_block(row_type, d.type) for d in dstates}
            A[label] = {str(s): _zero_block(row_type, s.type) for s in states}
            B[label] = {str(u): _zero_block(row_type, u.type) for u in input_perts}

            dE = expand(full_simplify(eq.delta()))
            if record is not None:
                record(f"dE[{label}]", dE)
            if _is_zero_block(dE):
                continue
            buckets = _partition_terms(dE, priority, label, param_names)

            # Route by position in ``priority`` (identity), never by printed name.
            for index, terms in enumerate(buckets):
                if not terms:
                    continue
                sym = priority[index]
                part = _combine(terms, row_type)
                key = str(sym)
                if record is not None:
                    record(f"bucket[{label}][{key}]", part)
                raw = extract_linear_coeff(part, sym)
                if record is not None:
                    record(f"raw[{label}][{key}]", raw)
                coeff = full_simplify(raw)
                _check_coefficient_free(coeff, all_perts, part, label)
                is_mass = index < n_mass
                if is_mass:
                    coeff = full_simplify(_negate(coeff))
                if row_type == ExprType.SCALAR and sym.type == ExprType.VECTOR:
                    coeff = Transpose(coeff)
                if record is not None:
                    record(f"block[{label}][{key}]", coeff)
                if debug:
                    check_tree(coeff)
                if is_mass:
                    M[label][key] = coeff
                elif index < n_state:
                    A[label][key] = coeff
                else:
                    B[label][key] = coeff

        return LinearizedSystem(
            rows=rows,
            states=states,
            dstates=dstates,
            inputs=input_perts,
            M=M,
            A=A,
            B=B,
            row_types=row_types,
            col_types=col_types,
        )
    except Exception as err:
        if trace is not None:
            trace._fail(_failed_row_stage(trace, stage), err)
        raise


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def _check_inputs(inputs, variables, state_symbols):
    """Inputs must be distinct, non-constant ``Scalar``/``Vector`` leaves that are neither
    configuration variables nor tangent vectors, and whose variation is not already a
    state or d/dt-of-state symbol."""
    if inputs is None:
        return []
    inputs = list(inputs)
    config = {str(v) for v in variables.scalars + variables.vectors + variables.matrices}
    states = {str(s) for s in state_symbols}
    seen = set()
    for u in inputs:
        if not isinstance(u, (Scalar, Vector)) or isinstance(u, S2):
            raise TypeError(
                f"linearize: input '{u}' must be a Scalar or Vector leaf, got {type(u).__name__}"
            )
        if isinstance(u, (TS2, TSO3)):
            raise ValueError(
                f"linearize: input '{u}' is a tangent or variation vector of a configuration"
                " variable; its variation is part of the perturbation state, not an input"
            )
        if u.is_constant or u.is_numeric or u.is_zero:
            raise ValueError(
                f"linearize: input '{u}' is constant: its variation vanishes,"
                " so it cannot be an input"
            )
        if str(u) in config:
            raise ValueError(f"linearize: input '{u}' is also a configuration variable")
        if str(u.delta()) in states:
            raise ValueError(
                f"linearize: the variation '{u.delta()}' of input '{u}' is already a"
                f" perturbation state (is '{u}' a tangent vector of a configuration"
                " variable?)"
            )
        if str(u) in seen:
            raise ValueError(f"linearize: input '{u}' is listed more than once")
        seen.add(str(u))
    return inputs


def _check_eom(eom, variables):
    if not isinstance(eom, dict):
        raise TypeError(
            f"linearize: eom must be the dict returned by compute_eom, got {type(eom)}"
        )
    expected = []
    for s in variables.scalars:
        expected.append(str(s.delta()))
    for v in variables.vectors:
        expected.append(str(v.get_variation_vector()))
    for m in variables.matrices:
        expected.append(str(m.get_variation_vector()))
    # compare the list length too: a repeated variable would collapse into one key
    if len(eom) != len(expected) or set(eom) != set(expected):
        raise ValueError(
            f"linearize: eom keys {sorted(eom)} do not match the variation vectors"
            f" {sorted(expected)} of the given SystemVariables"
        )
    for key, value in eom.items():
        if not (isinstance(value, tuple) and len(value) == 2):
            raise TypeError(
                f"linearize: eom['{key}'] must be a (variation_vector, expression) pair"
            )
        var_vec, eq = value
        if str(var_vec) != key or not isinstance(eq, Expr):
            raise TypeError(
                f"linearize: eom['{key}'] must be a (variation_vector, expression) pair"
            )
        if eq.type not in (ExprType.SCALAR, ExprType.VECTOR):
            raise TypeError(f"linearize: eom['{key}'] must be a scalar or vector expression")
        if eq.type != var_vec.type:
            raise TypeError(
                f"linearize: eom['{key}'] has type {eq.type.name}, but its variation vector"
                f" '{var_vec}' has type {var_vec.type.name}; the equation must have the"
                " type of the variation it was extracted against"
            )


def _perturbation_variables(variables):
    """Configuration variables with their configuration and velocity perturbations."""
    cfg_vars, cfg_perts, vel_perts = [], [], []
    for s in variables.scalars:
        cfg_vars.append(s)
        cfg_perts.append(s.delta())
        vel_perts.append(s.t_diff().delta())
    for v in variables.vectors:
        cfg_vars.append(v)
        if isinstance(v, S2):
            cfg_perts.append(v.get_variation_vector())
            vel_perts.append(Variation(v.get_tangent_vector()))
        else:
            cfg_perts.append(v.delta())
            vel_perts.append(v.t_diff().delta())
    for m in variables.matrices:
        if not isinstance(m, SO3):
            raise TypeError(
                f"linearize: matrix configuration variable '{m}' must be an SO3;"
                " general matrix configurations have no perturbation vector"
            )
        cfg_vars.append(m)
        cfg_perts.append(m.get_variation_vector())
        vel_perts.append(Variation(m.get_tangent_vector()))
    return cfg_vars, cfg_perts, vel_perts


# ---------------------------------------------------------------------------
# Term handling
# ---------------------------------------------------------------------------


def _flatten_addends(expr):
    if isinstance(expr, (Add, VAdd, MAdd)):
        return list(expr.nodes)
    return [expr]


def _combine(terms, row_type):
    if len(terms) == 1:
        return terms[0]
    if row_type == ExprType.SCALAR:
        return Add(*terms)
    return VAdd(*terms)


def _param_names(params):
    """Printed names of the declared prescribed parameters (``params`` of compute_eom)."""
    if isinstance(params, (Expr, str)) or not hasattr(params, "__iter__"):
        hint = f"; write params=[{params}]" if isinstance(params, Expr) else ""
        raise TypeError(
            f"linearize: params must be a sequence of leaves, got {type(params).__name__}"
            f" {str(params)!r}{hint}"
        )
    return {str(p) for p in params}


def _variation_symbols(term):
    """``{printed δs: printed s}`` for every ``Variation`` node in *term* (outermost)."""
    found = {}
    stack = [term]
    while stack:
        node = stack.pop()
        if isinstance(node, Variation):
            inner = node.expr
            while isinstance(inner, TimeDerivative):
                inner = inner.expr
            found[str(node)] = str(inner)
            continue
        stack.extend(getattr(node, "nodes", None) or ())
    return found


def _no_perturbation_hint(term, param_names):
    """The remedy part of the "no perturbation symbol" message for *term*."""
    symbols = _variation_symbols(term)
    params = sorted({s for s in symbols.values() if s in param_names})
    if params:
        names = ", ".join(f"'{p}'" for p in params)
        return (
            f" It carries the variation of the prescribed parameter(s) {names} declared in"
            " params: linearize does not support prescribed parameters yet (the variation"
            " of a prescribed parameter vanishes along the reference, but its rate does not"
            " in general)."
        )
    if symbols:
        listed = ", ".join(f"'{d}' (of '{s}')" for d, s in sorted(symbols.items()))
        return (
            f" It carries {listed}, which is neither a configuration, velocity nor input"
            " variation. If the symbol is an input, list it in `inputs`; if it is a"
            " constant, give it attr=['Constant']; if it is a prescribed parameter"
            " (compute_eom(..., params=[...])), linearize does not support prescribed"
            " parameters yet (pass params=[...] to linearize to have it say so)."
        )
    return " Is an input missing from `inputs`, or a symbol missing the 'Constant' attribute?"


def _partition_terms(dE, priority, label, param_names=frozenset()):
    """Assign every additive term of *dE* to exactly one perturbation symbol.

    Returns one list of terms per entry of *priority* (same order). A term that
    carries no perturbation symbol, one only inside a calculus operator, or more
    than one distinct symbol raises ``AlgebraicError``: the first two mean the
    linearization is incomplete, the last that the term is nonlinear or was not
    fully distributed by ``expand`` and would be extracted wrongly.
    """
    buckets = [[] for _ in priority]
    for term in _flatten_addends(dE):
        if _is_zero_block(term):
            continue
        owner = None
        for index, sym in enumerate(priority):
            try:
                found = _has_algebraic(term, sym)
            except AlgebraicError as err:
                raise AlgebraicError(
                    f"linearize: in row '{label}', the term '{term}' of the variation"
                    f" carries '{sym}' only inside a calculus operator; the"
                    f" linearization is incomplete ({err})"
                ) from None
            if found:
                owner = index
                break
        if owner is None:
            raise AlgebraicError(
                f"linearize: in row '{label}', the term '{term}' of the variation"
                " contains no perturbation symbol (configuration, velocity or input"
                " variation); the linearization is incomplete."
                + _no_perturbation_hint(term, param_names)
            )
        for other in priority[owner + 1 :]:
            if _occurs_outside(term, other, priority[owner]):
                raise AlgebraicError(
                    f"linearize: in row '{label}', the term '{term}' of the variation"
                    f" carries both '{priority[owner]}' and '{other}'; a term must be"
                    " linear in exactly one perturbation symbol"
                )
        buckets[owner].append(term)
    return buckets


def _occurs_outside(term, sym, owner):
    """Does *sym* occur in *term* anywhere outside the subtrees equal to *owner*?

    Subtrees identical to *owner* are skipped, so ``δω`` inside ``d/dt(δω)`` is not
    reported when ``d/dt(δω)`` is the owner. A calculus node other than the owner
    that contains *sym* counts as an occurrence (such a term has no linear
    coefficient either way).
    """
    if _same(term, owner):
        return False
    if _same(term, sym):
        return True
    if isinstance(term, (TimeDerivative, TimeIntegral, Variation)):
        return term.has(sym)
    return any(_occurs_outside(n, sym, owner) for n in getattr(term, "nodes", None) or [])


def _check_coefficient_free(coeff, all_perts, part, label):
    """A block must not depend on any perturbation symbol."""
    for other in all_perts:
        try:
            leftover = _has_algebraic(coeff, other)
            inner = None
        except AlgebraicError as err:
            leftover, inner = True, err
        if leftover:
            detail = f" ({inner})" if inner is not None else ""
            raise AlgebraicError(
                f"linearize: in row '{label}', the coefficient '{coeff}' extracted from"
                f" '{part}' still contains the perturbation '{other}'{detail}"
            ) from inner


def _negate(expr):
    match getattr(expr, "type", None):
        case ExprType.SCALAR:
            return Mul(expr, -1)
        case ExprType.VECTOR:
            return expr * -1
        case ExprType.MATRIX:
            return SMMul(expr, -1)
        case _:
            raise TypeError(f"linearize: cannot negate '{expr}'")


def _zero_block(row_type, col_type):
    """Zero constant of the block shape (row type × column type)."""
    match (row_type, col_type):
        case (ExprType.SCALAR, ExprType.SCALAR):
            return Zero
        case (ExprType.SCALAR, ExprType.VECTOR):
            return Transpose(ZeroVector)
        case (ExprType.VECTOR, ExprType.SCALAR):
            return ZeroVector
        case (ExprType.VECTOR, ExprType.VECTOR):
            return ZeroMatrix
        case _:
            raise TypeError(f"linearize: unsupported block shape {row_type} × {col_type}")


def _is_zero_block(expr) -> bool:
    return bool(getattr(expr, "is_zero", False))
