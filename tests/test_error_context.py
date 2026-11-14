"""Phase 4.2: engine errors re-raised with pipeline context, and the row-form kinetic
energy that ``collect`` now handles (phase 3 review F1).

Each re-raised error keeps the engine's exception type and chains the original
(``__cause__``), so ``except AlgebraicError`` code and ``pytest.raises(..., match=<engine
text>)`` keep working.
"""

from __future__ import annotations

import warnings

import pytest

from geomech import (
    SO3,
    AlgebraicError,
    Dot,
    Hat,
    Matrix,
    MMMul,
    NonLinearError,
    Scalar,
    SystemVariables,
    Transpose,
    Vector,
    Vee,
    VVMul,
    compute_eom,
    full_simplify,
    getScalars,
    linearize,
    separate_variations,
    to_standard_form,
)
from geomech.dynamics.vbl import _check_coefficient_free
from geomech.numeric import assert_eom
from tests.phase1_systems import so3_outer_product_on_velocity
from tests.test_eom_oracle import row_kinetic_energy, row_kinetic_energy_so3
from tests.vbl_systems import HALF, free_rigid_body

# ---------------------------------------------------------------------------
# compute_eom: a d/dt(η) that survived integration by parts (audit C3)
# ---------------------------------------------------------------------------


def _vee_lagrangian():
    R = SO3("R")
    Om = R.get_tangent_vector()
    A = Matrix("A", attr=["Constant"])
    w = Vector("w", attr=["Constant"])
    return Dot(w, Vee(MMMul(Hat(Om), A))), R


def test_surviving_derivative_names_the_term_and_collect():
    L, R = _vee_lagrangian()
    with pytest.raises(AlgebraicError) as info:
        compute_eom(L, None, SystemVariables(matrices=[R]))
    msg = str(info.value)
    assert msg.startswith(
        "compute_eom: while extracting the equation for '\\eta_{R}': the time derivative"
        " '\\frac{d}{dt}(\\eta_{R})' of the variation vector '\\eta_{R}' survived"
        " integration by parts in the term 'Dot(w,Vee("
    )
    assert "collect could not expose it as the left operand of a Dot" in msg
    assert "ibp.collect[\\frac{d}{dt}(\\eta_{R})]" in msg
    # the engine message is kept, and chained
    assert "(engine: extract_linear_coeff: '\\eta_{R}' appears only inside" in msg
    assert isinstance(info.value.__cause__, AlgebraicError)
    assert "jerk" not in msg


def test_row_product_hint_when_an_inner_product_kept_the_derivative():
    """``separate_variations`` on a δS that still holds an inner product with d/dt(δx)."""
    x, f = Vector("x"), Vector("f")
    dx = x.delta()
    dS = VVMul(Transpose(dx.t_diff()), f)
    with pytest.raises(AlgebraicError) as info:
        separate_variations(dS, [dx])
    msg = str(info.value)
    assert "survived integration by parts in the term '(\\frac{d}{dt}(\\delta{x}))'f'" in msg
    assert "write it as Dot(v, A * w)" in msg


def test_other_extraction_errors_keep_the_engine_text_after_the_context():
    m, x = Scalar("m", attr=["Constant"]), Scalar("x")
    dS = x.delta() * x.delta() * m  # quadratic in δx
    with pytest.raises(NonLinearError) as info:
        separate_variations(dS, [x.delta()])
    msg = str(info.value)
    assert msg.startswith("compute_eom: while extracting the equation for '\\delta{x}' from δS")
    assert "nonlinear" in msg
    assert isinstance(info.value.__cause__, NonLinearError)


def test_collect_error_during_ibp_carries_the_stage():
    s = so3_outer_product_on_velocity()
    with pytest.raises(AlgebraicError) as info:
        compute_eom(s.L, s.dW, s.variables)
    msg = str(info.value)
    assert msg.startswith(
        "compute_eom: while integrating by parts (stage ibp, targets"
        " \\frac{d}{dt}(\\eta_{R})): collect: the target '\\frac{d}{dt}(\\eta_{R})' occurs"
        " inside the matrix factor"
    )
    assert isinstance(info.value.__cause__, AlgebraicError)


# ---------------------------------------------------------------------------
# to_standard_form: the jerk case is told apart from the IBP case
# ---------------------------------------------------------------------------


def test_jerk_in_standard_form_names_equation_and_term():
    x = Scalar("x")
    (m,) = getScalars("m", attr=["Constant"])
    v = x.t_diff()
    variables = SystemVariables(scalars=[x])
    eom = compute_eom(HALF * m * v * v, x.delta() * v.t_diff().t_diff(), variables)
    with pytest.raises(AlgebraicError) as info:
        to_standard_form(eom, variables)
    msg = str(info.value)
    assert msg.startswith(
        "to_standard_form: in equation '\\delta{x}', term '\\frac{d}{dt}(\\frac{d}{dt}"
        "(\\frac{d}{dt}(x)))': the acceleration '\\frac{d}{dt}(\\frac{d}{dt}(x))' occurs only"
        " inside '\\frac{d}{dt}(\\frac{d}{dt}(\\frac{d}{dt}(x)))', a third or higher time"
        " derivative of the configuration (a jerk): the equation is not second order"
    )
    assert "survived integration by parts" not in msg
    assert isinstance(info.value.__cause__, AlgebraicError)


# ---------------------------------------------------------------------------
# linearize
# ---------------------------------------------------------------------------


def _pendulum_with_length():
    from geomech import S2

    m, g = getScalars("m g", attr=["Constant"])
    l = Scalar("l")
    e3 = Vector("e3", attr=["Constant"])
    q, f = S2("q"), Vector("f")
    om, xi = q.get_tangent_vector(), q.get_variation_vector()
    L = HALF * m * l * l * Dot(om, om) - m * g * l * Dot(q, e3)
    return L, Dot(xi, f), q, l, f


def test_linearize_says_params_are_not_supported_when_declared():
    """Phase 2 F-1: no 'list it in inputs' / 'mark it Constant' for a declared param."""
    L, dW, q, l, f = _pendulum_with_length()
    variables = SystemVariables(vectors=[q])
    eom = compute_eom(L, dW, variables, params=[l])
    with pytest.raises(AlgebraicError) as info:
        linearize(eom, variables, [f], params=[l])
    msg = str(info.value)
    assert "contains no perturbation symbol" in msg
    assert "prescribed parameter(s) 'l' declared in params" in msg
    assert "linearize does not support prescribed parameters yet" in msg
    assert "inputs" not in msg.split("incomplete.", 1)[1]
    assert "Constant" not in msg.split("incomplete.", 1)[1]


def test_linearize_hint_without_params_lists_every_remedy():
    L, dW, q, l, f = _pendulum_with_length()
    variables = SystemVariables(vectors=[q])
    eom = compute_eom(L, dW, variables, params=[l])
    with pytest.raises(AlgebraicError) as info:
        linearize(eom, variables, [f])
    msg = str(info.value)
    assert "It carries '\\delta{l}' (of 'l')" in msg
    assert "list it in `inputs`" in msg and "attr=['Constant']" in msg
    assert "linearize does not support prescribed parameters yet" in msg


def test_linearize_params_must_be_a_sequence():
    L, dW, q, l, f = _pendulum_with_length()
    variables = SystemVariables(vectors=[q])
    eom = compute_eom(L, dW, variables, params=[l])
    with pytest.raises(TypeError, match=r"params must be a sequence of leaves.*params=\[l\]"):
        linearize(eom, variables, [f], params=l)


def test_check_coefficient_free_keeps_the_inner_message():
    x = Vector("x")
    dv = x.t_diff().delta()
    coeff = Dot(dv.t_diff(), Vector("a"))  # d/dt(δv) inside, δv a perturbation
    with pytest.raises(AlgebraicError) as info:
        _check_coefficient_free(coeff, [dv], coeff, "\\delta{x}")
    msg = str(info.value)
    assert "still contains the perturbation" in msg
    assert "appears only inside the calculus operator" in msg
    assert isinstance(info.value.__cause__, AlgebraicError)


# ---------------------------------------------------------------------------
# Kinetic energy written with rows (phase 3 review F1)
# ---------------------------------------------------------------------------


def test_row_kinetic_energy_derives():
    """½ ẋᵀAẋ − ½ xᵀKx: the coefficient of δx is −½(A + Aᵀ)ẍ − Kx + f."""
    s = row_kinetic_energy()
    eom = compute_eom(s.L, s.dW, s.variables)
    (eq,) = [e for _, e in eom.values()]
    assert str(eq) == (
        "(A\\frac{d}{dt}(\\frac{d}{dt}(x))(-0.5)+(A)'\\frac{d}{dt}(\\frac{d}{dt}(x))(-0.5)"
        "+Kx(-1)+f)"
    )
    assert assert_eom(s.L, s.dW, s.variables, eom, inputs=s.inputs) < 1e-6


def test_row_kinetic_energy_free_rigid_body_matches_the_dot_form():
    """½ ΩᵀJΩ gives Lee et al. Eq. 6.16, the same string as ½ Ω·JΩ."""
    s_row, s_dot = row_kinetic_energy_so3(), free_rigid_body()
    eom_row = compute_eom(s_row.L, s_row.dW, s_row.variables)
    eom_dot = compute_eom(s_dot.L, s_dot.dW, s_dot.variables)
    assert {k: str(e) for k, (_, e) in eom_row.items()} == {
        k: str(e) for k, (_, e) in eom_dot.items()
    }
    assert assert_eom(s_row.L, s_row.dW, s_row.variables, eom_row, inputs=s_row.inputs) < 1e-6


def test_collect_turns_an_inner_product_with_the_target_into_a_dot():
    """(ẋ_δᵀA) ẋ → Dot(Aᵀẋ_δ, ẋ) → peel → Dot(ẋ_δ, A ẋ); ẋᵀA ẋ_δ → Dot(Aᵀẋ, ẋ_δ) → flip."""
    from geomech import collect

    x = Vector("x")
    A = Matrix("A", attr=["Constant"])
    target = x.delta().t_diff()
    v = x.t_diff()
    assert str(collect(VVMul(Transpose(target) * A, v), target)) == (
        "Dot(\\frac{d}{dt}(\\delta{x}),A\\frac{d}{dt}(x))"
    )
    assert str(collect(VVMul(Transpose(v) * A, target), target)) == (
        "Dot(\\frac{d}{dt}(\\delta{x}),(A)'\\frac{d}{dt}(x))"
    )
    # an inner product without the target is left alone
    other = VVMul(Transpose(v), v)
    assert collect(other, target) is other


# ---------------------------------------------------------------------------
# full_simplify cap warning wording (phase 2 F-6)
# ---------------------------------------------------------------------------


def test_cap_warning_says_unconfirmed():
    from geomech import S2, Cross, MVMul

    q = S2("q")
    om = q.get_tangent_vector()
    expr = Dot(MVMul(Hat(q), Cross(om, q)), q)
    with pytest.warns(RuntimeWarning, match="no confirmed fixed point after 1 passes"):
        full_simplify(expr, max_iter=1)
    with pytest.warns(
        RuntimeWarning, match=r"after 0 passes \(max_iter=0: no simplify pass ran\)"
    ):
        full_simplify(Scalar("a"), max_iter=0)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        full_simplify(expr)
