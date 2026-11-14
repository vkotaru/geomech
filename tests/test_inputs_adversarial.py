"""Adversarial review of the shared input validation (``check_inputs``, ``check_symbol_clash``).

Passing
tests record what held; each confirmed finding is a strict xfail whose reason starts with
its number (``F<n>``), so that fixing it turns the test into an XPASS failure.
"""

import time

import numpy as np
import pytest

from geomech import (
    Dot,
    Matrix,
    Scalar,
    SystemVariables,
    Vector,
    ZeroVector,
    compute_eom,
    linearize,
    to_standard_form,
)
from geomech.dynamics.eom import check_inputs, check_symbol_clash
from geomech.utils.errors import AlgebraicError, NonLinearError
from tests import vbl_systems
from tests.codegen_systems import ALL_SYSTEMS

HALF = Scalar("0.5", value=0.5)
CALLERS = {
    "to_standard_form": lambda eom, variables, inputs: to_standard_form(eom, variables, inputs),
    "linearize": lambda eom, variables, inputs: linearize(eom, variables, inputs),
}


def _mass(force, extra_L=None, params=()):
    """``L = ½ m ẋ·ẋ (+ extra_L)``, ``δW = δx·force``: ``(eom, variables)``."""
    m, x = Scalar("m", attr=["Constant"]), Vector("x")
    L = HALF * m * Dot(x.t_diff(), x.t_diff())
    if extra_L is not None:
        L = L + extra_L
    variables = SystemVariables(vectors=[x])
    return compute_eom(L, Dot(x.delta(), force), variables, params=params), variables


def _damped_mass():
    """``m ẍ = F − c ẋ``: the equation holds the node ``d/dt(x)``."""
    x, F, c = Vector("x"), Vector("F"), Scalar("c", attr=["Constant"])
    m = Scalar("m", attr=["Constant"])
    variables = SystemVariables(vectors=[x])
    L = HALF * m * Dot(x.t_diff(), x.t_diff())
    eom = compute_eom(L, Dot(x.delta(), F - c * x.t_diff()), variables)
    return eom, variables, x, F


def _message(caller, eom, variables, inputs, exc=ValueError):
    with pytest.raises(exc) as err:
        CALLERS[caller](eom, variables, inputs)
    return str(err.value)


# ---------------------------------------------------------------------------
# (1) False positives: legitimate systems still accepted
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("build", ALL_SYSTEMS, ids=lambda f: f.__name__)
def test_every_system_builder_passes_with_its_inputs(build):
    """compute_eom → to_standard_form → linearize with the builder's own inputs: no input
    validation error. linearize may still raise AlgebraicError for a system it does not
    support (a non-constant symbol that is not a configuration variable)."""
    s = build()
    eom = compute_eom(s.L, s.dW, s.variables)
    sf = to_standard_form(eom, s.variables, s.inputs)
    g_names = set().union(*(set(eq.G) for eq in sf.values())) if sf else set()
    assert g_names <= {str(u) for u in s.inputs}
    try:
        lin = linearize(eom, s.variables, s.inputs)
    except AlgebraicError:
        return
    assert [str(u) for u in lin.inputs] == [str(u.delta()) for u in s.inputs]


def test_unit_norm_input_used_as_the_same_symbol():
    F = Vector("F", attr=["UnitNorm"])
    eom, variables = _mass(F)
    assert set(to_standard_form(eom, variables, [F])["\\delta{x}"].G) == {"F"}
    assert [str(u) for u in linearize(eom, variables, [F]).inputs] == ["\\delta{F}"]


def test_attr_as_tuple_and_flag_order_are_the_same_symbol():
    F = Vector("F", attr=["UnitNorm"])
    eom, variables = _mass(F)
    rebuilt = Vector("F", attr=("UnitNorm",))
    assert set(to_standard_form(eom, variables, [rebuilt])["\\delta{x}"].G) == {"F"}
    a1 = Vector("a", attr=["UnitNorm", "Constant"])
    a2 = Vector("a", attr=["Constant", "UnitNorm"])
    check_symbol_clash([a2], [Dot(a1, Vector("y"))], "t")  # no raise


def test_input_in_the_lagrangian_and_the_virtual_work():
    """An input that is also a prescribed param of L (``u e₁·x``) and a force in δW."""
    u, F, e1 = Scalar("u"), Vector("F"), Vector("e1", attr=["Constant"])
    x = Vector("x")
    eom, variables = _mass(F, extra_L=-u * Dot(e1, x), params=[u])
    eq = to_standard_form(eom, variables, [F, u])["\\delta{x}"]
    assert set(eq.G) == {"F", "u"}
    assert "e1" in str(eq.G["u"])
    lin = linearize(eom, variables, [F, u])
    assert "e1" in str(lin.B["\\delta{x}"]["\\delta{u}"])


def test_param_used_only_inside_a_time_derivative():
    l, x, e1 = Scalar("l"), Vector("x"), Vector("e1", attr=["Constant"])
    eom, _ = _mass(Vector("F"), extra_L=Dot(x, e1) * l.t_diff(), params=[l])
    assert "\\frac{d}{dt}(l)" in str(eom["\\delta{x}"][1])
    # a rebuilt param (the same symbol) from a generator is accepted too
    eom2, _ = _mass(
        Vector("F"), extra_L=Dot(x, e1) * l.t_diff(), params=(p for p in [Scalar("l")])
    )
    assert str(eom2["\\delta{x}"][1]) == str(eom["\\delta{x}"][1])


def test_param_of_another_class_inside_a_time_derivative_is_found():
    """The walk enters calculus nodes: a Vector param against a Scalar under d/dt."""
    l, x, e1 = Scalar("l"), Vector("x"), Vector("e1", attr=["Constant"])
    with pytest.raises(
        ValueError, match="params entry 'l' \\(Vector\\) has the name of the Scalar"
    ):
        _mass(Vector("F"), extra_L=Dot(x, e1) * l.t_diff(), params=[Vector("l")])


def test_input_that_is_also_a_param_fails_in_the_engine():
    """An input that is a prescribed param of compute_eom, with its rate in the equations:
    since the fix log of the review, check_inputs raises the AlgebraicError (naming
    d/dt(l) and the param) in both functions, instead of the engine."""
    l, F = Scalar("l"), Vector("F")
    m, x = Scalar("m", attr=["Constant"]), Vector("x")
    variables = SystemVariables(vectors=[x])
    L = HALF * m * l * l * Dot(x.t_diff(), x.t_diff())
    eom = compute_eom(L, Dot(x.delta(), F), variables, params=[l])
    for caller in CALLERS:
        with pytest.raises(AlgebraicError, match="calculus operator"):
            CALLERS[caller](eom, variables, [F, l])


# ---------------------------------------------------------------------------
# (2) False negatives that held
# ---------------------------------------------------------------------------


def test_nonlinear_input():
    u, e1, x = Scalar("u"), Vector("e1", attr=["Constant"]), Vector("x")
    m = Scalar("m", attr=["Constant"])
    variables = SystemVariables(vectors=[x])
    L = HALF * m * Dot(x.t_diff(), x.t_diff())
    eom = compute_eom(L, Dot(x.delta(), e1) * u * u, variables)
    with pytest.raises(NonLinearError, match="'u' appears in both sides"):
        to_standard_form(eom, variables, [u])
    b = linearize(eom, variables, [u]).B["\\delta{x}"]["\\delta{u}"]
    assert "u" in str(b) and "(2)" in str(b)


@pytest.mark.parametrize("caller", sorted(CALLERS))
@pytest.mark.parametrize(
    "inputs, exc, fragment",
    [
        (lambda F: [F.T()], TypeError, "got Transpose"),
        (lambda F: [Matrix("F")], TypeError, "got Matrix"),
        (lambda F: ["F"], TypeError, "got str"),
        (lambda F: "F", TypeError, "inputs must be a sequence"),
        (lambda F: {"F": 1}, TypeError, "got str"),
        (lambda F: [Scalar("F")], ValueError, "(Scalar) has the name of the Vector 'F'"),
        (lambda F: [Vector("F", attr=["UnitNorm"])], ValueError, "a different symbol"),
    ],
    ids=["transpose", "matrix", "list-of-str", "str", "dict-of-str", "class", "flags"],
)
def test_misuse_rejected_by_both(caller, inputs, exc, fragment):
    F = Vector("F")
    eom, variables = _mass(F)
    msg = _message(caller, eom, variables, inputs(F), exc)
    assert msg.startswith(f"{caller}: ") and fragment in msg


def test_input_named_like_a_number_literal_of_the_equations():
    x, k = Vector("x"), Scalar("k", attr=["Constant"])
    eom, variables = _mass(Vector("F"), extra_L=-k * Dot(x, x))
    assert "(-2)" in str(eom["\\delta{x}"][1])
    for caller in CALLERS:
        msg = _message(caller, eom, variables, [Vector("F"), Scalar("(-2)")])
        assert "has the name of the constant '(-2)'" in msg


def test_input_named_like_the_zero_vector():
    check_symbol_clash([Vector("w")], [ZeroVector], "t")
    with pytest.raises(ValueError, match="constant '0v'"):
        check_symbol_clash([Vector("0v")], [Dot(ZeroVector, Vector("y"))], "t")


@pytest.mark.parametrize(
    "inputs",
    [
        lambda F: [F, Scalar("F")],
        lambda F: [Vector("F", size=(6,))],
        lambda F: [Scalar("F")],
    ],
    ids=["dup-other-class", "other-size", "other-class"],
)
def test_new_cases_give_the_same_message_in_both(inputs):
    F = Vector("F")
    eom, variables = _mass(F)
    msgs = {c: _message(c, eom, variables, inputs(F)).removeprefix(f"{c}: ") for c in CALLERS}
    assert msgs["to_standard_form"] == msgs["linearize"]


# ---------------------------------------------------------------------------
# (5) Performance
# ---------------------------------------------------------------------------


def _best(fn, n):
    best = float("inf")
    for _ in range(n):
        t = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - t)
    return best


@pytest.mark.parametrize(
    "build",
    [vbl_systems.quadrotor_with_load, vbl_systems.double_rigid_pendulum],
    ids=lambda f: f.__name__,
)
def test_check_inputs_is_cheap_next_to_to_standard_form(build):
    """Measured on review: 0.06-0.10 ms for check_inputs against 1.1-2.1 ms for
    to_standard_form (4-5 %); the bound leaves a wide margin for a loaded machine."""
    s = build()
    eom = compute_eom(s.L, s.dW, s.variables)
    eqs = [eq for _, eq in eom.values()]
    t_check = _best(lambda: check_inputs(s.inputs, s.variables, eqs, "t"), 20)
    t_sf = _best(lambda: to_standard_form(eom, s.variables, s.inputs), 3)
    assert t_check < 0.25 * t_sf


# ---------------------------------------------------------------------------
# Confirmed findings (strict xfail)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("size", [[3], 3], ids=["list", "int"])
def test_f1_vector_size_spelling_is_the_same_symbol(size):
    F = Vector("F")
    eom, variables = _mass(F)
    rebuilt = Vector("F", size=size)
    assert set(to_standard_form(eom, variables, [rebuilt])["\\delta{x}"].G) == {"F"}
    assert [str(u) for u in linearize(eom, variables, [rebuilt]).inputs] == ["\\delta{F}"]


def test_f1_matrix_param_size_spelling_is_the_same_symbol():
    x = Vector("x")
    J = Matrix("J", attr=["SymmetricMatrix"])
    rebuilt = Matrix("J", size=[3, 3], attr=["SymmetricMatrix"])
    eom, _ = _mass(Vector("F"), extra_L=Dot(x, J * x), params=[rebuilt])
    assert "J" in str(eom["\\delta{x}"][1])


def test_f2_size_clash_message_names_the_sizes():
    eom, variables = _mass(Vector("F"))
    msg = _message("linearize", eom, variables, [Vector("F", size=(6,))])
    assert "(6,)" in msg and "(3,)" in msg


def test_f3_numeric_leaf_hint_names_value():
    x, k = Vector("x"), Scalar("k", value=2.0)
    eom, variables = _mass(Vector("F"), extra_L=-HALF * k * Dot(x, x))
    msg = _message("to_standard_form", eom, variables, [Vector("F"), Scalar("k")])
    assert "attr=['Constant']" not in msg and "value" in msg


@pytest.mark.parametrize("caller", sorted(CALLERS))
def test_f4_input_printed_like_the_velocity_is_rejected(caller):
    eom, variables, x, F = _damped_mass()
    with pytest.raises(ValueError, match="d}{dt}"):
        CALLERS[caller](eom, variables, [F, Vector(str(x.t_diff()))])


@pytest.mark.parametrize("kind", ["set", "dict"])
@pytest.mark.parametrize("caller", sorted(CALLERS))
def test_f5_unordered_inputs_rejected(caller, kind):
    s = vbl_systems.quadrotor()
    eom = compute_eom(s.L, s.dW, s.variables)
    inputs = set(s.inputs) if kind == "set" else dict.fromkeys(s.inputs)
    with pytest.raises(TypeError):
        CALLERS[caller](eom, s.variables, inputs)


def test_value_on_a_vector_leaf_is_ignored_by_the_key():
    """Vector value= sets no flag (numeric evaluation ignores it, follow-up 2), so a
    rebuilt input with a value is the same symbol; Scalars with different values are
    both numeric and never reach the check as inputs or params."""
    F = Vector("F")
    eom, variables = _mass(F)
    rebuilt = Vector("F", value=np.zeros(3))
    assert set(to_standard_form(eom, variables, [rebuilt])["\\delta{x}"].G) == {"F"}
