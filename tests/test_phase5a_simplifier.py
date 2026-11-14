"""Phase 5a: simplifier completeness (fix plan 5.1–5.4, decision D1, phase-5 follow-ups).

Each rewrite is pinned by its printed result and checked numerically: the input and the
output evaluate to the same value at random points (on S² with |q| = 1 and ω ⊥ q, the
hypothesis of the tangency rules).
"""

from __future__ import annotations

import logging

import numpy as np
import pytest

from geomech import (
    S2,
    SO3,
    Cross,
    Dot,
    Hat,
    IdentityMatrix,
    MAdd,
    Matrix,
    MMMul,
    Mul,
    MVMul,
    Number,
    Scalar,
    SkewSymmMatrix,
    SMMul,
    SVMul,
    SystemVariables,
    TimeDerivative,
    Transpose,
    Vector,
    Vee,
    VVMul,
    expand,
    full_simplify,
    set_trace,
    to_latex,
)
from geomech.core.base.expressions import is_number_literal
from geomech.numeric import Evaluator, sample_state

q = S2("q")
w = q.get_tangent_vector()
a, b, v = Vector("a"), Vector("b"), Vector("v")
s, k = Scalar("s"), Scalar("k")
J, A = Matrix("J"), Matrix("A")
R = SO3("R")


def _same_value(before, after, trials=4):
    for seed in range(trials):
        rng = np.random.default_rng(seed)
        env = sample_state(SystemVariables(vectors=[q], matrices=[R]), [], rng)
        ev = Evaluator(env=env, rng=rng)
        x, y = ev(before), ev(after)
        assert np.allclose(np.asarray(x, float), np.asarray(y, float), atol=1e-12), (x, y)


# ---------------------------------------------------------------------------
# 5.1  vector_rules and pull reach every node class (engine audit F8, fuzz finding 11)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "expr, expected",
    [
        (Dot(v, Cross(q, Cross(w, q))), "Dot(v,\\omega_{q})"),
        (Dot(Cross(q, Cross(w, q)), v), "Dot(\\omega_{q},v)"),
        (Hat(Cross(a, a)), "0"),
        (MVMul(J, Cross(a, a)), "0v"),
        (TimeDerivative(Cross(a, a)), "0v"),
        (Hat(SVMul(v, s)), "Hat(v)s"),
        (Hat(Cross(SVMul(a, s), a)), "0"),
        (SMMul(J, Dot(q, q)), "J"),
        (MVMul(J, Cross(q, Cross(w, q))), "J\\omega_{q}"),
        (Vee(SMMul(J, s)), "Vee(J)s"),
        (MVMul(Transpose(SVMul(a, s)), J), "(a)'Js"),
    ],
    ids=lambda x: str(x)[:30],
)
def test_rules_fire_below_every_node(expr, expected):
    out = full_simplify(expr)
    assert str(out) == expected
    if "d/dt" not in str(expr) and "frac" not in str(expr):
        _same_value(expr, out)


# ---------------------------------------------------------------------------
# 5.2  expand distributes Vee over MAdd (fuzz finding 10)
# ---------------------------------------------------------------------------


def test_vee_distributes_over_madd():
    expr = Vee(MAdd(Hat(a), Hat(b)))
    assert str(expand(expr)) == "(Vee(Hat(a))+Vee(Hat(b)))"
    assert str(full_simplify(expr)) == "(a+b)"
    _same_value(expr, full_simplify(expr))


# ---------------------------------------------------------------------------
# 5.3  identity, skew transpose, MMMul associativity (engine audit F9)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "expr, expected",
    [
        (MVMul(Transpose(IdentityMatrix), v), "v"),
        (Transpose(IdentityMatrix), "I"),
        (Transpose(SkewSymmMatrix("S")), "S(-1)"),
        (Transpose(MAdd(SkewSymmMatrix("S"), Hat(a))), "(S+Hat(a))(-1)"),
        (MMMul(MMMul(A, Transpose(R)), R), "A"),
        (MMMul(Transpose(R), MMMul(R, A)), "A"),
        (MVMul(MMMul(A, Transpose(R)), MVMul(R, v)), "Av"),
    ],
    ids=lambda x: str(x)[:30],
)
def test_identity_skew_and_associativity(expr, expected):
    out = full_simplify(expr)
    assert str(out) == expected
    _same_value(expr, out)


def test_identity_is_flagged_symmetric():
    assert IdentityMatrix.is_symmetric


# ---------------------------------------------------------------------------
# 5.4  like terms on flattened products (engine audit F10, phase 3 follow-ups F3/F4)
# ---------------------------------------------------------------------------

x, y = Scalar("x"), Scalar("y")


@pytest.mark.parametrize(
    "expr, expected",
    [
        ((x * 2) * y + (x * 3) * y, "xy(5)"),
        (x * (y * 2) + x * (y * 3), "xy(5)"),
        (SVMul(v, Mul(s, 2)) + SVMul(v, s), "vs(3)"),
        (a * k - a * k, "0v"),
        (k * (a + b) - k * a, "bk"),
        (J * k - J * k, "0"),
        (VVMul(SVMul(Transpose(v), s), b) + Mul(s, VVMul(Transpose(v), b)), "(v)'bs(2)"),
        (Mul(Number(-0.5), Mul(x, 2)), "x(-1)"),
        (x * 0.5 * 2, "x"),
        (Mul(Number(-1), x), "x(-1)"),
    ],
    ids=lambda e: str(e)[:30],
)
def test_like_terms_and_numbers(expr, expected):
    out = full_simplify(expr)
    assert str(out) == expected
    _same_value(expr, out)


def test_integer_valued_float_prints_as_int():
    out = full_simplify(x * Number(0.5) * Number(-2))
    assert str(out) == "x(-1)"
    assert isinstance(out.right.value, int)


def test_non_integral_float_is_kept():
    assert str(full_simplify(x * Number(0.1) * Number(3))) == "x(0.30000000000000004)"


def test_so3_so3_relative_rate_terms_merge():
    """The three k Ω₁×QΩ₂ terms of the SO(3)×SO(3) coupling and k Q(Ω₂×Ω₂) cancel
    (phase 1 follow-ups 5.1 and 5.4)."""
    from tests import phase1_systems as ps

    s_, eom, _ = ps.derive(ps.so3_so3_relative_rate_coupling)
    eq = str(ps.equation(eom, s_.symbols["R1"]))
    assert "Cross(\\Omega_{R2},\\Omega_{R2})" not in eq
    assert eq.count("Cross(") == 2


# ---------------------------------------------------------------------------
# D1  named scalars with a value stay symbolic (engine audit F12)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "expr, expected",
    [
        (Scalar("m", value=1) * v, "vm"),
        (Scalar("m", value=2) + Scalar("l", value=3), "(m+l)"),
        (Scalar("m", value=2) * Scalar("l", value=3), "ml"),
        (Scalar("g", value=9.81) + Scalar("g", value=9.81), "g(2)"),
    ],
    ids=lambda e: str(e)[:20],
)
def test_named_scalar_with_value_is_a_symbol(expr, expected):
    out = full_simplify(expr)
    assert str(out) == expected
    _same_value(expr, out)


def test_named_scalar_prints_by_name_in_latex():
    assert to_latex(Vector("x") * Scalar("k", value=1e7)) == "k x"


@pytest.mark.parametrize(
    "leaf, literal",
    [
        (Number(2), True),
        (Number(-1.0), True),
        (Scalar("0.5", value=0.5), True),
        (Scalar("(-1)", value=-1), True),
        (Scalar("m", value=2), False),
        (Scalar("m"), False),
        (Scalar("(-1)"), False),
    ],
    ids=str,
)
def test_is_number_literal(leaf, literal):
    assert is_number_literal(leaf) is literal


# ---------------------------------------------------------------------------
# Follow-ups: repr of nested negation, rule log lines, distinct rule names
# ---------------------------------------------------------------------------


def test_nested_negation_repr_is_parenthesised():
    rho = Vector("rho", attr=["Constant"])
    assert repr(Hat(rho) * -1 * -1) == "-(-hat(v:rho*))"
    assert repr(SMMul(Hat(rho), Mul(Number(-1), Number(-1)))) == "hat(v:rho*)·(-(-1))"


@pytest.fixture
def rule_lines():
    lines = []

    class _H(logging.Handler):
        def emit(self, record):
            lines.append(record.getMessage())

    handler = _H()
    log = logging.getLogger("geomech.rules")
    log.addHandler(handler)
    set_trace(True)
    yield lines
    set_trace(False)
    log.removeHandler(handler)


def test_no_logged_rule_line_has_two_equal_sides(rule_lines):
    from geomech import compute_eom
    from tests import vbl_systems

    for build in (vbl_systems.rigid_pendulum_pivot, vbl_systems.quadrotor_with_load):
        sy = build()
        compute_eom(sy.L, sy.dW, sy.variables)
    assert rule_lines
    for line in rule_lines:
        name, rest = line.split(": ", 1)
        before, after = rest.split(" ==> ")
        assert before != after, line


def test_flip_and_product_rule_names(rule_lines):
    from geomech import compute_eom

    sy_q = S2("q")
    l, m = Scalar("l"), Scalar("m", attr=["Constant"])
    om, xi = sy_q.get_tangent_vector(), sy_q.get_variation_vector()
    L = Mul(Mul(m, l), Dot(om, om)) * 0.5
    compute_eom(L, Dot(xi, Vector("f")), SystemVariables(vectors=[sy_q]), params=[l])
    names = {line.split(":")[0] for line in rule_lines}
    assert "ibp.product_rule" in names
    assert "collect.flip" not in names
