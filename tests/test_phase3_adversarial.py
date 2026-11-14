"""Adversarial review of phase 3 (`specs/fix-plan.md` 3.1–3.6): operators, printing.

Confirmed findings are strict ``xfail`` tests whose reason starts with ``F<n>``; the
reproducers and diagnosis are in ``specs/reviews/phase3-review.md``. Everything else here
pins what held: operator trees equal to the explicit constructors, ``Transpose`` collapse,
numeric operands, error messages, precedence-aware LaTeX on random trees, ``repr`` on
every node, and derived equations for operator-built Lagrangians checked by ``check_eom``.
"""

from __future__ import annotations

import copy
import json
import pathlib
import pickle
import re

import numpy as np
import pytest

from geomech import (
    S2,
    Cross,
    Dot,
    ExpressionMismatchError,
    Hat,
    Matrix,
    MVMul,
    Number,
    Scalar,
    SVMul,
    SystemVariables,
    TimeDerivative,
    TimeIntegral,
    Transpose,
    Variation,
    Vector,
    VVMul,
    ZeroVector,
    check_eom,
    compute_eom,
    expand,
    full_simplify,
    to_latex,
    tree_str,
)
from geomech.core.operations.geometry import is_row_vector

ROOT = pathlib.Path(__file__).resolve().parent.parent

x, y, z = Scalar("x"), Scalar("y"), Scalar("z")
a, b, c = Vector("a"), Vector("b"), Vector("c")
A, K, R = Matrix("A"), Matrix("K"), Matrix("R")


# ---------------------------------------------------------------------------
# Findings
# ---------------------------------------------------------------------------


def test_f1_scalar_times_row_keeps_old_tree():
    # At 1e10349 ScalarExpr.__mul__ built SVMul(Transpose(a), x): one transpose, a row.
    assert str(x * a.T()) == "(a)'*x"
    assert x * a.T() == SVMul(Transpose(a), x)
    assert str(Number(2) * a.T()) == "(a)'*(2)"


def test_f1_row_forms_cancel():
    assert str(full_simplify((x * a.T()) * b - x * (a.T() * b))) == "0"


def test_f2_quadratic_form_row_times_matrix():
    m = Scalar("m", attr=["Constant"])
    Ks = Matrix("K", attr=["Constant", "SymmetricMatrix"])
    q, f = Vector("q"), Vector("f")
    v = q.t_diff()
    L = m * Dot(v, v) / 2 - q.T() * Ks * q / 2
    dW = Dot(q.delta(), f)
    variables = SystemVariables(vectors=[q])
    eom = compute_eom(L, dW, variables)
    assert check_eom(L, dW, variables, eom, inputs=[f]) < 1e-6


def test_f2_expand_row_times_matrix_sum():
    J = Matrix("J")
    out = expand((a.T() * (K + J)) * b)
    assert str(out) in ("((a)'*K*b+(a)'*J*b)", "(a)'Kb+(a)'Jb")


def test_f3_negated_transpose_in_product_latex():
    # x·(−aᵀb) must not read as x − aᵀb; A·(−K)ᵀ must not read as A − Kᵀ
    for e in (x * (-a.T() * b), A * (-K).T(), c * -a.T()):
        s = to_latex(e)
        assert not re.search(r"[^\s(]\s-", s), s


def test_f3_no_double_minus():
    for e in (-((-K).T()), VVMul(-a.T(), b) * -1):
        assert "--" not in to_latex(e), to_latex(e)


def test_f4_builtin_sum_of_vectors():
    s = sum([a, b, c])
    assert s.type == a.type


def test_f5_error_names_the_operand_written():
    big = Dot(a, a) + Dot(b, a) * x
    with pytest.raises(ExpressionMismatchError) as ei:
        big + a
    # the user wrote `big + a`; the message quotes only 'Dot(a,a)'
    assert "Dot(b,a)" in str(ei.value)


def test_f5_vector_minus_number_message():
    with pytest.raises(ExpressionMismatchError) as ei:
        a - 1
    assert "(-1)" not in str(ei.value)


def test_f6_tree_is_deterministic():
    # two equal points (phase 2: S2("q") == S2("q")) built in two runs of the same script
    q1, q2 = S2("q"), S2("q")
    assert q1 == q2
    t1 = tree_str(Dot(q1.get_tangent_vector(), q1), style="indent")
    t2 = tree_str(Dot(q2.get_tangent_vector(), q2), style="indent")
    assert t1 == t2


def test_f6_notebooks_have_no_addresses():
    for name in ("rigid_pendulum", "spherical_pendulum"):
        text = (ROOT / "notebooks" / f"{name}.ipynb").read_text()
        assert "@0x" not in text, name


def test_f7_quadrotor_notebook_outputs_regenerated():
    nb = json.loads((ROOT / "notebooks" / "quadrotor.ipynb").read_text())
    text = json.dumps([cell.get("outputs", []) for cell in nb["cells"]])
    # −Ω × (JΩ) printed without parentheses, the ambiguity 3.5 removed
    assert "\\\\Omega_{R} \\\\times J \\\\Omega_{R}" not in text


# ---------------------------------------------------------------------------
# What held: operators build the old trees
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "lhs, rhs",
    [
        (lambda: x * -1, lambda: -x),
        (lambda: x - y, lambda: x + y * -1),
        (lambda: R.T() * a, lambda: MVMul(Transpose(R), a)),
        (lambda: x / 2, lambda: x * 0.5),
        (lambda: x - 1, lambda: x + (-1)),
        (lambda: np.float64(2) * x, lambda: 2.0 * x),
        (lambda: x * np.int64(2), lambda: x * 2),
        (lambda: a * x, lambda: x * a),
        (lambda: Number(2) * a, lambda: SVMul(a, Number(2))),
        (lambda: -R, lambda: R * -1),
        (lambda: -Dot(a, b), lambda: Dot(a, b) * -1),
        (lambda: a.T() * x, lambda: Transpose(SVMul(a, x))),
    ],
)
def test_operator_tree_equals_constructor(lhs, rhs):
    l, r = lhs(), rhs()
    assert str(l) == str(r)
    assert l == r


def test_numeric_folding_after_numeric_operands():
    assert str(full_simplify(x - 1 + 1)) == "x"
    assert str(full_simplify(1 - Number(1) + x)) == "x"
    assert str(full_simplify(x * 2 - x * 2.0)) == "0"
    assert str(full_simplify(x / 2 + x / 2)) == "x"
    assert str(full_simplify(-(-x))) == "x"


def test_numpy_operands():
    with pytest.raises(TypeError, match="ndarray"):
        np.array([1, 2, 3]) * x
    with pytest.raises(TypeError, match="ndarray"):
        x * np.array([1, 2, 3])
    with pytest.raises(ZeroDivisionError):
        x / 0.0
    with pytest.raises(ZeroDivisionError):
        x / np.float64(0)
    with pytest.raises(TypeError, match="bool"):
        np.bool_(True) * x
    assert str(x / np.float64(2)) == "x*(0.5)"


# ---------------------------------------------------------------------------
# What held: Transpose collapse and row orientation
# ---------------------------------------------------------------------------


def test_transpose_collapse_and_copy():
    assert Transpose(Transpose(a)) is a
    row = a.T() * R
    assert Transpose(Transpose(row)) is row
    assert is_row_vector(row) and is_row_vector(row * x) and is_row_vector(x * a.T())
    t = Transpose(a)
    for clone in (copy.copy(t), copy.deepcopy(t), pickle.loads(pickle.dumps(t))):
        assert type(clone) is Transpose and clone == t
    assert str(copy.deepcopy(VVMul(Transpose(a), b))) == "(a)'*b"


def test_row_shapes_multiply():
    assert (a.T() * R * b).type.name == "SCALAR"
    assert ((x * a.T()) * b).type.name == "SCALAR"
    assert (b * (x * a.T())).type.name == "MATRIX"
    assert (TimeDerivative(a.T()) * b).type.name == "SCALAR"
    with pytest.raises(ExpressionMismatchError, match="row vector"):
        R * (x * a.T())


# ---------------------------------------------------------------------------
# What held: operator-built Lagrangians and virtual work derive correctly
# ---------------------------------------------------------------------------


def _ke(v):
    return Scalar("m", attr=["Constant"]) * Dot(v, v) / 2


@pytest.mark.parametrize(
    "build",
    [
        lambda q, v, f, k, A_: (_ke(v) - q.T() * (A_ * q) / 2, Dot(q.delta(), f)),
        lambda q, v, f, k, A_: (_ke(v) - (k * q).T() * q / 2, Dot(q.delta(), f)),
        lambda q, v, f, k, A_: (_ke(v) - k * q.T() * q / 2, Dot(q.delta(), f)),
        lambda q, v, f, k, A_: (_ke(v) - (A_ * q).T() * q, Dot(q.delta(), f)),
        lambda q, v, f, k, A_: (_ke(v), f.T() * q.delta()),
        lambda q, v, f, k, A_: (_ke(v), 2 * f.T() * q.delta() / 2),
        lambda q, v, f, k, A_: (_ke(v), (-k * v).T() * q.delta() + Dot(q.delta(), f)),
        lambda q, v, f, k, A_: (_ke(v), q.delta().T() * (-k * v) + Dot(q.delta(), f)),
        lambda q, v, f, k, A_: (1 - k * Dot(q, q) / 2 + _ke(v), Dot(q.delta(), f)),
        lambda q, v, f, k, A_: (-(-_ke(v)) - k * Dot(q, q), Dot(q.delta(), f)),
    ],
)
def test_operator_built_systems_pass_oracle(build):
    q, f = Vector("q"), Vector("f")
    k = Scalar("k", attr=["Constant"])
    A_ = Matrix("A", attr=["Constant"])
    L, dW = build(q, q.t_diff(), f, k, A_)
    variables = SystemVariables(vectors=[q])
    eom = compute_eom(L, dW, variables)
    assert check_eom(L, dW, variables, eom, inputs=[f]) < 1e-6


# ---------------------------------------------------------------------------
# What held: error messages
# ---------------------------------------------------------------------------


def test_error_operands_truncated():
    big = Dot(a, a)
    for i in range(40):
        big = big + Dot(Vector(f"a{i}"), a) * x
    for make in (lambda: Dot(big, a), lambda: Hat(big), lambda: big / a, lambda: big * "q"):
        with pytest.raises((ExpressionMismatchError, TypeError)) as ei:
            make()
        assert len(str(ei.value)) < 300


def test_exception_classes_unchanged():
    assert ExpressionMismatchError.__mro__[1] is Exception
    e = ExpressionMismatchError("Dot", x.type, a.type)
    assert str(e) == "Dot cannot be performed between SCALAR and VECTOR"


# ---------------------------------------------------------------------------
# What held: LaTeX and repr on random trees
# ---------------------------------------------------------------------------


def _random_trees(n, seed):
    rng = np.random.default_rng(seed)
    nums = [2, -1, 0.5, -0.5, 3]

    def pick(seq):
        return seq[rng.integers(len(seq))]

    def col(d):
        if d <= 0 or rng.random() < 0.2:
            return pick([a, b, c, a.t_diff(), Variation(b)])
        r = rng.random()
        if r < 0.2:
            return col(d - 1) + col(d - 1)
        if r < 0.3:
            return pick(nums) * col(d - 1)
        if r < 0.45:
            return sc(d - 1) * col(d - 1)
        if r < 0.6:
            return mat(d - 1) * col(d - 1)
        if r < 0.75:
            return Cross(col(d - 1), col(d - 1))
        if r < 0.85:
            return TimeIntegral(col(d - 1))
        return -col(d - 1)

    def sc(d):
        if d <= 0 or rng.random() < 0.2:
            return pick([x, y, Number(2), x.t_diff()])
        r = rng.random()
        if r < 0.25:
            return sc(d - 1) - pick(nums)
        if r < 0.5:
            return sc(d - 1) * sc(d - 1)
        if r < 0.7:
            return Dot(col(d - 1), col(d - 1))
        if r < 0.85:
            return col(d - 1).T() * col(d - 1)
        return 1 - sc(d - 1)

    def mat(d):
        if d <= 0 or rng.random() < 0.3:
            return pick([A, K, R])
        r = rng.random()
        if r < 0.3:
            return mat(d - 1) * mat(d - 1)
        if r < 0.5:
            return mat(d - 1).T()
        if r < 0.7:
            return Hat(col(d - 1))
        return col(d - 1) * col(d - 1).T()

    return [pick([sc, col, mat])(int(rng.integers(1, 5))) for _ in range(n)]


def _balanced(s):
    t = s.replace("\\left(", "(").replace("\\right)", ")")
    stack = []
    for ch in t:
        if ch in "({":
            stack.append(ch)
        elif ch in ")}":
            if not stack or stack.pop() != {")": "(", "}": "{"}[ch]:
                return False
    return not stack


def test_latex_and_repr_on_random_trees():
    for e in _random_trees(400, seed=3):
        s = to_latex(e)
        assert _balanced(s), (str(e), s)
        assert "^T^T" not in s and not re.search(r"\^T\s*\^", s), (str(e), s)
        repr(e)
        tree_str(e)
        tree_str(e, style="indent")


def test_latex_round_trip_cases():
    # structurally different trees that read the same before 3.5
    pairs = [
        (Variation(x * y), Variation(x) * y),
        (TimeIntegral(Cross(a, b)), Cross(TimeIntegral(a), b)),
        (TimeDerivative(Variation(a)), Variation(TimeDerivative(a))),
        (Transpose(SVMul(Transpose(a), x)), SVMul(Transpose(a), x)),
    ]
    for p, q in pairs:
        assert to_latex(p) != to_latex(q), (str(p), str(q))
    # a transposed scaled row is the column again: (s aᵀ)ᵀ = s a
    assert to_latex(Transpose(SVMul(Transpose(a), x))) == "x a"


def test_zero_vector_start_for_sum():
    assert str(sum([a, b, c], ZeroVector)) == "(0v+a+b+c)"
