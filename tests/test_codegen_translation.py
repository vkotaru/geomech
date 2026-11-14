"""Every row of the translation table, on hand-built trees.

The pipeline's output never holds some node shapes (a row times a matrix, an inner
``VVMul``, ``Vee``), so each node class is checked here directly: ``lambdify`` (both targets)
against ``evaluate`` at random values, with non-symmetric matrices so that a missing
transpose shows.
"""

from __future__ import annotations

import numpy as np
import pytest

from geomech import (
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
    One,
    Scalar,
    SMMul,
    SVMul,
    SystemVariables,
    Transpose,
    VAdd,
    Vector,
    Vee,
    VVMul,
    Zero,
    ZeroMatrix,
    ZeroVector,
)
from geomech.codegen import StateLayout, lambdify, to_source
from tests.codegen_systems import assert_close, reference

x = Vector("x")
R = SO3("R")
s = Scalar("s")
A = Matrix("A", attr=["Constant"])
B = Matrix("B", attr=["Constant"])
b = Vector("b", attr=["Constant"])
c = Vector("c", attr=["Constant"])
k = Scalar("k", attr=["Constant"])
F = Vector("F")
W = R.get_tangent_vector()
xd = x.t_diff()
xdd = xd.t_diff()

LAYOUT = StateLayout(SystemVariables(scalars=[s], vectors=[x], matrices=[R]), [F], [A, B, b, c, k])

CASES = {
    "add": k + s + Mul(k, s),
    "vadd": VAdd(x, b, SVMul(F, k)),
    "madd": MAdd(A, B, R),
    "mul": Mul(k, s),
    "svmul": SVMul(x, k),
    "smmul": SMMul(A, s),
    "mvmul": MVMul(A, x),
    "mvmul_row": VVMul(MVMul(Transpose(b), A), c),
    "mvmul_row_vvmul": VVMul(MVMul(Transpose(x), A), b),
    "mmmul": MMMul(A, MMMul(B, R)),
    "vvmul_inner": VVMul(Transpose(b), x),
    "vvmul_outer": MVMul(VVMul(b, Transpose(x)), c),
    "dot": Dot(x, MVMul(A, b)),
    "cross": Cross(x, b),
    "hat": MVMul(Hat(W), b),
    "vee": Vee(MAdd(A, SMMul(Transpose(A), -1))),
    "transpose_matrix": MVMul(Transpose(A), b),
    "transpose_row_scaled": VVMul(SVMul(Transpose(b), k), x),
    "row_sum": VVMul(VAdd(Transpose(b), Transpose(c)), x),
    "velocity": Dot(xd, W) + Mul(s.t_diff(), k),
    "acceleration": VAdd(xdd, SVMul(W.t_diff(), s.t_diff().t_diff())),
    "input": Dot(F, MVMul(R, b)),
    "zeros": VAdd(ZeroVector, MVMul(ZeroMatrix, x), SVMul(x, Zero)),
    "identity": MVMul(MAdd(IdentityMatrix, A), b),
    "numbers": Mul(Mul(k, Number(-2)), Number(0.25)) + One,
}


@pytest.mark.parametrize("name", sorted(CASES))
@pytest.mark.parametrize("target", ["numpy", "casadi"])
def test_node_translation(name, target):
    if target == "casadi":
        pytest.importorskip("casadi", reason="the CasADi target needs `pip install casadi`")
    expr = CASES[name]
    fn = lambdify(expr, LAYOUT, target)
    rng = np.random.default_rng(3)
    for _ in range(5):
        env = {str(A): rng.normal(size=(3, 3)), str(B): rng.normal(size=(3, 3))}
        state = rng.normal(size=LAYOUT.n_state)
        R_val = np.linalg.qr(rng.normal(size=(3, 3)))[0]
        state[LAYOUT.slot(R).slice] = R_val.reshape(9)
        a, u = rng.normal(size=LAYOUT.n_accel), rng.normal(size=LAYOUT.n_inputs)
        p = {"A": env["A"], "B": env["B"], "b": rng.normal(size=3), "c": rng.normal(size=3)}
        p["k"] = float(rng.normal())
        env = LAYOUT.env_from_arrays(state, a, u, p)
        got = np.array(fn(state, u, p, a), dtype=float)
        assert_close(got, reference(expr, env), f"{name} ({target})")


def test_row_mvmul_is_transposed_matrix_product():
    src = to_source(CASES["mvmul_row"], LAYOUT)
    assert "A.T @ b" in src


def test_literals_inlined():
    src = to_source(CASES["numbers"], LAYOUT)
    assert "-2" in src and "0.25" in src and "p[" in src
    assert "np.zeros(3)" in to_source(CASES["zeros"], LAYOUT)
    assert "np.eye(3)" in to_source(CASES["identity"], LAYOUT)


def test_common_subexpressions_computed_once():
    shared = Cross(x, MVMul(A, b))
    expr = Dot(shared, shared) + Dot(shared, c)
    src = to_source(expr, LAYOUT)
    assert src.count("_cross(") == 1
    assert src.count("A @ b") == 1


def test_long_expression_split_into_short_lines():
    expr = sum((Dot(MVMul(A, x), Cross(b, c)) * Number(i) for i in range(1, 30)), Zero)
    src = to_source(expr, LAYOUT)
    assert max(len(line) for line in src.splitlines()) <= 99
    fn = lambdify(expr, LAYOUT)
    rng = np.random.default_rng(0)
    state, u = rng.normal(size=LAYOUT.n_state), rng.normal(size=LAYOUT.n_inputs)
    p = {"A": rng.normal(size=(3, 3)), "B": np.eye(3), "b": rng.normal(size=3)}
    p.update(c=rng.normal(size=3), k=1.0)
    assert_close(
        fn(state, u, p), reference(expr, LAYOUT.env_from_arrays(state, None, u, p)), "sum"
    )


def test_returns_scalar_vector_matrix_shapes():
    rng = np.random.default_rng(0)
    state, u = rng.normal(size=LAYOUT.n_state), rng.normal(size=LAYOUT.n_inputs)
    p = {"A": np.eye(3), "B": np.eye(3), "b": np.ones(3), "c": np.ones(3), "k": 2.0}
    assert np.shape(lambdify(CASES["mul"], LAYOUT)(state, u, p)) == ()
    assert np.shape(lambdify(CASES["cross"], LAYOUT)(state, u, p)) == (3,)
    assert np.shape(lambdify(CASES["madd"], LAYOUT)(state, u, p)) == (3, 3)
