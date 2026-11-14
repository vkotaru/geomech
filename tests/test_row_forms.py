"""Row vectors built with the operators: ``vᵀ M``, ``s vᵀ``, ``aᵀ + bᵀ`` (phase 3 review, F2).

``x.T() * K`` builds the row ``MVMul(xᵀ, K)``. These tests check that such rows go through
``compute_eom`` (the derived equation is compared with the Euler-Lagrange equation derived
by hand, and with the finite-difference oracle), that ``expand`` distributes them in row
form, and that ``extract_linear_coeff`` takes their column structurally.
"""

import numpy as np
import pytest

from geomech import (
    Matrix,
    Scalar,
    SVMul,
    TimeDerivative,
    Transpose,
    Vector,
    VVMul,
    assert_eom,
    compute_eom,
    expand,
    extract_linear_coeff,
)
from geomech.core.operations.geometry import column_of, is_row_vector
from geomech.numeric import evaluate
from tests.test_eom_oracle import (
    row_damping_work,
    row_quadratic_form,
    row_shifted_spring,
    row_spring,
)

RNG = np.random.default_rng(7)
_S = RNG.standard_normal((3, 3))
ENV = {
    "m": 1.7,
    "x": RNG.standard_normal((3, 1)),
    "\\frac{d}{dt}(x)": RNG.standard_normal((3, 1)),
    "\\frac{d}{dt}(\\frac{d}{dt}(x))": RNG.standard_normal((3, 1)),
    "f": RNG.standard_normal((3, 1)),
    "e3": np.array([[0.0], [0.0], [1.0]]),
    "K": _S + _S.T,
    "A": RNG.standard_normal((3, 3)),
}


def _expected(name):
    """Hand-derived ∂L/∂x − d/dt ∂L/∂ẋ + f, the form ``compute_eom`` returns."""
    x, v, a, f = (
        ENV[k] for k in ("x", "\\frac{d}{dt}(x)", "\\frac{d}{dt}(\\frac{d}{dt}(x))", "f")
    )
    m, K, A, e3 = ENV["m"], ENV["K"], ENV["A"], ENV["e3"]
    return {
        # L = ½ m ẋ·ẋ − ½ xᵀKx
        "row_spring": -m * a - K @ x + f,
        # L = ½ m ẋ·ẋ − xᵀAx: ∂/∂x (xᵀAx) = (A + Aᵀ) x
        "row_quadratic_form": -m * a - (A + A.T) @ x + f,
        # L = ½ m ẋ·ẋ − (x − e₃)ᵀK(x − e₃)
        "row_shifted_spring": -m * a - 2 * K @ (x - e3) + f,
        # δW = −δxᵀ A ẋ + δx·f
        "row_damping_work": -m * a - A @ v + f,
    }[name]


@pytest.mark.parametrize(
    "build", [row_spring, row_quadratic_form, row_shifted_spring, row_damping_work]
)
def test_row_times_matrix_derives(build):
    s = build()
    eom = compute_eom(s.L, s.dW, s.variables)
    assert list(eom) == ["\\delta{x}"]
    _, eq = eom["\\delta{x}"]
    np.testing.assert_allclose(evaluate(eq, ENV), _expected(build.__name__), atol=1e-12)
    assert assert_eom(s.L, s.dW, s.variables, eom, inputs=s.inputs) < 1e-6


a, b, c = Vector("a"), Vector("b"), Vector("c")
K, J = Matrix("K"), Matrix("J")
x = Scalar("x")
VALS = {
    "a": RNG.standard_normal((3, 1)),
    "b": RNG.standard_normal((3, 1)),
    "c": RNG.standard_normal((3, 1)),
    "K": RNG.standard_normal((3, 3)),
    "J": RNG.standard_normal((3, 3)),
    "x": 0.7,
}


@pytest.mark.parametrize(
    "build, text",
    [
        (lambda: (a.T() * (K + J)) * b, "((a)'Kb+(a)'Jb)"),
        (lambda: ((a + b).T() * K) * c, "((a)'Kc+(b)'Kc)"),
        (lambda: VVMul(SVMul(Transpose(a + b), x), c), "((a)'xc+(b)'xc)"),
        (lambda: (a.T() + b.T()) * c, "((a)'c+(b)'c)"),
        (lambda: c * ((a + b).T() * K), "(c(a)'K+c(b)'K)"),
        (lambda: (a.T() * K) * (x + Scalar("y")), "((a)'Kx+(a)'Ky)"),
    ],
)
def test_expand_distributes_row_shapes(build, text):
    expr = build()
    out = expand(expr)
    assert str(out) == text
    env = dict(VALS, y=-1.3)
    np.testing.assert_allclose(evaluate(out, env), evaluate(expr, env), atol=1e-12)


def test_column_of_each_row_shape():
    assert column_of(Transpose(a)) is a
    assert str(column_of(a.T() * K)) == "(K)'a"
    assert str(column_of(SVMul(Transpose(a), x))) == "ax"
    assert str(column_of(a.T() + b.T())) == "(a+b)"
    assert str(column_of(TimeDerivative(Transpose(a)))) == "\\frac{d}{dt}(a)"
    for row in (a.T() * K, SVMul(Transpose(a), x), a.T() + b.T(), (a.T() * K) * J):
        assert is_row_vector(row) and not is_row_vector(column_of(row))
        np.testing.assert_allclose(evaluate(column_of(row), VALS).T, evaluate(row, VALS))


@pytest.mark.parametrize(
    "build, query, coeff",
    [
        # aᵀK b: coefficient of a is K b, of b is Kᵀa
        (lambda: (a.T() * K) * b, a, lambda v: v["K"] @ v["b"]),
        (lambda: (a.T() * K) * b, b, lambda v: v["K"].T @ v["a"]),
        # x aᵀ b with the scalar inside the row
        (lambda: VVMul(SVMul(Transpose(a), x), b), a, lambda v: v["x"] * v["b"]),
    ],
)
def test_extract_through_row_shapes(build, query, coeff):
    out = extract_linear_coeff(build(), query)
    np.testing.assert_allclose(evaluate(out, VALS), coeff(VALS), atol=1e-12)
