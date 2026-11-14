"""Self-tests of ``geomech.numeric``: ``evaluate`` on known trees and the ``check_eom`` oracle.

The oracle tests assert two things: a correct equation scores at finite-difference noise
(~1e-6) and a wrong one scores far above it. The second half is what makes the first half
evidence rather than a tautology.
"""

import numpy as np
import pytest

from geomech import (
    S2,
    SO3,
    Cross,
    Dot,
    Hat,
    I,
    Matrix,
    Scalar,
    SystemVariables,
    TimeDerivative,
    Transpose,
    Variation,
    Vector,
    Vee,
    Zero,
    ZeroMatrix,
    ZeroVector,
    assert_eom,
    check_eom,
    compute_eom,
    evaluate,
    getScalars,
)
from geomech.numeric import Evaluator, expm_so3, hat, rand_rot, sample_state, vee
from tests import vbl_systems

E1 = np.array([[1.0], [0.0], [0.0]])
E2 = np.array([[0.0], [1.0], [0.0]])
E3 = np.array([[0.0], [0.0], [1.0]])
HALF = Scalar("0.5", value=0.5, attr=["Constant"])


# ---------------------------------------------------------------------------
# SO(3) primitives
# ---------------------------------------------------------------------------


class TestPrimitives:
    def test_hat_is_cross_product(self):
        a, b = np.array([1.0, 2.0, 3.0]), np.array([-2.0, 0.5, 4.0])
        assert np.allclose(hat(a) @ b, np.cross(a, b))

    def test_vee_inverts_hat(self):
        v = np.array([[0.3], [-1.2], [2.0]])
        assert np.allclose(vee(hat(v)), v)

    def test_rand_rot_is_a_rotation(self):
        R = rand_rot(np.random.default_rng(1))
        assert np.allclose(R.T @ R, np.eye(3))
        assert np.isclose(np.linalg.det(R), 1.0)

    def test_expm_so3_matches_series(self):
        v = np.array([0.01, -0.02, 0.005])  # |v| ≈ 0.023: the 4th-order series is exact to 1e-9
        K = hat(v)
        series = np.eye(3) + K + K @ K / 2 + K @ K @ K / 6 + K @ K @ K @ K / 24
        assert np.allclose(expm_so3(v), series, atol=1e-9)
        E = expm_so3(np.array([1.0, 2.0, -0.5]))
        assert np.allclose(E.T @ E, np.eye(3)) and np.isclose(np.linalg.det(E), 1.0)


# ---------------------------------------------------------------------------
# evaluate on known trees
# ---------------------------------------------------------------------------


class TestEvaluate:
    def setup_method(self):
        self.x, self.y = Vector("x"), Vector("y")
        self.A = Matrix("A")
        self.env = {"x": E1, "y": E2, "A": np.diag([1.0, 2.0, 3.0])}

    def test_dot(self):
        assert evaluate(Dot(self.x, self.x), self.env) == 1.0
        assert evaluate(Dot(self.x, self.y), self.env) == 0.0

    def test_cross(self):
        assert np.allclose(evaluate(Cross(self.x, self.y), self.env), E3)

    def test_hat_and_vee(self):
        assert np.allclose(evaluate(Hat(self.x), self.env), hat(E1))
        assert np.allclose(evaluate(Vee(Hat(self.y)), self.env), E2)

    def test_transpose_and_matrix_vector(self):
        assert np.allclose(evaluate(Transpose(self.A), self.env), np.diag([1.0, 2.0, 3.0]))
        assert np.allclose(evaluate(self.A * self.y, self.env), 2 * E2)
        assert np.allclose(evaluate(self.A * self.A, self.env), np.diag([1.0, 4.0, 9.0]))

    def test_scalar_arithmetic_and_numbers(self):
        a, b = getScalars("a b")
        assert evaluate(a * b + a, {"a": 2.0, "b": 3.0}) == 8.0
        assert evaluate(a * 0.5, {"a": 2.0}) == 1.0
        assert evaluate(self.x * a, {"x": E1, "a": 3.0})[0, 0] == 3.0

    def test_time_derivative_keys(self):
        v = self.x.t_diff()
        env = {str(v): 2 * E3, str(v.t_diff()): E1, "x": E1}
        assert np.allclose(evaluate(v, env), 2 * E3)
        assert np.allclose(evaluate(v.t_diff(), env), E1)
        assert np.allclose(evaluate(Cross(self.x, v), env), 2 * Cross_numeric(E1, E3))

    def test_manifold_tangent_keys(self):
        q, R = S2("q"), SO3("R")
        omega, Om = q.get_tangent_vector(), R.get_tangent_vector()
        env = {"q": E3, str(omega): E1, "R": np.eye(3), str(Om): E2}
        assert np.allclose(evaluate(q.t_diff(), env), Cross_numeric(E1, E3))  # ω × q
        assert np.allclose(evaluate(R.t_diff(), env), hat(E2))  # R Ω̂

    def test_constants_resolved_by_type(self):
        assert evaluate(Zero, {}) == 0.0
        assert np.allclose(evaluate(ZeroVector, {}), np.zeros((3, 1)))
        assert np.allclose(evaluate(ZeroMatrix, {}), np.zeros((3, 3)))
        assert np.allclose(evaluate(I, {}), np.eye(3))
        assert np.allclose(evaluate(I * self.y + ZeroVector, self.env), E2)

    def test_scalar_value_is_used(self):
        assert evaluate(HALF * Scalar("m"), {"m": 4.0}) == 2.0

    def test_unlisted_leaves_are_deterministic_in_seed(self):
        expr = Dot(self.x, Matrix("B") * self.y)
        a, b = evaluate(expr, {}, seed=3), evaluate(expr, {}, seed=3)
        assert a == b
        assert evaluate(expr, {}, seed=4) != a

    def test_env_is_not_mutated(self):
        env = {"x": E1}
        evaluate(Dot(self.x, self.y), env)
        assert set(env) == {"x"}

    def test_flagged_leaves_respect_their_structure(self):
        e = Vector("e", attr=["Constant", "UnitNorm"])
        q, R = S2("q"), SO3("R")
        S = Matrix("S", attr=["SymmetricMatrix"])
        K = Matrix("K", attr=["SkewSymmetricMatrix"])
        Q = Matrix("Q", attr=["OrthogonalMatrix"])
        assert np.isclose(np.linalg.norm(evaluate(e, {})), 1.0)
        assert np.isclose(np.linalg.norm(evaluate(q, {})), 1.0)
        for M in (R, Q):
            val = evaluate(M, {})
            assert np.allclose(val.T @ val, np.eye(3)) and np.isclose(np.linalg.det(val), 1.0)
        assert np.allclose(evaluate(S, {}), evaluate(S, {}).T)
        assert np.allclose(evaluate(K, {}), -evaluate(K, {}).T)

    def test_variation_needs_a_value(self):
        with pytest.raises(KeyError):
            evaluate(Variation(self.x), {})
        with pytest.raises(KeyError):
            evaluate(TimeDerivative(Variation(self.x)), {})
        assert np.allclose(evaluate(Variation(self.x), {str(Variation(self.x)): E2}), E2)

    def test_evaluator_shares_drawn_leaves(self):
        ev = Evaluator(env={}, rng=np.random.default_rng(0))
        first = ev(Dot(self.x, self.y))
        assert "x" in ev.env and "y" in ev.env
        assert ev(Dot(self.x, self.y)) == first


def Cross_numeric(a, b):
    return np.cross(a.reshape(3), b.reshape(3)).reshape(3, 1)


# ---------------------------------------------------------------------------
# check_eom: correct equations score at noise level
# ---------------------------------------------------------------------------

ORACLE_SYSTEMS = [
    vbl_systems.spherical_pendulum_s2,
    vbl_systems.free_rigid_body,
    vbl_systems.double_rigid_pendulum,
    vbl_systems.quadrotor_with_load,
]


@pytest.mark.parametrize("build", ORACLE_SYSTEMS, ids=lambda f: f.__name__)
def test_check_eom_is_small_on_pinned_systems(build):
    s = build()
    eom = compute_eom(s.L, s.dW, s.variables)
    err = check_eom(s.L, s.dW, s.variables, eom, inputs=s.inputs)
    assert err < 1e-5, err


def test_check_eom_with_a_prescribed_parameter():
    """Spherical pendulum with a prescribed cable length l(t) (not a configuration variable).

    L = ½ m l² ω·ω − m g l q·e₃ gives E = −m l² ω̇ − 2 m l l̇ ω − m g l (q × e₃) + f. The
    l̇ term is only right when l is moved along its own rate, so ``params=[l]`` is needed.
    """
    m, g = getScalars("m g", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    q, l, f = S2("q"), Scalar("l"), Vector("f")
    omega, xi = q.get_tangent_vector(), q.get_variation_vector()
    L = HALF * m * l * l * Dot(omega, omega) - m * g * l * Dot(q, e3)
    variables = SystemVariables(vectors=[q])
    eom = compute_eom(L, Dot(xi, f), variables, params=[l])
    assert check_eom(L, Dot(xi, f), variables, eom, inputs=[f], params=[l]) < 1e-5
    assert check_eom(L, Dot(xi, f), variables, eom, inputs=[f]) > 1e-2  # l frozen: wrong


def test_check_eom_without_forcing():
    (m,) = getScalars("m", attr=["Constant"])
    x = Vector("x")
    v = x.t_diff()
    L = m * Dot(v, v) * HALF
    variables = SystemVariables(vectors=[x])
    dW = Scalar("0", value=0, attr=["Constant", "Zero"])
    eom = compute_eom(L, dW, variables)
    assert check_eom(L, dW, variables, eom) < 1e-5
    assert check_eom(L, None, variables, eom) < 1e-5


def test_check_eom_returns_zero_without_variables():
    assert check_eom(Scalar("a"), None, SystemVariables(), {}) == 0.0


# ---------------------------------------------------------------------------
# check_eom: wrong equations are detected
# ---------------------------------------------------------------------------


def _flip_term(eom, key, index):
    """Copy of ``eom`` with the sign of one additive term of equation ``key`` flipped."""
    vec, eq = eom[key]
    term = eq.nodes[index]
    out = dict(eom)
    out[key] = (vec, eq - term * 2)
    return out


def _normal_to_sphere(term, q, variables, inputs):
    """True if *term* has no tangent component at random points of the manifold.

    The coefficient of ξ ⊥ q is only defined on the tangent plane, so the sign of a term
    like ω × (q × ẋ) = q (ω·ẋ) carries no information and the oracle must not see it.
    """
    rng = np.random.default_rng(11)
    for _ in range(3):
        env = sample_state(variables, inputs, rng)
        val = np.asarray(Evaluator(env=env, rng=rng)(term)).reshape(3, 1)
        qv = env[str(q)]
        if np.max(np.abs(val - qv * float((qv.T @ val)[0, 0]))) > 1e-9:
            return False
    return True


@pytest.mark.parametrize("build", ORACLE_SYSTEMS, ids=lambda f: f.__name__)
def test_check_eom_detects_a_flipped_term(build):
    s = build()
    eom = compute_eom(s.L, s.dW, s.variables)
    spheres = {str(q.get_variation_vector()): q for q in s.variables.vectors if isinstance(q, S2)}
    checked = 0
    for key, (_, eq) in eom.items():
        for index, term in enumerate(eq.nodes):
            if key in spheres and _normal_to_sphere(term, spheres[key], s.variables, s.inputs):
                continue
            wrong = _flip_term(eom, key, index)
            err = check_eom(s.L, s.dW, s.variables, wrong, inputs=s.inputs)
            assert err > 1e-2, f"{key} term {index} ({term}) flipped undetected"
            checked += 1
    assert checked >= 3


def test_check_eom_detects_a_wrong_coriolis_term():
    """Rigid body with the sign of the gyroscopic term flipped by hand."""
    s = vbl_systems.free_rigid_body()
    R, J = s.symbols["R"], s.symbols["J"]
    Om = R.get_tangent_vector()
    key = str(R.get_variation_vector())
    eom = compute_eom(s.L, s.dW, s.variables)
    # Correct: −J Ω̇ − Ω × JΩ + M.  Wrong: −J Ω̇ + Ω × JΩ + M.
    right = (J * Om.t_diff()) * -1 - Cross(Om, J * Om) + s.symbols["M"]
    wrong = (J * Om.t_diff()) * -1 + Cross(Om, J * Om) + s.symbols["M"]
    assert check_eom(s.L, s.dW, s.variables, {key: (eom[key][0], right)}, inputs=s.inputs) < 1e-5
    assert check_eom(s.L, s.dW, s.variables, {key: (eom[key][0], wrong)}, inputs=s.inputs) > 1e-1


def test_assert_eom_names_the_variable():
    s = vbl_systems.quadrotor()
    eom = compute_eom(s.L, s.dW, s.variables)
    assert assert_eom(s.L, s.dW, s.variables, eom, inputs=s.inputs) < 1e-5
    key = str(s.symbols["R"].get_variation_vector())
    wrong = _flip_term(eom, key, 0)
    with pytest.raises(AssertionError, match=r"\\eta_\{R\}.*relative error"):
        assert_eom(s.L, s.dW, s.variables, wrong, inputs=s.inputs)


def test_check_eom_rejects_missing_equation():
    s = vbl_systems.spring_mass_r()
    with pytest.raises(ValueError, match="no equation"):
        check_eom(s.L, s.dW, s.variables, {}, inputs=s.inputs)
