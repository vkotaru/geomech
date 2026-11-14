"""Adversarial review tests for variation-based linearization (``geomech.dynamics.vbl``).

Every expected value below was derived by hand (the derivation is in the docstring) before
the code was run. The numeric checks here deliberately do **not** go through
``tests.numeric_helpers.finite_difference_check``: ``_ev`` is an independent evaluator with
an explicit environment (no random draws), the kinematic rows are compared against *exact*
perturbed trajectory families rather than against the formula they were written from, and
the elastic spherical pendulum is compared block by block against its analytic
linearization. ``finite_difference_check`` is used only for the "held" attack cases, where
the question is whether a hand-built equation linearizes consistently at all.

Convention (``vbl`` module docstring): ``E`` is the coefficient of the variation after IBP
(the negative of the Euler-Lagrange expression) and ``δE = −M δv̇ + A δx + B δu``.

Section 8 covers outer and inner products with a sum as a factor, tangent vectors offered
as inputs, duplicate inputs, and equation/variation type mismatches.
"""

from __future__ import annotations

import numpy as np
import pytest

from geomech import (
    S2,
    SO3,
    Add,
    Cross,
    Dot,
    Hat,
    MAdd,
    Matrix,
    MMMul,
    Mul,
    MVMul,
    Scalar,
    SMMul,
    SVMul,
    SystemVariables,
    TimeDerivative,
    Transpose,
    VAdd,
    Variation,
    Vector,
    Vee,
    VVMul,
    compute_eom,
    getScalars,
    linearize,
)
from geomech.core.base.types import ExprType
from geomech.utils.errors import AlgebraicError
from tests import vbl_systems
from tests.numeric_helpers import finite_difference_check

HALF = Scalar("0.5", value=0.5, attr=["Constant"])

# ---------------------------------------------------------------------------
# Independent numerics (no shared code with tests/numeric_helpers.py)
# ---------------------------------------------------------------------------


def _hat(v):
    v = np.asarray(v, dtype=float).reshape(3)
    return np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])


def _vee(M):
    return np.array([[M[2, 1]], [M[0, 2]], [M[1, 0]]])


def _cross(a, b):
    return np.cross(np.asarray(a).reshape(3), np.asarray(b).reshape(3)).reshape(3, 1)


def _expm(v):
    """Rodrigues: exp(Hat(v))."""
    v = np.asarray(v, dtype=float).reshape(3)
    th = np.linalg.norm(v)
    K = _hat(v)
    if th < 1e-14:
        return np.eye(3) + K
    return np.eye(3) + np.sin(th) / th * K + (1.0 - np.cos(th)) / th**2 * (K @ K)


def _unit(rng):
    v = rng.normal(size=(3, 1))
    return v / np.linalg.norm(v)


def _perp(q, rng):
    v = rng.normal(size=(3, 1))
    return v - q * (q.T @ v).item()


def _ev(node, env):
    """Evaluate an expression tree; every non-numeric leaf must be present in *env*.

    Calculus nodes (``TimeDerivative``, ``Variation``) are looked up by their printed name,
    like leaves. Missing keys raise ``KeyError`` instead of drawing a random value.
    """
    if isinstance(node, (Add, VAdd, MAdd)):
        return sum(_ev(n, env) for n in node.nodes)
    if isinstance(node, (Mul, SVMul, SMMul)):
        return _ev(node.left, env) * _ev(node.right, env)
    if isinstance(node, (MVMul, MMMul, VVMul)):
        out = _ev(node.left, env) @ _ev(node.right, env)
        return float(out[0, 0]) if out.shape == (1, 1) else out
    if isinstance(node, Dot):
        return float((_ev(node.left, env).T @ _ev(node.right, env))[0, 0])
    if isinstance(node, Cross):
        return _cross(_ev(node.left, env), _ev(node.right, env))
    if isinstance(node, Hat):
        return _hat(_ev(node.expr, env))
    if isinstance(node, Vee):
        return _vee(_ev(node.expr, env))
    if isinstance(node, Transpose):
        return np.atleast_2d(_ev(node.expr, env)).T
    if isinstance(node, (TimeDerivative, Variation)):
        return env[str(node)]
    if node.is_zero:
        return {ExprType.MATRIX: np.zeros((3, 3)), ExprType.VECTOR: np.zeros((3, 1))}.get(
            node.type, 0.0
        )
    if getattr(getattr(node, "flags", None), "is_identity", False):
        return np.eye(3)
    if isinstance(node, Scalar) and node.value is not None:
        return float(node.value)
    return env[str(node)]


def _m2(x):
    return np.atleast_2d(np.asarray(x, dtype=float))


# ---------------------------------------------------------------------------
# 1. Kinematic rows against exact perturbed trajectory families
# ---------------------------------------------------------------------------


class TestKinematicRowsIndependent:
    """The kinematic rows are not covered by ``finite_difference_check`` (it only differences
    the dynamic rows, and ``_dstate_values`` re-uses the very formulas under test), so they
    are checked here against trajectory families whose ε- and t-derivatives are exact.

    S²: q(ε, t) = exp(t Hat(a + ε b)) exp(ε Hat(c)) q₀, so q̇ = (a + ε b) × q exactly and
    ω = q × q̇ (⊥ q, ω × q = q̇). With δq = ∂_ε q, the tangent perturbation is ξ = q × δq
    (ξ ⊥ q and ξ × q = δq). Then d/dt(δq) = ∂_ε q̇ = b × q + a × δq exactly, and
    ξ̇ = q̇ × δq + q × d/dt(δq). The row must give ξ̇ = A[ξ] ξ + A[δω] δω.

    SO(3): R(ε, t) = exp(t Hat(a + ε b)) R₀ exp(ε Hat(c)), so Ṙ = Hat(a + ε b) R and
    Ω = Rᵀ(a + ε b) exactly. With δR = ∂_ε R: η = Vee(Rᵀ δR), δΩ = δRᵀ a + Rᵀ b,
    d/dt(δR) = Hat(b) R + Hat(a) δR, η̇ = Vee(Ṝ̇ᵀ δR + Rᵀ d/dt(δR)). The row must give
    η̇ = A[η] η + A[δΩ] δΩ.
    """

    @pytest.mark.parametrize("seed", range(4))
    def test_s2_row_reproduces_xi_dot_of_a_true_perturbed_curve(self, seed):
        rng = np.random.default_rng(seed)
        s = vbl_systems.spherical_pendulum_s2()
        eom = compute_eom(s.L, s.dW, s.variables)
        lin = linearize(eom, s.variables, s.inputs)
        q_sym = s.symbols["q"]
        om_sym, xi_sym = q_sym.get_tangent_vector(), q_sym.get_variation_vector()

        a, b, c = (rng.normal(size=(3, 1)) for _ in range(3))
        q0, t0 = _unit(rng), float(rng.uniform(0.2, 1.0))

        def q_of(eps):
            return _expm(t0 * (a + eps * b)) @ _expm(eps * c) @ q0

        def qdot_of(eps):
            return _cross(a + eps * b, q_of(eps))

        def omega_of(eps):
            return _cross(q_of(eps), qdot_of(eps))

        eps = 1e-5
        q, qd, w = q_of(0.0), qdot_of(0.0), omega_of(0.0)
        dq = (q_of(eps) - q_of(-eps)) / (2 * eps)
        dw = (omega_of(eps) - omega_of(-eps)) / (2 * eps)
        xi = _cross(q, dq)
        ddq = _cross(b, q) + _cross(a, dq)
        xidot = _cross(qd, dq) + _cross(q, ddq)

        # sanity of the family: ξ ⊥ q, ξ × q = δq, q·δω = q·(ξ × ω) (constraint q·ω = 0)
        assert abs((q.T @ xi).item()) < 1e-9
        assert np.allclose(_cross(xi, q), dq, atol=1e-9)
        assert abs((q.T @ dw).item() - (q.T @ _cross(xi, w)).item()) < 1e-7

        env = {"q": q, str(om_sym): w}
        A_xi = _ev(lin.A["kin:q"][str(xi_sym)], env)
        A_dw = _ev(lin.A["kin:q"][str(Variation(om_sym))], env)
        M_xd = _ev(lin.M["kin:q"][str(xi_sym.t_diff())], env)
        assert np.allclose(M_xd, np.eye(3))
        assert np.allclose(A_xi @ xi + A_dw @ dw, xidot, atol=1e-7)
        # and the closed form I derived independently of the spec
        P = q @ q.T
        assert np.allclose(P @ _hat(w) @ xi + (np.eye(3) - P) @ dw, xidot, atol=1e-7)

    @pytest.mark.parametrize("seed", range(4))
    def test_so3_row_reproduces_eta_dot_of_a_true_perturbed_curve(self, seed):
        rng = np.random.default_rng(10 + seed)
        s = vbl_systems.free_rigid_body()
        eom = compute_eom(s.L, s.dW, s.variables)
        lin = linearize(eom, s.variables, s.inputs)
        R_sym = s.symbols["R"]
        Om_sym, eta_sym = R_sym.get_tangent_vector(), R_sym.get_variation_vector()

        a, b, c = (rng.normal(size=(3, 1)) for _ in range(3))
        R0, t0 = _expm(rng.normal(size=(3, 1))), float(rng.uniform(0.2, 1.0))

        def R_of(eps):
            return _expm(t0 * (a + eps * b)) @ R0 @ _expm(eps * c)

        eps = 1e-5
        R = R_of(0.0)
        Om = R.T @ a
        dR = (R_of(eps) - R_of(-eps)) / (2 * eps)
        eta = _vee(R.T @ dR)
        dOm = dR.T @ a + R.T @ b
        ddR = _hat(b) @ R + _hat(a) @ dR
        Rdot = _hat(a) @ R
        etadot = _vee(Rdot.T @ dR + R.T @ ddR)

        assert np.allclose(R.T @ dR, -(R.T @ dR).T, atol=1e-7)  # δR = R Hat(η) is skew

        env = {"R": R, str(Om_sym): Om}
        A_eta = _ev(lin.A["kin:R"][str(eta_sym)], env)
        A_dOm = _ev(lin.A["kin:R"][str(Variation(Om_sym))], env)
        assert np.allclose(_ev(lin.M["kin:R"][str(eta_sym.t_diff())], env), np.eye(3))
        assert np.allclose(A_eta @ eta + A_dOm @ dOm, etadot, atol=1e-7)
        assert np.allclose(-_hat(Om) @ eta + dOm, etadot, atol=1e-7)


# ---------------------------------------------------------------------------
# 2. Analytic linearization of a mixed R × S² system, block by block
# ---------------------------------------------------------------------------


class TestElasticSphericalPendulumAnalytic:
    """Pinned equations (``tests/test_scalar_adversarial.py``):

        E_l = −m l̈ + m l ω·ω − m g q·e₃ − k (l − l₀) + f_l
        E_ξ = −m l² ω̇ − 2 m l l̇ ω − m g l q × e₃ + f_q

    Variation with δq = ξ × q, (ξ × q)·e₃ = ξ·(q × e₃), (ξ × q) × e₃ = Hat(e₃) Hat(q) ξ:

        δE_l = −m δl̈ + (m ω·ω − k) δl + 2 m l ωᵀ δω − m g (q × e₃)ᵀ ξ + δf_l
        δE_ξ = −m l² δω̇ + (−2 m l ω̇ − 2 m l̇ ω − m g q × e₃) δl − 2 m l ω δl̇
               − 2 m l l̇ δω − m g l Hat(e₃) Hat(q) ξ + δf_q

    so with M = −coeff(δv̇): M[l][δl̈] = m, M[ξ][δω̇] = m l² I, every other block as listed in
    ``_expected`` below (including the zero blocks, which pins the shapes 1×1, 1×3, 3×1, 3×3).
    """

    def _expected(self, env, sym):
        m, g, k, l = env["m"], env["g"], env["k"], env["l"]
        ld, ldd = env[str(sym["l"].t_diff())], env[str(sym["l"].t_diff().t_diff())]
        del ldd
        q, e3 = env["q"], env["e3"]
        w, wd = env[str(sym["om"])], env[str(sym["om"].t_diff())]
        I3, Z33, Z31, Z13 = np.eye(3), np.zeros((3, 3)), np.zeros((3, 1)), np.zeros((1, 3))
        qe3 = _cross(q, e3)
        l_row, x_row = str(sym["dl"]), str(sym["xi"])
        dl, xi, dld, dw = (str(x) for x in (sym["dl"], sym["xi"], sym["dld"], sym["dw"]))
        dldd, dwd, dxid = (str(x) for x in (sym["dldd"], sym["dwd"], sym["xi"].t_diff()))
        dfl, dfq = str(sym["f_l"].delta()), str(sym["f_q"].delta())
        return {
            ("M", l_row, dld): 0.0,
            ("M", l_row, dxid): Z13,
            ("M", l_row, dldd): m,
            ("M", l_row, dwd): Z13,
            ("A", l_row, dl): m * (w.T @ w).item() - k,
            ("A", l_row, xi): (-m * g * qe3).T,
            ("A", l_row, dld): 0.0,
            ("A", l_row, dw): (2 * m * l * w).T,
            ("B", l_row, dfl): 1.0,
            ("B", l_row, dfq): Z13,
            ("M", x_row, dld): Z31,
            ("M", x_row, dxid): Z33,
            ("M", x_row, dldd): Z31,
            ("M", x_row, dwd): m * l * l * I3,
            ("A", x_row, dl): -2 * m * l * wd - 2 * m * ld * w - m * g * qe3,
            ("A", x_row, xi): -m * g * l * _hat(e3) @ _hat(q),
            ("A", x_row, dld): -2 * m * l * w,
            ("A", x_row, dw): -2 * m * l * ld * I3,
            ("B", x_row, dfl): Z31,
            ("B", x_row, dfq): I3,
        }

    @pytest.mark.parametrize("seed", range(3))
    def test_every_block_matches_the_hand_linearization(self, seed):
        rng = np.random.default_rng(100 + seed)
        s = vbl_systems.elastic_spherical_pendulum()
        eom = compute_eom(s.L, s.dW, s.variables)
        lin = linearize(eom, s.variables, s.inputs)
        l, q = s.symbols["l"], s.symbols["q"]
        om, xi = q.get_tangent_vector(), q.get_variation_vector()
        sym = dict(
            l=l,
            om=om,
            xi=xi,
            dl=l.delta(),
            dld=l.delta().t_diff(),
            dldd=l.delta().t_diff().t_diff(),
            dw=Variation(om),
            dwd=Variation(om).t_diff(),
            f_l=s.symbols["f_l"],
            f_q=s.symbols["f_q"],
        )
        qv = _unit(rng)
        env = {
            "m": float(rng.uniform(0.5, 2)),
            "g": 9.81,
            "k": float(rng.uniform(0.5, 2)),
            "l0": float(rng.uniform(0.5, 2)),
            "l": float(rng.uniform(0.5, 2)),
            str(l.t_diff()): float(rng.normal()),
            str(l.t_diff().t_diff()): float(rng.normal()),
            "q": qv,
            str(om): _perp(qv, rng),
            str(om.t_diff()): _perp(qv, rng),
            "e3": rng.normal(size=(3, 1)),
            "f_l": float(rng.normal()),
            "f_q": rng.normal(size=(3, 1)),
        }
        expected = self._expected(env, sym)
        tables = {"M": lin.M, "A": lin.A, "B": lin.B}
        seen = set()
        for row in (str(sym["dl"]), str(sym["xi"])):
            for name, table in tables.items():
                for col, block in table[row].items():
                    key = (name, row, col)
                    assert key in expected, f"unexpected column {key}"
                    seen.add(key)
                    got, want = _m2(_ev(block, env)), _m2(expected[key])
                    assert got.shape == want.shape, (key, block, got.shape, want.shape)
                    assert np.allclose(got, want, atol=1e-12), (key, block, got, want)
        assert seen == set(expected)


# ---------------------------------------------------------------------------
# 3. S² mass matrix depending on q (q̇-form spherical pendulum)
# ---------------------------------------------------------------------------


class TestSphericalPendulumQdotForm:
    """L = ½ m l² q̇·q̇ − m g l q·e₃ with the eager q̇ = ω × q gives (observed, then checked)

        E = −m l² q × (ω̇ × q) − m g l q × e₃ + f = m l² Hat(q)Hat(q) ω̇ − m g l q × e₃ + f.

    The mass block depends on q, so M[δω̇] = −m l² Hat(q)Hat(q) = m l² (I − q qᵀ) and the
    variation of the mass term lands in A[ξ]:

        δ_ξ E = −m l² [ (ξ × q) × (ω̇ × q) + q × (ω̇ × (ξ × q)) ] − m g l (ξ × q) × e₃.
    """

    def setup_method(self):
        m, g, l = getScalars("m g l", attr=["Constant"])
        e3 = Vector("e3", attr=["Constant"])
        q, f = S2("q"), Vector("f")
        self.q, self.f = q, f
        self.om, self.xi = q.get_tangent_vector(), q.get_variation_vector()
        qd = q.t_diff()
        L = HALF * m * l * l * Dot(qd, qd) - m * g * l * Dot(q, e3)
        self.variables = SystemVariables(vectors=[q])
        self.eom = compute_eom(L, Dot(self.xi, f), self.variables)
        self.lin = linearize(self.eom, self.variables, [f])

    def test_mass_block_is_projected_inertia(self):
        rng = np.random.default_rng(5)
        qv = _unit(rng)
        env = {"q": qv, "m": 1.7, "l": 0.6}
        M = _ev(self.lin.M[str(self.xi)][str(Variation(self.om).t_diff())], env)
        assert np.allclose(M, 1.7 * 0.36 * (np.eye(3) - qv @ qv.T), atol=1e-12)
        assert self.lin.M[str(self.xi)][str(self.xi.t_diff())].is_zero

    @pytest.mark.parametrize("seed", range(3))
    def test_configuration_block_contains_the_mass_variation(self, seed):
        rng = np.random.default_rng(seed)
        qv = _unit(rng)
        wd, xi = _perp(qv, rng), _perp(qv, rng)
        e3 = rng.normal(size=(3, 1))
        m, g, l = 1.3, 9.81, 0.7
        env = {"q": qv, str(self.om.t_diff()): wd, "e3": e3, "m": m, "g": g, "l": l}
        A_xi = _ev(self.lin.A[str(self.xi)][str(self.xi)], env)
        dq = _cross(xi, qv)
        expected = -m * l * l * (_cross(dq, _cross(wd, qv)) + _cross(qv, _cross(wd, dq)))
        expected = expected - m * g * l * _cross(dq, e3)
        assert np.allclose(A_xi @ xi, expected, atol=1e-12)

    def test_finite_difference_holds(self):
        finite_difference_check(self.lin, self.eom, self.variables, [self.f])


# ---------------------------------------------------------------------------
# 4. δω is not a tangent vector: no S² tangent rule may fire on it
# ---------------------------------------------------------------------------


class TestVariationOfOmegaIsNotTangent:
    """E = c q × (ω × q) equals c ω on the constraint set, but δE must be computed on the
    unconstrained δω: δE = c[(ξ × q) × (ω × q) + q × (δω × q) + q × (ω × (ξ × q))].
    If the rule q × (δω × q) → δω fired, the row would be off by c q qᵀ δω = c ξ × ω.
    On a tangent direction (ξ ⊥ q, δω = ξ × ω + w, w ⊥ q) the correct row sum is c δω.
    """

    @pytest.mark.parametrize("seed", range(3))
    def test_row_sum_equals_c_delta_omega_on_tangent_directions(self, seed):
        rng = np.random.default_rng(20 + seed)
        c = Scalar("c", attr=["Constant"])
        q = S2("q")
        om, xi = q.get_tangent_vector(), q.get_variation_vector()
        E = Cross(q, Cross(om, q)) * c
        var = SystemVariables(vectors=[q])
        lin = linearize({str(xi): (xi, E)}, var)
        qv = _unit(rng)
        w, xiv, wv = _perp(qv, rng), _perp(qv, rng), _perp(qv, rng)
        dw = _cross(xiv, w) + wv
        env = {"q": qv, str(om): w, "c": 2.5}
        row = str(xi)
        total = _ev(lin.A[row][str(xi)], env) @ xiv + _ev(lin.A[row][str(Variation(om))], env) @ dw
        assert np.allclose(total, 2.5 * dw, atol=1e-12)
        # the δω block must be the projector, not the identity
        assert np.allclose(
            _ev(lin.A[row][str(Variation(om))], env), 2.5 * (np.eye(3) - qv @ qv.T), atol=1e-12
        )


# ---------------------------------------------------------------------------
# 5. Smaller exact pins: damping, cubic input, scalar/vector mixed shapes
# ---------------------------------------------------------------------------


class TestExactPins:
    def test_scalar_damping_goes_to_velocity_block(self):
        """E = −m s̈ − k s + f − c ṡ: A[δṡ] = −c, M = m, A[δs] = −k."""
        m, k, c = getScalars("m k c", attr=["Constant"])
        x, f = Scalar("x"), Scalar("f")
        v = x.t_diff()
        L = HALF * m * v * v - HALF * k * x * x
        var = SystemVariables(scalars=[x])
        eom = compute_eom(L, x.delta() * (f - c * v), var)
        lin = linearize(eom, var, [f])
        row = str(x.delta())
        assert str(lin.A[row][str(x.delta().t_diff())]) == "c(-1)"
        assert str(lin.A[row][str(x.delta())]) == "k(-1)"
        assert str(lin.M[row][str(x.delta().t_diff().t_diff())]) == "m"
        assert str(lin.B[row][str(f.delta())]) == "1"

    def test_vector_damping_goes_to_velocity_block(self):
        """E = −m ẍ + F − c ẋ: A[δẋ] = −c I."""
        m, c = getScalars("m c", attr=["Constant"])
        x, F = Vector("x"), Vector("F")
        v = x.t_diff()
        var = SystemVariables(vectors=[x])
        eom = compute_eom(HALF * m * Dot(v, v), Dot(x.delta(), F - v * c), var)
        lin = linearize(eom, var, [F])
        row = str(x.delta())
        assert str(lin.A[row][str(x.delta().t_diff())]) == "Ic(-1)"
        assert lin.A[row][str(x.delta())].is_zero
        assert str(lin.M[row][str(x.delta().t_diff().t_diff())]) == "Im"

    def test_cubic_input_block(self):
        """E = −m ẍ + F (F·F): B[δF] = (F·F) I + 2 F Fᵀ."""
        (m,) = getScalars("m", attr=["Constant"])
        x, F = Vector("x"), Vector("F")
        v = x.t_diff()
        var = SystemVariables(vectors=[x])
        eom = compute_eom(HALF * m * Dot(v, v), Dot(x.delta(), F * Dot(F, F)), var)
        lin = linearize(eom, var, [F])
        Fv = np.array([[0.3], [-1.2], [0.8]])
        B = _ev(lin.B[str(x.delta())][str(F.delta())], {"F": Fv, "m": 1.0})
        assert np.allclose(B, (Fv.T @ Fv).item() * np.eye(3) + 2 * Fv @ Fv.T, atol=1e-12)

    def test_cubic_velocity_term(self):
        """E = Ω (Ω·Ω) + Ω × JΩ: A[δΩ] = (Ω·Ω) I + 2 Ω Ωᵀ − Hat(JΩ) + Hat(Ω) J."""
        J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
        R = SO3("R")
        Om, eta = R.get_tangent_vector(), R.get_variation_vector()
        E = Om * Dot(Om, Om) + Cross(Om, J * Om)
        lin = linearize({str(eta): (eta, E)}, SystemVariables(matrices=[R]))
        rng = np.random.default_rng(1)
        Jv = rng.normal(size=(3, 3))
        Jv = Jv + Jv.T
        W = rng.normal(size=(3, 1))
        A = _ev(lin.A[str(eta)][str(Variation(Om))], {"J": Jv, str(Om): W})
        expected = (W.T @ W).item() * np.eye(3) + 2 * W @ W.T - _hat(Jv @ W) + _hat(W) @ Jv
        assert np.allclose(A, expected, atol=1e-12)
        assert lin.A[str(eta)][str(eta)].is_zero

    def test_quadrotor_with_load_thrust_and_coupling_blocks(self):
        """From E_{x_L} = … + m_Q l ω̇ × q − m_Q l (ω·ω) q + f R e₃ and
        E_ξ = … − l f q × R e₃ (CDC 2013 load system):
        M[x_L][δω̇] = −coeff = −m_Q l Hat(q)ᵀ… = m_Q l Hat(q) (since ω̇ × q = −Hat(q) ω̇),
        B[x_L][δf] = R e₃, A[x_L][η] = −f R Hat(e₃), B[ξ][δf] = −l q × R e₃,
        A[ξ][η] = −l f Hat(q) δ(R e₃) coefficient = l f Hat(q) R Hat(e₃).
        """
        s = vbl_systems.quadrotor_with_load()
        eom = compute_eom(s.L, s.dW, s.variables)
        lin = linearize(eom, s.variables, s.inputs)
        xL, q, R = s.symbols["xL"], s.symbols["q"], s.symbols["R"]
        om, xi, eta = q.get_tangent_vector(), q.get_variation_vector(), R.get_variation_vector()
        rng = np.random.default_rng(9)
        qv = _unit(rng)
        Rv = _expm(rng.normal(size=(3, 1)))
        e3 = np.array([[0.0], [0.0], [1.0]])
        env = {"q": qv, "R": Rv, "e3": e3, "f": 2.0, "l": 0.5, "m_Q": 1.5}
        rowx, rowq = str(xL.delta()), str(xi)
        assert np.allclose(
            _ev(lin.M[rowx][str(Variation(om).t_diff())], env), 1.5 * 0.5 * _hat(qv), atol=1e-12
        )
        assert np.allclose(_ev(lin.B[rowx][str(s.symbols["f"].delta())], env), Rv @ e3)
        assert np.allclose(_ev(lin.A[rowx][str(eta)], env), -2.0 * Rv @ _hat(e3), atol=1e-12)
        assert np.allclose(
            _ev(lin.B[rowq][str(s.symbols["f"].delta())], env), -0.5 * _cross(qv, Rv @ e3)
        )
        assert np.allclose(
            _ev(lin.A[rowq][str(eta)], env), 0.5 * 2.0 * _hat(qv) @ Rv @ _hat(e3), atol=1e-12
        )


# ---------------------------------------------------------------------------
# 6. Sums inside products: the attack surface that held
# ---------------------------------------------------------------------------

_q = S2("q")
_om, _xi = _q.get_tangent_vector(), _q.get_variation_vector()
_e3 = Vector("e3", attr=["Constant"])
_c = Scalar("c", attr=["Constant"])
_s = Scalar("s")

HELD_VECTOR_EQUATIONS = {
    "svmul_q_add": _q * (Dot(_om, _om) + Dot(_q, _e3)),
    "svmul_vadd_scalar": (_q + _om) * Dot(_om, _e3),
    "hat_vadd": Hat(_q + _om) * _e3,
    "madd_times_vadd": (Hat(_q) + Hat(_om)) * (_q + _om),
    "cross_vadd_vadd": Cross(_q + _om, _om + _e3),
    "dot_vadd_vadd_q": _q * Dot(_q + _om, _om + _e3),
    "nested_scalars": (_q * Dot(_om, _om)) * Dot(_q, _e3),
    "smmul_madd_dot": SMMul(Hat(_q) + Hat(_om), Dot(_q, _e3)) * _om,
    "mmmul_hat_madd": MMMul(Hat(_q), Hat(_q) + Hat(_om)) * _om,
    "hat_hat_scaled": SMMul(MMMul(Hat(_om), Hat(_om)), _c) * _q,
    "transpose_hat_vadd": MVMul(Transpose(Hat(_q + _om)), _e3),
    "dot_q_omega_times_q": _q * (Dot(_q, _om) * _c) + _om * _c,
    "outer_sum_times_hat": MVMul(MMMul(VVMul(_q + _om, Transpose(_e3)), Hat(_e3)), _om),
}

HELD_SCALAR_EQUATIONS = {
    "dot_vadd_times_add": Dot(_q + _om, _om) * (Dot(_q, _e3) + Dot(_om, _om)),
    "scalar_config_times_dot": _s * Dot(_q + _om, _e3) + _s.t_diff() * _s,
}


@pytest.mark.parametrize("name", sorted(HELD_VECTOR_EQUATIONS))
def test_sums_inside_products_vector_row_hold(name):
    """``expand`` distributes these; each term carries one perturbation symbol."""
    E = HELD_VECTOR_EQUATIONS[name]
    var = SystemVariables(vectors=[_q])
    eom = {str(_xi): (_xi, E)}
    lin = linearize(eom, var)
    assert finite_difference_check(lin, eom, var, []) < 1e-6


@pytest.mark.parametrize("name", sorted(HELD_SCALAR_EQUATIONS))
def test_sums_inside_products_scalar_row_hold(name):
    E = HELD_SCALAR_EQUATIONS[name]
    var = SystemVariables(scalars=[_s], vectors=[_q])
    eom = {str(_s.delta()): (_s.delta(), E), str(_xi): (_xi, _q * E)}
    lin = linearize(eom, var)
    assert finite_difference_check(lin, eom, var, []) < 1e-6


# ---------------------------------------------------------------------------
# 7. Error paths that held
# ---------------------------------------------------------------------------


class TestErrorPathsHeld:
    def test_stray_variation_vector_in_equation_raises(self):
        R = SO3("R")
        Om, eta = R.get_tangent_vector(), R.get_variation_vector()
        with pytest.raises(ValueError, match="variation vector"):
            linearize({str(eta): (eta, Cross(eta, Om) + Om)}, SystemVariables(matrices=[R]))

    def test_stray_variation_node_in_equation_raises(self):
        R = SO3("R")
        Om, eta = R.get_tangent_vector(), R.get_variation_vector()
        with pytest.raises(ValueError, match="variation of a variation"):
            linearize({str(eta): (eta, Cross(Variation(Om), Om))}, SystemVariables(matrices=[R]))

    def test_variation_vector_as_input_raises(self):
        s = vbl_systems.spherical_pendulum_s2()
        eom = compute_eom(s.L, s.dW, s.variables)
        with pytest.raises(ValueError):
            linearize(eom, s.variables, [s.symbols["f"], s.symbols["q"].get_variation_vector()])

    def test_time_derivative_of_input_raises_algebraic_error(self):
        (m,) = getScalars("m", attr=["Constant"])
        x, F = Vector("x"), Vector("F")
        v = x.t_diff()
        var = SystemVariables(vectors=[x])
        eom = compute_eom(HALF * m * Dot(v, v), Dot(x.delta(), F.t_diff()), var)
        with pytest.raises(AlgebraicError, match="calculus operator"):
            linearize(eom, var, [F])

    @pytest.mark.parametrize(
        "E",
        [
            MVMul(Transpose(VVMul(_q, Transpose(_om))), _om),
            Vee(Hat(_q) * Hat(_om) + Hat(_om) * Hat(_q)),
        ],
        ids=["transpose_outer", "vee_madd"],
    )
    def test_unsupported_structures_raise_not_drop(self, E):
        """Transpose of a sum/outer product and Vee are opaque to extraction; must raise."""
        with pytest.raises(AlgebraicError):
            linearize({str(_xi): (_xi, E)}, SystemVariables(vectors=[_q]))

    @pytest.mark.parametrize(
        "E",
        [
            MVMul(Transpose(Hat(_q) + Hat(_om)), _om),
            Vee(Hat(_q + _om)),
        ],
        ids=["transpose_madd", "vee_hat_vadd"],
    )
    def test_formerly_unsupported_structures_linearize_correctly(self, E):
        """Phase 5a: (Hat(a) + Hat(b))ᵀ simplifies to −(Hat(a) + Hat(b)) (a sum of skew
        matrices is skew) and Vee distributes over the sum, so these two used to raise and
        now linearize. Not dropped: the blocks match finite differences."""
        eom, variables = {str(_xi): (_xi, E)}, SystemVariables(vectors=[_q])
        lin = linearize(eom, variables)
        assert finite_difference_check(lin, eom, variables, []) < 1e-6

    def test_input_named_like_a_configuration_variable_raises(self):
        s = vbl_systems.free_rigid_body()
        eom = compute_eom(s.L, s.dW, s.variables)
        with pytest.raises(ValueError, match="configuration variable"):
            linearize(eom, s.variables, [Vector("R")])

    def test_exports(self):
        import geomech
        import geomech.dynamics.vbl as vbl

        assert geomech.linearize is vbl.linearize
        assert geomech.LinearizedSystem is vbl.LinearizedSystem
        assert "linearize" in geomech.dynamics.__all__


# ---------------------------------------------------------------------------
# 8. Outer products of sums, tangent-vector inputs, duplicate inputs, type mismatches
# ---------------------------------------------------------------------------


class TestOuterProductOfSumsLinearizes:
    """Outer and inner products with a sum as a factor. After δ, a factor (q + ω) becomes
    (ξ × q + δω), so the term (ξ × q + δω) e₃ᵀ e₃ carries two perturbation symbols until
    ``expand`` distributes the product. ``linearize`` must either split it and return the
    correct block, or refuse with ``AlgebraicError``; it must never return a zero block for
    the dropped half (the completeness rule: exactly one perturbation symbol per term).
    """

    def test_vector_row_left_factor(self):
        """E = (q + ω) e₃ᵀ e₃ = (e₃·e₃)(q + ω): A[ξ] = −(e₃·e₃) Hat(q), A[δω] = (e₃·e₃) I."""
        E = MVMul(VVMul(_q + _om, Transpose(_e3)), _e3)
        var = SystemVariables(vectors=[_q])
        try:
            lin = linearize({str(_xi): (_xi, E)}, var)
        except AlgebraicError:
            return
        qv, e3 = _unit(np.random.default_rng(0)), np.array([[0.2], [-0.5], [1.0]])
        env = {"q": qv, "e3": e3}
        assert np.allclose(_ev(lin.A[str(_xi)][str(_xi)], env), -(e3.T @ e3).item() * _hat(qv))

    def test_vector_row_right_factor(self):
        """E = e₃ (q + ω)ᵀ e₃: δ_ξ E = e₃ (ξ × q)·e₃ = e₃ (q × e₃)ᵀ ξ, so A[ξ] = e₃ (q × e₃)ᵀ."""
        E = MVMul(VVMul(_e3, Transpose(_q + _om)), _e3)
        var = SystemVariables(vectors=[_q])
        try:
            lin = linearize({str(_xi): (_xi, E)}, var)
        except AlgebraicError:
            return
        qv, e3 = _unit(np.random.default_rng(1)), np.array([[0.2], [-0.5], [1.0]])
        env = {"q": qv, "e3": e3}
        assert np.allclose(_ev(lin.A[str(_xi)][str(_xi)], env), e3 @ _cross(qv, e3).T)

    def test_scalar_row_inner_product(self):
        """E_s = (q + ω)ᵀ e₃: δE_s = (ξ × q)·e₃ + δω·e₃ = (q × e₃)ᵀ ξ + e₃ᵀ δω."""
        E = VVMul(Transpose(_q + _om), _e3)
        var = SystemVariables(scalars=[_s], vectors=[_q])
        eom = {str(_s.delta()): (_s.delta(), E), str(_xi): (_xi, _q * _c)}
        try:
            lin = linearize(eom, var)
        except AlgebraicError:
            return
        qv, e3 = _unit(np.random.default_rng(2)), np.array([[0.2], [-0.5], [1.0]])
        env = {"q": qv, "e3": e3}
        got = _m2(_ev(lin.A[str(_s.delta())][str(_xi)], env))
        assert got.shape == (1, 3)
        assert np.allclose(got, _cross(qv, e3).T)


class TestTangentVectorAsInputRejected:
    """A tangent vector Ω (or ω) is not an input: its variation δΩ is already a
    perturbation state, so accepting it would make δΩ both a state and an input column and
    move the A[δΩ] block into B. ``linearize`` must raise ``ValueError``.
    """

    def test_so3_tangent_vector_as_input_raises(self):
        s = vbl_systems.free_rigid_body()
        eom = compute_eom(s.L, s.dW, s.variables)
        Om = s.symbols["R"].get_tangent_vector()
        with pytest.raises(ValueError):
            linearize(eom, s.variables, [s.symbols["M"], Om])

    def test_s2_tangent_vector_as_input_raises(self):
        s = vbl_systems.spherical_pendulum_s2()
        eom = compute_eom(s.L, s.dW, s.variables)
        om = s.symbols["q"].get_tangent_vector()
        with pytest.raises(ValueError):
            linearize(eom, s.variables, [s.symbols["f"], om])


class TestDuplicateInputsRejected:
    """``inputs=[M, M]`` would give ``lin.inputs`` two entries while every ``B[row]`` has
    one column, so ``B`` would not be dense over ``inputs``. ``linearize`` must raise.
    """

    def test_duplicate_input_raises(self):
        s = vbl_systems.free_rigid_body()
        eom = compute_eom(s.L, s.dW, s.variables)
        with pytest.raises(ValueError):
            linearize(eom, s.variables, [s.symbols["M"], s.symbols["M"]])


class TestEquationTypeMismatchRejected:
    """An equation must have the type of the variation it was extracted against. A scalar
    equation under a vector key would yield a 3×1 block in a slot ``block_shape`` reports as
    3×3, so ``linearize`` must raise.
    """

    def test_scalar_equation_for_vector_variation_raises(self):
        R = SO3("R")
        Om, eta = R.get_tangent_vector(), R.get_variation_vector()
        with pytest.raises((TypeError, ValueError)):
            linearize({str(eta): (eta, Dot(Om, Om))}, SystemVariables(matrices=[R]))
