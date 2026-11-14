"""Adversarial review of phase 1 (`specs/fix-plan.md`): S² rule, collect, Hat(Vee), extract,
oracle.

The confirmed findings (``F<n>``, details and the fix log in
``specs/reviews/phase1-review.md``) were first written as strict ``xfail`` tests and are
now all passing: F1 (``collect`` loop) and F3 (default step) were fixed in the engine and
the oracle. F2 (a gravity term next to a stiff ``k = 1e7`` term) is a resolution limit of
the oracle's normalisation, not a bug: the flipped term moves the relative error by about
5e-7 at any step size, far below the 1e-5 gate, so no single call can both accept the
correct equation and reject the flipped one. The resolution is that ``assert_eom`` names
the unresolved terms instead of passing silently; the F2 tests assert that behaviour.

Independent checks used here, none of which go through the library's variational rules:

- the S² variation formula is derived numerically from a family of curves on the sphere
  (``TestS2VariationFormula``);
- the equations of three S² systems are integrated with RK4 and compared with the R³
  constrained Lagrangian with a multiplier whose force is written by hand in the docstring
  (``TestTrajectoryAgreement``);
- ``check_eom`` is run on Lagrangians whose ``∂L/∂ω`` has a normal component and explicit
  q-dependence, on products with R³/SO(3), on two S² variables sharing a parameter and on
  virtual work with ξ under a matrix (``TestS2RuleOnAdversarialLagrangians``).
"""

from __future__ import annotations

import numpy as np
import pytest

from geomech import (
    SO3,
    Cross,
    Dot,
    Hat,
    I,
    Matrix,
    MMMul,
    MVMul,
    Scalar,
    SkewSymmMatrix,
    SMMul,
    SVMul,
    SystemVariables,
    Transpose,
    VAdd,
    Vector,
    Vee,
    VVMul,
    compute_eom,
    extract_linear_coeff,
    full_simplify,
    linearize,
    to_standard_form,
)
from geomech.core.math.collect import collect
from geomech.core.operations import MAdd
from geomech.numeric import (
    Evaluator,
    assert_eom,
    check_eom,
    evaluate,
    finite_difference_check,
    sample_state,
)
from geomech.utils.errors import AlgebraicError, NonLinearError
from tests import phase1_systems as ps
from tests import test_eom_oracle
from tests.phase1_systems import (
    ADVERSARIAL_SYSTEMS,
    HALF,
    TOL,
    adv_q_dependent_gyro,
    adv_r3_s2_pendulum_with_gyro,
    s2_symbols,
)

# h=1e-4 keeps the finite-difference noise two orders below the gate; it was pinned here
# while the default was 1e-5 (F3) and is kept explicit because the engine, not the
# oracle's default step, is what is under test.
H = 1e-4


def _cross(a, b):
    return np.cross(a.reshape(3), b.reshape(3)).reshape(3, 1)


# ---------------------------------------------------------------------------
# 1. The S² variation formula, derived without the library
# ---------------------------------------------------------------------------


class TestS2VariationFormula:
    """δω = ξ̇ − 2 ω × ξ on a family of curves q(t, ε) = (p(t) + ε r(t)) / |p(t) + ε r(t)|.

    ω = q × q̇ and ξ = q × ∂q/∂ε are computed from the curve analytically (q̇ and ∂q/∂ε of
    the normalised polynomial are closed form), so only one level of central differences is
    used: ∂ω/∂ε in ε and dξ/dt in t. Lee et al. eq. (5.13) in this notation reads
    δω = −S(ω)ξ + (I − qqᵀ)ξ̇; the test checks the engine's closed form against both.
    """

    @staticmethod
    def _curve(rng):
        P = rng.normal(size=(3, 4))  # p(t) = Σ P[:, k] t^k
        Rr = rng.normal(size=(3, 4))

        def poly(C, t):
            return sum(C[:, [k]] * t**k for k in range(4))

        def dpoly(C, t):
            return sum(k * C[:, [k]] * t ** (k - 1) for k in range(1, 4))

        def q_qdot(t, eps):
            p, pd = poly(P, t) + eps * poly(Rr, t), dpoly(P, t) + eps * dpoly(Rr, t)
            n = np.linalg.norm(p)
            q = p / n
            qd = (pd - q * float((q.T @ pd)[0, 0])) / n
            return q, qd

        def q_eps(t):
            """∂q/∂ε at ε = 0 (closed form)."""
            p, r = poly(P, t), poly(Rr, t)
            n = np.linalg.norm(p)
            q = p / n
            return (r - q * float((q.T @ r)[0, 0])) / n

        return q_qdot, q_eps

    def test_delta_omega_has_coefficient_minus_two(self):
        rng = np.random.default_rng(3)
        h = 1e-4
        for _ in range(10):
            q_qdot, q_eps = self._curve(rng)
            t = float(rng.uniform(-1, 1))
            q, qd = q_qdot(t, 0.0)
            omega = _cross(q, qd)

            def xi(tt):
                qq, _ = q_qdot(tt, 0.0)
                return _cross(qq, q_eps(tt))

            def omega_eps(eps):
                qq, qqd = q_qdot(t, eps)
                return _cross(qq, qqd)

            xi0 = xi(t)
            xi_dot = (xi(t + h) - xi(t - h)) / (2 * h)
            d_omega = (omega_eps(h) - omega_eps(-h)) / (2 * h)
            # constraints of the setting
            assert abs(float((q.T @ omega)[0, 0])) < 1e-12
            assert abs(float((q.T @ xi0)[0, 0])) < 1e-12
            engine_rule = xi_dot - 2 * _cross(omega, xi0)
            book_513 = -_cross(omega, xi0) + (np.eye(3) - q @ q.T) @ xi_dot
            assert np.max(np.abs(d_omega - engine_rule)) < 1e-6
            assert np.max(np.abs(d_omega - book_513)) < 1e-6
            # the pre-fix rule is only the tangential projection: it is wrong along q
            old_rule = xi_dot - _cross(omega, xi0)
            P = np.eye(3) - q @ q.T
            assert np.max(np.abs(P @ (d_omega - old_rule))) < 1e-6
            assert abs(float((q.T @ (d_omega - old_rule))[0, 0])) > 1e-3


# ---------------------------------------------------------------------------
# 2. Adversarial Lagrangians through compute_eom + check_eom
# ---------------------------------------------------------------------------


class TestS2RuleOnAdversarialLagrangians:
    @pytest.mark.parametrize("build", ADVERSARIAL_SYSTEMS, ids=lambda f: f.__name__)
    def test_matches_euler_lagrange_residual(self, build):
        s = build()
        eom = compute_eom(s.L, s.dW, s.variables)
        for seed in (0, 1):
            err = check_eom(s.L, s.dW, s.variables, eom, inputs=s.inputs, seed=seed, h=H)
            assert err < TOL, f"{build.__name__} seed {seed}: {err:.2e}"

    def test_prescribed_parameter_multiplying_a_gyroscopic_term(self):
        """L = ½ ω·ω + p(t) c·ω − g q·e₃ needs the ṗ term:

        E = −ω̇ − ṗ c + 2 p ω × c − g q × e₃ + f.
        """
        m, g, l, k, e3, a, c, d, A, S, J, q, f = s2_symbols()
        om, xi = q.get_tangent_vector(), q.get_variation_vector()
        p = Scalar("p")
        L = HALF * Dot(om, om) + p * Dot(c, om) - g * Dot(q, e3)
        variables = SystemVariables(vectors=[q])
        eom = compute_eom(L, Dot(xi, f), variables, params=[p])
        assert "c\\frac{d}{dt}(p)(-1)" in str(eom[str(xi)][1])
        assert check_eom(L, Dot(xi, f), variables, eom, inputs=[f], params=[p], h=H) < TOL


# ---------------------------------------------------------------------------
# 3. Trajectory agreement with a hand-written R³ constrained Lagrangian
# ---------------------------------------------------------------------------


def _rk4(fun, y, dt):
    k1 = fun(y)
    k2 = fun(y + 0.5 * dt * k1)
    k3 = fun(y + 0.5 * dt * k2)
    k4 = fun(y + dt * k3)
    return y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def _engine_rhs(eq, consts, q):
    """ω̇ from the engine's E = 0 solved on the tangent plane (E is affine in ω̇)."""
    kq, kw = str(q), str(q.get_tangent_vector())
    kwd = str(q.get_tangent_vector().t_diff())

    def rhs(y):
        qv, w = y[:3].reshape(3, 1), y[3:].reshape(3, 1)
        qv = qv / np.linalg.norm(qv)
        P = np.eye(3) - qv @ qv.T
        w = P @ w

        def E(wd):
            env = {**consts, kq: qv, kw: w, kwd: wd}
            return np.asarray(Evaluator(env=env, rng=np.random.default_rng(0))(eq)).reshape(3, 1)

        e0 = E(np.zeros((3, 1)))
        A = np.zeros((3, 3))
        for i in range(3):
            ei = np.zeros((3, 1))
            ei[i] = 1.0
            A[:, [i]] = E(ei) - e0
        U = np.linalg.svd(P)[0][:, :2]  # basis of the tangent plane
        wd = U @ np.linalg.solve(U.T @ A @ U, -(U.T @ e0))
        return np.concatenate([_cross(w, qv), wd]).reshape(-1)

    return rhs


def _r3_rhs(force, m, l):
    """m l² q̈ = F(q, q̇) + λ q with λ from q·q̈ = −|q̇|² (|q| = 1)."""

    def rhs(y):
        qv, qd = y[:3].reshape(3, 1), y[3:].reshape(3, 1)
        Fq = force(qv, qd)
        lam = -m * l * l * (qd.T @ qd).item() - (qv.T @ Fq).item()
        qdd = (Fq + lam * qv) / (m * l * l)
        return np.concatenate([qd, qdd]).reshape(-1)

    return rhs


def _integrate_both(build, force, T=0.5, dt=2e-3):
    s = build()
    eom = compute_eom(s.L, s.dW, s.variables)
    q = s.symbols["q"]
    eq = eom[str(q.get_variation_vector())][1]
    rng = np.random.default_rng(42)
    ev = Evaluator(env={}, rng=rng)
    ev(s.L)
    ev(eq)
    consts = dict(ev.env)
    for key in (str(q), str(q.get_tangent_vector()), str(q.get_tangent_vector().t_diff())):
        consts.pop(key, None)
    consts["m"], consts["l"], consts["g"] = 1.3, 0.8, 9.81
    consts["e3"] = np.array([[0.0], [0.0], [1.0]])
    q0 = rng.normal(size=(3, 1))
    q0 /= np.linalg.norm(q0)
    w0 = rng.normal(size=(3, 1))
    w0 -= q0 * (q0.T @ w0).item()
    y_e = np.concatenate([q0, w0]).reshape(-1)
    y_r = np.concatenate([q0, _cross(w0, q0)]).reshape(-1)
    f_e = _engine_rhs(eq, consts, q)
    f_r = _r3_rhs(lambda qq, qd: force(qq, qd, consts), consts["m"], consts["l"])
    worst = 0.0
    for _ in range(int(T / dt)):
        y_e, y_r = _rk4(f_e, y_e, dt), _rk4(f_r, y_r, dt)
        worst = max(worst, float(np.max(np.abs(y_e[:3] - y_r[:3]))))
    return worst


def _force_gyroscopic(qv, qd, c):
    """s2_gyroscopic: L̃ = ½ m l² |q̇|² + c·(q × q̇) − m g l q·e₃.

    ∂L̃/∂q̇ = m l² q̇ + c × q, ∂L̃/∂q = q̇ × c − m g l e₃, so the Euler-Lagrange equation is
    m l² q̈ + c × q̇ − q̇ × c + m g l e₃ = λ q + f × q, i.e.
    m l² q̈ = −2 c × q̇ − m g l e₃ + f × q + λ q (δq = ξ × q turns ξ·f into the R³ force f × q).
    """
    return -2 * _cross(c["c"], qd) - c["m"] * c["g"] * c["l"] * c["e3"] + _cross(c["f"], qv)


def _force_gyroscopic_half(qv, qd, c):
    """What the pre-fix engine integrates: the gyroscopic force at half strength."""
    return -1 * _cross(c["c"], qd) - c["m"] * c["g"] * c["l"] * c["e3"] + _cross(c["f"], qv)


def _force_rotating_frame(qv, qd, c):
    """s2_rotating_frame: L̃ = ½ m l² (|q̇|² + 2 W₀·(q × q̇) + |W₀|² − (W₀·q)²) − m g l q·e₃.

    ∂L̃/∂q̇ = m l² (q̇ + W₀ × q), ∂L̃/∂q = m l² (q̇ × W₀ − (W₀·q) W₀) − m g l e₃, hence
    m l² q̈ = −2 m l² W₀ × q̇ − m l² (W₀·q) W₀ − m g l e₃ + f × q + λ q.
    """
    m, l, g, W0 = c["m"], c["l"], c["g"], c["W0"]
    return (
        -2 * m * l * l * _cross(W0, qd)
        - m * l * l * (W0.T @ qv).item() * W0
        - m * g * l * c["e3"]
        + _cross(c["f"], qv)
    )


def _force_q_dependent_gyro(qv, qd, c):
    """adv_q_dependent_gyro: L̃ = ½ m l² |q̇|² + (q·a) q̇·(c × q) − m g l q·e₃.

    ∂L̃/∂q̇ = m l² q̇ + (q·a)(c × q); d/dt of it is m l² q̈ + (q̇·a)(c × q) + (q·a)(c × q̇);
    ∂L̃/∂q = a (c·(q × q̇)) + (q·a)(q̇ × c) − m g l e₃. Hence
    m l² q̈ = −(q̇·a)(c × q) − 2 (q·a)(c × q̇) + a (c·(q × q̇)) − m g l e₃ + f × q + λ q.
    """
    m, l, g, a, cc = c["m"], c["l"], c["g"], c["a"], c["c"]
    qa, qda = (qv.T @ a).item(), (qd.T @ a).item()
    return (
        -qda * _cross(cc, qv)
        - 2 * qa * _cross(cc, qd)
        + a * (cc.T @ _cross(qv, qd)).item()
        - m * g * l * c["e3"]
        + _cross(c["f"], qv)
    )


class TestTrajectoryAgreement:
    """RK4 integration of the engine's equation vs the hand-written R³ multiplier form."""

    def test_gyroscopic_pendulum(self):
        assert _integrate_both(ps.s2_gyroscopic, _force_gyroscopic) < 1e-8

    def test_gyroscopic_pendulum_rejects_half_strength_rule(self):
        assert _integrate_both(ps.s2_gyroscopic, _force_gyroscopic_half) > 1e-2

    def test_rotating_frame_pendulum(self):
        assert _integrate_both(ps.s2_rotating_frame, _force_rotating_frame) < 1e-8

    def test_q_dependent_gyroscopic_coupling(self):
        assert _integrate_both(adv_q_dependent_gyro, _force_q_dependent_gyro) < 1e-8


# ---------------------------------------------------------------------------
# 4. collect
# ---------------------------------------------------------------------------


def _collect_preserves_value(expr, vec, seeds=3):
    out = collect(expr, vec)
    rng = np.random.default_rng(0)
    for _ in range(seeds):
        ev = Evaluator(env={}, rng=rng)
        x, y = ev(expr), ev(out)
        assert abs(x - y) <= 1e-9 * max(1.0, abs(x)), f"{expr} -> {out}: {x} vs {y}"
    return out


class TestCollectPeel:
    def setup_method(self):
        self.v, self.a, self.b, self.c = Vector("v"), Vector("a"), Vector("b"), Vector("c")
        self.M, self.N, self.R, self.s = Matrix("M"), Matrix("N"), SO3("R"), Scalar("s")

    def _ok(self, expr, exposed=True):
        out = _collect_preserves_value(expr, self.v)
        if exposed:
            assert isinstance(out, Dot) and out.left == self.v, str(out)
        return out

    def test_two_matrices(self):
        v, a, b, M, N = self.v, self.a, self.b, self.M, self.N
        self._ok(Dot(MVMul(M, MVMul(N, v)), b))
        self._ok(Dot(MVMul(M, MVMul(N, Cross(v, a))), b))
        self._ok(Dot(b, MVMul(M, MVMul(N, Cross(a, v)))))

    def test_under_transpose_R_and_smmul(self):
        v, a, b, R, M, s = self.v, self.a, self.b, self.R, self.M, self.s
        self._ok(Dot(MVMul(Transpose(R), Cross(v, a)), b))
        self._ok(Dot(MVMul(SMMul(M, s), Cross(v, a)), b))

    def test_both_sides_matrices_target_on_right(self):
        v, a, b, M, N, R = self.v, self.a, self.b, self.M, self.N, self.R
        self._ok(Dot(MVMul(M, a), MVMul(N, v)))
        self._ok(Dot(MVMul(M, a), MVMul(N, Cross(b, v))))
        out = self._ok(Dot(MVMul(R, a), MVMul(R, Cross(b, v))))
        assert "R" not in str(out)  # common rotation stripped

    def test_cross_in_mvmul_in_cross(self):
        v, a, b, c, M = self.v, self.a, self.b, self.c, self.M
        self._ok(Dot(a, Cross(b, MVMul(M, Cross(v, c)))))
        self._ok(Dot(Cross(b, MVMul(M, Cross(v, c))), a))
        self._ok(Dot(MVMul(M, Cross(v, a)), Cross(b, c)))
        self._ok(Dot(MVMul(M, Cross(a, MVMul(self.N, v))), b))

    def test_unexposable_shapes_terminate_unchanged_in_value(self):
        v, a, b, M, s = self.v, self.a, self.b, self.M, self.s
        self._ok(Dot(MVMul(M, SVMul(v, s)), b), exposed=False)
        self._ok(Dot(MVMul(M, VAdd(v, a)), b), exposed=False)
        self._ok(Dot(a, Cross(SVMul(v, s), b)), exposed=False)

    def test_target_on_both_sides_raises(self):
        v, a, b, M = self.v, self.a, self.b, self.M
        with pytest.raises(NonLinearError, match="both sides"):
            collect(Dot(MVMul(M, v), Cross(a, v)), v)
        with pytest.raises(NonLinearError):
            collect(Dot(Cross(a, Cross(b, v)), MVMul(M, Cross(v, a))), v)

    def test_random_dots_with_target_in_a_vector_slot(self):
        """300 random Dots with the target in a vector slot: terminates, value preserved."""
        rng = np.random.default_rng(20261113)
        v, s, leaves, mats = self.v, self.s, [self.a, self.b, self.c], [self.M, self.N, self.R]

        def vec(depth, with_v):
            if depth == 0 or rng.random() < 0.25:
                return v if with_v else rng.choice(leaves)
            kind = rng.integers(0, 4)
            lh = rng.random() < 0.5
            if kind == 0:
                return Cross(vec(depth - 1, with_v and lh), vec(depth - 1, with_v and not lh))
            if kind == 1:
                return MVMul(rng.choice(mats), vec(depth - 1, with_v))
            if kind == 2:
                return SVMul(vec(depth - 1, with_v), s)
            return VAdd(vec(depth - 1, with_v and lh), vec(depth - 1, with_v and not lh))

        exposed = 0
        for _ in range(300):
            lh = rng.random() < 0.5
            expr = Dot(vec(3, lh), vec(3, not lh))
            out = _collect_preserves_value(expr, v, seeds=2)
            exposed += isinstance(out, Dot) and out.left == v
        assert exposed > 100

    def test_target_inside_matrix_factor_terminates(self):
        """F1 (fixed): the peel is guarded on the vector operand; the matrix-factor case
        raises ``AlgebraicError`` instead of looping (``tests/test_collect_fix.py``)."""
        v, a, b, M, c = self.v, self.a, self.b, self.M, self.c
        for expr in (
            Dot(MVMul(Hat(v), a), b),
            Dot(b, MVMul(Hat(v), a)),
            Dot(MVMul(VVMul(v, Transpose(c)), a), b),
            Dot(MVMul(MMMul(M, Hat(v)), a), b),
            Dot(MVMul(Hat(v), a), MVMul(M, b)),
        ):
            try:
                _collect_preserves_value(expr, v)
            except (NonLinearError, AlgebraicError):
                pass  # a typed error is acceptable; an infinite loop is not


# ---------------------------------------------------------------------------
# 5. Hat(Vee) and skew propagation
# ---------------------------------------------------------------------------


class TestSkewPropagation:
    def setup_method(self):
        self.K, self.K2 = SkewSymmMatrix("K"), SkewSymmMatrix("K2")
        self.J = Matrix("J", attr=["SymmetricMatrix"])
        self.M, self.R, self.v, self.s = Matrix("M"), SO3("R"), Vector("v"), Scalar("s")

    def test_leaf_flag_is_set(self):
        assert self.K.is_skew_symmetric and self.K.flags.is_skew_symmetric
        assert SkewSymmMatrix("K3", attr=["Constant"]).is_constant

    def test_is_skew_symmetric_is_never_wrongly_true(self):
        K, K2, J, M, v, s = self.K, self.K2, self.J, self.M, self.v, self.s
        assert not MAdd(K, J).is_skew_symmetric
        assert not MMMul(K, K2).is_skew_symmetric
        assert not MMMul(MMMul(self.R, K), Transpose(self.R)).is_skew_symmetric  # conservative
        assert not VVMul(v, Transpose(v)).is_skew_symmetric
        assert not Hat(v).is_symmetric
        assert not MAdd(J, M).is_symmetric
        assert not SMMul(M, s).is_skew_symmetric
        assert MAdd(K, K2).is_skew_symmetric and SMMul(K, s).is_skew_symmetric
        assert Transpose(K).is_skew_symmetric and Hat(v).is_skew_symmetric

    @pytest.mark.parametrize(
        "name",
        ["I", "J", "K", "Hat", "MAdd_KK", "MAdd_KJ", "SMMul", "Transpose", "MMMul_KK", "RKRt"],
    )
    def test_hat_vee_simplification_preserves_value(self, name):
        K, K2, J, R, v, s = self.K, self.K2, self.J, self.R, self.v, self.s
        X = {
            "I": I,
            "J": J,
            "K": K,
            "Hat": Hat(v),
            "MAdd_KK": MAdd(K, K2),
            "MAdd_KJ": MAdd(K, J),
            "SMMul": SMMul(K, s),
            "Transpose": Transpose(K),
            "MMMul_KK": MMMul(K, K2),
            "RKRt": MMMul(MMMul(R, K), Transpose(R)),
        }[name]
        out = full_simplify(Hat(Vee(X)))
        for seed in range(3):
            assert np.allclose(evaluate(Hat(Vee(X)), seed=seed), evaluate(out, seed=seed))
        cancelled = name in ("K", "Hat", "MAdd_KK", "SMMul", "Transpose")
        assert (str(out) == str(X)) == cancelled


# ---------------------------------------------------------------------------
# 6. extract on everything the pipeline produces
# ---------------------------------------------------------------------------


class TestExtractOnPipelineOutput:
    @pytest.mark.parametrize("build", test_eom_oracle.SYSTEMS, ids=lambda f: f.__name__)
    def test_standard_form_and_linearize_do_not_raise_and_agree_numerically(self, build):
        s = build()
        eom = compute_eom(s.L, s.dW, s.variables)
        sf = to_standard_form(eom, s.variables, inputs=s.inputs)
        rng = np.random.default_rng(3)
        for _ in range(2):
            env = sample_state(s.variables, s.inputs, rng)
            ev = Evaluator(env=env, rng=rng)
            for key, (_, eq) in eom.items():
                e_val = np.atleast_2d(np.asarray(ev(eq), dtype=float)).reshape(-1, 1)
                total = np.zeros_like(e_val)
                if sf[key].f is not None:
                    total = total + np.atleast_2d(np.asarray(ev(sf[key].f), dtype=float)).reshape(
                        -1, 1
                    )
                for blocks in (sf[key].M, sf[key].G):
                    for col, coeff in blocks.items():
                        cval = np.atleast_2d(np.asarray(ev(coeff), dtype=float))
                        xval = np.atleast_2d(np.asarray(ev.env[col], dtype=float)).reshape(-1, 1)
                        if cval.shape == (3, 1) and xval.shape == (3, 1):
                            cval = cval.T  # Dot(coeff, vector accel) in a scalar equation
                        total = total + cval @ xval
                assert np.max(np.abs(total - e_val)) < 1e-9
        lin = linearize(eom, s.variables, inputs=s.inputs)
        assert finite_difference_check(lin, eom, s.variables, s.inputs, trials=2) < 1e-6

    def test_both_sides_guard_in_cross_keeps_identically_zero_cases(self):
        q, a, s = Vector("q"), Vector("a"), Scalar("s")
        assert extract_linear_coeff(Cross(q, q), q).is_zero
        assert extract_linear_coeff(full_simplify(Cross(q, SVMul(q, s))), q).is_zero
        with pytest.raises(NonLinearError):
            extract_linear_coeff(Cross(q, Cross(a, q)), q)
        with pytest.raises(NonLinearError):
            extract_linear_coeff(Cross(Cross(a, q), q), q)

    def test_zero_or_raise_returns_zero_only_for_query_free_terms(self):
        q, a, b, s, M = Vector("q"), Vector("a"), Vector("b"), Scalar("s"), Matrix("M")
        assert extract_linear_coeff(SVMul(a, Dot(b, MVMul(M, b))), q).is_zero
        assert extract_linear_coeff(Hat(a), q).is_zero
        assert extract_linear_coeff(Transpose(SMMul(VVMul(a, Transpose(b)), s)), q).is_zero
        with pytest.raises(AlgebraicError):
            extract_linear_coeff(Transpose(SMMul(VVMul(a, Transpose(b)), Dot(q, b))), q)


# ---------------------------------------------------------------------------
# 7. The oracle itself
# ---------------------------------------------------------------------------


def _stiff_point_mass(k_value):
    s = ps.stiff_point_mass(k_value)
    eom = compute_eom(s.L, s.dW, s.variables)
    return s.L, s.dW, s.variables, eom, str(s.symbols["x"].delta()), s.inputs


def _with_term_changed(eom, key, pred, factor):
    vec, eq = eom[key]
    (term,) = [t for t in eq.nodes if pred(str(t))]
    return {key: (vec, eq - term * factor)}


class TestOracle:
    def test_scalar_with_value_is_folded_by_engine_and_used_by_evaluator(self):
        m2, g = Scalar("m", value=2, attr=["Constant"]), Scalar("g", attr=["Constant"])
        e3, x, f = Vector("e3", attr=["Constant"]), Vector("x"), Vector("f")
        v = x.t_diff()
        L = HALF * m2 * Dot(v, v) - m2 * g * Dot(x, e3)
        variables = SystemVariables(vectors=[x])
        eom = compute_eom(L, Dot(x.delta(), f), variables)
        assert "m" not in str(eom[str(x.delta())][1])  # folded to the literal 2
        assert check_eom(L, Dot(x.delta(), f), variables, eom, inputs=[f]) < TOL

    def test_flipped_terms_are_detected_with_random_parameters(self):
        L, dW, variables, eom, key, inputs = _stiff_point_mass(None)
        wrong = _with_term_changed(eom, key, lambda t: "e3" in t, 2)
        assert check_eom(L, dW, variables, wrong, inputs=inputs) > 1e-1

    def test_flipped_terms_are_detected_at_h_1e_4(self):
        from tests import vbl_systems

        s = vbl_systems.free_rigid_body()
        eom = compute_eom(s.L, s.dW, s.variables)
        key = next(iter(eom))
        vec, eq = eom[key]
        for term in eq.nodes:
            wrong = {key: (vec, eq - term * 2)}
            assert check_eom(s.L, s.dW, s.variables, wrong, inputs=s.inputs, h=H) > 1e-1

    def test_stiff_term_hides_gravity_and_the_oracle_says_so(self):
        """F2 (resolved by declaring the blind spot): next to k = 1e7 the gravity term
        (tangent size of order one) is below the resolution rtol * max(1, |E_true|) of
        ``check_eom`` (|E_true| is of order 1e7), so a flipped or dropped gravity term
        moves the relative error by only ~5e-7 and cannot be detected by it alone. Instead of
        silently passing, ``assert_eom`` reports every unresolved term of the correct
        equation, naming the gravity term.
        """
        L, dW, variables, eom, key, inputs = _stiff_point_mass(1e7)
        assert check_eom(L, dW, variables, eom, inputs=inputs) < TOL  # the correct one passes
        with pytest.raises(AssertionError, match="below the resolution") as info:
            assert_eom(L, dW, variables, eom, inputs=inputs)
        assert "e3" in str(info.value)  # the gravity term is one of the unresolved ones
        # The same equation is accepted on the relative error alone when asked.
        assert assert_eom(L, dW, variables, eom, inputs=inputs, require_resolved=False) < TOL

    def test_dropped_gravity_is_detected_once_the_parameters_are_comparable(self):
        """The resolution message tells the user to scale the parameters; with random
        parameters of order one the dropped and flipped gravity terms are detected."""
        L, dW, variables, eom, key, inputs = _stiff_point_mass(None)
        dropped = _with_term_changed(eom, key, lambda t: "e3" in t, 1)
        assert check_eom(L, dW, variables, dropped, inputs=inputs) > 1e-2
        assert_eom(L, dW, variables, eom, inputs=inputs)  # correct and fully resolved

    def test_correct_equation_scores_well_below_the_gate_at_default_h(self):
        """F3 (fixed): the default step is 1e-4; at 1e-5 the nested central differences
        were roundoff-dominated and this correct equation scored 1.4e-5 at seed 2."""
        s = adv_r3_s2_pendulum_with_gyro()
        eom = compute_eom(s.L, s.dW, s.variables)
        assert check_eom(s.L, s.dW, s.variables, eom, inputs=s.inputs, h=H) < 1e-6  # it is right
        worst = max(
            check_eom(s.L, s.dW, s.variables, eom, inputs=s.inputs, seed=seed) for seed in range(3)
        )
        assert worst < 5e-6, worst  # a correct equation must sit well below the 1e-5 gate

    def test_default_step_is_near_the_noise_minimum(self):
        """F3 (fixed): the default step is at the noise minimum of the nested scheme."""
        s = adv_r3_s2_pendulum_with_gyro()
        eom = compute_eom(s.L, s.dW, s.variables)
        e_default = max(
            check_eom(s.L, s.dW, s.variables, eom, inputs=s.inputs, seed=sd) for sd in range(3)
        )
        e_coarser = max(
            check_eom(s.L, s.dW, s.variables, eom, inputs=s.inputs, seed=sd, h=1e-4)
            for sd in range(3)
        )
        assert e_default <= 3 * e_coarser


# ---------------------------------------------------------------------------
# 8. The collect loop is reachable from compute_eom
# ---------------------------------------------------------------------------


class TestPipelineReachesCollectLoop:
    """F1 (fixed): L = ½ Ω·JΩ + a·((Ω cᵀ) b) is a legal Lagrangian ((a·Ω)(c·b)); the outer
    product survives simplify, so collect sees the target inside a matrix factor. It used to
    recurse forever; it now raises a typed ``AlgebraicError`` naming the term (the message
    is pinned in ``tests/test_collect_fix.py``)."""

    def test_outer_product_on_the_velocity_so3(self):
        s = ps.so3_outer_product_on_velocity()
        try:
            eom = compute_eom(s.L, s.dW, s.variables)
        except AlgebraicError:
            return  # a typed error naming the term would be acceptable
        assert check_eom(s.L, s.dW, s.variables, eom, inputs=s.inputs, h=H) < TOL

    def test_outer_product_on_the_velocity_s2(self):
        s = ps.s2_outer_product_on_velocity()
        try:
            eom = compute_eom(s.L, s.dW, s.variables)
        except AlgebraicError:
            return
        assert check_eom(s.L, s.dW, s.variables, eom, inputs=s.inputs, h=H) < TOL
