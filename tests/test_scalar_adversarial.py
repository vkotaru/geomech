"""Adversarial review tests for scalar configuration variables.

Every expected value below was derived by hand from the Lagrangian before
the code was run; the derivation (or the reason the value is what it is)
is in each docstring.

Convention (same as the S2 / SO3 pipeline tests): the pipeline returns the
coefficient of the variation in the action integrand after IBP, so each
equation reads  ``-(d/dt ∂L/∂q̇ - ∂L/∂q) + f = 0``, and ``to_standard_form``
reports ``Σ M[a] a + f + Σ G[u] u = 0`` with a *negative* mass coefficient.

Tests marked ``xfail(strict=True)`` record bugs; the reason string says
which.
"""

import pytest

from geomech import (
    S2,
    SO3,
    Add,
    Cross,
    Dot,
    Hat,
    Matrix,
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
    to_standard_form,
)
from geomech.core.base.expressions import One
from geomech.core.math.extract import extract_linear_coeff
from geomech.core.math.ibp import integrate_by_parts
from geomech.utils.errors import AlgebraicError, NonLinearError

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

HALF = Scalar("0.5", value=0.5, attr=["Constant"])

DDX = "\\frac{d}{dt}(\\frac{d}{dt}(x))"
DDL = "\\frac{d}{dt}(\\frac{d}{dt}(l))"
DL = "\\frac{d}{dt}(l)"
DW = "\\frac{d}{dt}(\\omega_{q})"


def _strs(d):
    """str() every value of a {name: Expr} dict for exact comparison."""
    return {k: str(v) for k, v in d.items()}


def _numeric_values(expr, out=None):
    """Collect the Python values of every numeric leaf in *expr*."""
    if out is None:
        out = []
    nodes = getattr(expr, "nodes", None)
    if nodes:
        for n in nodes:
            _numeric_values(n, out)
    elif isinstance(expr, Scalar) and expr.value is not None:
        out.append(expr.value)
    return out


def _spring_mass():
    """L = ½ m ẋ² − ½ k x²,  δW = f δx."""
    m, k = getScalars("m k", attr=["Constant"])
    x, f = Scalar("x"), Scalar("f")
    v = x.t_diff()
    L = HALF * m * v * v - HALF * k * x * x
    variables = SystemVariables(scalars=[x])
    eom = compute_eom(L, x.delta() * f, variables)
    return eom, to_standard_form(eom, variables, [f]), x, f


def _spring_mass_damper():
    """L = ½ m ẋ² − ½ k x²,  δW = (f − c ẋ) δx."""
    m, k, c = getScalars("m k c", attr=["Constant"])
    x, f = Scalar("x"), Scalar("f")
    v = x.t_diff()
    L = HALF * m * v * v - HALF * k * x * x
    dW = x.delta() * f - c * v * x.delta()
    variables = SystemVariables(scalars=[x])
    eom = compute_eom(L, dW, variables)
    return eom, to_standard_form(eom, variables, [f]), x, f


def _double_spring_mass():
    """L = ½ m₁ ẋ₁² + ½ m₂ ẋ₂² − ½ k₁ x₁² − ½ k₂ (x₂ − x₁)²,  δW = f₁ δx₁ + f₂ δx₂."""
    m1, m2, k1, k2 = getScalars("m1 m2 k1 k2", attr=["Constant"])
    x1, x2, f1, f2 = Scalar("x1"), Scalar("x2"), Scalar("f1"), Scalar("f2")
    v1, v2 = x1.t_diff(), x2.t_diff()
    KE = HALF * m1 * v1 * v1 + HALF * m2 * v2 * v2
    PE = HALF * k1 * x1 * x1 + HALF * k2 * (x2 - x1) * (x2 - x1)
    dW = x1.delta() * f1 + x2.delta() * f2
    variables = SystemVariables(scalars=[x1, x2])
    eom = compute_eom(KE - PE, dW, variables)
    return eom, to_standard_form(eom, variables, [f1, f2]), x1, x2


def _elastic_pendulum(qdot_form=False):
    """Elastic spherical pendulum on R × S² (spec table, row 4).

    x = l q, KE = ½ m (l̇² + l² ω·ω) (or ½ m ẋ·ẋ when ``qdot_form``),
    PE = m g l (q·e₃) + ½ k (l − l₀)², δW = f_l δl + ξ·f_q.
    """
    m, g, k, l0 = getScalars("m g k l0", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    l, q = Scalar("l"), S2("q")
    omega, xi = q.get_tangent_vector(), q.get_variation_vector()
    f_l, f_q = Scalar("f_l"), Vector("f_q")
    ld = l.t_diff()
    if qdot_form:
        v = (l * q).t_diff()
        KE = HALF * m * Dot(v, v)
    else:
        KE = HALF * m * (ld * ld + l * l * Dot(omega, omega))
    PE = m * g * l * Dot(q, e3) + HALF * k * (l - l0) * (l - l0)
    dW = l.delta() * f_l + Dot(xi, f_q)
    variables = SystemVariables(scalars=[l], vectors=[q])
    eom = compute_eom(KE - PE, dW, variables)
    sf = to_standard_form(eom, variables, [f_l, f_q])
    return eom, sf, l, q


def _r3_x_r_coupled():
    """Point mass on a prismatic slide: R³ × R with an inertial cross term.

    v = ẋ + l̇ e₁ (‖e₁‖ = 1),  KE = ½ m v·v = ½ m (ẋ·ẋ + 2 l̇ ẋ·e₁ + l̇²),
    δW = δx·F + f_l δl.
    """
    (m,) = getScalars("m", attr=["Constant"])
    e1 = Vector("e1", attr=["Constant", "UnitNorm"])
    x, l = Vector("x"), Scalar("l")
    F, f_l = Vector("F"), Scalar("f_l")
    v = x.t_diff() + e1 * l.t_diff()
    L = HALF * m * Dot(v, v)
    dW = Dot(x.delta(), F) + l.delta() * f_l
    variables = SystemVariables(scalars=[l], vectors=[x])
    eom = compute_eom(L, dW, variables)
    return eom, to_standard_form(eom, variables, [F, f_l]), x, l


# ---------------------------------------------------------------------------
# 1. Exact equations for the three pure-scalar systems
# ---------------------------------------------------------------------------


class TestSpringMassExact:
    """δL = m ẋ δẋ − k x δx;  IBP: m ẋ δẋ → −m ẍ δx.

    Coefficient of δx:  −m ẍ − k x + f = 0.
    Standard form:      M[ẍ] = −m,  f = −k x,  G[f] = 1.
    """

    def setup_method(self):
        self.eom, self.sf, self.x, self.f = _spring_mass()
        self.key = str(Variation(self.x))

    def test_eom_exact(self):
        _, eqn = self.eom[self.key]
        assert str(eqn) == f"(m*{DDX}*(-1)+k*x*(-1)+f)"

    def test_M_exact(self):
        assert _strs(self.sf[self.key].M) == {DDX: "m*(-1)"}

    def test_f_exact(self):
        assert str(self.sf[self.key].f) == "k*x*(-1)"

    def test_G_exact(self):
        assert _strs(self.sf[self.key].G) == {"f": "1"}

    def test_sign_convention_matches_s2_pendulum(self):
        """Both equations are −(inertia·accel) + … + input, so M is negative
        for the scalar system exactly as it is for the S² pendulum
        (``Imll(-1)`` in test_manifold_pipelines)."""
        m, g, l = getScalars("m g l", attr=["Constant"])
        e3 = Vector("e3", attr=["Constant"])
        q, f = S2("q"), Vector("f")
        omega, xi = q.get_tangent_vector(), q.get_variation_vector()
        L = HALF * m * l * l * Dot(omega, omega) - m * g * l * Dot(q, e3)
        variables = SystemVariables(vectors=[q])
        s2 = to_standard_form(compute_eom(L, Dot(xi, f), variables), variables, [f])
        assert str(s2[str(xi)].M[DW]) == "I*m*l*l*(-1)"
        assert str(self.sf[self.key].M[DDX]) == "m*(-1)"
        assert str(s2[str(xi)].G["f"]) == "I"
        assert str(self.sf[self.key].G["f"]) == "1"


class TestSpringMassDamperExact:
    """δW = f δx − c ẋ δx adds −c ẋ to the δx coefficient:

    −m ẍ − k x + f − c ẋ = 0.
    Standard form:  M[ẍ] = −m,  f = −k x − c ẋ,  G[f] = 1.
    """

    def setup_method(self):
        self.eom, self.sf, self.x, self.f = _spring_mass_damper()
        self.key = str(Variation(self.x))

    def test_eom_exact(self):
        _, eqn = self.eom[self.key]
        assert str(eqn) == f"(m*{DDX}*(-1)+k*x*(-1)+f+c*\\frac{{d}}{{dt}}(x)*(-1))"

    def test_M_exact(self):
        assert _strs(self.sf[self.key].M) == {DDX: "m*(-1)"}

    def test_f_exact(self):
        assert str(self.sf[self.key].f) == "(k*x*(-1)+c*\\frac{d}{dt}(x)*(-1))"

    def test_G_exact(self):
        assert _strs(self.sf[self.key].G) == {"f": "1"}

    def test_damping_term_not_classified_as_acceleration(self):
        """ẋ must not be mistaken for ẍ by the ``has`` partition."""
        eq = self.sf[self.key]
        assert not eq.f.has(TimeDerivative(TimeDerivative(self.x)))
        assert eq.f.has(TimeDerivative(self.x))

    def test_numeric_coefficients_have_one_type(self):
        """All numeric constants in one equation come from the same
        arithmetic (−1 × …); they should not print as both ``(-1)`` and
        ``(-1.0)``. The S² pipeline strings use ``(-1)`` throughout.
        Fixed in phase 5a (a float with an integer value folds to an int)."""
        _, eqn = self.eom[self.key]
        types = {type(v) for v in _numeric_values(eqn)}
        assert len(types) == 1, str(eqn)


class TestDoubleSpringMassExact:
    """∂PE/∂x₁ = k₁ x₁ − k₂ (x₂ − x₁),  ∂PE/∂x₂ = k₂ (x₂ − x₁).

    δx₁:  −m₁ ẍ₁ − k₁ x₁ + k₂ x₂ − k₂ x₁ + f₁ = 0
    δx₂:  −m₂ ẍ₂ − k₂ x₂ + k₂ x₁ + f₂ = 0
    Each M has exactly one key (no cross inertia), each G exactly one input.
    """

    def setup_method(self):
        self.eom, self.sf, self.x1, self.x2 = _double_spring_mass()
        self.k1, self.k2 = str(Variation(self.x1)), str(Variation(self.x2))

    def test_x1_eom_exact(self):
        _, eqn = self.eom[self.k1]
        ddx1 = "\\frac{d}{dt}(\\frac{d}{dt}(x1))"
        assert str(eqn) == f"(m1*{ddx1}*(-1)+k1*x1*(-1)+k2*x2+k2*x1*(-1)+f1)"

    def test_x2_eom_exact(self):
        _, eqn = self.eom[self.k2]
        ddx2 = "\\frac{d}{dt}(\\frac{d}{dt}(x2))"
        assert str(eqn) == f"(m2*{ddx2}*(-1)+k2*x2*(-1)+k2*x1+f2)"

    def test_x1_standard_form_exact(self):
        eq = self.sf[self.k1]
        assert _strs(eq.M) == {"\\frac{d}{dt}(\\frac{d}{dt}(x1))": "m1*(-1)"}
        assert str(eq.f) == "(k1*x1*(-1)+k2*x2+k2*x1*(-1))"
        assert _strs(eq.G) == {"f1": "1"}

    def test_x2_standard_form_exact(self):
        eq = self.sf[self.k2]
        assert _strs(eq.M) == {"\\frac{d}{dt}(\\frac{d}{dt}(x2))": "m2*(-1)"}
        assert str(eq.f) == "(k2*x2*(-1)+k2*x1)"
        assert _strs(eq.G) == {"f2": "1"}

    def test_no_cross_contamination(self):
        for key, (_, eqn) in self.eom.items():
            assert not eqn.has(Variation(self.x1)), key
            assert not eqn.has(Variation(self.x2)), key


class TestCoupledScalarAccelerations:
    """L = ½ m (ẋ₁ + ẋ₂)²:  ∂L/∂ẋᵢ = m (ẋ₁ + ẋ₂), so both equations read
    −m ẍ₁ − m ẍ₂ + fᵢ = 0 and M has two keys, both −m."""

    def test_both_accelerations_in_M(self):
        (m,) = getScalars("m", attr=["Constant"])
        x1, x2, f1, f2 = Scalar("x1"), Scalar("x2"), Scalar("f1"), Scalar("f2")
        s = x1.t_diff() + x2.t_diff()
        variables = SystemVariables(scalars=[x1, x2])
        eom = compute_eom(HALF * m * s * s, x1.delta() * f1 + x2.delta() * f2, variables)
        sf = to_standard_form(eom, variables, [f1, f2])
        ddx1 = "\\frac{d}{dt}(\\frac{d}{dt}(x1))"
        ddx2 = "\\frac{d}{dt}(\\frac{d}{dt}(x2))"
        eq1 = sf[str(Variation(x1))]
        eq2 = sf[str(Variation(x2))]
        assert str(eom[str(Variation(x1))][1]) == f"(m*{ddx1}*(-1)+m*{ddx2}*(-1)+f1)"
        assert _strs(eq1.M) == {ddx1: "m*(-1)", ddx2: "m*(-1)"}
        assert _strs(eq2.M) == {ddx1: "m*(-1)", ddx2: "m*(-1)"}
        assert str(eq1.f) == "0" and str(eq2.f) == "0"
        assert _strs(eq1.G) == {"f1": "1"} and _strs(eq2.G) == {"f2": "1"}


class TestVaryingInertia:
    """L = ½ (m₀ + m₁ x) ẋ²  (the "varying inertia" case named in the spec goal).

    ∂L/∂x = ½ m₁ ẋ²,  d/dt ∂L/∂ẋ = (m₀ + m₁ x) ẍ + m₁ ẋ².
    Coefficient of δx:  −(m₀ + m₁ x) ẍ − m₁ ẋ² + ½ m₁ ẋ² + f
                      = −(m₀ + m₁ x) ẍ − ½ m₁ ẋ² + f = 0.
    """

    def setup_method(self):
        m0, m1 = getScalars("m0 m1", attr=["Constant"])
        self.x, self.f = Scalar("x"), Scalar("f")
        v = self.x.t_diff()
        L = HALF * (m0 + m1 * self.x) * v * v
        variables = SystemVariables(scalars=[self.x])
        self.eom = compute_eom(L, self.x.delta() * self.f, variables)
        self.sf = to_standard_form(self.eom, variables, [self.f])
        self.eq = self.sf[str(Variation(self.x))]

    def test_mass_matrix_exact(self):
        assert _strs(self.eq.M) == {DDX: "(m0*(-1)+m1*x*(-1))"}
        assert _strs(self.eq.G) == {"f": "1"}

    def test_f_is_mathematically_correct_but_uncombined(self):
        """The two ẋ² terms are present with coefficients +0.5 and −1.0.
        The new printed form is the same expression simplified further (phase 5a merges
        the two terms into −0.5 m₁ ẋ²)."""
        assert str(self.eq.f) == ("m1*\\frac{d}{dt}(x)*\\frac{d}{dt}(x)*(-0.5)")

    def test_f_combines_like_terms(self):
        """½ m₁ ẋ² − m₁ ẋ² must collapse to a single term −½ m₁ ẋ².
        Fixed in phase 5a (like terms are compared on flattened Mul chains)."""
        assert not isinstance(self.eq.f, Add), str(self.eq.f)

    def test_x_squared_inertia_exact(self):
        """L = ½ m x² ẋ²:  ∂L/∂x = m x ẋ²,  d/dt ∂L/∂ẋ = 2 m x ẋ² + m x² ẍ,
        so the δx coefficient is −m x ẋ² − m x² ẍ + f."""
        (m,) = getScalars("m", attr=["Constant"])
        x, f = Scalar("x"), Scalar("f")
        v = x.t_diff()
        variables = SystemVariables(scalars=[x])
        eom = compute_eom(HALF * m * x * x * v * v, x.delta() * f, variables)
        _, eqn = eom[str(Variation(x))]
        dx = "\\frac{d}{dt}(x)"
        assert str(eqn) == f"(m*x*{dx}*{dx}*(-1)+m*x*x*{DDX}*(-1)+f)"
        eq = to_standard_form(eom, variables, [f])[str(Variation(x))]
        assert _strs(eq.M) == {DDX: "m*x*x*(-1)"}


# ---------------------------------------------------------------------------
# 2. Mixed R × S²: elastic spherical pendulum
# ---------------------------------------------------------------------------


class TestElasticSphericalPendulum:
    """Hand derivation (S² convention q̇ = ω × q, δq = ξ × q, δω = ξ̇ − 2 ω × ξ).

    δKE = m l̇ δl̇ + m l (ω·ω) δl + m l² ω·ξ̇          (ω·(ω×ξ) = 0)
    δPE = m g (q·e₃) δl + m g l ξ·(q×e₃) + k (l−l₀) δl  ((ξ×q)·e₃ = ξ·(q×e₃))
    IBP: m l̇ δl̇ → −m l̈ δl ;  m l² ω·ξ̇ → −ξ·d/dt(m l² ω) = −ξ·(2 m l l̇ ω + m l² ω̇)

    δl:  −m l̈ + m l ω·ω − m g q·e₃ − k l + k l₀ + f_l = 0
    ξ :  −m l² ω̇ − 2 m l l̇ ω − m g l (q×e₃) + f_q = 0

    Standard form:
      δl:  M[l̈] = −m,          f = m l ω·ω − m g q·e₃ − k l + k l₀,     G[f_l] = 1
      ξ :  M[ω̇] = −m l² I,     f = −2 m l l̇ ω − m g l q×e₃,           G[f_q] = I
    """

    def setup_method(self):
        self.eom, self.sf, self.l, self.q = _elastic_pendulum()
        self.omega = self.q.get_tangent_vector()
        self.xi = self.q.get_variation_vector()
        self.kl, self.kx = str(Variation(self.l)), str(self.xi)

    def test_two_equations(self):
        assert set(self.eom) == {self.kl, self.kx}

    def test_l_eom_exact(self):
        _, eqn = self.eom[self.kl]
        assert str(eqn) == (
            f"(m*{DDL}*(-1)+m*l*Dot(\\omega_{{q}},\\omega_{{q}})+m*g*Dot(q,e3)*(-1)+k*l*(-1)+k*l0+f_l)"
        )

    def test_l_standard_form_exact(self):
        eq = self.sf[self.kl]
        assert _strs(eq.M) == {DDL: "m*(-1)"}
        assert str(eq.f) == "(m*l*Dot(\\omega_{q},\\omega_{q})+m*g*Dot(q,e3)*(-1)+k*l*(-1)+k*l0)"
        assert _strs(eq.G) == {"f_l": "1"}

    def test_xi_mass_and_input_exact(self):
        eq = self.sf[self.kx]
        assert _strs(eq.M) == {DW: "I*m*l*l*(-1)"}
        assert _strs(eq.G) == {"f_q": "I"}

    def test_no_cross_contamination(self):
        _, eq_l = self.eom[self.kl]
        _, eq_x = self.eom[self.kx]
        assert not eq_l.has(self.xi) and not eq_l.has(Variation(self.l))
        assert not eq_x.has(self.xi) and not eq_x.has(Variation(self.l))
        assert not eq_x.has(TimeDerivative(self.xi))

    def test_xi_eom_has_coriolis_term(self):
        """The ξ equation must contain l̇ (from d/dt(m l² ω) = 2 m l l̇ ω + m l² ω̇).
        Expected: (d/dt(ω)mll(-1) + ω l̇ l m (-2) + Cross(q,e3)mgl(-1) + f_q), in
        some term order. The code returns the l̇-free spherical-pendulum equation."""
        _, eqn = self.eom[self.kx]
        assert eqn.has(TimeDerivative(self.l)), str(eqn)

    def test_xi_f_has_coriolis_term(self):
        """f = −2 m l l̇ ω − m g l q×e₃: two addends, one containing l̇.

        Fixed by F1. (Arbitration: the original assertion asked for ``Add``; a
        vector sum is a ``VAdd`` in this library.)"""
        f = self.sf[self.kx].f
        assert f.has(TimeDerivative(self.l)), str(f)
        assert isinstance(f, VAdd) and len(f.nodes) == 2, str(f)


class TestElasticSphericalPendulumQdotForm:
    """Same system with KE = ½ m ẋ·ẋ, x = l q, ẋ = l̇ q + l (ω × q).

    ẋ·ẋ = l̇² + l² (ω×q)·(ω×q) + 2 l l̇ q·(ω×q) = l̇² + l² ω·ω, so the physics is
    identical to the ω form. The ξ̇ terms are m l² ω·ξ̇ (IBP → −2 m l l̇ ω − m l² ω̇);
    the remaining ξ terms, m ẋ·[l̇ (ξ×q) − l (ω×ξ)×q + l ω×(ξ×q)], sum to
    −m l l̇ ω·ξ + m l l̇ ω·ξ = 0.

    ξ:  −m l² q×(ω̇×q) − 2 m l l̇ ω − m g l (q×e₃) + f_q = 0
    f = −2 m l l̇ ω − m g l q×e₃  (two addends).
    """

    def setup_method(self):
        self.eom, self.sf, self.l, self.q = _elastic_pendulum(qdot_form=True)
        self.xi = self.q.get_variation_vector()

    def test_l_eom_exact(self):
        """δl: −m l̈ + m l (ω×q)·(ω×q) − m g q·e₃ − k l + k l₀ + f_l = 0 —
        correct, although (ω×q)·(ω×q) is left unreduced (it equals ω·ω)."""
        _, eqn = self.eom[str(Variation(self.l))]
        assert str(eqn) == (
            "(m*l*Dot(Cross(\\omega_{q},q),Cross(\\omega_{q},q))"
            f"+m*{DDL}*(-1)+m*g*Dot(q,e3)*(-1)+k*l*(-1)+k*l0+f_l)"
        )

    def test_xi_mass_matrix_exact(self):
        """M[ω̇] = m l² q̂ q̂ = −m l² (I − q qᵀ): the tangent projection of −m l² I."""
        assert _strs(self.sf[str(self.xi)].M) == {DW: "Hat(q)*Hat(q)*m*l*l"}

    def test_xi_f_exact_structure(self):
        """f must be −2 m l l̇ ω − m g l q×e₃ (two addends).
        Fixed in phase 5a (the pair +m l l̇ ω − m l l̇ ω now cancels)."""
        f = self.sf[str(self.xi)].f
        assert isinstance(f, VAdd) and len(f.nodes) == 2, str(f)


# ---------------------------------------------------------------------------
# 3. Cross-coupling in M for R³ × R (the cleanest exhibit of the IBP bug)
# ---------------------------------------------------------------------------


class TestR3xRCoupledMassMatrix:
    """KE = ½ m (ẋ·ẋ + 2 l̇ ẋ·e₁ + l̇²).

    δKE = m ẋ·δẋ + m l̇ e₁·δẋ + m (ẋ·e₁) δl̇ + m l̇ δl̇
    IBP: → −δx·(m ẍ + m l̈ e₁) − δl (m ẍ·e₁ + m l̈)

    δx:  −m ẍ − m l̈ e₁ + F = 0        M[ẍ] = −m I,  M[l̈] = −m e₁
    δl:  −m e₁·ẍ − m l̈ + f_l = 0      M[ẍ] = −m e₁ᵀ, M[l̈] = −m
    The mass matrix is symmetric: M_x[l̈] = (M_l[ẍ])ᵀ.
    """

    def setup_method(self):
        self.eom, self.sf, self.x, self.l = _r3_x_r_coupled()
        self.kx, self.kl = str(Variation(self.x)), str(Variation(self.l))

    def test_l_equation_exact(self):
        _, eqn = self.eom[self.kl]
        assert str(eqn) == f"(m*Dot(e1,{DDX})*(-1)+m*{DDL}*(-1)+f_l)"
        eq = self.sf[self.kl]
        assert _strs(eq.M) == {DDL: "m*(-1)", DDX: "e1*m*(-1)"}
        assert str(eq.f) == "0"
        assert _strs(eq.G) == {"f_l": "1"}

    def test_x_equation_has_l_ddot(self):
        """Expected (d/dt(d/dt(x))(-1)m + e1 d/dt(d/dt(l)) m (-1) + F); the code
        returns (d/dt(d/dt(x))(-1)m + F)."""
        _, eqn = self.eom[self.kx]
        assert eqn.has(TimeDerivative(TimeDerivative(self.l))), str(eqn)


# ---------------------------------------------------------------------------
# 4. Root cause at unit level: integrate_by_parts
# ---------------------------------------------------------------------------


class TestIntegrateByPartsScalarFactors:
    """∫ a Dot(ξ̇, ω) dt = −∫ Dot(ξ, d/dt(a ω)) dt = −∫ (ȧ Dot(ξ, ω) + a Dot(ξ, ω̇)) dt
    for a non-constant scalar a. The scalar-target path (new in this feature)
    gets this right; the vector-target path (pre-existing) does not."""

    def test_scalar_target_inside_dot_after_pull(self):
        """Dot(q δl̇, w) → pull → δl̇ Dot(q, w) → IBP → −δl (Dot(q̇, w) + Dot(q, ẇ))."""
        l = Scalar("l")
        dld = l.delta().t_diff()
        q, w = Vector("q"), Vector("w")
        out = integrate_by_parts(Dot(SVMul(q, dld), w), [dld])
        assert str(out) == (
            "(\\delta{l}*Dot(\\frac{d}{dt}(q),w)*(-1)+\\delta{l}*Dot(q,\\frac{d}{dt}(w))*(-1))"
        )

    def test_scalar_target_with_nonconstant_prefactor(self):
        """a Dot(q δl̇, w) → −δl (ȧ Dot(q,w) + a Dot(q̇,w) + a Dot(q,ẇ)): ȧ is kept."""
        l, a = Scalar("l"), Scalar("a")
        dld = l.delta().t_diff()
        q, w = Vector("q"), Vector("w")
        out = integrate_by_parts(Mul(a, Dot(SVMul(q, dld), w)), [dld])
        assert out.has(TimeDerivative(a))
        assert str(out) == (
            "(\\delta{l}*\\frac{d}{dt}(a)*Dot(q,w)*(-1)+\\delta{l}*a*Dot(\\frac{d}{dt}(q),w)*(-1)"
            "+\\delta{l}*a*Dot(q,\\frac{d}{dt}(w))*(-1))"
        )

    def test_vector_target_with_nonconstant_prefactor(self):
        """a Dot(ξ̇, ω) → −ȧ Dot(ξ, ω) − a Dot(ξ, ω̇). The code returns only the second term."""
        a = Scalar("a")
        xi, omega = Vector("xi"), Vector("omega")
        target = TimeDerivative(xi)
        out = integrate_by_parts(Mul(a, Dot(target, omega)), [target])
        assert out.has(TimeDerivative(a)), str(out)

    def test_s2_pendulum_with_time_varying_length(self):
        """Same bug without SystemVariables.scalars: a prescribed l(t) gives
        m l² ω̇ + 2 m l l̇ ω + m g l q×e₃ = f, so l̇ must appear."""
        m, g = getScalars("m g", attr=["Constant"])
        l = Scalar("l")  # not constant, not a config variable
        e3 = Vector("e3", attr=["Constant"])
        q, f = S2("q"), Vector("f")
        omega, xi = q.get_tangent_vector(), q.get_variation_vector()
        L = HALF * m * l * l * Dot(omega, omega) - m * g * l * Dot(q, e3)
        eom = compute_eom(L, Dot(xi, f), SystemVariables(vectors=[q]), params=[l])
        _, eqn = eom[str(xi)]
        assert eqn.has(TimeDerivative(l)), str(eqn)


# ---------------------------------------------------------------------------
# 5. Scalar config variable inside vector expressions (extract)
# ---------------------------------------------------------------------------


class TestScalarQueryInsideVectorOps:
    """Coefficient of l in each expression, by linearity of the vector ops."""

    def setup_method(self):
        self.l = Scalar("l")
        self.a, self.b, self.c = Vector("a"), Vector("b"), Vector("c")
        self.M = Matrix("M")

    def test_cross_left(self):
        """(a l) × b = (a × b) l = −b̂ a l."""
        out = extract_linear_coeff(Cross(SVMul(self.a, self.l), self.b), self.l)
        assert str(out) == "Hat(b)*(-1)*a"

    def test_cross_right(self):
        """a × (b l) = â b l."""
        out = extract_linear_coeff(Cross(self.a, SVMul(self.b, self.l)), self.l)
        assert str(out) == "Hat(a)*b"

    def test_dot_of_cross(self):
        """((a l) × b)·c = (−b̂ a · c) l."""
        out = extract_linear_coeff(Dot(Cross(SVMul(self.a, self.l), self.b), self.c), self.l)
        assert str(out) == "Dot(Hat(b)*(-1)*a,c)"

    def test_mvmul_scalar_in_matrix(self):
        """(M l) a = (M a) l."""
        out = extract_linear_coeff(MVMul(SMMul(self.M, self.l), self.a), self.l)
        assert str(out) == "M*a"

    def test_vvmul_inner_product_both_sides(self):
        """aᵀ(b l) = (aᵀ b) l and (a l)ᵀ b = (aᵀ b) l."""
        out1 = extract_linear_coeff(VVMul(Transpose(self.a), SVMul(self.b, self.l)), self.l)
        out2 = extract_linear_coeff(VVMul(Transpose(SVMul(self.a, self.l)), self.b), self.l)
        assert str(out1) == "(a)'*b"
        assert str(out2) == "(a)'*b"

    def test_quadratic_in_svmul_raises(self):
        with pytest.raises(NonLinearError):
            extract_linear_coeff(SVMul(SVMul(self.a, self.l), self.l), self.l)
        with pytest.raises(NonLinearError):
            extract_linear_coeff(Cross(SVMul(self.a, self.l), SVMul(self.b, self.l)), self.l)

    def test_dot_both_sides_raises(self):
        with pytest.raises(NonLinearError):
            extract_linear_coeff(Dot(SVMul(self.a, self.l), SVMul(self.b, self.l)), self.l)

    @pytest.mark.parametrize(
        "label, build, expected",
        [
            ("Hat", lambda a, l, M: Hat(SVMul(a, l)), "Hat(a)"),
            ("Hat*vec", lambda a, l, M: MVMul(Hat(SVMul(a, l)), Vector("w")), "Hat(a)*w"),
            ("Transpose", lambda a, l, M: Transpose(SVMul(a, l)), "(a)'"),
            ("Vee", lambda a, l, M: Vee(SMMul(M, l)), "Vee(M)"),
        ],
    )
    def test_silent_zero_branches(self, label, build, expected):
        """Hat, Transpose and Vee are linear, so the coefficient of l is the
        operator applied to the remaining factor — never zero. The spec requires
        unsupported patterns to raise a typed error, not return Zero."""
        out = extract_linear_coeff(build(self.a, self.l, self.M), self.l)
        assert str(out) == expected, f"{label}: {out}"


# ---------------------------------------------------------------------------
# 6. Nonlinear and unsupported acceleration patterns
# ---------------------------------------------------------------------------


class TestNonLinearAcceleration:
    def setup_method(self):
        (self.m,) = getScalars("m", attr=["Constant"])
        self.x = Scalar("x")
        self.xdd = self.x.t_diff().t_diff()
        self.variables = SystemVariables(scalars=[self.x])
        self.key = str(Variation(self.x))

    def _sf(self, eqn):
        return to_standard_form({self.key: (Variation(self.x), eqn)}, self.variables, [])

    def test_xdd_squared_raises(self):
        with pytest.raises(NonLinearError):
            self._sf(Mul(self.xdd, self.xdd))

    def test_m_xdd_xdd_raises(self):
        with pytest.raises(NonLinearError):
            self._sf(Mul(Mul(self.m, self.xdd), self.xdd))

    def test_xdd_squared_plus_linear_raises(self):
        """A partial coefficient must never be returned for ẍ² + m x."""
        with pytest.raises(NonLinearError):
            self._sf(Add(Mul(self.xdd, self.xdd), Mul(self.m, self.x)))

    def test_xdd_inside_vector_quadratic_raises(self):
        a = Vector("a")
        with pytest.raises(NonLinearError):
            extract_linear_coeff(SVMul(a, Mul(self.xdd, self.xdd)), self.xdd)

    def test_matrix_query_raises_typed(self):
        """(SCALAR expr, MATRIX query) is outside the dispatch grid: AlgebraicError
        (was NotImplementedError before fix plan 1.4)."""
        with pytest.raises(AlgebraicError, match="no left coefficient"):
            extract_linear_coeff(Mul(self.m, self.x), Matrix("A"))

    def test_jerk_term_is_not_silently_dropped(self):
        """m x⃛ + k x: there is no ẍ, so either raise or keep m x⃛ in f."""
        (k,) = getScalars("k", attr=["Constant"])
        jerk = TimeDerivative(self.xdd)
        eqn = Add(Mul(self.m, jerk), Mul(k, self.x))
        with pytest.raises((NonLinearError, AlgebraicError, NotImplementedError)):
            sf = self._sf(eqn)[self.key]
            assert sf.f.has(jerk), (_strs(sf.M), str(sf.f))

    def test_higher_derivative_lagrangian_raises(self):
        """L = ½ m ẍ² gives δL = m ẍ δẍ; δẍ is not an IBP target, so the term
        survives IBP and is then opaque to extraction. Correct EOM would be
        −m d⁴x/dt⁴ + f = 0 (two IBPs); at minimum a typed error is required."""
        f = Scalar("f")
        L = HALF * self.m * self.xdd * self.xdd
        with pytest.raises((NonLinearError, AlgebraicError, NotImplementedError)):
            eom = compute_eom(L, self.x.delta() * f, self.variables)
            assert str(eom[self.key][1]) != "f", str(eom[self.key][1])


class TestAccelerationPartition:
    """Terms holding both ẋ and ẍ (damping × something) must land in M with
    the ẋ factor kept as the coefficient."""

    def test_c_xdot_xddot_coefficient(self):
        m, c, k = getScalars("m c k", attr=["Constant"])
        x = Scalar("x")
        xd, xdd = x.t_diff(), x.t_diff().t_diff()
        eqn = Add(Mul(Mul(c, xd), xdd), Mul(k, x), Mul(m, xdd))
        key = str(Variation(x))
        sf = to_standard_form({key: (Variation(x), eqn)}, SystemVariables(scalars=[x]), [])
        assert _strs(sf[key].M) == {DDX: "(c*\\frac{d}{dt}(x)+m)"}
        assert str(sf[key].f) == "k*x"

    def test_xdot_term_never_matches_xddot(self):
        """has() recurses into TimeDerivative, so ẍ.has(ẋ) is True; the reverse
        must be False or damping would be classified as inertia."""
        x = Scalar("x")
        xd, xdd = x.t_diff(), x.t_diff().t_diff()
        assert xdd.has(xd)
        assert not xd.has(xdd)
        assert not Mul(Scalar("c"), xd).has(xdd)


# ---------------------------------------------------------------------------
# 7. Names that are substrings of one another
# ---------------------------------------------------------------------------


class TestNameCollisions:
    def test_x_and_x1_coupled_exact(self):
        """L = ½ m ẋ² + ½ m ẋ₁² − k x x₁:
        δx:  −m ẍ − k x₁ + f = 0;   δx₁:  −m ẍ₁ − k x + f₁ = 0."""
        m, k = getScalars("m k", attr=["Constant"])
        x, x1, f, f1 = Scalar("x"), Scalar("x1"), Scalar("f"), Scalar("f1")
        L = HALF * m * x.t_diff() * x.t_diff() + HALF * m * x1.t_diff() * x1.t_diff()
        L = L - k * x * x1
        variables = SystemVariables(scalars=[x, x1])
        eom = compute_eom(L, x.delta() * f + x1.delta() * f1, variables)
        sf = to_standard_form(eom, variables, [f, f1])
        ddx1 = "\\frac{d}{dt}(\\frac{d}{dt}(x1))"
        assert str(eom[str(Variation(x))][1]) == f"(m*{DDX}*(-1)+k*x1*(-1)+f)"
        assert str(eom[str(Variation(x1))][1]) == f"(m*{ddx1}*(-1)+k*x*(-1)+f1)"
        assert _strs(sf[str(Variation(x))].M) == {DDX: "m*(-1)"}
        assert _strs(sf[str(Variation(x1))].M) == {ddx1: "m*(-1)"}
        assert _strs(sf[str(Variation(x))].G) == {"f": "1"}
        assert _strs(sf[str(Variation(x1))].G) == {"f1": "1"}

    def test_has_is_whole_string_and_typed(self):
        k, x = Scalar("k", attr=["Constant"]), Scalar("x")
        assert not Mul(k, x).has(Scalar("kx"))
        assert not Scalar("x1").has(x)
        assert not Scalar("l0").has(Scalar("l"))

    def test_product_string_does_not_alias_a_leaf(self):
        """k·x contains no symbol named kx, so the coefficient of kx is 0."""
        k, x = Scalar("k", attr=["Constant"]), Scalar("x")
        assert str(extract_linear_coeff(Mul(k, x), Scalar("kx"))) == "0"

    def test_x_times_one_does_not_alias_x1(self):
        x = Scalar("x")
        assert str(extract_linear_coeff(Mul(x, One), Scalar("x1"))) == "0"


# ---------------------------------------------------------------------------
# 8. Input validation and standard-form hygiene
# ---------------------------------------------------------------------------


class TestInputValidation:
    def test_constant_scalar_in_scalars_raises(self):
        m, k = getScalars("m k", attr=["Constant"])
        x, f = Scalar("x"), Scalar("f")
        xc = Scalar("xc", attr=["Constant"])
        L = HALF * m * x.t_diff() * x.t_diff() - HALF * k * x * x
        with pytest.raises((ValueError, TypeError, NotImplementedError)):
            eom = compute_eom(L, x.delta() * f, SystemVariables(scalars=[x, xc]))
            assert "0" not in eom, list(eom)


class TestStandardFormHygiene:
    """Sliding mass on a rotating body, SO(3) × R, x = R ρ l (‖ρ‖ = 1):

    ẋ = R (Ω × ρ) l + R ρ l̇,  ẋ·ẋ = l² ‖Ω × ρ‖² + l̇²  ((Ω×ρ)·ρ = 0)
    L = ½ m (l̇² + l² ‖Ω×ρ‖²) + ½ Ω·JΩ
    δl:  −m l̈ + m l ‖Ω×ρ‖² + f_l = 0   — no Ω̇ in this equation.
    """

    def setup_method(self):
        (self.m,) = getScalars("m", attr=["Constant"])
        J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
        rho = Vector("rho", attr=["Constant", "UnitNorm"])
        R = SO3("R")
        Om, eta = R.get_tangent_vector(), R.get_variation_vector()
        self.l = Scalar("l")
        f_l, Mt = Scalar("f_l"), Vector("M")
        v = (R * (rho * self.l)).t_diff()
        L = HALF * self.m * Dot(v, v) + HALF * Dot(Om, J * Om)
        dW = self.l.delta() * f_l + Dot(eta, Mt)
        variables = SystemVariables(scalars=[self.l], matrices=[R])
        self.eom = compute_eom(L, dW, variables)
        self.sf = to_standard_form(self.eom, variables, [f_l, Mt])
        self.eq_l = self.sf[str(Variation(self.l))]
        self.eta = eta

    def test_l_equation_mass_and_input(self):
        assert str(self.eq_l.M[DDL]) == "m*Dot(R*rho,R*rho)*(-1)"
        assert _strs(self.eq_l.G) == {"f_l": "1"}

    @pytest.mark.xfail(
        strict=True,
        reason="QUALITY: M keeps an identically-zero coefficient (M[Ω̇] = 0v) because"
        " ρ·(Ω̇×ρ) is only recognised as zero after extraction",
    )
    def test_zero_coefficients_are_pruned_from_M(self):
        assert set(self.eq_l.M) == {DDL}, _strs(self.eq_l.M)

    @pytest.mark.xfail(
        strict=True,
        reason="F6 (deferred): the IBP term J̇_eff Ω is present since F1, but the two"
        " cancelling non-IBP l̇ terms are not combined, so f has three l̇ addends.",
    )
    def test_eta_equation_has_inertia_rate_term(self):
        """J_eff = J − m l² ρ̂ρ̂ is time-varying: Euler's equation is
        J_eff Ω̇ + J̇_eff Ω + Ω × J_eff Ω = M with J̇_eff Ω = −2 m l l̇ ρ×(ρ×Ω).
        The code's η equation has two l̇ terms that cancel, so the l̇ content of
        f must be a single non-zero term after combining like terms."""
        f = self.sf[str(self.eta)].f
        l_dot_terms = [n for n in f.nodes if n.has(TimeDerivative(self.l))]
        assert len(l_dot_terms) == 1, [str(n) for n in l_dot_terms]


# ---------------------------------------------------------------------------
# 9. Fix coverage (added with the fixes for F1–F5; see the review's "Fix log")
# ---------------------------------------------------------------------------


class TestFixF1PipelineExact:
    """Exact equations that the F1 fix (IBP under a non-constant prefactor) unlocks."""

    def test_elastic_pendulum_xi_equation_exact(self):
        """ξ:  −m l² ω̇ − 2 m l l̇ ω − m g l (q×e₃) + f_q = 0  (hand derivation in
        TestElasticSphericalPendulum). The Coriolis term is a single addend
        with coefficient −2 because d/dt(m l l) is formed as 2 m l l̇."""
        eom, sf, l, q = _elastic_pendulum()
        xi = q.get_variation_vector()
        _, eqn = eom[str(xi)]
        assert str(eqn) == (
            f"({DW}*m*l*l*(-1)+\\omega_{{q}}*m*l*{DL}*(-2)+Cross(q,e3)*m*g*l*(-1)+f_q)"
        )
        f = sf[str(xi)].f
        assert str(f) == f"(\\omega_{{q}}*m*l*{DL}*(-2)+Cross(q,e3)*m*g*l*(-1))"

    def test_prescribed_length_s2_pendulum_exact(self):
        """m l² ω̇ + 2 m l l̇ ω + m g l q×e₃ = f with l(t) prescribed (no scalars)."""
        m, g = getScalars("m g", attr=["Constant"])
        l = Scalar("l")
        e3 = Vector("e3", attr=["Constant"])
        q, f = S2("q"), Vector("f")
        omega, xi = q.get_tangent_vector(), q.get_variation_vector()
        L = HALF * m * l * l * Dot(omega, omega) - m * g * l * Dot(q, e3)
        eom = compute_eom(L, Dot(xi, f), SystemVariables(vectors=[q]), params=[l])
        _, eqn = eom[str(xi)]
        assert (
            str(eqn) == f"({DW}*m*l*l*(-1)+\\omega_{{q}}*m*l*{DL}*(-2)+Cross(q,e3)*m*g*l*(-1)+f)"
        )

    def test_r3_x_r_cross_inertia_exact(self):
        """δx:  −m ẍ − m l̈ e₁ + F = 0, so M_x[l̈] = −m e₁ = M_l[ẍ] (same column vector,
        the library stores the row −m e₁ᵀ of the δl equation as that vector)."""
        eom, sf, x, l = _r3_x_r_coupled()
        kx, kl = str(Variation(x)), str(Variation(l))
        _, eqn = eom[kx]
        assert str(eqn) == f"({DDX}*m*(-1)+e1*m*{DDL}*(-1)+F)"
        assert _strs(sf[kx].M) == {DDL: "e1*m*(-1)", DDX: "I*m*(-1)"}
        assert str(sf[kx].M[DDL]) == str(sf[kl].M[DDX])


class TestCollectTwoRotations:
    """Double rigid pendulum on SO(3)×SO(3) (test_manifold_pipelines._double_rigid_pendulum):
    the coupling KE_c = m₂ (R₁(Ω₁×l₁))·(R₂(Ω₂×ρ₂)) used to survive IBP as
    Dot(R₁(η̇₁×l₁), R₂(Ω₂×ρ₂)) — collect could not rotate η̇₁ out of a Dot of two
    different rotations — and extraction then dropped it silently. With the F3
    rule such a term raises, so collect now peels R₁ off: Dot(R₁ f, w) = Dot(f, R₁ᵀ w).

    Hand derivation: KE_c = m₂ Ω₁·(l₁ × R₁ᵀR₂(Ω₂×ρ₂)) = −m₂ Ω₁·(l̂₁ R₁ᵀ R₂ ρ̂₂ Ω₂), so
    ∂KE_c/∂Ω₁ = −m₂ l̂₁ R₁ᵀ R₂ ρ̂₂ Ω₂ and the η₁ equation (−d/dt ∂L/∂Ω₁ + …) has
    M₁[Ω̇₂] = m₂ l̂₁ R₁ᵀ R₂ ρ̂₂; symmetrically M₂[Ω̇₁] = m₂ ρ̂₂ R₂ᵀ R₁ l̂₁ = M₁[Ω̇₂]ᵀ.
    """

    def test_cross_inertia_blocks_exact(self):
        J1 = Matrix("J1", attr=["Constant", "SymmetricMatrix"])
        J2 = Matrix("J2", attr=["Constant", "SymmetricMatrix"])
        rho1 = Vector("\\rho_1", attr=["Constant"])
        rho2 = Vector("\\rho_2", attr=["Constant"])
        l1 = Vector("l_1", attr=["Constant"])
        m1, m2, g = getScalars("m1 m2 g", attr=["Constant"])
        e3 = Vector("e3", attr=["Constant"])
        M1, M2 = Vector("M1"), Vector("M2")
        R1, R2 = SO3("R1"), SO3("R2")
        Om1, Om2 = R1.get_tangent_vector(), R2.get_tangent_vector()
        eta1, eta2 = R1.get_variation_vector(), R2.get_variation_vector()
        x1 = R1 * rho1
        x2 = R1 * l1 + R2 * rho2
        v1, v2 = x1.t_diff(), x2.t_diff()
        KE = (
            Dot(Om1, J1 * Om1) * HALF
            + Dot(Om2, J2 * Om2) * HALF
            + Dot(v1, v1) * m1 * HALF
            + Dot(v2, v2) * m2 * HALF
        )
        PE = m1 * g * Dot(x1, e3) + m2 * g * Dot(x2, e3)
        variables = SystemVariables(matrices=[R1, R2])
        eom = compute_eom(KE - PE, Dot(eta1, M1) + Dot(eta2, M2), variables)
        sf = to_standard_form(eom, variables, [M1, M2])
        dOm1, dOm2 = str(TimeDerivative(Om1)), str(TimeDerivative(Om2))
        assert str(sf[str(eta1)].M[dOm2]) == "Hat(l_1)*(R1)'*R2*Hat(\\rho_2)*m2"
        assert str(sf[str(eta2)].M[dOm1]) == "Hat(\\rho_2)*(R2)'*R1*Hat(l_1)*m2"
        for key, (_, eqn) in eom.items():
            assert not eqn.has(TimeDerivative(eta1)) and not eqn.has(TimeDerivative(eta2)), key


class TestFixF2LinearUnaryOps:
    """Linearity of the remaining unary/outer-product nodes for a scalar query."""

    def setup_method(self):
        self.l = Scalar("l")
        self.a, self.b, self.w = Vector("a"), Vector("b"), Vector("w")
        self.M = Matrix("M")

    def test_row_vector_times_matrix(self):
        """(a l)ᵀ M = l (aᵀ M): coefficient is aᵀ M (was an ExpressionMismatchError)."""
        out = extract_linear_coeff(MVMul(Transpose(SVMul(self.a, self.l)), self.M), self.l)
        assert str(out) == "(a)'*M"

    def test_outer_product_both_sides(self):
        """(a l) bᵀ = l (a bᵀ) and a (b l)ᵀ = l (a bᵀ)."""
        out1 = extract_linear_coeff(VVMul(SVMul(self.a, self.l), Transpose(self.b)), self.l)
        out2 = extract_linear_coeff(VVMul(self.a, Transpose(SVMul(self.b, self.l))), self.l)
        assert str(out1) == "a*(b)'"
        assert str(out2) == "a*(b)'"

    def test_matrix_transpose(self):
        """(M l)ᵀ = l Mᵀ."""
        assert str(extract_linear_coeff(Transpose(SMMul(self.M, self.l)), self.l)) == "(M)'"

    def test_hat_times_vector_quadratic_raises(self):
        """Hat(a l)(w l) is quadratic in l."""
        with pytest.raises(NonLinearError):
            extract_linear_coeff(MVMul(Hat(SVMul(self.a, self.l)), SVMul(self.w, self.l)), self.l)


class TestFixF3HigherDerivativeRule:
    """A query hidden under a calculus operator that is not the query itself raises
    AlgebraicError naming the higher-order derivative, in extraction and in the
    standard-form partition alike."""

    def setup_method(self):
        (self.m,) = getScalars("m", attr=["Constant"])
        (self.k,) = getScalars("k", attr=["Constant"])
        self.x = Scalar("x")
        self.xdd = self.x.t_diff().t_diff()
        self.jerk_eqn = Add(Mul(self.m, TimeDerivative(self.xdd)), Mul(self.k, self.x))

    def test_extract_raises_with_message(self):
        with pytest.raises(AlgebraicError, match="higher order"):
            extract_linear_coeff(self.jerk_eqn, self.xdd)

    def test_standard_form_partition_raises(self):
        key = str(Variation(self.x))
        with pytest.raises(AlgebraicError, match="higher order"):
            to_standard_form(
                {key: (Variation(self.x), self.jerk_eqn)}, SystemVariables(scalars=[self.x]), []
            )

    def test_lower_derivative_is_still_opaque_not_an_error(self):
        """ẋ contains no ẍ: c ẋ + m ẍ partitions normally (M = m, f = c ẋ)."""
        (c,) = getScalars("c", attr=["Constant"])
        key = str(Variation(self.x))
        eqn = Add(Mul(c, self.x.t_diff()), Mul(self.m, self.xdd))
        sf = to_standard_form(
            {key: (Variation(self.x), eqn)}, SystemVariables(scalars=[self.x]), []
        )
        assert _strs(sf[key].M) == {DDX: "m"}
        assert str(sf[key].f) == "c*\\frac{d}{dt}(x)"

    def test_variation_under_derivative_raises(self):
        """δẍ = d/dt(d/dt(δx)) survives IBP (target is d/dt(δx)); extraction must raise."""
        f = Scalar("f")
        L = HALF * self.m * self.xdd * self.xdd
        with pytest.raises(AlgebraicError, match="higher order"):
            compute_eom(L, self.x.delta() * f, SystemVariables(scalars=[self.x]))


class TestFixF4TypedIdentity:
    def test_leaf_named_like_a_product_is_found_only_as_a_leaf(self):
        """k x + kx: the coefficient of the leaf kx is 1, from the leaf term only."""
        k, x = Scalar("k", attr=["Constant"]), Scalar("x")
        assert str(extract_linear_coeff(Add(Mul(k, x), Scalar("kx")), Scalar("kx"))) == "1"

    def test_scalar_and_vector_with_same_name_do_not_alias(self):
        """A vector named x is not the scalar x."""
        x = Scalar("x")
        assert str(extract_linear_coeff(Dot(Vector("x"), Vector("w")), x)) == "0"


class TestFixF5Validation:
    def setup_method(self):
        (self.m,) = getScalars("m", attr=["Constant"])
        self.x = Vector("x")
        self.L = HALF * self.m * Dot(self.x.t_diff(), self.x.t_diff())
        self.dW = Dot(self.x.delta(), Vector("F"))

    def test_constant_vector_raises(self):
        xc = Vector("xc", attr=["Constant"])
        with pytest.raises(ValueError, match="constant"):
            compute_eom(self.L, self.dW, SystemVariables(vectors=[self.x, xc]))

    def test_constant_matrix_raises(self):
        Jc = Matrix("Jc", attr=["Constant"])
        with pytest.raises(ValueError, match="constant"):
            compute_eom(self.L, self.dW, SystemVariables(vectors=[self.x], matrices=[Jc]))

    def test_wrong_type_raises(self):
        with pytest.raises(TypeError):
            compute_eom(self.L, self.dW, SystemVariables(scalars=[self.x]))
        with pytest.raises(TypeError):
            compute_eom(self.L, self.dW, SystemVariables(vectors=[Scalar("s")]))
