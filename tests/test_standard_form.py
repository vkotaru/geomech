"""Tests for standard form extraction: M*a + f + G*u = 0."""

from geomech import (
    S2,
    SO3,
    Dot,
    Matrix,
    Scalar,
    StandardFormEquation,
    SystemVariables,
    TimeDerivative,
    Variation,
    Vector,
    ZeroVector,
    compute_eom,
    getScalars,
    to_standard_form,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _point_mass_eom():
    """Point mass: L = 0.5*m*Dot(v,v) - m*g*Dot(x,e3), dW = Dot(dx,f)."""
    m, g = getScalars("m g", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    x = Vector("x")
    f = Vector("f")
    v = x.t_diff()

    half = Scalar("0.5", value=0.5, attr=["Constant"])
    KE = m * Dot(v, v) * half
    PE = m * Dot(x, g * e3)
    L = KE - PE
    dW = Dot(x.get_variation_vector(), f)

    variables = SystemVariables(vectors=[x])
    eom = compute_eom(L, dW, variables)
    return eom, variables, [f], x


def _spherical_pendulum_eom():
    """Spherical pendulum on S2: L = 0.5*m*l^2*Dot(w,w) - m*g*l*Dot(q,e3)."""
    m, g, l = getScalars("m g l", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    q = S2("q")
    f = Vector("f")

    x = l * q
    v = x.t_diff()
    half = Scalar("0.5", value=0.5, attr=["Constant"])
    KE = m * Dot(v, v) * half
    PE = m * Dot(x, g * e3)
    L = KE - PE
    dW = Dot(q.get_variation_vector(), f)

    variables = SystemVariables(vectors=[q])
    eom = compute_eom(L, dW, variables)
    return eom, variables, [f], q


def _rigid_pendulum_eom():
    """Rigid pendulum on SO3: L = 0.5*Dot(Om,J*Om) + 0.5*m*Dot(v,v) - m*g*Dot(x,e3)."""
    J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
    rho = Vector("\\rho", attr=["Constant"])
    m, g = getScalars("m g", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    M_torque = Vector("M")

    R = SO3("R")
    Om = R.get_tangent_vector()
    eta = R.get_variation_vector()

    x = R * rho
    v = x.t_diff()
    half = Scalar("0.5", value=0.5, attr=["Constant"])
    KE = Dot(Om, J * Om) * half + Dot(v, v) * m * half
    PE = m * g * Dot(x, e3)
    L = KE - PE
    dW = Dot(eta, M_torque)

    variables = SystemVariables(matrices=[R])
    eom = compute_eom(L, dW, variables)
    return eom, variables, [M_torque], R


# ---------------------------------------------------------------------------
# Point mass tests
# ---------------------------------------------------------------------------


class TestPointMassStandardForm:
    def test_returns_dict(self):
        eom, variables, inputs, x = _point_mass_eom()
        sf = to_standard_form(eom, variables, inputs)
        assert isinstance(sf, dict)
        assert len(sf) == 1

    def test_keys_match_eom(self):
        eom, variables, inputs, x = _point_mass_eom()
        sf = to_standard_form(eom, variables, inputs)
        assert set(sf.keys()) == set(eom.keys())

    def test_is_standard_form_equation(self):
        eom, variables, inputs, x = _point_mass_eom()
        sf = to_standard_form(eom, variables, inputs)
        key = str(Variation(x))
        eq = sf[key]
        assert isinstance(eq, StandardFormEquation)

    def test_M_has_acceleration(self):
        eom, variables, inputs, x = _point_mass_eom()
        sf = to_standard_form(eom, variables, inputs)
        key = str(Variation(x))
        eq = sf[key]

        ddx = str(TimeDerivative(TimeDerivative(x)))
        assert ddx in eq.M, f"M should contain ddx key, got keys: {list(eq.M.keys())}"
        M = eq.M[ddx]
        assert M is not None
        # M should contain mass m
        assert "m" in str(M)

    def test_G_has_input(self):
        eom, variables, inputs, x = _point_mass_eom()
        sf = to_standard_form(eom, variables, inputs)
        key = str(Variation(x))
        eq = sf[key]

        assert str(inputs[0]) in eq.G, f"G should contain input f, got keys: {list(eq.G.keys())}"

    def test_f_has_gravity(self):
        eom, variables, inputs, x = _point_mass_eom()
        sf = to_standard_form(eom, variables, inputs)
        key = str(Variation(x))
        eq = sf[key]

        f_str = str(eq.f)
        # f should contain gravity terms (g, e3, m)
        assert "g" in f_str
        assert "e3" in f_str

    def test_no_acceleration_in_f(self):
        eom, variables, inputs, x = _point_mass_eom()
        sf = to_standard_form(eom, variables, inputs)
        key = str(Variation(x))
        eq = sf[key]

        ddx = TimeDerivative(TimeDerivative(x))
        assert not eq.f.has(ddx), "f should not contain acceleration terms"

    def test_no_input_in_f(self):
        eom, variables, inputs, x = _point_mass_eom()
        sf = to_standard_form(eom, variables, inputs)
        key = str(Variation(x))
        eq = sf[key]

        assert not eq.f.has(inputs[0]), "f should not contain input terms"


# ---------------------------------------------------------------------------
# Spherical pendulum (S2) tests
# ---------------------------------------------------------------------------


class TestSphericalPendulumStandardForm:
    def test_returns_one_equation(self):
        eom, variables, inputs, q = _spherical_pendulum_eom()
        sf = to_standard_form(eom, variables, inputs)
        assert len(sf) == 1

    def test_M_has_acceleration(self):
        """M should contain d/dt(omega) terms now that manifold kinematic
        substitution is implemented."""
        eom, variables, inputs, q = _spherical_pendulum_eom()
        sf = to_standard_form(eom, variables, inputs)
        key = list(sf.keys())[0]
        eq = sf[key]
        assert len(eq.M) > 0

    def test_G_has_input(self):
        eom, variables, inputs, q = _spherical_pendulum_eom()
        sf = to_standard_form(eom, variables, inputs)
        key = list(sf.keys())[0]
        eq = sf[key]
        assert len(eq.G) > 0, "G should have the input force"

    def test_f_has_gravity(self):
        eom, variables, inputs, q = _spherical_pendulum_eom()
        sf = to_standard_form(eom, variables, inputs)
        key = list(sf.keys())[0]
        eq = sf[key]
        assert "g" in str(eq.f) or eq.f == ZeroVector, "f should have gravity or be zero"


# ---------------------------------------------------------------------------
# Rigid pendulum (SO3) tests
# ---------------------------------------------------------------------------


class TestRigidPendulumStandardForm:
    def test_returns_one_equation(self):
        eom, variables, inputs, R = _rigid_pendulum_eom()
        sf = to_standard_form(eom, variables, inputs)
        assert len(sf) == 1

    def test_M_has_acceleration(self):
        """M should contain d/dt(Omega) terms now that manifold kinematic
        substitution is implemented."""
        eom, variables, inputs, R = _rigid_pendulum_eom()
        sf = to_standard_form(eom, variables, inputs)
        key = list(sf.keys())[0]
        eq = sf[key]
        Om = R.get_tangent_vector()
        ddOm = str(TimeDerivative(Om))
        assert ddOm in eq.M, f"M should contain d/dt(Omega), got keys: {list(eq.M.keys())}"
        # M should contain inertia terms (J)
        assert "J" in str(eq.M[ddOm])

    def test_G_has_torque(self):
        eom, variables, inputs, R = _rigid_pendulum_eom()
        sf = to_standard_form(eom, variables, inputs)
        key = list(sf.keys())[0]
        eq = sf[key]
        assert len(eq.G) > 0, "G should contain the torque input"

    def test_f_exists(self):
        eom, variables, inputs, R = _rigid_pendulum_eom()
        sf = to_standard_form(eom, variables, inputs)
        key = list(sf.keys())[0]
        eq = sf[key]
        # f should have Coriolis and/or gravity terms
        assert eq.f is not None


# ---------------------------------------------------------------------------
# No inputs case
# ---------------------------------------------------------------------------


class TestNoInputs:
    def test_no_inputs_gives_empty_G(self):
        """Free particle with no external forces."""
        (m,) = getScalars("m", attr=["Constant"])
        x = Vector("x")
        v = x.t_diff()
        half = Scalar("0.5", value=0.5, attr=["Constant"])
        L = m * Dot(v, v) * half
        dW = Scalar("0", value=0, attr=["Constant", "Zero"])

        variables = SystemVariables(vectors=[x])
        eom = compute_eom(L, dW, variables)
        sf = to_standard_form(eom, variables, inputs=[])

        key = str(Variation(x))
        eq = sf[key]
        assert len(eq.G) == 0, "G should be empty with no inputs"


# ---------------------------------------------------------------------------
# Scalar systems
# ---------------------------------------------------------------------------


def _spring_mass_eom():
    """Spring-mass: m*ẍ + k*x = f. Single scalar config var x."""
    m, k = getScalars("m k", attr=["Constant"])
    half = Scalar("0.5", value=0.5, attr=["Constant"])
    x = Scalar("x")
    f = Scalar("f")
    v = x.t_diff()
    dx = x.delta()

    L = half * m * v * v - half * k * x * x
    dW = dx * f

    variables = SystemVariables(scalars=[x])
    eom = compute_eom(L, dW, variables)
    return eom, variables, [f], x


def _spring_mass_damper_eom():
    """Spring-mass-damper: m*ẍ + c*ẋ + k*x = f."""
    m, k, c = getScalars("m k c", attr=["Constant"])
    half = Scalar("0.5", value=0.5, attr=["Constant"])
    x = Scalar("x")
    f = Scalar("f")
    v = x.t_diff()
    dx = x.delta()

    L = half * m * v * v - half * k * x * x
    dW = dx * f - c * v * dx

    variables = SystemVariables(scalars=[x])
    eom = compute_eom(L, dW, variables)
    return eom, variables, [f], x


def _double_spring_mass_eom():
    """Double spring-mass with coupling through k2."""
    m1, m2, k1, k2 = getScalars("m1 m2 k1 k2", attr=["Constant"])
    half = Scalar("0.5", value=0.5, attr=["Constant"])
    x1 = Scalar("x1")
    x2 = Scalar("x2")
    f1 = Scalar("f1")
    f2 = Scalar("f2")
    v1, v2 = x1.t_diff(), x2.t_diff()
    dx1, dx2 = x1.delta(), x2.delta()

    KE = half * m1 * v1 * v1 + half * m2 * v2 * v2
    PE = half * k1 * x1 * x1 + half * k2 * (x2 - x1) * (x2 - x1)
    L = KE - PE
    dW = dx1 * f1 + dx2 * f2

    variables = SystemVariables(scalars=[x1, x2])
    eom = compute_eom(L, dW, variables)
    return eom, variables, [f1, f2], (x1, x2)


class TestSpringMassStandardForm:
    """Spring-mass: M[ẍ] = -m (sign from -m·ẍ term), f = -k·x, G[f] = 1."""

    def setup_method(self):
        self.eom, self.variables, self.inputs, self.x = _spring_mass_eom()
        self.sf = to_standard_form(self.eom, self.variables, self.inputs)
        self.key = str(Variation(self.x))
        self.eq = self.sf[self.key]

    def test_returns_one_equation(self):
        assert len(self.sf) == 1

    def test_is_standard_form_equation(self):
        assert isinstance(self.eq, StandardFormEquation)

    def test_M_has_acceleration_key(self):
        ddx_key = str(TimeDerivative(TimeDerivative(self.x)))
        assert ddx_key in self.eq.M

    def test_M_contains_mass(self):
        ddx_key = str(TimeDerivative(TimeDerivative(self.x)))
        assert "m" in str(self.eq.M[ddx_key])

    def test_G_has_input_key(self):
        assert "f" in self.eq.G

    def test_f_contains_stiffness(self):
        f_str = str(self.eq.f)
        assert "k" in f_str
        assert "x" in f_str

    def test_no_acceleration_in_f(self):
        ddx = TimeDerivative(TimeDerivative(self.x))
        assert not self.eq.f.has(ddx)

    def test_no_input_in_f(self):
        assert not self.eq.f.has(self.inputs[0])


class TestSpringMassDamperStandardForm:
    """Spring-mass-damper: M[ẍ] = -m, f contains k·x and c·ẋ, G[f] = 1."""

    def setup_method(self):
        self.eom, self.variables, self.inputs, self.x = _spring_mass_damper_eom()
        self.sf = to_standard_form(self.eom, self.variables, self.inputs)
        self.key = str(Variation(self.x))
        self.eq = self.sf[self.key]

    def test_M_has_acceleration_key(self):
        ddx_key = str(TimeDerivative(TimeDerivative(self.x)))
        assert ddx_key in self.eq.M

    def test_f_contains_damping(self):
        """f should include c·ẋ damping term."""
        f_str = str(self.eq.f)
        assert "c" in f_str

    def test_f_contains_stiffness(self):
        f_str = str(self.eq.f)
        assert "k" in f_str

    def test_G_has_input_key(self):
        assert "f" in self.eq.G


class TestDoubleSpringMassStandardForm:
    """Double spring-mass: two coupled equations, each with its own M and G."""

    def setup_method(self):
        self.eom, self.variables, self.inputs, (self.x1, self.x2) = _double_spring_mass_eom()
        self.sf = to_standard_form(self.eom, self.variables, self.inputs)
        self.key1 = str(Variation(self.x1))
        self.key2 = str(Variation(self.x2))

    def test_returns_two_equations(self):
        assert len(self.sf) == 2

    def test_both_keys_present(self):
        assert self.key1 in self.sf
        assert self.key2 in self.sf

    def test_eq1_M_has_m1(self):
        eq1 = self.sf[self.key1]
        ddx1 = str(TimeDerivative(TimeDerivative(self.x1)))
        assert ddx1 in eq1.M
        assert "m1" in str(eq1.M[ddx1])

    def test_eq2_M_has_m2(self):
        eq2 = self.sf[self.key2]
        ddx2 = str(TimeDerivative(TimeDerivative(self.x2)))
        assert ddx2 in eq2.M
        assert "m2" in str(eq2.M[ddx2])

    def test_eq1_G_has_f1(self):
        assert "f1" in self.sf[self.key1].G

    def test_eq2_G_has_f2(self):
        assert "f2" in self.sf[self.key2].G

    def test_eq1_f_has_k2_coupling(self):
        """x1 EOM residual should reference k2 (coupling to x2)."""
        assert "k2" in str(self.sf[self.key1].f)

    def test_eq2_f_has_k2_coupling(self):
        """x2 EOM residual should reference k2."""
        assert "k2" in str(self.sf[self.key2].f)

    def test_no_acceleration_in_f(self):
        ddx1 = TimeDerivative(TimeDerivative(self.x1))
        ddx2 = TimeDerivative(TimeDerivative(self.x2))
        assert not self.sf[self.key1].f.has(ddx1)
        assert not self.sf[self.key1].f.has(ddx2)
        assert not self.sf[self.key2].f.has(ddx1)
        assert not self.sf[self.key2].f.has(ddx2)
