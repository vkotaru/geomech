"""Integration by parts of a vector target under a scalar prefactor (review finding F1).

Identity (boundary terms dropped), for a scalar a(t) and vectors x(t), r(t):

    ∫ a Dot(ẋ, r) dt = −∫ Dot(x, d/dt(a r)) dt = −∫ [ȧ Dot(x, r) + a Dot(x, ṙ)] dt

Before the fix, ``_apply_ibp`` differentiated only ``r`` and treated any scalar
prefactor hoisted outside the Dot as a constant, dropping ``ȧ Dot(x, r)``.
Every expected string below was derived by hand from this identity; the
``a Dot(x, ṙ)`` term keeps the pre-fix tree shape so constant-prefactor outputs
are unchanged.
"""

from geomech import Add, Dot, Mul, Scalar, SVMul, TimeDerivative, Vector
from geomech.core.math.ibp import _apply_ibp, _t_diff_product, integrate_by_parts

DA = "\\frac{d}{dt}(a)"
DB = "\\frac{d}{dt}(b)"
DL = "\\frac{d}{dt}(l)"
DOM = "\\frac{d}{dt}(omega)"


def _setup():
    a, b, l = Scalar("a"), Scalar("b"), Scalar("l")
    c = Scalar("c", attr=["Constant"])
    xi, omega = Vector("xi"), Vector("omega")
    target = TimeDerivative(xi)
    return a, b, l, c, xi, omega, target


class TestConstantPrefactorUnchanged:
    """A constant prefactor c has ċ = 0, so the result is the pre-fix one."""

    def test_full_pipeline(self):
        """∫ c Dot(ξ̇, ω) dt = −∫ c Dot(ξ, ω̇) dt."""
        _, _, _, c, _, omega, target = _setup()
        out = integrate_by_parts(Mul(c, Dot(target, omega)), [target])
        assert str(out) == f"c*Dot(xi,{DOM})*(-1)"

    def test_single_pass_keeps_tree_shape(self):
        """_apply_ibp replaces the Dot in place: Mul(Mul(c, Dot), b) → Mul(Mul(c, Dot'), b)
        plus the ḃ term for the non-constant factor b."""
        _, b, _, c, _, omega, target = _setup()
        out = _apply_ibp(Mul(Mul(c, Dot(target, omega)), b), target)
        assert isinstance(out, Add) and len(out.nodes) == 2
        first = out.nodes[0]
        assert isinstance(first, Mul) and isinstance(first.left, Mul)
        assert str(first) == f"c*Dot(xi*(-1),{DOM})*b"
        assert str(out.nodes[1]) == f"Dot(xi*(-1),omega*c*{DB})"


class TestNonConstantPrefactor:
    def test_single_factor(self):
        """∫ a Dot(ξ̇, ω) dt = −∫ [a Dot(ξ, ω̇) + ȧ Dot(ξ, ω)] dt."""
        a, _, _, _, _, omega, target = _setup()
        out = integrate_by_parts(Mul(a, Dot(target, omega)), [target])
        assert str(out) == f"(a*Dot(xi,{DOM})*(-1)+{DA}*Dot(xi,omega)*(-1))"

    def test_prefactor_inside_dot_is_pulled_first(self):
        """Dot(a ξ̇, ω) is the same integrand: pull hoists a, then the identity applies."""
        a, _, _, _, _, omega, target = _setup()
        out = integrate_by_parts(Dot(SVMul(target, a), omega), [target])
        assert str(out) == f"(a*Dot(xi,{DOM})*(-1)+{DA}*Dot(xi,omega)*(-1))"

    def test_two_factors_product_rule(self):
        """d/dt(a b) = ȧ b + a ḃ:
        ∫ a b Dot(ξ̇, ω) dt = −∫ [a b Dot(ξ, ω̇) + ȧ b Dot(ξ, ω) + a ḃ Dot(ξ, ω)] dt."""
        a, b, _, _, _, omega, target = _setup()
        out = integrate_by_parts(Mul(Mul(a, b), Dot(target, omega)), [target])
        assert str(out) == (
            f"(a*b*Dot(xi,{DOM})*(-1)+{DA}*b*Dot(xi,omega)*(-1)+a*{DB}*Dot(xi,omega)*(-1))"
        )

    def test_repeated_factor_is_one_term(self):
        """d/dt(l²) = 2 l l̇:  ∫ l² Dot(ξ̇, ω) dt = −∫ [l² Dot(ξ, ω̇) + 2 l l̇ Dot(ξ, ω)] dt.
        The Coriolis-type term must be a single addend with coefficient 2, not
        two uncombined copies of l l̇."""
        _, _, l, _, _, omega, target = _setup()
        out = integrate_by_parts(Mul(Mul(l, l), Dot(target, omega)), [target])
        assert str(out) == f"(l*l*Dot(xi,{DOM})*(-1)+l*{DL}*Dot(xi,omega)*(-2))"

    def test_constant_rhs_keeps_only_prefactor_rate(self):
        """ė = 0 for a constant vector e:  ∫ a Dot(ξ̇, e) dt = −∫ ȧ Dot(ξ, e) dt."""
        a, _, _, _, _, _, target = _setup()
        e = Vector("e", attr=["Constant"])
        out = integrate_by_parts(Mul(a, Dot(target, e)), [target])
        assert str(out) == f"{DA}*Dot(xi,e)*(-1)"

    def test_dot_in_middle_of_chain_not_double_counted(self):
        """Mul(Mul(a, Dot(ξ̇, ω)), b) has prefactor a b; the ȧ b term must appear once:
        −[a b Dot(ξ, ω̇) + ȧ b Dot(ξ, ω) + a ḃ Dot(ξ, ω)]."""
        a, b, _, _, _, omega, target = _setup()
        out = integrate_by_parts(Mul(Mul(a, Dot(target, omega)), b), [target])
        assert str(out) == (
            f"(a*Dot(xi,{DOM})*b*(-1)+{DA}*b*Dot(xi,omega)*(-1)+a*{DB}*Dot(xi,omega)*(-1))"
        )

    def test_mixed_sum(self):
        """Terms are treated independently: a constant prefactor on one Dot and a
        non-constant one on another."""
        a, _, _, c, _, omega, target = _setup()
        w = Vector("w")
        out = integrate_by_parts(Add(Mul(c, Dot(target, omega)), Mul(a, Dot(target, w))), [target])
        assert str(out) == (
            f"(c*Dot(xi,{DOM})*(-1)+a*Dot(xi,\\frac{{d}}{{dt}}(w))*(-1)+{DA}*Dot(xi,w)*(-1))"
        )


class TestProductDerivative:
    """_t_diff_product: d/dt(∏ fᵢ^{nᵢ}) = Σᵢ nᵢ fᵢ^{nᵢ−1} ḟᵢ ∏_{j≠i} fⱼ^{nⱼ}."""

    def test_repeated_factor(self):
        """d/dt(c l l) = 2 c l l̇ (c constant)."""
        _, _, l, c, _, _, _ = _setup()
        assert str(_t_diff_product(Mul(Mul(c, l), l))) == f"c*l*{DL}*(2)"

    def test_distinct_factors(self):
        """d/dt(a b) = ȧ b + a ḃ."""
        a, b, _, _, _, _, _ = _setup()
        assert str(_t_diff_product(Mul(a, b))) == f"({DA}*b+a*{DB})"

    def test_all_constant_is_zero(self):
        """d/dt(c · 2) = 0."""
        _, _, _, c, _, _, _ = _setup()
        two = Scalar("2", value=2, attr=["Constant"])
        assert str(_t_diff_product(Mul(c, two))) == "0"
