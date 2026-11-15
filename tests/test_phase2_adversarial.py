"""Adversarial review of phase 2 (2.1–2.12): error instead of silence.

Every new guard is attacked two ways: a legitimate system that must still derive (false
positive) and a misuse that must raise (false negative). Confirmed findings are strict
``xfail`` tests whose reason starts with ``F<n>``. Everything else here pins what held.

The numerical check (``check_eom``) is used wherever a derived equation is claimed correct,
so the tests do not depend on the printed form alone.
"""

from __future__ import annotations

import re
import warnings

import pytest

from geomech import (
    S2,
    SO3,
    AlgebraicError,
    Cross,
    Dot,
    Hat,
    Matrix,
    MVMul,
    NonLinearError,
    Number,
    Scalar,
    SystemVariables,
    TimeDerivative,
    Transpose,
    Variation,
    Vector,
    VVMul,
    Zero,
    collect,
    compute_eom,
    extract_linear_coeff,
    full_simplify,
    linearize,
    to_standard_form,
)
from geomech.core.base.expressions import TS2
from geomech.core.base.types import VALID_ATTRS
from geomech.numeric import check_eom
from tests import vbl_systems

HALF = Number(0.5)
TOL = 1e-5


def point_mass():
    """L = ½ m ẋ·ẋ − m g e₃·x, δW = δx·f; returns (L, dW, x, f, m, g, e3)."""
    m, g = Scalar("m", attr=["Constant"]), Scalar("g", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    x, f = Vector("x"), Vector("f")
    L = HALF * m * Dot(x.t_diff(), x.t_diff()) - m * g * Dot(e3, x)
    return L, Dot(x.delta(), f), x, f, m, g, e3


def spherical_pendulum(q=None):
    """L = ½ m l² ω·ω − m g l q·e₃, δW = ξ·f; returns (L, q, ω, ξ, f)."""
    m, g, l = (Scalar(s, attr=["Constant"]) for s in "mgl")
    e3 = Vector("e3", attr=["Constant"])
    q = S2("q") if q is None else q
    f = Vector("f")
    om, xi = q.get_tangent_vector(), q.get_variation_vector()
    L = HALF * m * l * l * Dot(om, om) - m * g * l * Dot(q, e3)
    return L, q, om, xi, f


def rigid_body():
    """L = ½ Ω·JΩ, δW = η·M; returns (L, R, Ω, η, M, J)."""
    J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
    R, M = SO3("R"), Vector("M")
    Om, eta = R.get_tangent_vector(), R.get_variation_vector()
    return HALF * Dot(Om, J * Om), R, Om, eta, M, J


STRAY = r"neither entries of SystemVariables nor declared in params"
REMEDIES = "list it in SystemVariables, declare it in params=..., or mark it 'Constant'"


# ---------------------------------------------------------------------------
# F1 — dbg-A25: a time derivative of a variation vector in δW is accepted silently
# ---------------------------------------------------------------------------


class TestF1DerivativeOfVariationInVirtualWork:
    """``dW = Dot(ẋ.delta(), f)`` (dbg-A25) builds ``d/dt(δx)·f``; IBP turns it into
    ``−δx·ḟ`` and the equation ends in ``−ḟ`` instead of ``+f``. Virtual work is a
    one-form on configuration space: a derivative of a variation vector in δW is never
    legitimate and ``check_eom`` cannot even evaluate it (``KeyError`` on ``d/dt(δx)``).
    The fix plan lists A25 under 2.2; at c18fb55 all four variable kinds pass silently."""

    def test_r3_variation_of_velocity_raises(self):
        L, _, x, f, *_ = point_mass()
        with pytest.raises(ValueError, match=r"inf_work.*\\delta\{x\}"):
            compute_eom(L, Dot(x.t_diff().delta(), f), SystemVariables(vectors=[x]))

    def test_s2_derivative_of_variation_vector_raises(self):
        L, q, om, xi, f = spherical_pendulum()
        with pytest.raises(ValueError, match=r"inf_work.*\\xi_\{q\}"):
            compute_eom(L, Dot(xi.t_diff(), f), SystemVariables(vectors=[q]))

    def test_s2_variation_of_velocity_raises(self):
        L, q, om, xi, f = spherical_pendulum()
        with pytest.raises(ValueError, match=r"inf_work.*\\xi_\{q\}"):
            compute_eom(L, Dot(q.t_diff().delta(), f), SystemVariables(vectors=[q]))

    def test_so3_derivative_of_variation_vector_raises(self):
        L, R, Om, eta, M, J = rigid_body()
        with pytest.raises(ValueError, match=r"inf_work.*\\eta_\{R\}"):
            compute_eom(L, Dot(eta.t_diff(), M), SystemVariables(matrices=[R]))

    def test_scalar_variation_of_velocity_raises(self):
        m, k = Scalar("m", attr=["Constant"]), Scalar("k", attr=["Constant"])
        x, f = Scalar("x"), Scalar("f")
        L = HALF * m * x.t_diff() * x.t_diff() - HALF * k * x * x
        with pytest.raises(ValueError, match=r"inf_work.*\\delta\{x\}"):
            compute_eom(L, x.t_diff().delta() * f, SystemVariables(scalars=[x]))

    def test_message_names_the_derivative_and_the_fix(self):
        """Before the fix this derived ``… + d/dt(f)(-1)``; the correct δW still ends in +f."""
        L, dW, x, f, *_ = point_mass()
        bad_term = str(TimeDerivative(x.delta()))
        with pytest.raises(ValueError, match=re.escape(f"'{bad_term}'")) as info:
            compute_eom(L, Dot(x.t_diff().delta(), f), SystemVariables(vectors=[x]))
        assert "not Dot(x.t_diff().delta(), f)" in str(info.value)
        good = compute_eom(L, dW, SystemVariables(vectors=[x]))
        assert str(good[str(x.delta())][1]).endswith("+f)")

    def test_velocity_dependent_force_is_legitimate(self):
        """δx·(−c ẋ) (damping) has ẋ in the force, not in the variation: must derive."""
        L, _, x, f, m, *_ = point_mass()
        c = Scalar("c", attr=["Constant"])
        dW = Dot(x.delta(), x.t_diff() * c * -1)
        eom = compute_eom(L, dW, SystemVariables(vectors=[x]))
        assert check_eom(L, dW, SystemVariables(vectors=[x]), eom) < TOL


# ---------------------------------------------------------------------------
# F2 — a symbol printed like a listed variable's tangent vector passes the name checks
# ---------------------------------------------------------------------------


class TestF2TangentVectorNameClash:
    """``_check_variables`` rejects an entry printed like another entry's *variation*
    vector (``ξ_q``, ``η_R``) and ``_check_params`` rejects a param that is a configuration
    variable, but neither looks at the *tangent* vector names (``ω_q``, ``Ω_R``). A plain
    ``Vector("\\omega_{q}")`` next to ``S2("q")`` is then two different leaves with one
    printed name; every str-keyed consumer (printing, ``check_eom``, ``linearize``) merges
    them and the equation is wrong (``check_eom`` error ≈ 2 for the params case)."""

    def test_param_named_like_the_tangent_vector_raises(self):
        L, q, om, xi, f = spherical_pendulum()
        w = Vector(q.manifold_info.tangent_vector_name)
        with pytest.raises(ValueError, match=r"params entry '\\omega_\{q\}'.*tangent"):
            compute_eom(L + Dot(w, om), Dot(xi, f), SystemVariables(vectors=[q]), params=[w])

    def test_variable_named_like_the_tangent_vector_raises(self):
        L, q, om, xi, f = spherical_pendulum()
        w = Vector(q.manifold_info.tangent_vector_name)
        with pytest.raises(ValueError, match=r"'\\omega_\{q\}'.*tangent vector of 'q'"):
            compute_eom(
                L + HALF * Dot(w.t_diff(), w.t_diff()),
                Dot(xi, f) + Dot(w.delta(), f),
                SystemVariables(vectors=[q, w]),
            )

    def test_variable_named_like_so3_tangent_vector_raises(self):
        L, R, Om, eta, M, J = rigid_body()
        w = Vector(R.manifold_info.tangent_vector_name)
        with pytest.raises(ValueError, match=r"'\\Omega_\{R\}'.*tangent vector of 'R'"):
            compute_eom(
                L + HALF * Dot(w.t_diff(), w.t_diff()),
                Dot(eta, M) + Dot(w.delta(), Vector("F")),
                SystemVariables(vectors=[w], matrices=[R]),
            )

    def test_renamed_param_passes_the_oracle(self):
        """Before the fix the clashing name derived an equation with ``check_eom`` error ≈ 2; the
        same Lagrangian with ``w`` renamed passes, so the rename the message asks for works."""
        L, q, om, xi, f = spherical_pendulum()
        variables = SystemVariables(vectors=[q])
        w_good = Vector("w")
        L_good = L + Dot(w_good, om)
        eom_good = compute_eom(L_good, Dot(xi, f), variables, params=[w_good])
        assert (
            check_eom(L_good, Dot(xi, f), variables, eom_good, inputs=[f], params=[w_good]) < TOL
        )

    def test_variation_vector_name_clash_is_rejected(self):
        """The existing check (ξ_q, η_R) holds in both groups and for params."""
        L, q, om, xi, f = spherical_pendulum()
        with pytest.raises(ValueError, match="printed name of the variation vector of 'q'"):
            compute_eom(L, Dot(xi, f), SystemVariables(vectors=[q, Vector(str(xi))]))
        R = SO3("R")
        eta_name = R.manifold_info.variation_vector_name
        with pytest.raises(
            ValueError, match=re.escape(f"'R' has the variation vector '{eta_name}'")
        ):
            compute_eom(Zero, Zero, SystemVariables(vectors=[Vector(eta_name)], matrices=[R]))


# ---------------------------------------------------------------------------
# F3 — input-shape mistakes with Python's own message
# ---------------------------------------------------------------------------


class TestF3InputShape:
    def test_bare_leaf_for_params_names_the_argument(self):
        L, dW, x, f, *_ = point_mass()
        l = Scalar("l")
        with pytest.raises(TypeError, match="params"):
            compute_eom(L + l * Dot(x, x), dW, SystemVariables(vectors=[x]), params=l)

    def test_bare_leaf_for_vectors_names_the_group(self):
        L, dW, x, f, *_ = point_mass()
        with pytest.raises(TypeError, match="SystemVariables.vectors"):
            compute_eom(L, dW, SystemVariables(vectors=x))

    def test_no_configuration_variables_raises(self):
        m = Scalar("m", attr=["Constant"])
        with pytest.raises(ValueError, match="no configuration variables"):
            compute_eom(m * m, None, SystemVariables())

    def test_empty_variables_with_a_real_lagrangian_raises_stray(self):
        L, dW, x, *_ = point_mass()
        with pytest.raises(ValueError, match=r"'\\delta\{x\}' of the symbols 'x'"):
            compute_eom(L, dW, SystemVariables())

    def test_sequence_forms_for_params_are_equivalent(self):
        L, dW, x, f, m, g, e3 = point_mass()
        l = Scalar("l")
        L2 = L + l * Dot(x, e3)
        variables = SystemVariables(vectors=[x])
        key = str(x.delta())
        ref = str(compute_eom(L2, dW, variables, params=[l])[key][1])
        assert str(compute_eom(L2, dW, variables, params=(l,))[key][1]) == ref
        assert str(compute_eom(L2, dW, variables, params=(p for p in [l]))[key][1]) == ref
        assert str(compute_eom(L2, dW, variables, params={l})[key][1]) == ref
        assert str(TimeDerivative(l)) not in ref and "e3*l" in ref


# ---------------------------------------------------------------------------
# F4 — a hand-built Variation(TimeDerivative(x)) is reported as a stray "d/dt(x)"
# ---------------------------------------------------------------------------


class TestF4VariationOfDerivativeMessage:
    """The engine always builds ``d/dt(δx)``; a user who writes ``Variation(x.t_diff())``
    gets "variation symbols 'δ(d/dt(x))' of the symbols 'd/dt(x)' … list it in
    SystemVariables", which ``_check_variables`` would then reject. The message should
    name ``x`` and say that δ(ẋ) is written ``x.t_diff().delta()`` (= d/dt(δx))."""

    def test_message_names_the_leaf_and_the_commuted_form(self):
        L, dW, x, f, *_ = point_mass()
        with pytest.raises(ValueError, match=r"of the symbols 'x'") as info:
            compute_eom(L, dW + Dot(Variation(x.t_diff()), f), SystemVariables(vectors=[x]))
        assert "d/dt" in str(info.value) or "t_diff" in str(info.value)

    def test_it_does_raise(self):
        L, dW, x, f, *_ = point_mass()
        with pytest.raises(ValueError, match=STRAY):
            compute_eom(L, dW + Dot(Variation(x.t_diff()), f), SystemVariables(vectors=[x]))


# ---------------------------------------------------------------------------
# F5 — attr given as a string raises about the attribute 'C'
# ---------------------------------------------------------------------------


class TestF5AttrAsString:
    def test_string_attr_is_named_as_the_cause(self):
        with pytest.raises((TypeError, ValueError), match=r"attr=\['Constant'\]"):
            Scalar("m", attr="Constant")

    def test_string_attr_does_raise(self):
        """It never silently accepts, and the message names the string, not its first
        character (before the fix: ``unknown attribute 'C'``)."""
        for a in sorted(VALID_ATTRS):
            with pytest.raises(TypeError, match="attr must be a list of strings, got str"):
                Vector("v", attr=a)

    def test_non_string_entries_raise(self):
        for bad in (None, 1, "", " Constant", "constant", "CONSTANT"):
            with pytest.raises(ValueError, match="unknown attribute"):
                Scalar("m", attr=[bad])


# ---------------------------------------------------------------------------
# 2.2 held: legitimate systems with prescribed parameters derive and pass the numerical check
# ---------------------------------------------------------------------------


class TestParamsLegitimateSystems:
    def test_input_in_both_l_and_dw(self):
        L, dW, x, f, *_ = point_mass()
        variables = SystemVariables(vectors=[x])
        L2 = L + Dot(x, f)
        with pytest.raises(ValueError, match=r"'\\delta\{f\}' of the symbols 'f'"):
            compute_eom(L2, dW, variables)
        eom = compute_eom(L2, dW, variables, params=[f])
        assert str(eom[str(x.delta())][1]).endswith("+f*(2))")
        assert check_eom(L2, dW, variables, eom, inputs=[f], params=[f]) < TOL

    def test_vector_param_only_under_time_derivative(self):
        """δ(ḃ) is built as d/dt(δb), so the walk sees Variation(b) with b declared."""
        L, dW, x, f, m, *_ = point_mass()
        b = Vector("b")
        L2 = L + m * Dot(x.t_diff(), b.t_diff())
        variables = SystemVariables(vectors=[x])
        eom = compute_eom(L2, dW, variables, params=[b])
        eq = str(eom[str(x.delta())][1])
        assert str(TimeDerivative(TimeDerivative(b))) in eq
        assert check_eom(L2, dW, variables, eom, inputs=[f], params=[b]) < TOL

    def test_scalar_param_only_under_time_derivative(self):
        L, dW, x, f, m, g, e3 = point_mass()
        l = Scalar("l")
        L2 = L + m * l.t_diff() * Dot(x.t_diff(), e3)
        variables = SystemVariables(vectors=[x])
        eom = compute_eom(L2, dW, variables, params=[l])
        assert check_eom(L2, dW, variables, eom, inputs=[f], params=[l]) < TOL

    def test_prescribed_inertia_matrix(self):
        J = Matrix("J")
        R, M = SO3("R"), Vector("M")
        Om, eta = R.get_tangent_vector(), R.get_variation_vector()
        L = HALF * Dot(Om, J * Om)
        variables = SystemVariables(matrices=[R])
        with pytest.raises(ValueError, match=r"'\\delta\{J\}' of the symbols 'J'"):
            compute_eom(L, Dot(eta, M), variables)
        eom = compute_eom(L, Dot(eta, M), variables, params=[J])
        eq = str(eom[str(eta)][1])
        assert str(TimeDerivative(J)) in eq
        assert check_eom(L, Dot(eta, M), variables, eom, inputs=[M], params=[J]) < TOL
        row = to_standard_form(eom, variables, inputs=[M])[str(eta)]
        assert "J" in str(row.M)

    def test_param_shared_by_two_bodies(self):
        m, g = Scalar("m", attr=["Constant"]), Scalar("g", attr=["Constant"])
        e3 = Vector("e3", attr=["Constant"])
        l = Scalar("l")
        q1, q2 = S2("q1"), S2("q2")
        f1, f2 = Vector("f1"), Vector("f2")
        o1, x1 = q1.get_tangent_vector(), q1.get_variation_vector()
        o2, x2 = q2.get_tangent_vector(), q2.get_variation_vector()
        L = (
            HALF * m * l * l * Dot(o1, o1)
            + HALF * m * l * l * Dot(o2, o2)
            - m * g * l * Dot(q1, e3)
            - m * g * l * Dot(q2, e3)
        )
        dW = Dot(x1, f1) + Dot(x2, f2)
        variables = SystemVariables(vectors=[q1, q2])
        eom = compute_eom(L, dW, variables, params=[l])
        assert set(eom) == {str(x1), str(x2)}
        assert check_eom(L, dW, variables, eom, inputs=[f1, f2], params=[l]) < TOL

    def test_param_under_hat_and_inside_an_outer_product(self):
        L, R, Om, eta, M, J = rigid_body()
        b, c = Vector("b"), Vector("c", attr=["Constant"])
        variables = SystemVariables(matrices=[R])
        for L2 in (L + Dot(Om, Hat(b) * Om), L + Dot(c, MVMul(VVMul(b, Transpose(c)), Om))):
            eom = compute_eom(L2, Dot(eta, M), variables, params=[b])
            assert check_eom(L2, Dot(eta, M), variables, eom, inputs=[M], params=[b]) < TOL

    def test_moving_attachment_point(self):
        """Pendulum hung from a prescribed point b(t): x = b + l q."""
        m, g, l = (Scalar(s, attr=["Constant"]) for s in "mgl")
        e3 = Vector("e3", attr=["Constant"])
        q, f, b = S2("q"), Vector("f"), Vector("b")
        xi = q.get_variation_vector()
        x = b + q * l
        L = HALF * m * Dot(x.t_diff(), x.t_diff()) - m * g * Dot(x, e3)
        variables = SystemVariables(vectors=[q])
        eom = compute_eom(L, Dot(xi, f), variables, params=[b])
        assert str(TimeDerivative(TimeDerivative(b))) in str(eom[str(xi)][1])
        assert check_eom(L, Dot(xi, f), variables, eom, inputs=[f], params=[b]) < TOL
        to_standard_form(eom, variables, inputs=[f])
        with pytest.raises(AlgebraicError, match=r"\\delta\{b\}.*contains no perturbation symbol"):
            linearize(eom, variables, [f])

    def test_non_constant_force_only_in_dw_needs_no_params(self):
        L, dW, x, f, *_ = point_mass()
        b = Vector("b")
        eom = compute_eom(L, Dot(x.delta(), f + b), SystemVariables(vectors=[x]))
        assert str(eom[str(x.delta())][1]).endswith("+f+b)")

    def test_input_through_rotation_in_dw(self):
        """Quadrotor thrust δx·(R e₃ u): u and R appear only in δW, never varied."""
        m, g = Scalar("m", attr=["Constant"]), Scalar("g", attr=["Constant"])
        e3 = Vector("e3", attr=["Constant"])
        x, u, M = Vector("x"), Scalar("u"), Vector("M")
        L_rb, R, Om, eta, _, J = rigid_body()
        L = HALF * m * Dot(x.t_diff(), x.t_diff()) + L_rb - m * g * Dot(x, e3)
        dW = Dot(x.delta(), R * e3 * u) + Dot(eta, M)
        variables = SystemVariables(vectors=[x], matrices=[R])
        eom = compute_eom(L, dW, variables)
        assert check_eom(L, dW, variables, eom, inputs=[u, M]) < TOL
        to_standard_form(eom, variables, inputs=[u, M])
        linearize(eom, variables, [u, M])

    @pytest.mark.parametrize(
        "build",
        [
            vbl_systems.point_mass_r3,
            vbl_systems.spring_mass_r,
            vbl_systems.free_rigid_body,
            vbl_systems.rigid_pendulum_pivot,
            vbl_systems.spherical_pendulum_s2,
            vbl_systems.quadrotor,
            vbl_systems.double_rigid_pendulum,
            vbl_systems.elastic_spherical_pendulum,
            vbl_systems.quadrotor_with_load,
        ],
    )
    def test_pinned_systems_run_the_whole_pipeline(self, build):
        s = build()
        eom = compute_eom(s.L, s.dW, s.variables)
        to_standard_form(eom, s.variables, inputs=s.inputs)
        linearize(eom, s.variables, s.inputs)


# ---------------------------------------------------------------------------
# 2.2 held: misuse that raises, with the symbol and the three remedies in the message
# ---------------------------------------------------------------------------


class TestStrayVariationMisuse:
    def _assert_stray(self, info, symbols, names):
        msg = str(info.value)
        for s in symbols:
            assert f"'{s}'" in msg, msg
        assert f"of the symbols {', '.join(repr(n) for n in names)}" in msg, msg
        assert REMEDIES in msg

    def test_unlisted_vector_variable(self):
        L, dW, x, f, m, *_ = point_mass()
        y = Vector("y")
        with pytest.raises(ValueError) as info:
            compute_eom(
                L + HALF * m * Dot(y.t_diff(), y.t_diff()) - Dot(x, y),
                dW,
                SystemVariables(vectors=[x]),
            )
        self._assert_stray(info, ["\\delta{y}"], ["y"])

    def test_unlisted_scalar_variable(self):
        L, dW, x, f, m, g, e3 = point_mass()
        s = Scalar("s")
        with pytest.raises(ValueError) as info:
            compute_eom(
                L + HALF * m * s.t_diff() * s.t_diff() - s * Dot(x, e3),
                dW,
                SystemVariables(vectors=[x]),
            )
        self._assert_stray(info, ["\\delta{s}"], ["s"])

    def test_unlisted_s2_point_appearing_only_statically(self):
        L, dW, x, *_ = point_mass()
        p = S2("p")
        with pytest.raises(ValueError) as info:
            compute_eom(L - Dot(x, p), dW, SystemVariables(vectors=[x]))
        self._assert_stray(info, ["\\xi_{p}"], ["p"])

    def test_unlisted_so3_point_through_a_body_fixed_vector(self):
        L, dW, x, f, m, g, e3 = point_mass()
        R, rho = SO3("R"), Vector("rho", attr=["Constant"])
        L2 = L - m * g * Dot(R * rho, e3) + Dot(x.t_diff(), (R * rho).t_diff())
        with pytest.raises(ValueError) as info:
            compute_eom(L2, dW, SystemVariables(vectors=[x]))
        self._assert_stray(info, ["\\eta_{R}", "\\delta{\\Omega_{R}}"], ["R"])

    def test_s2_in_lagrangian_but_vector_listed_and_vice_versa(self):
        L, q, om, xi, f = spherical_pendulum()
        with pytest.raises(ValueError) as info:
            compute_eom(L, Dot(xi, f), SystemVariables(vectors=[Vector("q")]))
        self._assert_stray(info, ["\\xi_{q}", "\\delta{\\omega_{q}}"], ["q"])
        m, g = Scalar("m", attr=["Constant"]), Scalar("g", attr=["Constant"])
        qv = Vector("q")
        Lv = HALF * m * Dot(qv.t_diff(), qv.t_diff()) - m * g * Dot(
            qv, Vector("e3", attr=["Constant"])
        )
        with pytest.raises(ValueError) as info:
            compute_eom(Lv, Dot(qv.delta(), f), SystemVariables(vectors=[S2("q")]))
        self._assert_stray(info, ["\\delta{q}"], ["q"])

    def test_several_strays_are_all_named(self):
        L, dW, x, *_ = point_mass()
        y, p, l = Vector("y"), S2("p"), Scalar("l")
        with pytest.raises(ValueError) as info:
            compute_eom(L + Dot(x, y) * l + Dot(p, x), dW, SystemVariables(vectors=[x]))
        self._assert_stray(info, ["\\delta{l}", "\\delta{y}", "\\xi_{p}"], ["l", "p", "y"])

    def test_param_that_is_a_variable_or_a_constant(self):
        L, dW, x, f, m, g, *_ = point_mass()
        with pytest.raises(ValueError, match="params entry 'x' is also a configuration variable"):
            compute_eom(L, dW, SystemVariables(vectors=[x]), params=[x])
        with pytest.raises(ValueError, match="params entry 'g' is constant"):
            compute_eom(L, dW, SystemVariables(vectors=[x]), params=[g])

    def test_cancelling_stray_terms_are_not_reported(self):
        """Dot(x, y) − Dot(y, x) vanishes before the check: nothing survives to report."""
        L, dW, x, *_ = point_mass()
        y = Vector("y")
        eom = compute_eom(L + Dot(x, y) - Dot(y, x), dW, SystemVariables(vectors=[x]))
        assert "y" not in str(eom[str(x.delta())][1])


# ---------------------------------------------------------------------------
# 2.3 held
# ---------------------------------------------------------------------------


class TestSystemVariablesMisuse:
    def test_wrong_group_for_a_manifold_point(self):
        with pytest.raises(TypeError, match="vectors entry 'R' has type ExprType.MATRIX"):
            compute_eom(Zero, Zero, SystemVariables(vectors=[SO3("R")]))
        with pytest.raises(TypeError, match="matrices entry 'q' has type ExprType.VECTOR"):
            compute_eom(Zero, Zero, SystemVariables(matrices=[S2("q")]))

    def test_non_expression_entries(self):
        for bad in (None, "x", 3):
            with pytest.raises(TypeError, match="has type None"):
                compute_eom(Zero, Zero, SystemVariables(vectors=[bad]))

    def test_orthogonal_matrix_that_is_not_so3(self):
        with pytest.raises(TypeError, match="'R' is a general Matrix; only SO3"):
            compute_eom(
                Zero, Zero, SystemVariables(matrices=[Matrix("R", attr=["OrthogonalMatrix"])])
            )

    def test_two_distinct_objects_with_one_name(self):
        L, dW, x, f, *_ = point_mass()
        with pytest.raises(ValueError, match="'x' is listed more than once"):
            compute_eom(L, dW, SystemVariables(vectors=[x, Vector("x")]))
        with pytest.raises(ValueError, match="'q' is listed more than once"):
            compute_eom(Zero, Zero, SystemVariables(vectors=[Vector("q"), S2("q")]))

    def test_linearize_runs_the_same_checks(self):
        L, dW, x, f, *_ = point_mass()
        eom = compute_eom(L, dW, SystemVariables(vectors=[x]))
        with pytest.raises(ValueError, match="'x' is listed more than once"):
            linearize(eom, SystemVariables(vectors=[x, Vector("x")]), [f])
        with pytest.raises(TypeError, match="tangent or variation vector"):
            linearize(eom, SystemVariables(vectors=[S2("q").get_tangent_vector()]), [f])


# ---------------------------------------------------------------------------
# 2.5 held: the tangent rules compare the parent by name
# ---------------------------------------------------------------------------


class TestTangentParentByName:
    def test_two_objects_mixed_in_one_lagrangian_pass_the_oracle(self):
        m, g, l = (Scalar(s, attr=["Constant"]) for s in "mgl")
        e3 = Vector("e3", attr=["Constant"])
        q1, q2, f = S2("q"), S2("q"), Vector("f")
        om, xi = q1.get_tangent_vector(), q1.get_variation_vector()
        L = HALF * m * l * l * Dot(om, om) - m * g * l * Dot(q2, e3)
        variables = SystemVariables(vectors=[q2])
        eom = compute_eom(L, Dot(xi, f), variables)
        ref = compute_eom(
            HALF * m * l * l * Dot(om, om) - m * g * l * Dot(q1, e3),
            Dot(xi, f),
            SystemVariables(vectors=[q1]),
        )
        assert str(eom[str(xi)][1]) == str(ref[str(xi)][1])
        assert check_eom(L, Dot(xi, f), variables, eom, inputs=[f]) < TOL

    def test_q_and_q1_never_cross_fire(self):
        q, q1 = S2("q"), S2("q1")
        om = q.get_tangent_vector()
        assert str(full_simplify(Dot(om, q1))) == "Dot(\\omega_{q},q1)"
        assert str(full_simplify(Cross(q1, Cross(om, q1)))) == "Cross(q1,Cross(\\omega_{q},q1))"

    def test_parent_with_different_attrs_is_the_same_point(self):
        q, q_const = S2("q"), S2("q", attr=["Constant"])
        om = q.get_tangent_vector()
        assert full_simplify(Dot(om, q_const)).is_zero
        assert str(full_simplify(Cross(q_const, Cross(om, q_const)))) == str(om)

    def test_plain_vector_with_the_point_name_is_not_reduced(self):
        q = S2("q")
        om = q.get_tangent_vector()
        assert not full_simplify(Dot(om, Vector("q"))).is_zero

    def test_parentless_ts2_is_not_reduced(self):
        assert not full_simplify(Dot(TS2("w"), S2("q"))).is_zero


# ---------------------------------------------------------------------------
# 2.8 held: every legal calculus form still builds
# ---------------------------------------------------------------------------


class TestCalculusGuardsScope:
    def test_legal_forms(self):
        q, x = S2("q"), Vector("x")
        om, xi = q.get_tangent_vector(), q.get_variation_vector()
        assert isinstance(TimeDerivative(om), TimeDerivative)
        assert isinstance(Variation(om), Variation)
        assert isinstance(TimeDerivative(xi), TimeDerivative)
        assert isinstance(TimeDerivative(TimeDerivative(x)), TimeDerivative)
        assert isinstance(Variation(TimeDerivative(x)), Variation)
        assert isinstance(TimeDerivative(Variation(x)), TimeDerivative)
        assert str((q * Scalar("l")).t_diff()) == "(Cross(\\omega_{q},q)*l+q*\\frac{d}{dt}(l))"
        assert str((Matrix("A", attr=["Constant"]) * q).t_diff()) == "A*Cross(\\omega_{q},q)"
        assert str(Hat(q).t_diff()) == "Hat(Cross(\\omega_{q},q))"
        assert str(q.T()) == "(q)'"

    def test_manifold_leaf_under_calculus_raises_for_every_entry_point(self):
        q, R = S2("q"), SO3("R")
        with pytest.raises(ValueError, match=r"use q\.t_diff\(\)"):
            TimeDerivative(q)
        with pytest.raises(ValueError, match=r"use R\.delta\(\)"):
            Variation(R)


# ---------------------------------------------------------------------------
# 2.9 held: collect only rejects what is not linear
# ---------------------------------------------------------------------------


class TestCollectOrder:
    def test_dot_of_target_with_itself_is_now_rejected(self):
        """Before 2.9 the ``l == vec`` shortcut returned ``Dot(v, v)`` as "canonical"."""
        v = Vector("v")
        with pytest.raises(NonLinearError, match="'v' appears on both sides"):
            collect(Dot(v, v), v)

    def test_single_occurrence_anywhere_is_collected_or_typed_error(self):
        v, a, b, M = Vector("v"), Vector("a"), Vector("b"), Matrix("M")
        shapes = [
            Dot(a, v),
            Dot(v, a),
            Dot(a, Cross(b, v)),
            Dot(Cross(v, a), b),
            Dot(a, M * v),
            Dot(M * v, a),
            Dot(a, M * Cross(b, v)),
            Dot(a, b) * Dot(v, a),
        ]
        for s in shapes:
            out = collect(s, v)
            assert out.has(v)
            assert "NonLinear" not in type(out).__name__

    def test_target_inside_a_matrix_factor_is_a_typed_error(self):
        v, a, b = Vector("v"), Vector("a"), Vector("b")
        with pytest.raises(AlgebraicError, match="inside the matrix factor"):
            collect(Dot(MVMul(Hat(v), a), b), v)


# ---------------------------------------------------------------------------
# 2.10 / 2.11 held, with the edges pinned
# ---------------------------------------------------------------------------


class TestCapWarningAndNumbers:
    def _two_pass(self):
        q = S2("q")
        om = q.get_tangent_vector()
        return Dot(MVMul(Hat(q), Cross(om, q)), q)

    def test_warning_fires_only_when_the_loop_exhausts(self):
        """max_iter=2 reaches 0 on the second pass but has no confirming pass: it warns,
        as documented ("may be unsimplified"); max_iter=3 confirms and is silent."""
        expr = self._two_pass()
        with pytest.warns(RuntimeWarning, match="after 1 passes"):
            assert not full_simplify(expr, max_iter=1).is_zero
        with pytest.warns(RuntimeWarning, match="after 2 passes"):
            assert full_simplify(expr, max_iter=2).is_zero
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            assert full_simplify(expr, max_iter=3).is_zero

    def test_zero_flag_does_not_break_like_term_merging(self):
        a = Scalar("a")
        assert str(full_simplify(a + Number(0) * a + a)) == "a*(2)"
        assert str(full_simplify(Number(0) + a + Number(0))) == "a"
        assert str(full_simplify(Number(2) + Number(-2) + a)) == "a"

    def test_ones_is_a_genuine_coefficient_where_it_matters(self):
        a, v, w = Scalar("a"), Vector("v"), Vector("w")
        assert str(full_simplify(Number(1) * a + a)) == "a*(2)"
        assert str(full_simplify(v * Number(1) + v)) == "v*(2)"
        assert str(full_simplify(Number(1) + Number(1))) == "(2)"
        assert str(full_simplify(Number(1) * Number(1))) == "(1)"
        assert str(extract_linear_coeff(Dot(v, w) * Number(1), v)) == "w*(1)"
        assert str(extract_linear_coeff(Dot(v, w * Number(0)), v)) == "w*(0)"

    def test_float_and_signed_zero_literals(self):
        a = Scalar("a")
        assert full_simplify(a * -0.0).is_zero
        assert full_simplify(a * 0.0).is_zero
        assert str(full_simplify(a * 1.0)) == "a"
        assert Number(1e-300).is_zero is False

    def test_bool_literals_are_rejected(self):
        """bool is an int; before the fix ``Number(True)`` printed ``(True)`` with the Ones
        flag and ``compute_eom(L, False, …)`` meant δW = 0. Both now raise."""
        with pytest.raises(TypeError, match="got bool True"):
            Number(True)
        with pytest.raises(TypeError, match="got bool False"):
            Number(False)
        with pytest.raises(TypeError, match="got bool True"):
            Scalar("x") * True
        L, _, x, *_ = point_mass()
        with pytest.raises(TypeError, match="inf_work .*got bool False"):
            compute_eom(L, False, SystemVariables(vectors=[x]))

    def test_named_numeric_scalar_keeps_its_name(self):
        """D1 is open: Scalar('m', value=0) is not folded to the Zero leaf."""
        m = Scalar("m", value=0)
        assert not m.is_zero and str(full_simplify(m * Scalar("a"))) == "m*a"
