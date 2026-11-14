"""Phase 2 of ``specs/fix-plan.md``: ill-posed input raises instead of returning silently.

One class per item of the phase 2 table. Every test asserts the exception type and that
the message names the offending symbol, in the style of ``tests/test_vbl_validation.py``.
The legitimate inputs next to each guard are pinned too, so the guards cannot widen.
"""

import pathlib
import re
import warnings

import pytest

from geomech import (
    S2,
    SO3,
    AlgebraicError,
    Cross,
    Dot,
    ExpressionMismatchError,
    Hat,
    Matrix,
    MVMul,
    NonLinearError,
    Number,
    Scalar,
    SizeMismatchError,
    SkewSymmMatrix,
    SystemVariables,
    TimeDerivative,
    Transpose,
    Variation,
    Vector,
    VectorExpr,
    VVMul,
    Zero,
    collect,
    compute_eom,
    full_simplify,
    getScalars,
    getVectors,
    linearize,
    to_standard_form,
)
from geomech.core.base.types import VALID_ATTRS
from geomech.core.transformations.manifold_rules import apply_manifold_rules
from geomech.dynamics.vbl import _check_eom
from geomech.utils import errors

HALF = Number(0.5)


def _point_mass():
    """L = ½ m ẋ·ẋ − m g e₃·x, δW = δx·f."""
    m, g = getScalars("m g", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    x, f = Vector("x"), Vector("f")
    L = HALF * m * Dot(x.t_diff(), x.t_diff()) - m * g * Dot(e3, x)
    return L, Dot(x.delta(), f), x, f


def _spherical_pendulum(q):
    """L = ½ m l² ω·ω − m g l q·e₃, δW = ξ·f for the given S2 point."""
    m, g, l = getScalars("m g l", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    f = Vector("f")
    om, xi = q.get_tangent_vector(), q.get_variation_vector()
    L = HALF * m * l * l * Dot(om, om) - m * g * l * Dot(q, e3)
    return L, Dot(xi, f), f


# ---------------------------------------------------------------------------
# 2.1 unknown attr strings
# ---------------------------------------------------------------------------


class TestUnknownAttr:
    def test_lowercase_constant_raises_naming_symbol_and_attribute(self):
        with pytest.raises(ValueError, match="'constant'.*symbol 'm'") as info:
            Scalar("m", attr=["constant"])
        assert "valid attributes" in str(info.value)
        assert "'Constant'" in str(info.value)

    def test_vector_matrix_and_manifold_leaves(self):
        with pytest.raises(ValueError, match="'Const'.*'e3'"):
            Vector("e3", attr=["Const"])
        with pytest.raises(ValueError, match="'Symmetric'.*'J'") as info:
            Matrix("J", attr=["Symmetric"])
        assert "'SymmetricMatrix'" in str(info.value)
        with pytest.raises(ValueError, match="'unitnorm'.*'q'"):
            S2("q", attr=["unitnorm"])
        with pytest.raises(ValueError, match="'Orthogonal'.*'R'"):
            SO3("R", attr=["Orthogonal"])
        with pytest.raises(ValueError, match="'skew'.*'K'"):
            SkewSymmMatrix("K", attr=["skew"])

    def test_helpers_propagate(self):
        with pytest.raises(ValueError, match="'constant'.*'m'"):
            getScalars("m g", attr=["constant"])
        with pytest.raises(ValueError, match="'constant'.*'a'"):
            getVectors("a b", attr=["constant"])

    def test_every_valid_attribute_is_accepted(self):
        for a in sorted(VALID_ATTRS):
            Scalar("s", attr=[a])
            Vector("v", attr=[a])
            Matrix("M", attr=[a])

    def test_tangent_vector_marker_accepted_and_sets_no_flag(self):
        w = Vector("w", attr=["TangentVector"])
        assert not w.is_constant and not w.is_unit_norm and not w.is_manifold
        q, R = S2("q"), SO3("R")
        assert "TangentVector" in q.get_tangent_vector().attr
        assert "TangentVector" in R.get_variation_vector().attr

    def test_valid_combination_unchanged(self):
        J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
        assert J.is_constant and J.is_symmetric


# ---------------------------------------------------------------------------
# 2.2 variation symbols of an unlisted manifold point
# ---------------------------------------------------------------------------


class TestStrayManifoldVariation:
    def test_unlisted_s2_point_raises_naming_symbols_and_point(self):
        """x is listed, the S2 point q of the second term is not: its ξ and δω survive."""
        L, dW, x, f = _point_mass()
        q = S2("q")
        om = q.get_tangent_vector()
        e3 = Vector("e3", attr=["Constant"])
        L2 = L + HALF * Dot(om, om) - Dot(q, e3)
        with pytest.raises(ValueError, match=r"of the symbols 'q'") as info:
            compute_eom(L2, dW, SystemVariables(vectors=[x]))
        msg = str(info.value)
        assert "\\xi_{q}" in msg and "\\delta{\\omega_{q}}" in msg
        assert "list it in SystemVariables, declare it in params=..., or mark it 'Constant'" in msg

    def test_vector_listed_for_an_s2_point_raises(self):
        """Building L with S2("q") but listing Vector("q") used to return 0v silently."""
        L, dW, _ = _spherical_pendulum(S2("q"))
        with pytest.raises(ValueError, match="'q'"):
            compute_eom(L, dW, SystemVariables(vectors=[Vector("q")]))

    def test_unlisted_so3_point_raises(self):
        L, dW, x, f = _point_mass()
        R = SO3("R")
        J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
        Om = R.get_tangent_vector()
        L2 = L + HALF * Dot(Om, J * Om)
        with pytest.raises(ValueError, match=r"of the symbols 'R'") as info:
            compute_eom(L2, dW, SystemVariables(vectors=[x]))
        assert "\\delta{\\Omega_{R}}" in str(info.value)

    def test_stray_variation_vector_in_dw_raises(self):
        L, dW, x, f = _point_mass()
        p = S2("p")
        with pytest.raises(ValueError, match=r"'\\xi_\{p\}'.*of the symbols 'p'"):
            compute_eom(L, dW + Dot(p.get_variation_vector(), f), SystemVariables(vectors=[x]))

    def test_force_in_dw_is_not_a_variation(self):
        """A non-constant F in δW = δx·F never appears under a Variation: accepted."""
        L, dW, x, f = _point_mass()
        assert not f.is_constant
        eom = compute_eom(L, dW, SystemVariables(vectors=[x]))
        assert str(eom[str(x.delta())][1]).endswith("+f)")

    def test_unlisted_vector_degree_of_freedom_raises(self):
        """dbg-A24: y is used in L but not listed; its equation and its terms used to vanish."""
        L, dW, x, f = _point_mass()
        y = Vector("y")
        L2 = L + HALF * Dot(y.t_diff(), y.t_diff()) - Dot(x, y)
        with pytest.raises(ValueError, match=r"'\\delta\{y\}' of the symbols 'y'") as info:
            compute_eom(L2, dW, SystemVariables(vectors=[x]))
        assert "declare it in params=..." in str(info.value)

    def test_variation_in_dw_of_unlisted_vector_raises(self):
        """dyn-F5: δW = δx·F + δp·F with p unlisted used to return the δx equation alone."""
        L, dW, x, f = _point_mass()
        p = Vector("p")
        with pytest.raises(ValueError, match=r"'\\delta\{p\}' of the symbols 'p'"):
            compute_eom(L, dW + Dot(p.delta(), f), SystemVariables(vectors=[x]))

    def test_unflagged_parameter_raises_and_names_the_fixes(self):
        """dyn-F5: e3 without 'Constant' and without params is neither a variable nor a
        parameter; before, the equation was the same as for a constant e3 and the user
        could not tell which was meant."""
        m, g = getScalars("m g", attr=["Constant"])
        x, f = Vector("x"), Vector("f")
        e3 = Vector("e3")
        L = HALF * m * Dot(x.t_diff(), x.t_diff()) - m * g * Dot(e3, x)
        with pytest.raises(ValueError, match=r"'\\delta\{e3\}' of the symbols 'e3'") as info:
            compute_eom(L, Dot(x.delta(), f), SystemVariables(vectors=[x]))
        assert "mark it 'Constant'" in str(info.value)
        ok = compute_eom(L, Dot(x.delta(), f), SystemVariables(vectors=[x]), params=[e3])
        assert "e3" in str(ok[str(x.delta())][1])

    def test_prescribed_scalar_parameter_needs_params(self):
        """l(t) declared in params: its rate l̇ is derived and δl is dropped; undeclared,
        the surviving δl is an error naming l."""
        m, g = getScalars("m g", attr=["Constant"])
        l = Scalar("l")
        e3 = Vector("e3", attr=["Constant"])
        q, f = S2("q"), Vector("f")
        om, xi = q.get_tangent_vector(), q.get_variation_vector()
        L = HALF * m * l * l * Dot(om, om) - m * g * l * Dot(q, e3)
        with pytest.raises(ValueError, match=r"'\\delta\{l\}' of the symbols 'l'"):
            compute_eom(L, Dot(xi, f), SystemVariables(vectors=[q]))
        eom = compute_eom(L, Dot(xi, f), SystemVariables(vectors=[q]), params=[l])
        assert eom[str(xi)][1].has(TimeDerivative(l))
        assert "\\delta{l}" not in str(eom[str(xi)][1])

    def test_listed_manifold_points_pass(self):
        q, R = S2("q"), SO3("R")
        L, dW, f = _spherical_pendulum(q)
        J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
        Om, eta = R.get_tangent_vector(), R.get_variation_vector()
        M = Vector("M")
        eom = compute_eom(
            L + HALF * Dot(Om, J * Om),
            dW + Dot(eta, M),
            SystemVariables(vectors=[q], matrices=[R]),
        )
        assert set(eom) == {str(q.get_variation_vector()), str(eta)}


class TestParamsValidation:
    def _pendulum_with_length(self):
        m, g = getScalars("m g", attr=["Constant"])
        l = Scalar("l")
        e3 = Vector("e3", attr=["Constant"])
        q, f = S2("q"), Vector("f")
        om, xi = q.get_tangent_vector(), q.get_variation_vector()
        L = HALF * m * l * l * Dot(om, om) - m * g * l * Dot(q, e3)
        return L, Dot(xi, f), q, l, f

    def test_tangent_vector_rejected(self):
        L, dW, q, l, f = self._pendulum_with_length()
        with pytest.raises(TypeError, match=r"params entry '\\omega_\{q\}' is a tangent"):
            compute_eom(L, dW, SystemVariables(vectors=[q]), params=[q.get_tangent_vector()])

    def test_manifold_point_rejected(self):
        L, dW, q, l, f = self._pendulum_with_length()
        with pytest.raises(TypeError, match="params entry 'p' is a manifold point"):
            compute_eom(L, dW, SystemVariables(vectors=[q]), params=[l, S2("p")])

    def test_non_leaf_rejected(self):
        L, dW, q, l, f = self._pendulum_with_length()
        with pytest.raises(
            TypeError, match="must be a plain Scalar, Vector or Matrix leaf, got Mul"
        ):
            compute_eom(L, dW, SystemVariables(vectors=[q]), params=[l * l])

    def test_constant_rejected(self):
        L, dW, q, l, f = self._pendulum_with_length()
        with pytest.raises(ValueError, match="params entry 'g' is constant"):
            compute_eom(
                L, dW, SystemVariables(vectors=[q]), params=[l, Scalar("g", attr=["Constant"])]
            )

    def test_configuration_variable_rejected(self):
        L, dW, x, f = _point_mass()
        with pytest.raises(ValueError, match="params entry 'x' is also a configuration variable"):
            compute_eom(L, dW, SystemVariables(vectors=[x]), params=[x])
        with pytest.raises(ValueError, match="params entry 'x' is also a configuration variable"):
            compute_eom(Zero, Zero, SystemVariables(scalars=[Scalar("x")]), params=[Vector("x")])

    def test_duplicate_rejected(self):
        L, dW, q, l, f = self._pendulum_with_length()
        with pytest.raises(ValueError, match="params entry 'l' is listed more than once"):
            compute_eom(L, dW, SystemVariables(vectors=[q]), params=[l, Scalar("l")])

    def test_unused_param_is_harmless(self):
        L, dW, q, l, f = self._pendulum_with_length()
        ref = compute_eom(L, dW, SystemVariables(vectors=[q]), params=[l])
        out = compute_eom(L, dW, SystemVariables(vectors=[q]), params=[l, Vector("unused")])
        key = str(q.get_variation_vector())
        assert str(out[key][1]) == str(ref[key][1])

    def test_standard_form_and_linearize_on_a_params_equation(self):
        """The prescribed l appears as a time-varying leaf: to_standard_form splits it, and
        linearize keeps raising its clear AlgebraicError on the δl term."""
        L, dW, q, l, f = self._pendulum_with_length()
        variables = SystemVariables(vectors=[q])
        eom = compute_eom(L, dW, variables, params=[l])
        sf = to_standard_form(eom, variables, inputs=[f])
        row = sf[str(q.get_variation_vector())]
        assert "l" in str(row.M) and str(TimeDerivative(l)) in str(row.f)
        with pytest.raises(AlgebraicError, match=r"\\delta\{l\}.*contains no perturbation symbol"):
            linearize(eom, variables, [f])


# ---------------------------------------------------------------------------
# 2.3 SystemVariables validation
# ---------------------------------------------------------------------------


class TestSystemVariablesValidation:
    def test_same_vector_listed_twice(self):
        L, dW, x, _ = _point_mass()
        with pytest.raises(ValueError, match="'x' is listed more than once"):
            compute_eom(L, dW, SystemVariables(vectors=[x, x]))

    def test_two_s2_objects_with_the_same_name(self):
        q1, q2 = S2("q"), S2("q")
        L, dW, _ = _spherical_pendulum(q1)
        with pytest.raises(ValueError, match="'q' is listed more than once"):
            compute_eom(L, dW, SystemVariables(vectors=[q1, q2]))

    def test_scalar_and_vector_with_the_same_name(self):
        L, dW, x, _ = _point_mass()
        with pytest.raises(ValueError, match="'x' is listed more than once.*scalars and vectors"):
            compute_eom(L, dW, SystemVariables(scalars=[Scalar("x")], vectors=[x]))

    def test_entry_named_like_another_entrys_variation_vector(self):
        """A Vector named η_R next to R: η_R would be both a variable and R's variation."""
        R = SO3("R")
        eta_name = R.manifold_info.variation_vector_name
        v = Vector(eta_name)
        with pytest.raises(
            ValueError,
            match=re.escape(f"entry 'R' has the variation vector '{eta_name}'")
            + ".*listed in vectors",
        ):
            compute_eom(Zero, Zero, SystemVariables(vectors=[v], matrices=[R]))
        q = S2("q")
        xi_name = q.manifold_info.variation_vector_name
        with pytest.raises(
            ValueError,
            match=re.escape(
                f"entry '{xi_name}' has the printed name of the variation vector of 'q'"
            ),
        ):
            compute_eom(Zero, Zero, SystemVariables(vectors=[q, Vector(xi_name)]))

    def test_tangent_vector_as_configuration_variable(self):
        q = S2("q")
        L, dW, _ = _spherical_pendulum(q)
        om = q.get_tangent_vector()
        with pytest.raises(TypeError, match=r"'\\omega_\{q\}' is a tangent or variation vector"):
            compute_eom(L, dW, SystemVariables(vectors=[om]))

    def test_variation_vector_as_configuration_variable(self):
        R = SO3("R")
        eta = R.get_variation_vector()
        with pytest.raises(TypeError, match=r"'\\eta_\{R\}'.*list the manifold point 'R'"):
            compute_eom(Zero, Zero, SystemVariables(vectors=[eta]))

    def test_time_derivative_as_configuration_variable(self):
        L, dW, x, _ = _point_mass()
        with pytest.raises(TypeError, match="TimeDerivative node.*underlying variable 'x'"):
            compute_eom(L, dW, SystemVariables(vectors=[x.t_diff()]))

    def test_variation_node_as_configuration_variable(self):
        L, dW, x, _ = _point_mass()
        with pytest.raises(TypeError, match="Variation node.*underlying variable 'x'"):
            compute_eom(L, dW, SystemVariables(vectors=[x.delta()]))
        with pytest.raises(TypeError, match="Variation node.*underlying variable 's'"):
            compute_eom(L, dW, SystemVariables(scalars=[Scalar("s").delta()]))

    def test_general_matrix_as_configuration_variable(self):
        L, dW, x, _ = _point_mass()
        with pytest.raises(TypeError, match="'J' is a general Matrix; only SO3"):
            compute_eom(L, dW, SystemVariables(vectors=[x], matrices=[Matrix("J")]))
        with pytest.raises(TypeError, match="'K' is a general SkewSymmMatrix"):
            compute_eom(L, dW, SystemVariables(matrices=[SkewSymmMatrix("K")]))

    def test_constant_still_reported_as_constant(self):
        """The existing constant check keeps its message (and precedes the SO3 check)."""
        with pytest.raises(ValueError, match="'K' is constant"):
            compute_eom(Zero, Zero, SystemVariables(matrices=[Matrix("K", attr=["Constant"])]))

    def test_linearize_rejects_duplicates(self):
        L, dW, x, f = _point_mass()
        eom = compute_eom(L, dW, SystemVariables(vectors=[x]))
        with pytest.raises(ValueError, match="'x' is listed more than once"):
            linearize(eom, SystemVariables(vectors=[x, x]), [f])

    def test_check_eom_compares_key_count(self):
        """Direct call: one key for two (duplicate) variables is a length mismatch."""
        L, dW, x, f = _point_mass()
        eom = compute_eom(L, dW, SystemVariables(vectors=[x]))
        with pytest.raises(ValueError, match="do not match the variation vectors"):
            _check_eom(eom, SystemVariables(vectors=[x, x]))
        _check_eom(eom, SystemVariables(vectors=[x]))


# ---------------------------------------------------------------------------
# 2.4 constructors copy the attr list (done in phase 1; pinned here)
# ---------------------------------------------------------------------------


class TestAttrListCopied:
    @pytest.mark.parametrize(
        "build",
        [
            lambda attr: Scalar("s", attr=attr),
            lambda attr: Vector("v", attr=attr),
            lambda attr: Matrix("M", attr=attr),
            lambda attr: S2("q", attr=attr),
            lambda attr: SO3("R", attr=attr),
            lambda attr: SkewSymmMatrix("K", attr=attr),
        ],
    )
    def test_callers_list_is_not_aliased(self, build):
        shared = ["Constant"]
        leaf = build(shared)
        assert shared == ["Constant"]
        assert leaf.is_constant
        shared.append("Zero")
        assert not leaf.is_zero


# ---------------------------------------------------------------------------
# 2.5 S2 tangent rules compare the parent by name
# ---------------------------------------------------------------------------


class TestTangentParentByName:
    def test_dot_with_a_second_s2_object_is_zero(self):
        q1, q2 = S2("q"), S2("q")
        om = q1.get_tangent_vector()
        assert full_simplify(Dot(om, q2)).is_zero
        assert full_simplify(Dot(q2, om)).is_zero

    def test_bac_cab_with_a_second_s2_object(self):
        q1, q2 = S2("q"), S2("q")
        om = q1.get_tangent_vector()
        assert str(full_simplify(Cross(q2, Cross(om, q2)))) == str(om)
        # (ω × q) × q = −ω
        assert str(full_simplify(Cross(Cross(om, q2), q2))) == str(om) + "*(-1)"

    def test_a_different_point_is_not_reduced(self):
        q, p = S2("q"), S2("p")
        om = q.get_tangent_vector()
        assert not full_simplify(Dot(om, p)).is_zero
        assert isinstance(full_simplify(Cross(p, Cross(om, p))), Cross)

    def test_pipeline_is_the_same_for_two_objects(self):
        q1, q2 = S2("q"), S2("q")
        L, dW, _ = _spherical_pendulum(q1)
        expected = compute_eom(L, dW, SystemVariables(vectors=[q1]))
        eom = compute_eom(L, dW, SystemVariables(vectors=[q2]))
        key = str(q1.get_variation_vector())
        assert str(eom[key][1]) == str(expected[key][1])
        assert "Dot(\\omega_{q},q)" not in str(eom[key][1])


# ---------------------------------------------------------------------------
# 2.6 inf_work None / 0 / wrong type
# ---------------------------------------------------------------------------


class TestInfWork:
    def test_none_and_zero_mean_no_work(self):
        L, _, x, _ = _point_mass()
        variables = SystemVariables(vectors=[x])
        ref = str(compute_eom(L, Zero, variables)[str(x.delta())][1])
        assert str(compute_eom(L, None, variables)[str(x.delta())][1]) == ref
        assert str(compute_eom(L, 0, variables)[str(x.delta())][1]) == ref
        assert str(compute_eom(L, 0.0, variables)[str(x.delta())][1]) == ref

    def test_vector_work_raises_naming_type(self):
        L, _, x, f = _point_mass()
        with pytest.raises(TypeError, match="inf_work must be a scalar expression.*VECTOR.*'f'"):
            compute_eom(L, f, SystemVariables(vectors=[x]))

    def test_non_expression_work_raises_naming_type(self):
        L, _, x, _ = _point_mass()
        with pytest.raises(TypeError, match="inf_work must be a scalar expression, got str"):
            compute_eom(L, "0", SystemVariables(vectors=[x]))
        with pytest.raises(TypeError, match="got int 3"):
            compute_eom(L, 3, SystemVariables(vectors=[x]))

    def test_non_scalar_lagrangian_raises(self):
        _, dW, x, _ = _point_mass()
        with pytest.raises(TypeError, match="lagrangian must be a scalar expression.*VECTOR"):
            compute_eom(x, dW, SystemVariables(vectors=[x]))
        with pytest.raises(
            TypeError, match="lagrangian must be a scalar expression, got NoneType"
        ):
            compute_eom(None, dW, SystemVariables(vectors=[x]))


# ---------------------------------------------------------------------------
# 2.7 manifold_rules._rebuild rejects unknown nodes
# ---------------------------------------------------------------------------


class _Opaque(VectorExpr):
    """A vector node class the rebuild table does not know."""

    def __init__(self, *nodes):
        self.nodes = list(nodes)

    def __str__(self):
        return "Opaque(" + ",".join(str(n) for n in self.nodes) + ")"


class TestRebuildUnknownNode:
    def test_unknown_node_with_a_substituted_child_raises(self):
        q = S2("q")
        om = q.get_tangent_vector()
        expr = _Opaque(Variation(om))
        with pytest.raises(TypeError, match="cannot rebuild a _Opaque node"):
            apply_manifold_rules(expr, SystemVariables(vectors=[q]))

    def test_unknown_node_without_substitution_is_returned(self):
        q = S2("q")
        expr = _Opaque(Vector("a"))
        assert apply_manifold_rules(expr, SystemVariables(vectors=[q])) is expr

    def test_known_nodes_are_rebuilt(self):
        q = S2("q")
        om, xi = q.get_tangent_vector(), q.get_variation_vector()
        out = apply_manifold_rules(Cross(Variation(om), q), SystemVariables(vectors=[q]))
        assert isinstance(out, Cross)
        assert str(TimeDerivative(xi)) in str(out)


# ---------------------------------------------------------------------------
# 2.8 calculus nodes over a manifold leaf
# ---------------------------------------------------------------------------


class TestCalculusOfManifoldLeaf:
    def test_time_derivative_of_s2_points_to_t_diff(self):
        with pytest.raises(ValueError, match=r"TimeDerivative\(q\).*use q\.t_diff\(\)"):
            TimeDerivative(S2("q"))

    def test_variation_of_s2_points_to_delta(self):
        with pytest.raises(ValueError, match=r"Variation\(q\).*use q\.delta\(\)"):
            Variation(S2("q"))

    def test_so3(self):
        with pytest.raises(ValueError, match=r"TimeDerivative\(R\).*use R\.t_diff\(\)"):
            TimeDerivative(SO3("R"))
        with pytest.raises(ValueError, match=r"Variation\(R\).*use R\.delta\(\)"):
            Variation(SO3("R"))

    def test_eager_forms_still_work(self):
        q, R = S2("q"), SO3("R")
        assert str(q.t_diff()) == "Cross(\\omega_{q},q)"
        assert str(q.delta()) == "Cross(\\xi_{q},q)"
        assert str(R.t_diff()) == "R*Hat(\\Omega_{R})"
        assert str(R.delta()) == "R*Hat(\\eta_{R})"

    def test_plain_leaves_and_tangent_vectors_unaffected(self):
        x, q = Vector("x"), S2("q")
        assert isinstance(TimeDerivative(x), TimeDerivative)
        assert isinstance(Variation(q.get_tangent_vector()), Variation)
        assert isinstance(TimeDerivative(q.get_variation_vector()), TimeDerivative)


# ---------------------------------------------------------------------------
# 2.9 collect: both-sides check first, scalar argument only
# ---------------------------------------------------------------------------


class TestCollectGuards:
    def test_non_scalar_argument_raises_naming_it(self):
        v = Vector("v")
        with pytest.raises(TypeError, match="expected a scalar expression, got Vector 'v'"):
            collect(v, v)
        with pytest.raises(TypeError, match="got Matrix 'M' of type MATRIX"):
            collect(Matrix("M"), v)

    def test_target_on_both_sides_with_the_flip_shortcut(self):
        """Dot(f(η), η): the r == vec shortcut used to flip it and report success."""
        eta, a, b = Vector("eta"), Vector("a"), Vector("b")
        quadratic = Dot(MVMul(VVMul(eta, Transpose(a)), b), eta)
        with pytest.raises(NonLinearError, match="'eta' appears on both sides"):
            collect(quadratic, eta)

    def test_target_on_both_sides_with_the_canonical_shortcut(self):
        v, w = Vector("v"), Vector("w")
        with pytest.raises(NonLinearError, match="'v' appears on both sides"):
            collect(Dot(v, Cross(v, w)), v)

    def test_linear_dots_still_collected(self):
        v, a, b = Vector("v"), Vector("a"), Vector("b")
        out = collect(Dot(a, Cross(v, b)), v)
        assert isinstance(out, Dot) and out.left == v
        out = collect(Dot(a, v) + Dot(Matrix("M") * v, b), v)
        assert all(n.left == v for n in out.nodes)
        assert collect(Dot(a, b), v) == Dot(a, b)


# ---------------------------------------------------------------------------
# 2.10 full_simplify warns at the iteration cap
# ---------------------------------------------------------------------------


class TestFullSimplifyCapWarning:
    def _two_pass_expression(self):
        """Hat(q)·(ω×q) → q×(ω×q) in pass 1 (eliminate), → ω only in pass 2 (vector rules)."""
        q = S2("q")
        om = q.get_tangent_vector()
        return Dot(MVMul(Hat(q), Cross(om, q)), q)

    def test_warns_with_the_expression_when_the_cap_is_hit(self):
        expr = self._two_pass_expression()
        with pytest.warns(RuntimeWarning, match="no fixed point after 1 passes") as record:
            out = full_simplify(expr, max_iter=1)
        assert repr(out) in str(record[0].message)
        assert not out.is_zero

    def test_no_warning_when_a_fixed_point_is_reached(self):
        expr = self._two_pass_expression()
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            assert full_simplify(expr).is_zero
            assert full_simplify(expr, max_iter=3).is_zero

    def test_cap_zero_always_warns(self):
        with pytest.warns(RuntimeWarning, match="after 0 passes"):
            full_simplify(Scalar("a"), max_iter=0)


# ---------------------------------------------------------------------------
# 2.12 example symbol names
# ---------------------------------------------------------------------------


class TestExampleSymbolNames:
    def test_double_rigid_pendulum_inertias_are_distinct(self):
        src = pathlib.Path(__file__).resolve().parents[1] / "examples" / "rigid_pendulum.py"
        text = src.read_text()
        assert 'J1 = Matrix("J1"' in text
        assert 'J2 = Matrix("J2"' in text
        assert 'J2 = Matrix("J1"' not in text


# ---------------------------------------------------------------------------
# error classes are importable from the package
# ---------------------------------------------------------------------------


class TestErrorExports:
    def test_exported_classes_are_the_engine_classes(self):
        assert AlgebraicError is errors.AlgebraicError
        assert NonLinearError is errors.NonLinearError
        assert ExpressionMismatchError is errors.ExpressionMismatchError
        assert SizeMismatchError is errors.SizeMismatchError

    def test_catchable_from_the_public_names(self):
        a, m = Scalar("a"), Vector("m")
        with pytest.raises(ExpressionMismatchError):
            Dot(a, m)
        with pytest.raises(NonLinearError):
            collect(Dot(m, Cross(m, Vector("w"))), m)
