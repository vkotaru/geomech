"""An input (or a prescribed parameter) printed like a different leaf of the equations.

The engine matches leaves by printed name. Before this check, ``Scalar("g")`` listed as an
input of a system whose Lagrangian holds the constant ``Scalar("g", attr=['Constant'])``
was silently half-used: ``to_standard_form`` matched it by name and returned
``G[g] = −m e₃``, while ``linearize`` found no ``δg`` (the constant has no variation) and
returned a zero ``B[δg]`` column. Now ``linearize``, ``to_standard_form`` and
``compute_eom(..., params=...)`` raise ``ValueError`` naming both symbols; an input that is
the same symbol as its leaf in the equations (same class, flags and size) stays legal.
"""

import pytest

from geomech import (
    Dot,
    Scalar,
    SystemVariables,
    Vector,
    compute_eom,
    linearize,
    to_standard_form,
)

HALF = Scalar("0.5", value=0.5)

CONSTANT_G = (
    "input 'g' has the name of the constant 'g' in the equations; a symbol is either a"
    " constant parameter or an input (rename one, or drop attr=['Constant'] if 'g' is an"
    " input)"
)


def _falling_mass(g=None):
    """``L = ½ m ẋ·ẋ − m g e₃·x`` with the force ``f`` in δW; ``g`` constant by default."""
    m = Scalar("m", attr=["Constant"])
    g = Scalar("g", attr=["Constant"]) if g is None else g
    e3 = Vector("e3", attr=["Constant"])
    x, f = Vector("x"), Vector("f")
    L = HALF * m * Dot(x.t_diff(), x.t_diff()) - m * g * Dot(e3, x)
    variables = SystemVariables(vectors=[x])
    return compute_eom(L, Dot(x.delta(), f), variables), variables, f


def _forced_mass(F_eq):
    """``L = ½ m ẋ·ẋ`` with ``δW = δx·F_eq``."""
    m = Scalar("m", attr=["Constant"])
    x = Vector("x")
    variables = SystemVariables(vectors=[x])
    L = HALF * m * Dot(x.t_diff(), x.t_diff())
    return compute_eom(L, Dot(x.delta(), F_eq), variables), variables


class TestConstantNamedLikeAnInput:
    def test_linearize_rejects_the_reproducer(self):
        eom, variables, f = _falling_mass()
        with pytest.raises(ValueError) as err:
            linearize(eom, variables, [f, Scalar("g")])
        assert str(err.value) == f"linearize: {CONSTANT_G}"

    def test_to_standard_form_rejects_the_reproducer(self):
        eom, variables, _ = _falling_mass()
        with pytest.raises(ValueError) as err:
            to_standard_form(eom, variables, inputs=[Scalar("g")])
        assert str(err.value) == f"to_standard_form: {CONSTANT_G}"

    def test_constant_vector_named_like_an_input(self):
        eom, variables = _forced_mass(Vector("F", attr=["Constant"]))
        for call, name in (
            (lambda: linearize(eom, variables, [Vector("F")]), "linearize"),
            (lambda: to_standard_form(eom, variables, inputs=[Vector("F")]), "to_standard_form"),
        ):
            with pytest.raises(ValueError, match="has the name of the constant 'F'") as err:
                call()
            assert str(err.value).startswith(f"{name}: input 'F'")


class TestOtherSymbolNamedLikeAnInput:
    def test_scalar_input_named_like_a_constant_vector(self):
        eom, variables = _forced_mass(Vector("F", attr=["Constant"]))
        expected = (
            "input 'F' (Scalar) has the name of the constant Vector 'F' in the equations, a"
            " different symbol; leaves are matched by printed name, so the two would be"
            " confused (rename one, or build the input with the class and attr it has in"
            " the equations)"
        )
        with pytest.raises(ValueError) as err:
            linearize(eom, variables, [Scalar("F")])
        assert str(err.value) == f"linearize: {expected}"
        with pytest.raises(ValueError) as err:
            to_standard_form(eom, variables, inputs=[Scalar("F")])
        assert str(err.value) == f"to_standard_form: {expected}"

    def test_input_with_other_flags_than_its_leaf(self):
        eom, variables = _forced_mass(Vector("F", attr=["UnitNorm"]))
        with pytest.raises(ValueError, match=r"\(Vector\) has the name of the Vector with"):
            to_standard_form(eom, variables, inputs=[Vector("F")])


class TestLegalInputs:
    def test_same_symbol_as_in_the_equations(self):
        """A freshly built ``Vector("f")`` is the symbol ``f`` of the equations."""
        eom, variables, _ = _falling_mass()
        sf = to_standard_form(eom, variables, inputs=[Vector("f")])
        assert str(sf["\\delta{x}"].G["f"]) == "I"
        lin = linearize(eom, variables, [Vector("f")])
        assert str(lin.B["\\delta{x}"]["\\delta{f}"]) == "I"

    def test_gravity_as_an_input(self):
        """``g`` without ``Constant``, entering through the virtual work, is an input."""
        m, g = Scalar("m", attr=["Constant"]), Scalar("g")
        e3 = Vector("e3", attr=["Constant"])
        x = Vector("x")
        variables = SystemVariables(vectors=[x])
        dW = Dot(x.delta(), e3) * m * g * -1
        eom = compute_eom(HALF * m * Dot(x.t_diff(), x.t_diff()), dW, variables)
        sf = to_standard_form(eom, variables, inputs=[Scalar("g")])
        assert str(sf["\\delta{x}"].G["g"]) == "e3*m*(-1)"
        lin = linearize(eom, variables, [Scalar("g")])
        assert str(lin.B["\\delta{x}"]["\\delta{g}"]) == "e3*m*(-1)"

    def test_constant_parameter_and_a_separate_input(self):
        eom, variables, f = _falling_mass()
        sf = to_standard_form(eom, variables, inputs=[f])
        assert set(sf["\\delta{x}"].G) == {"f"}
        lin = linearize(eom, variables, [f])
        assert set(lin.B["\\delta{x}"]) == {"\\delta{f}"}

    def test_input_absent_from_the_equations(self):
        eom, variables, f = _falling_mass()
        sf = to_standard_form(eom, variables, inputs=[f, Vector("u")])
        assert set(sf["\\delta{x}"].G) == {"f"}


class TestParamNamedLikeAConstant:
    def test_compute_eom_rejects_a_param_named_like_a_constant(self):
        m, l = Scalar("m", attr=["Constant"]), Scalar("l", attr=["Constant"])
        x = Vector("x")
        L = HALF * m * l * l * Dot(x.t_diff(), x.t_diff())
        with pytest.raises(ValueError) as err:
            compute_eom(L, None, SystemVariables(vectors=[x]), params=[Scalar("l")])
        assert str(err.value) == (
            "compute_eom: params entry 'l' has the name of the constant 'l' in the Lagrangian"
            " or the virtual work; a symbol is either a constant parameter or a prescribed"
            " parameter (rename one, or drop attr=['Constant'] if 'l' is a prescribed"
            " parameter)"
        )

    def test_params_from_a_generator_still_checked(self):
        l = Scalar("l")
        x = Vector("x")
        L = HALF * l * Dot(x.t_diff(), x.t_diff())
        eom = compute_eom(L, None, SystemVariables(vectors=[x]), params=(p for p in [l]))
        assert "\\delta{x}" in eom
