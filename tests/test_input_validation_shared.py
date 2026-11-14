"""``to_standard_form`` and ``linearize`` accept and reject exactly the same ``inputs``.

Both call ``geomech.dynamics.eom.check_inputs``. Before, ``to_standard_form`` only ran the
name-clash check: a ``Constant`` input returned a ``G[g]`` block, a duplicate input was
silently merged, and a configuration variable, a tangent or variation vector, a manifold
point or an expression was either accepted or failed later with an engine error about
the term it occurred in. Each rejection case below raises the same exception type from
both functions, with the same message after the caller prefix, naming the input.
"""

import pytest

from geomech import (
    Dot,
    Matrix,
    Scalar,
    SystemVariables,
    Vector,
    compute_eom,
    linearize,
    to_standard_form,
)
from tests import vbl_systems

CALLERS = {
    "to_standard_form": lambda eom, variables, inputs: to_standard_form(eom, variables, inputs),
    "linearize": lambda eom, variables, inputs: linearize(eom, variables, inputs),
}


def _derive(build):
    s = build()
    return s, compute_eom(s.L, s.dW, s.variables)


# (system builder, inputs from the system's symbols, exception, printed name of the culprit)
REJECTED = {
    "constant input of L": (
        vbl_systems.rigid_pendulum_pivot,
        lambda s: [s.symbols["M"], s.symbols["g"]],
        ValueError,
        "g",
    ),
    "constant input absent from the equations": (
        vbl_systems.point_mass_r3,
        lambda s: [Vector("u", attr=["Constant"])],
        ValueError,
        "u",
    ),
    "numeric input": (
        vbl_systems.point_mass_r3,
        lambda s: [Scalar("2", value=2)],
        ValueError,
        "2",
    ),
    "configuration variable (Vector)": (
        vbl_systems.point_mass_r3,
        lambda s: [s.symbols["x"]],
        ValueError,
        "x",
    ),
    "configuration variable (Scalar)": (
        vbl_systems.spring_mass_r,
        lambda s: [s.symbols["f"], s.symbols["x"]],
        ValueError,
        "x",
    ),
    "configuration variable (S2)": (
        vbl_systems.spherical_pendulum_s2,
        lambda s: [s.symbols["q"]],
        TypeError,
        "q",
    ),
    "configuration variable (SO3)": (
        vbl_systems.rigid_pendulum_pivot,
        lambda s: [s.symbols["R"]],
        TypeError,
        "R",
    ),
    "tangent vector (TS2)": (
        vbl_systems.spherical_pendulum_s2,
        lambda s: [s.symbols["q"].get_tangent_vector()],
        ValueError,
        "\\omega_{q}",
    ),
    "variation vector (TS2)": (
        vbl_systems.spherical_pendulum_s2,
        lambda s: [s.symbols["q"].get_variation_vector()],
        ValueError,
        "\\xi_{q}",
    ),
    "tangent vector (TSO3)": (
        vbl_systems.rigid_pendulum_pivot,
        lambda s: [s.symbols["R"].get_tangent_vector()],
        ValueError,
        "\\Omega_{R}",
    ),
    "variation vector (TSO3)": (
        vbl_systems.rigid_pendulum_pivot,
        lambda s: [s.symbols["R"].get_variation_vector()],
        ValueError,
        "\\eta_{R}",
    ),
    "plain Vector printed like a tangent vector": (
        vbl_systems.rigid_pendulum_pivot,
        lambda s: [s.symbols["M"], Vector(str(s.symbols["R"].get_tangent_vector()))],
        ValueError,
        "\\Omega_{R}",
    ),
    "non-leaf: time derivative": (
        vbl_systems.point_mass_r3,
        lambda s: [s.symbols["x"].t_diff()],
        TypeError,
        "\\frac{d}{dt}(x)",
    ),
    "non-leaf: sum": (
        vbl_systems.point_mass_r3,
        lambda s: [s.symbols["F"] + Vector("u")],
        TypeError,
        "(F+u)",
    ),
    "matrix input": (
        vbl_systems.point_mass_r3,
        lambda s: [Matrix("K")],
        TypeError,
        "K",
    ),
    "duplicate": (
        vbl_systems.point_mass_r3,
        lambda s: [s.symbols["F"], s.symbols["F"]],
        ValueError,
        "F",
    ),
    "name clash with a constant of L": (
        vbl_systems.rigid_pendulum_pivot,
        lambda s: [s.symbols["M"], Scalar("g")],
        ValueError,
        "g",
    ),
    "name clash with another class": (
        vbl_systems.point_mass_r3,
        lambda s: [Scalar("F")],
        ValueError,
        "F",
    ),
    "one leaf instead of a sequence": (
        vbl_systems.point_mass_r3,
        lambda s: s.symbols["F"],
        TypeError,
        "F",
    ),
}


@pytest.mark.parametrize("caller", sorted(CALLERS))
@pytest.mark.parametrize("case", sorted(REJECTED))
def test_rejected_input(caller, case):
    build, make_inputs, exc, name = REJECTED[case]
    s, eom = _derive(build)
    with pytest.raises(exc) as err:
        CALLERS[caller](eom, s.variables, make_inputs(s))
    assert type(err.value) is exc
    msg = str(err.value)
    assert msg.startswith(f"{caller}: ")
    assert f"'{name}'" in msg


@pytest.mark.parametrize("case", sorted(REJECTED))
def test_both_functions_give_the_same_message(case):
    build, make_inputs, _, _ = REJECTED[case]
    s, eom = _derive(build)
    messages = {}
    for caller, call in CALLERS.items():
        with pytest.raises((TypeError, ValueError)) as err:
            call(eom, s.variables, make_inputs(s))
        messages[caller] = (type(err.value), str(err.value).removeprefix(f"{caller}: "))
    assert messages["to_standard_form"] == messages["linearize"]


def test_constant_input_message():
    s, eom = _derive(vbl_systems.rigid_pendulum_pivot)
    g = Scalar("g", attr=["Constant"])
    with pytest.raises(ValueError) as err:
        to_standard_form(eom, s.variables, inputs=[g])
    assert str(err.value) == (
        "to_standard_form: input 'g' is constant: its variation vanishes, so it cannot be an input"
    )


def test_not_a_sequence_message():
    s, eom = _derive(vbl_systems.point_mass_r3)
    with pytest.raises(TypeError) as err:
        linearize(eom, s.variables, s.symbols["F"])
    assert str(err.value) == (
        "linearize: inputs must be a sequence of Scalar/Vector leaves, got Vector 'F';"
        " write inputs=[F]"
    )


# ---------------------------------------------------------------------------
# Legal inputs: accepted by both, with the same input set
# ---------------------------------------------------------------------------


def _input_names(eom, variables, inputs_factory):
    sf = to_standard_form(eom, variables, inputs_factory())
    lin = linearize(eom, variables, inputs_factory())
    g_names = set().union(*(set(eq.G) for eq in sf.values()))
    b_names = {str(u) for u in lin.inputs}
    return g_names, b_names


@pytest.mark.parametrize(
    "kind",
    ["list", "tuple", "generator"],
)
def test_legal_inputs_of_the_quadrotor(kind):
    s, eom = _derive(vbl_systems.quadrotor)
    f, M = s.symbols["f"], s.symbols["M"]
    factory = {
        "list": lambda: [f, M],
        "tuple": lambda: (f, M),
        "generator": lambda: (u for u in [f, M]),
    }[kind]
    g_names, b_names = _input_names(eom, s.variables, factory)
    assert g_names == {"f", "M"}
    assert b_names == {"\\delta{f}", "\\delta{M}"}


@pytest.mark.parametrize(
    "build",
    [
        vbl_systems.point_mass_r3,
        vbl_systems.spring_mass_r,
        vbl_systems.free_rigid_body,
        vbl_systems.spherical_pendulum_s2,
        vbl_systems.quadrotor,
    ],
    ids=lambda f: f.__name__,
)
def test_system_inputs_are_legal_for_both(build):
    s, eom = _derive(build)
    g_names, b_names = _input_names(eom, s.variables, lambda: list(s.inputs))
    assert g_names == {str(u) for u in s.inputs}
    assert b_names == {str(u.delta()) for u in s.inputs}


def test_none_and_empty_list_mean_no_inputs():
    m = Scalar("m", attr=["Constant"])
    x = Vector("x")
    L = Scalar("0.5", value=0.5) * m * Dot(x.t_diff(), x.t_diff())
    variables = SystemVariables(vectors=[x])
    eom = compute_eom(L, None, variables)
    for inputs in (None, []):
        sf = to_standard_form(eom, variables, inputs)
        assert sf["\\delta{x}"].G == {}
        lin = linearize(eom, variables, inputs)
        assert lin.inputs == []
        assert lin.B["\\delta{x}"] == {}


def test_input_absent_from_the_equations_is_legal_for_both():
    s, eom = _derive(vbl_systems.point_mass_r3)
    u = Vector("u")
    sf = to_standard_form(eom, s.variables, [s.symbols["F"], u])
    assert set(sf["\\delta{x}"].G) == {"F"}
    lin = linearize(eom, s.variables, [s.symbols["F"], u])
    assert [str(v) for v in lin.inputs] == ["\\delta{F}", "\\delta{u}"]
    assert lin.B["\\delta{x}"]["\\delta{u}"].is_zero
