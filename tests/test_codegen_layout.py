"""``StateLayout`` slots and the errors of code generation (spec verification item 4)."""

from __future__ import annotations

import numpy as np
import pytest

from geomech import (
    S2,
    SO3,
    Cross,
    Dot,
    Matrix,
    Number,
    Scalar,
    SystemVariables,
    TimeDerivative,
    Variation,
    Vector,
)
from geomech.codegen import (
    StateLayout,
    collect_params,
    generate_module,
    lambdify,
    load_module,
    to_source,
)
from tests.codegen_systems import VBL_SYSTEMS, derived

s, x, q, R = Scalar("s"), Vector("x"), S2("q"), SO3("R")
f, M = Scalar("f"), Vector("M")
m = Scalar("m", attr=["Constant"])
J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
e3 = Vector("e3", attr=["Constant"])
VARS = SystemVariables(scalars=[s], vectors=[x, q], matrices=[R])


def layout():
    return StateLayout(VARS, inputs=[f, M], params=[m, J, e3])


def test_slot_order_and_sizes():
    lay = layout()
    cfg = [(sl.key, sl.start, sl.stop, sl.manifold) for sl in lay.configuration]
    assert cfg == [("s", 0, 1, "R"), ("x", 1, 4, "R3"), ("q", 4, 7, "S2"), ("R", 7, 16, "SO3")]
    vel = [(sl.key, sl.start, sl.stop) for sl in lay.velocity]
    assert vel == [
        (str(s.t_diff()), 16, 17),
        (str(x.t_diff()), 17, 20),
        (str(q.get_tangent_vector()), 20, 23),
        (str(R.get_tangent_vector()), 23, 26),
    ]
    acc = [(sl.key, sl.start, sl.stop) for sl in lay.acceleration]
    assert acc == [
        (str(s.t_diff().t_diff()), 0, 1),
        (str(x.t_diff().t_diff()), 1, 4),
        (str(q.get_tangent_vector().t_diff()), 4, 7),
        (str(R.get_tangent_vector().t_diff()), 7, 10),
    ]
    assert (lay.n_config, lay.n_state, lay.n_accel, lay.n_inputs) == (16, 26, 10, 4)
    assert lay.input_slices == {"f": slice(0, 1), "M": slice(1, 4)}
    assert lay.slices["R"] == slice(7, 16)
    assert lay.param_names == ["m", "J", "e3"]


def test_identifiers():
    lay = layout()
    names = {sl.key: sl.name for sl in lay}
    assert names[str(q.get_tangent_vector())] == "omega_q"
    assert names[str(q.get_tangent_vector().t_diff())] == "omega_q_dot"
    assert names[str(R.get_tangent_vector())] == "Omega_R"
    assert names[str(x.t_diff())] == "x_dot"
    assert names[str(s.t_diff().t_diff())] == "s_ddot"
    rho = Vector("\\rho_1", attr=["Constant"])
    p = Scalar("p", attr=["Constant"])  # reserved argument name
    lay = StateLayout(VARS, params=[rho, p])
    assert lay.slot(rho).name == "rho_1" and lay.slot(p).name == "p_"


def test_rotation_is_row_major():
    lay = StateLayout(SystemVariables(matrices=[R]))
    Rv = np.arange(9.0).reshape(3, 3)
    env = {"R": Rv, str(R.get_tangent_vector()): np.ones((3, 1))}
    state = lay.state_from_env(env)
    np.testing.assert_array_equal(state[:9], [0, 1, 2, 3, 4, 5, 6, 7, 8])
    fn = lambdify(R, lay)
    np.testing.assert_array_equal(fn(state, None, {}), Rv)
    np.testing.assert_array_equal(lay.env_from_arrays(state)["R"], Rv)


def test_to_dict_is_plain_data():
    data = layout().to_dict()
    assert data["n_state"] == 26
    assert data["configuration"][3] == {
        "symbol": "R",
        "name": "R",
        "kind": "configuration",
        "manifold": "SO3",
        "shape": [3, 3],
        "slice": [7, 16],
    }
    assert data["params"][1] == {"symbol": "J", "name": "J", "kind": "param", "shape": [3, 3]}
    assert [i["symbol"] for i in data["inputs"]] == ["f", "M"]


def test_collect_params_leaves_out_literals_and_state():
    d = derived(VBL_SYSTEMS[4])  # spherical pendulum
    exprs = [d.system.L] + [eq for _, eq in d.eom.values()]
    params = collect_params(exprs, d.system.variables, d.system.inputs)
    assert [str(p) for p in params] == ["e3", "g", "l", "m"]


# ---------------------------------------------------------------------------
# Errors, each naming the symbol
# ---------------------------------------------------------------------------


def _module():
    d = derived(VBL_SYSTEMS[4])
    return load_module(generate_module(d.sf, d.layout)), d


def test_missing_parameter():
    mod, d = _module()
    p = {"e3": np.array([0, 0, 1.0]), "g": 9.8, "l": 1.0}
    with pytest.raises(KeyError, match="missing parameter 'm'"):
        mod.mass(np.zeros(6), p)


def test_extra_parameter():
    mod, d = _module()
    p = {"e3": np.array([0, 0, 1.0]), "g": 9.8, "l": 1.0, "m": 1.0, "k": 2.0}
    with pytest.raises(ValueError, match="unknown parameter 'k'"):
        mod.bias(np.zeros(6), p)


def test_wrong_state_size():
    mod, d = _module()
    p = {"e3": np.array([0, 0, 1.0]), "g": 9.8, "l": 1.0, "m": 1.0}
    with pytest.raises(ValueError, match="state must have 6 entries, got 5"):
        mod.bias(np.zeros(5), p)


def test_stray_leaf():
    stray = Vector("w", attr=["Constant"])
    with pytest.raises(ValueError, match="symbol 'w'.*not a configuration"):
        to_source(Dot(x, stray), layout())


def test_stray_leaf_in_a_block():
    d = derived(VBL_SYSTEMS[4])
    lay = StateLayout(d.system.variables, d.system.inputs, [p for p in d.layout.params][:-1])
    with pytest.raises(ValueError, match="symbol 'm'"):
        generate_module(d.sf, lay)


def test_time_derivative_of_non_leaf():
    with pytest.raises(ValueError, match="TimeDerivative of the non-leaf expression 'Cross"):
        to_source(TimeDerivative(Cross(x, e3)), layout())


def test_time_derivative_of_unknown_leaf():
    l = Scalar("l")
    lay = StateLayout(VARS, params=[l])
    with pytest.raises(ValueError, match="time derivative of the leaf 'l'"):
        to_source(l.t_diff(), lay)
    lay = StateLayout(VARS, params=[l, l.t_diff()])  # a prescribed rate is a parameter
    assert "\\\\frac{d}{dt}(l)" in to_source(l.t_diff(), lay)


def test_variation_rejected():
    with pytest.raises(ValueError, match="Variation"):
        to_source(Dot(Variation(x), e3), layout())


def test_kind_restriction():
    # a mass block may not depend on an input
    from geomech.codegen.source import function_source, get_target

    with pytest.raises(ValueError, match="contains the input 'M'"):
        function_source(
            "mass",
            ("state", "p"),
            "doc",
            layout(),
            get_target("numpy"),
            [Dot(M, e3)],
            lambda em, codes: [f"return {codes[0].text}"],
            kinds={"configuration", "velocity", "param"},
            what="mass",
        )


def test_acceleration_needs_a():
    fn = lambdify(Dot(x.t_diff().t_diff(), e3), layout())
    with pytest.raises(ValueError, match="a must have 10 entries"):
        fn(np.zeros(26), np.zeros(4), {"m": 1.0, "J": np.eye(3), "e3": np.ones(3)})


@pytest.mark.parametrize(
    "kwargs, err, match",
    [
        (dict(params=[Number(2)]), ValueError, "numeric literal"),
        (dict(params=[Dot(x, e3)]), TypeError, "parameter 'Dot"),
        (dict(params=[Scalar("s", attr=["Constant"])]), ValueError, "'s'.*prints like"),
        (dict(inputs=[q]), TypeError, "input 'q'"),
        (dict(inputs=[f, f]), ValueError, "'f'.*prints like"),
        (dict(inputs=f), TypeError, "inputs=\\[f\\]"),
        (dict(params=[Vector("v5", size=(5,))]), ValueError, "'v5' has size"),
    ],
)
def test_layout_validation(kwargs, err, match):
    with pytest.raises(err, match=match):
        StateLayout(VARS, **kwargs)


def test_layout_rejects_bad_variables():
    with pytest.raises(TypeError, match="SystemVariables"):
        StateLayout([x])
    with pytest.raises(TypeError, match="StateLayout"):
        StateLayout(SystemVariables(matrices=[Matrix("A")]))


def test_unknown_target():
    with pytest.raises(ValueError, match="unknown target 'jax'"):
        to_source(x, layout(), target="jax")


def test_generate_module_rejects_other_objects():
    with pytest.raises(TypeError, match="generate_module: system must be"):
        generate_module({"a": 1}, layout())
