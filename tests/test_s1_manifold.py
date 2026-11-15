"""The one-sphere S¹ as a configuration manifold (Lee et al. ch. 4).

An ``S1`` point is the unit circle in the plane normal to a constant unit axis ``e``:
``q̇ = ω e × q`` and ``δq = γ e × q`` with scalar ω, γ, and ``δω = γ̇`` (eq. 4.11). In that
plane ``e × ·`` is the book's rotation ``S`` by π/2, so the derived equations are those of
eq. 4.12. The tests below cover, in order: the leaf and its validation, the manifold rule,
the S¹ vector identities (each against a numeric evaluation on the circle), the derived
equations of §4.6 (pinned strings, compared with the book's), the validation of
``compute_eom``, the standard form and VBL, the numeric oracle, code generation,
simulation, LQR, and the spacecraft with wheels on SO(3) × (S¹)³ against its R³ twin.
"""

from __future__ import annotations

import numpy as np
import pytest

from geomech import (
    S1,
    S2,
    TS1,
    Cross,
    Dot,
    EOMTrace,
    Number,
    Scalar,
    SVMul,
    SystemVariables,
    TimeDerivative,
    Variation,
    Vector,
    check_eom,
    compute_eom,
    linearize,
    to_standard_form,
)
from geomech.codegen import StateLayout, collect_params, load_module, module_source
from geomech.core.transformations import full_simplify
from geomech.core.transformations.manifold_rules import apply_manifold_rules
from geomech.core.transformations.vector_identities import vector_rules
from geomech.numeric import Evaluator, expm_so3, finite_difference_check, sample_state
from geomech.simulation import check_state, simulate
from geomech.systems import get
from geomech.utils.invariants import check_tree
from geomech.utils.printing import tree_str
from geomech.viz.latex import latex_names
from tests import s1_systems

N = Vector("n", attr=["Constant", "UnitNorm"])
Q = S1("q", axis=N)
P = S1("p", axis=N)  # a second point of the same circle
F = S1("f", axis=Vector("m", attr=["Constant", "UnitNorm"]))  # a point of another circle
A = Vector("a")


def _env(seed=0):
    """Values on the circle: ``q``, ``p`` ⊥ ``n``, unit; ``f`` ⊥ ``m``; random ``a``."""
    rng = np.random.default_rng(seed)
    variables = SystemVariables(vectors=[Q, P, F])
    return sample_state(variables, [], rng)


def _value(expr, env):
    return np.asarray(Evaluator(env=dict(env), rng=np.random.default_rng(1))(expr), float)


def _same_value(a, b, seeds=(0, 1, 2)):
    for seed in seeds:
        env = _env(seed)
        assert np.allclose(_value(a, env), _value(b, env), atol=1e-12), (a, b)


# ---------------------------------------------------------------------------
# The leaf
# ---------------------------------------------------------------------------


def test_kinematics_are_eager_and_scalar():
    w, g = Q.get_tangent_vector(), Q.get_variation_vector()
    assert isinstance(w, TS1) and isinstance(g, TS1) and w.S1 is Q and g.S1 is Q
    assert (str(w), str(g)) == ("\\omega_{q}", "\\gamma_{q}")
    assert Q.t_diff() == SVMul(Cross(N, Q), w)
    assert Q.delta() == SVMul(Cross(N, Q), g)
    assert Q.is_manifold and Q.is_unit_norm and Q.axis is N
    assert repr(Q) == "S1:q"
    assert "axis='n'" in tree_str(Q, style="indent")


@pytest.mark.parametrize(
    "axis, error, match",
    [
        (None, TypeError, "an axis is required"),
        (Scalar("n", attr=["Constant"]), TypeError, "must be a Vector leaf"),
        (S2("n"), TypeError, "must be a Vector leaf"),
        (Cross(N, A), TypeError, "must be a Vector leaf"),
        (Vector("n", attr=["UnitNorm"]), ValueError, "nonzero constant"),
        (Vector("n", attr=["Constant"]), ValueError, "UnitNorm"),
        (Vector("n", size=(2,), attr=["Constant", "UnitNorm"]), ValueError, "size"),
        (Vector("q", attr=["Constant", "UnitNorm"]), ValueError, "same name"),
    ],
)
def test_axis_must_be_a_constant_unit_vector_leaf(axis, error, match):
    with pytest.raises(error, match=match):
        S1("q", axis=axis)


def test_calculus_over_the_point_and_variation_of_gamma_are_rejected():
    with pytest.raises(ValueError, match=r"use q.t_diff\(\)"):
        TimeDerivative(Q)
    with pytest.raises(ValueError, match=r"use q.delta\(\)"):
        Variation(Q)
    with pytest.raises(ValueError, match="variation of variation vector"):
        Variation(Q.get_variation_vector())
    check_tree(Q.t_diff().t_diff())
    check_tree(Q.delta())


def test_latex_names_index_omega_and_gamma():
    q1 = S1("q1", axis=N)
    names = latex_names([q1])
    assert names["\\omega_{q1}"] == "\\omega_1"
    assert names["\\gamma_{q1}"] == "\\gamma_1"


# ---------------------------------------------------------------------------
# Manifold rule δω = γ̇ (Lee et al. eq. 4.11)
# ---------------------------------------------------------------------------


def test_variation_of_omega_is_gamma_dot():
    w, g = Q.get_tangent_vector(), Q.get_variation_vector()
    V = SystemVariables(vectors=[Q])
    out = apply_manifold_rules(Variation(w) * A.dot(A), V)
    assert str(out) == str(TimeDerivative(g) * A.dot(A))
    # an unlisted point keeps its δω
    assert apply_manifold_rules(Variation(P.get_tangent_vector()), V) == Variation(
        P.get_tangent_vector()
    )


def test_variation_of_omega_matches_the_variation_of_the_velocity():
    """δ(q̇) computed by the product rule, with δω → γ̇, equals d/dt(δq) (both on the
    circle): the commutation that eq. 4.11 encodes."""
    V = SystemVariables(vectors=[Q])
    lhs = full_simplify(apply_manifold_rules(Q.t_diff().delta(), V))
    rhs = full_simplify(Q.delta().t_diff())
    assert lhs != rhs  # different trees (factor order) ...
    _same_value(lhs, rhs)  # ... with one value


# ---------------------------------------------------------------------------
# S¹ vector identities
# ---------------------------------------------------------------------------

IDENTITIES = [
    (Dot(N, Q), "0"),
    (Dot(Q, N), "0"),
    (Cross(N, Cross(N, Q)), "q*(-1)"),
    (Cross(Cross(Q, N), N), "q*(-1)"),
    (Cross(N, Cross(Q, N)), "q"),
    (Cross(Cross(N, Q), N), "q"),
    (Cross(Q, Cross(N, Q)), "n"),
    (Cross(Cross(Q, N), Q), "n"),
    (Cross(Q, Cross(Q, N)), "n*(-1)"),
    (Cross(Cross(N, Q), Q), "n*(-1)"),
    (Dot(Cross(N, Q), Cross(N, Q)), "1"),
    (Dot(Cross(N, Q), Cross(N, P)), "Dot(q,p)"),
    (Dot(Cross(Q, N), Cross(N, P)), "Dot(q,p)*(-1)"),
    (Dot(Cross(Q, N), Cross(P, N)), "Dot(q,p)"),
    (Dot(Cross(N, Q), P), "Dot(p,Cross(n,q))"),
    (Dot(P, Cross(N, Q)), "Dot(p,Cross(n,q))"),
    (Dot(Q, Cross(N, P)), "Dot(p,Cross(n,q))*(-1)"),
]


@pytest.mark.parametrize("expr, expected", IDENTITIES, ids=[str(e) for e, _ in IDENTITIES])
def test_s1_identity(expr, expected):
    out = full_simplify(expr)
    assert str(out) == expected
    _same_value(expr, out)


NEAR_MISSES = [
    Dot(F, N),  # f lies on the circle normal to m, not n
    Cross(N, Cross(N, F)),
    Cross(N, Cross(N, A)),  # a generic vector
    Cross(Vector("u", attr=["UnitNorm"]), Cross(Vector("u", attr=["UnitNorm"]), Q)),
    Dot(Cross(N, Q), Cross(F.axis, F)),  # two circles
    Dot(Cross(N, A), Cross(N, Q)),
    Cross(Q, Cross(N, P)),
]


@pytest.mark.parametrize("expr", NEAR_MISSES, ids=[str(e) for e in NEAR_MISSES])
def test_s1_identities_need_their_hypotheses(expr):
    assert vector_rules(expr) == expr


# ---------------------------------------------------------------------------
# Derived equations of Lee et al. §4.6
# ---------------------------------------------------------------------------


def _derive(build):
    s = build()
    return s, compute_eom(s.L, s.dW, s.variables)


def _eq(eom, point):
    return eom[str(point.get_variation_vector())][1]


def test_planar_pendulum_qdot_form_eq_4_33():
    """``L = ½ m l² q̇·q̇ − m g l e₃·q``: the S¹ identities reduce ``q̇·q̇`` to ``ω²`` and the
    equation is that of the (q, ω) form, eq. 4.33: ``m l² ω̇ + m g l e₁ᵀq = 0`` with the
    book's ``e₁ = −S e₂ = −n × e₃`` (vertical e₃): ``E = −m l² ω̇ − m g l e₃·(n × q) + τ``."""
    s, eom = _derive(s1_systems.planar_pendulum_qdot)
    q = s.symbols["q"]
    assert str(_eq(eom, q)) == (
        "(m*l*l*\\frac{d}{dt}(\\omega_{q})*(-1)+m*g*l*Dot(e3,Cross(n,q))*(-1)+\\tau)"
    )


def test_rotating_hoop_eq_4_40():
    """Eq. 4.40: ``m L² ω̇ + m L² Ω² qᵀe₂ e₁ᵀq = 0`` with ``e₂ = S e₁ = n × e₁``; since
    ``qᵀ(n × e₁) = −e₁·(n × q)``, ``E = −m L² ω̇ + m L² Ω² e₁·(n × q) e₁·q + τ``."""
    s, eom = _derive(s1_systems.rotating_hoop)
    q = s.symbols["q"]
    eq = _eq(eom, q)
    assert str(eq) == (
        "(m*L*L*\\frac{d}{dt}(\\omega_{q})*(-1)+m*L*L*\\Omega*\\Omega*Dot(e1,Cross(n,q))"
        "*Dot(e1,q)+\\tau)"
    )
    rng = np.random.default_rng(3)
    for _ in range(5):
        env = sample_state(s.variables, s.inputs, rng)
        ev = Evaluator(env=env, rng=rng)
        val = float(ev(eq))
        m, L, Om = (ev(Scalar(k, attr=["Constant"])) for k in ("m", "L", "\\Omega"))
        qv, e1, n = (np.asarray(ev(x)).reshape(3) for x in (q, s.symbols["e1"], Vector("n")))
        e2 = np.cross(n, e1)
        book = m * L * L * ev(q.get_tangent_vector().t_diff()) + m * L * L * Om**2 * (qv @ e2) * (
            e1 @ qv
        )
        assert val == pytest.approx(-book + ev(Scalar("\\tau")), abs=1e-12)


def test_elastic_planar_pendulums_eqs_4_49_4_50():
    """Eqs. 4.49–4.50: ``m L² ω̇₁ + κ q₁ᵀSq₂ = 0``, ``m L² ω̇₂ + κ q₂ᵀSq₁ = 0``; with
    ``q₂ᵀSq₁ = −q₁ᵀSq₂`` and ``q₁ᵀSq₂ = q₁·(n × q₂)`` both equations print the one triple
    product ``Dot(q_1,Cross(n,q_2))`` (``vector.s1_planar_triple``)."""
    s, eom = _derive(s1_systems.elastic_planar_pendulums)
    q1, q2 = s.symbols["q1"], s.symbols["q2"]
    assert str(_eq(eom, q1)) == (
        "(m*L*L*\\frac{d}{dt}(\\omega_{q_1})*(-1)+\\kappa*Dot(q_1,Cross(n,q_2))*(-1)+\\tau_1)"
    )
    assert str(_eq(eom, q2)) == (
        "(m*L*L*\\frac{d}{dt}(\\omega_{q_2})*(-1)+\\kappa*Dot(q_1,Cross(n,q_2))+\\tau_2)"
    )


def test_closed_one_form_exerts_no_force():
    """``L = ½ J ω² + c ω a·q − k e₃·q``: the ``c`` terms cancel after integration by
    parts with a non-constant prefactor (``ibp.scalar_factor`` and the product rule)."""
    s, eom = _derive(s1_systems.gyroscopic_s1)
    eq = _eq(eom, s.symbols["q"])
    assert str(eq) == "(J*\\frac{d}{dt}(\\omega_{q})*(-1)+k*Dot(e3,Cross(n,q))*(-1)+\\tau)"


@pytest.mark.parametrize("build", s1_systems.SYSTEMS, ids=lambda f: f.__name__)
def test_standard_form_and_linearization_agree_with_the_equations(build):
    s = build()
    eom = compute_eom(s.L, s.dW, s.variables)
    assert check_eom(s.L, s.dW, s.variables, eom, inputs=s.inputs) < 1e-5
    sf = to_standard_form(eom, s.variables, s.inputs)
    rng = np.random.default_rng(5)
    for _ in range(3):
        env = sample_state(s.variables, s.inputs, rng)
        ev = Evaluator(env=env, rng=rng)
        for key, (_, eq) in eom.items():
            total = np.atleast_1d(np.asarray(ev(sf[key].f), float)).reshape(-1)
            for col, coeff in [*sf[key].M.items(), *sf[key].G.items()]:
                c = np.asarray(ev(coeff), float)
                x = np.asarray(ev.env[col], float).reshape(-1)
                if x.size == 1:  # a scalar column (an S¹ ω̇, a scalar input)
                    term = c.reshape(-1) * x[0]
                elif c.size == 3:  # a scalar row, a vector column: Dot(c, x)
                    term = np.atleast_1d(c.reshape(-1) @ x)
                else:
                    term = c @ x
                total = total + term
            assert np.allclose(total, np.asarray(ev(eq), float).reshape(-1), atol=1e-10)
    lin = linearize(eom, s.variables, s.inputs)
    assert finite_difference_check(lin, eom, s.variables, s.inputs, trials=2) < 1e-6


def test_planar_pendulum_standard_form_and_kinematic_row():
    s = get("planar_pendulum")
    q = s.symbols["q"]
    eom = s.derive()
    key = str(q.get_variation_vector())
    sf = to_standard_form(eom, s.variables, s.inputs)[key]
    acc = str(q.get_tangent_vector().t_diff())
    assert str(sf.M[acc]) == "m*l*l*(-1)"
    assert str(sf.G["\\tau"]) == "1"
    lin = linearize(eom, s.variables, s.inputs)
    gamma, dw = str(q.get_variation_vector()), str(Variation(q.get_tangent_vector()))
    assert [str(x) for x in lin.states] == [gamma, dw]
    row = "kin:q"
    assert str(lin.M[row][str(TimeDerivative(q.get_variation_vector()))]) == "1"
    assert str(lin.A[row][dw]) == "1" and str(lin.A[row][gamma]) == "0"
    assert str(lin.M[key][str(TimeDerivative(Variation(q.get_tangent_vector())))]) == "m*l*l"


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def test_compute_eom_rejects_misused_s1_symbols():
    w, g = Q.get_tangent_vector(), Q.get_variation_vector()
    L = Number(0.5) * w * w
    with pytest.raises(TypeError, match="expected ExprType.SCALAR"):
        compute_eom(L, None, SystemVariables(scalars=[Q]))
    with pytest.raises(TypeError, match="tangent or variation vector"):
        compute_eom(L, None, SystemVariables(scalars=[w]))
    with pytest.raises(TypeError, match="tangent or variation vector"):
        compute_eom(L, None, SystemVariables(vectors=[Q]), params=[w])
    with pytest.raises(TypeError, match="manifold point"):
        compute_eom(L, None, SystemVariables(vectors=[P]), params=[Q])
    # a forgotten S¹ point: its γ is left in δS and named
    wp = P.get_tangent_vector()
    L2 = Number(0.5) * wp * wp + Dot(Q, Vector("c", attr=["Constant"]))
    with pytest.raises(ValueError, match="gamma_\\{q\\}"):
        compute_eom(L2, P.get_variation_vector() * Scalar("u"), SystemVariables(vectors=[P]))
    assert g.type.name == "SCALAR"


def test_inputs_reject_s1_points_and_their_coordinates():
    s = get("planar_pendulum")
    eom = s.derive()
    q = s.symbols["q"]
    with pytest.raises(TypeError, match="Scalar or Vector leaf"):
        to_standard_form(eom, s.variables, [q])
    with pytest.raises(ValueError, match="tangent or variation vector"):
        to_standard_form(eom, s.variables, [q.get_tangent_vector()])


def test_trace_names_the_s1_stages():
    s = get("planar_pendulum")
    tr = EOMTrace()
    compute_eom(s.L, s.dW, s.variables, trace=tr)
    assert "ibp.collect[\\frac{d}{dt}(\\gamma_{q})]" in tr.stage_names()
    assert "\\frac{d}{dt}(\\gamma_{q})" in str(tr["manifold_rules"])


# ---------------------------------------------------------------------------
# Numeric oracle
# ---------------------------------------------------------------------------


def test_sample_state_puts_s1_points_on_their_circle():
    s = get("double_planar_pendulum")
    env = sample_state(s.variables, s.inputs, np.random.default_rng(0), env=dict(s.params))
    for q in (s.symbols["q1"], s.symbols["q2"]):
        qv = np.asarray(env[str(q)]).reshape(3)
        assert abs(np.linalg.norm(qv) - 1) < 1e-12 and abs(qv @ s.params["e2"]) < 1e-12
        assert np.isscalar(env[str(q.get_tangent_vector())])


def test_check_eom_detects_a_wrong_s1_equation():
    s = get("planar_pendulum")
    eom = s.derive()
    key, (vec, eq) = next(iter(eom.items()))
    q = s.symbols["q"]
    wrong = {
        key: (vec, eq + Dot(Cross(q.axis, q), s.symbols["e3"]) * Scalar("m", attr=["Constant"]))
    }
    assert check_eom(s.L, s.dW, s.variables, wrong, inputs=s.inputs) > 1e-3


# ---------------------------------------------------------------------------
# Code generation, simulation and LQR
# ---------------------------------------------------------------------------


def _load(system, layout, **kwargs):
    return load_module(module_source(system, layout, **kwargs))


def _modules(name):
    s = get(name)
    eom = s.derive()
    layout = s.layout()
    sf = _load(to_standard_form(eom, s.variables, s.inputs), layout, lagrangian=s.L)
    lin = _load(linearize(eom, s.variables, s.inputs), layout)
    return s, layout, sf, lin


def test_layout_has_a_scalar_velocity_and_the_axis_as_a_parameter():
    s, layout, sf, _ = _modules("planar_pendulum")
    q = layout.configuration[0]
    w = layout.velocity[0]
    assert (q.manifold, q.size, q.axis) == ("S1", 3, "e2")
    assert (w.name, w.size, layout.acceleration[0].name) == ("omega_q", 1, "omega_q_dot")
    assert layout.to_dict()["configuration"][0]["axis"] == "e2"
    assert sf.LAYOUT["configuration"][0]["axis"] == "e2"
    # a free rotor: L holds no axis, but simulate needs it
    rotor = S1("r", axis=N)
    w_r = rotor.get_tangent_vector()
    L = Number(0.5) * Scalar("J", attr=["Constant"]) * w_r * w_r
    params = collect_params([L], SystemVariables(vectors=[rotor]))
    assert [str(p) for p in params] == ["J", "n"]
    # kinematics(state, p) is ω e × q
    x = s.x0
    qv, wv = x[:3], x[3]
    assert np.allclose(sf.kinematics(x, dict(s.params)), wv * np.cross(s.params["e2"], qv))
    assert StateLayout(s.variables, s.inputs).n_state == 4


def _theta_reference(theta0, dtheta0, g_over_l, t):
    from scipy.integrate import solve_ivp

    sol = solve_ivp(
        lambda _, y: [y[1], -g_over_l * np.sin(y[0])],
        (0.0, t[-1]),
        [theta0, dtheta0],
        t_eval=t,
        rtol=1e-12,
        atol=1e-12,
        method="DOP853",
    )
    return sol.y[0]


@pytest.mark.parametrize("method", ["rkmk4", "solve_ivp"])
def test_planar_pendulum_simulates_like_the_angle_ode(method):
    """``q = (sin θ, 0, −cos θ)`` and ``ω = −θ̇`` about e₂: the S¹ simulation matches
    ``θ̈ = −(g/l) sin θ``; ``q`` stays a unit vector orthogonal to e₂ and the energy is
    conserved."""
    pytest.importorskip("scipy")
    s, _, sf, _ = _modules("planar_pendulum")
    p = dict(s.params)
    traj = simulate(sf, s.x0, 3.0, p, method=method, dt=1e-3)
    q = traj["q"]
    theta = np.arctan2(q[:, 0], -q[:, 2])
    theta0 = float(np.arctan2(s.x0[0], -s.x0[2]))
    ref = _theta_reference(theta0, -s.x0[3], p["g"] / p["l"], traj.t)
    assert np.max(np.abs(np.angle(np.exp(1j * (theta - ref))))) < 1e-8
    assert np.max(np.abs(np.linalg.norm(q, axis=1) - 1)) < 1e-12
    assert np.max(np.abs(q @ p["e2"])) < 1e-12
    E = traj.energy()
    assert np.max(np.abs(E - E[0])) < 1e-9 * max(1.0, abs(E[0]))
    if method == "solve_ivp":
        assert "q.e2" in traj.constraint_drift


def test_check_state_rejects_points_off_the_circle():
    s, _, sf, _ = _modules("planar_pendulum")
    x = np.array(s.x0, dtype=float)
    x[:3] = [np.sin(1.0), 0.0, -np.cos(1.0)] + np.array([0.0, 1e-3, 0.0])
    x[:3] /= np.linalg.norm(x[:3])
    check_state(sf, x)  # without p the axis is unknown: only |q| = 1 is checked
    with pytest.raises(ValueError, match="orthogonal to its axis 'e2'"):
        check_state(sf, x, p=dict(s.params))
    x[:3] *= 1.1
    with pytest.raises(ValueError, match="unit norm"):
        check_state(sf, x)


def test_lqr_stabilises_the_inverted_planar_pendulum():
    pytest.importorskip("scipy")
    from geomech.simulation import lqr

    s, _, sf, lin = _modules("planar_pendulum")
    p = dict(s.params)
    x_eq = np.r_[0.0, 0.0, 1.0, 0.0]  # upright: q = +e₃, ω = 0
    ctrl = lqr(lin, x_eq, [0.0], p, [100.0, 10.0], 1.0, module=sf)
    assert ctrl.K.shape == (1, 2) and ctrl.labels == ["gamma_q", "delta_omega_q"]
    assert np.max(np.linalg.eigvals(ctrl.A).real) > 1.0
    assert np.max(ctrl.closed_loop_eigenvalues.real) < 0
    tilt = 0.3
    x0 = np.r_[expm_so3(tilt * p["e2"]) @ [0.0, 0.0, 1.0], 0.0]
    assert ctrl.error(x0)[0] == pytest.approx(tilt)
    traj = simulate(sf, x0, 6.0, p, u=ctrl)
    assert abs(ctrl.error(traj.state[-1])[0]) < 1e-3


# ---------------------------------------------------------------------------
# Spacecraft with reaction wheels: SO(3) × (S¹)³ against SO(3) × R³
# ---------------------------------------------------------------------------


def _renamed(text):
    """The printed equations of the R³ model in the names of the S¹ model."""
    for i in (1, 2, 3):
        th = f"\\theta_{i}"
        text = text.replace(
            f"\\frac{{d}}{{dt}}(\\frac{{d}}{{dt}}({th}))", f"\\frac{{d}}{{dt}}(\\omega_{{q_{i}}})"
        )
        text = text.replace(f"\\frac{{d}}{{dt}}({th})", f"\\omega_{{q_{i}}}")
        text = text.replace(f"\\delta{{{th}}}", f"\\gamma_{{q_{i}}}")
    return text


def test_spacecraft_s1_equations_are_the_r3_equations():
    """With ``θ̇ᵢ → ωᵢ``, ``θ̈ᵢ → ω̇ᵢ`` and ``δθᵢ → γᵢ`` the two derivations print the same
    equations, and the standard forms agree numerically at random states."""
    r3, s1 = get("spacecraft_reaction_wheels"), get("spacecraft_reaction_wheels_s1")
    e_r3, e_s1 = r3.derive(), s1.derive()
    assert {_renamed(k): _renamed(str(e)) for k, (_, e) in e_r3.items()} == {
        k: str(e) for k, (_, e) in e_s1.items()
    }
    sf_r3 = to_standard_form(e_r3, r3.variables, r3.inputs)
    sf_s1 = to_standard_form(e_s1, s1.variables, s1.inputs)
    rng = np.random.default_rng(4)
    for _ in range(5):
        env = sample_state(
            s1.variables,
            s1.inputs,
            rng,
            env={
                k: np.asarray(v, float).reshape(3, 1) if np.ndim(v) == 1 else v
                for k, v in s1.params.items()
            },
        )
        for i, th in enumerate(r3.symbols["theta"], 1):
            env[str(th.t_diff())] = env[f"\\omega_{{q_{i}}}"]
            env[str(th.t_diff().t_diff())] = env[f"\\frac{{d}}{{dt}}(\\omega_{{q_{i}}})"]
        ev = Evaluator(env=env, rng=rng)
        for k_r3, k_s1 in zip(sf_r3, sf_s1):
            assert _renamed(k_r3) == k_s1
            a, b = sf_r3[k_r3], sf_s1[k_s1]
            assert np.allclose(ev(a.f), ev(b.f), atol=1e-12)
            for col in a.M:
                assert np.allclose(ev(a.M[col]), ev(b.M[_renamed(col)]), atol=1e-12)
            for col in a.G:
                assert np.allclose(ev(a.G[col]), ev(b.G[col]), atol=1e-12)


def test_spacecraft_s1_simulates_like_the_r3_model():
    """From the same attitude, body rate and wheel rates, with constant motor torques, the
    two models give the same R, Ω and wheel rates, and each S¹ wheel direction is its
    initial direction turned about bᵢ by θᵢ(t) − θᵢ(0) of the R³ model."""
    r3, s1 = get("spacecraft_reaction_wheels"), get("spacecraft_reaction_wheels_s1")
    mods = {}
    for s in (r3, s1):
        eom = s.derive()
        mods[s.name] = _load(to_standard_form(eom, s.variables, s.inputs), s.layout())
    p = dict(s1.params)
    u = [0.01, -0.02, 0.015]
    t1 = simulate(mods[r3.name], r3.x0, 2.0, p, u=u, dt=1e-3)
    t2 = simulate(mods[s1.name], s1.x0, 2.0, p, u=u, dt=1e-3)
    assert np.allclose(t1["R"], t2["R"], atol=1e-9)
    assert np.allclose(t1["Omega_R"], t2["Omega_R"], atol=1e-9)
    for i in (1, 2, 3):
        b = p[f"b{i}"]
        assert np.allclose(t1[f"theta_{i}_dot"], t2[f"omega_q_{i}"], atol=1e-8)
        th = t1[f"theta_{i}"] - t1[f"theta_{i}"][0]
        q0 = t2[f"q_{i}"][0]
        expected = np.array([expm_so3(a * b) @ q0 for a in th])
        assert np.allclose(t2[f"q_{i}"], expected, atol=1e-8)
        assert np.max(np.abs(t2[f"q_{i}"] @ b)) < 1e-12


@pytest.mark.parametrize(
    "name",
    [
        "planar_pendulum",
        "double_planar_pendulum",
        "spacecraft_reaction_wheels_s1",
        "spacecraft_solar_arrays",
    ],
)
def test_x0_lies_on_the_circles(name):
    s = get(name)
    vals = s.layout().unpack(s.x0)
    for slot in s.layout().configuration:
        if slot.manifold == "S1":
            q, e = vals[slot.name], s.params[slot.axis]
            assert abs(np.linalg.norm(q) - 1) < 1e-12 and abs(q @ e) < 1e-12
            assert abs(vals[f"omega_{slot.name}"]) > 0.1


# ---------------------------------------------------------------------------
# Generated code of every S¹ system against the symbolic blocks, and CasADi against NumPy
# ---------------------------------------------------------------------------


def _rows(d):
    """``[(equation key, row slice)]`` in layout order; an S¹ equation has one row."""
    out, pos = [], 0
    for var in d.layout.config_variables:
        key = str(var.get_variation_vector() if hasattr(var, "manifold_info") else var.delta())
        n = 1 if d.eom[key][1].type.name == "SCALAR" else 3
        out.append((key, slice(pos, pos + n)))
        pos += n
    return out


@pytest.mark.parametrize("build", s1_systems.SYSTEMS, ids=lambda f: f.__name__)
def test_generated_blocks_match_the_symbolic_ones(build):
    from tests.codegen_systems import assert_close, derived, modules, reference, sample_points

    d = derived(build)
    sf_mod, eom_mod, lin_mod = modules(build)
    acc = {s.key: s.slice for s in d.layout.acceleration}
    inp = {s.key: s.slice for s in d.layout.input_slots}
    for pt in sample_points(d, n=5):
        M, f = sf_mod.mass(pt.state, pt.p), sf_mod.bias(pt.state, pt.p)
        G = sf_mod.input_matrix(pt.state, pt.p)
        E = eom_mod.residual(pt.state, pt.a, pt.u, pt.p)
        for key, rows in _rows(d):
            assert_close(E[rows], reference(d.eom[key][1], pt.env), f"E[{key}]")
            eq = d.sf[key]
            assert_close(f[rows], reference(eq.f, pt.env), f"f[{key}]")
            for col, block in eq.M.items():
                assert_close(M[rows, acc[col]], reference(block, pt.env), f"M[{key}][{col}]")
            for col, block in eq.G.items():
                assert_close(G[rows, inp[col]], reference(block, pt.env), f"G[{key}][{col}]")
        kin = sf_mod.kinematics(pt.state, pt.p)
        for var, slot in zip(d.layout.config_variables, d.layout.configuration):
            assert_close(kin[slot.slice], reference(var.t_diff(), pt.env), f"d/dt {var}")
        assert lin_mod is not None
        assert lin_mod.M_lin(pt.state, pt.a, pt.u, pt.p).shape[0] == 2 * d.layout.n_accel


@pytest.mark.parametrize("build", s1_systems.SYSTEMS, ids=lambda f: f.__name__)
def test_casadi_target_equals_numpy(build):
    pytest.importorskip("casadi")
    from tests.codegen_systems import derived, modules, sample_points

    d = derived(build)
    np_mods, ca_mods = modules(build, "numpy"), modules(build, "casadi")
    for pt in sample_points(d, n=3):
        calls = {
            "mass": (pt.state, pt.p),
            "bias": (pt.state, pt.p),
            "input_matrix": (pt.state, pt.p),
            "kinematics": (pt.state, pt.p),
            "energy": (pt.state, pt.p),
            "residual": (pt.state, pt.a, pt.u, pt.p),
            "M_lin": (pt.state, pt.a, pt.u, pt.p),
            "A_lin": (pt.state, pt.a, pt.u, pt.p),
            "B_lin": (pt.state, pt.a, pt.u, pt.p),
        }
        for np_mod, ca_mod in zip(np_mods, ca_mods):
            for name, args in calls.items():
                if hasattr(np_mod, name):
                    ref = np.asarray(getattr(np_mod, name)(*args), float)
                    got = np.asarray(getattr(ca_mod, name)(*args), float).reshape(ref.shape)
                    scale = max(1.0, float(np.max(np.abs(ref)))) if ref.size else 1.0
                    assert np.max(np.abs(got - ref), initial=0.0) <= 1e-12 * scale, name
