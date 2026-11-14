"""Adversarial review of ``geomech.codegen`` (``specs/reviews/codegen-review.md``).

Passing tests record what held under attack; ``xfail(strict=True)`` tests are confirmed
findings, named ``F<n>`` after the review. A strict xfail turns into a failure once the
finding is fixed, so the marker has to be removed together with the fix.

1. Translation on trees the 64 derived systems never produce: the typed random trees of
   ``specs/reviews/scripts/engine_fuzz.py`` (rows, outer products, nested ``Transpose``,
   ``Vee``, literals) through ``lambdify`` for both targets against ``numeric.evaluate``.
2. CSE and identifiers: names that are keywords, builtins, helpers, temporaries, LaTeX,
   unicode; no textual common subexpression left in the checked-in modules.
3. ``StateLayout``: several S² and SO(3) variables, flatten/unflatten, kinematics against
   the exact flow, the ``LAYOUT`` slices, scalar/vector coupling blocks.
4. ``energy``: five systems not in the suite against hand-derived closed forms, and
   ``dE/dt`` equal to the power of nonzero inputs.
5. Generated modules: deterministic across ``PYTHONHASHSEED``; ``--check`` and stale files.
6. Performance: CSE keeps shared trees linear; the cost of ``np.cross``.
"""

from __future__ import annotations

import ast
import collections
import importlib.util
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from geomech import (
    S2,
    SO3,
    Cross,
    Dot,
    Matrix,
    MVMul,
    Number,
    One,
    Scalar,
    SVMul,
    SystemVariables,
    Transpose,
    VAdd,
    Vector,
    VVMul,
    compute_eom,
    getScalars,
    linearize,
    to_standard_form,
)
from geomech.codegen import (
    StateLayout,
    collect_params,
    lambdify,
    load_module,
    module_source,
    to_source,
)
from geomech.core.base.types import ExprType
from geomech.numeric import Evaluator, evaluate, expm_so3, rand_rot, sample_state
from geomech.simulation import accelerations

ROOT = Path(__file__).resolve().parents[1]
GENERATED = ROOT / "examples" / "generated"
SCRIPT = ROOT / "scripts" / "regenerate_generated.py"
HALF = Scalar("0.5", value=0.5, attr=["Constant"])


def _has_casadi() -> bool:
    try:
        import casadi  # noqa: F401
    except ImportError:
        return False
    return True


TARGETS = [
    "numpy",
    pytest.param("casadi", marks=pytest.mark.skipif(not _has_casadi(), reason="casadi")),
]


def _flat(x):
    return np.atleast_1d(np.asarray(x, dtype=float)).reshape(-1)


def _close(got, ref, tol=1e-9):
    got, ref = _flat(got), _flat(ref)
    if got.shape != ref.shape:
        return False
    scale = max(1.0, float(np.max(np.abs(ref)))) if ref.size else 1.0
    return bool(np.all(np.abs(got - ref) <= tol * scale))


def _load_script(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------------------
# 1. Fuzzed translation
# ---------------------------------------------------------------------------


def _fuzz_setup(seed):
    ef = _load_script("engine_fuzz", ROOT / "specs" / "reviews" / "scripts" / "engine_fuzz.py")
    L = ef.Leaves()
    # no variation symbols in generated code: map them to ordinary vectors
    L.xi, L.xip, L.eta = L.a, L.b, L.c
    f_in, tau = Vector("f_{in}"), Scalar("tau")
    L.vectors = [L.a, L.b, L.c, L.e3, L.u, L.q, L.w, ef.ZeroVector, f_in]
    L.scalars = L.scalars + [
        Number(-2.5),
        Number(-3),
        Number(1e-20),
        Number(-1e-20),
        Scalar("0.5", value=0.5),
        tau,
    ]
    L.variation_leaves = set()

    class Gen(ef.Gen):
        def wrap_calc(self, e):  # product rules over whole subtrees instead of δ
            return e.t_diff()

    rng = np.random.default_rng(seed)
    gen = Gen(rng, L, calc=True, variation=False)
    variables = SystemVariables(scalars=[L.s1, L.s2], vectors=[L.a, L.b, L.q, L.p], matrices=[L.R])
    return rng, gen, variables, [tau, f_in]


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_fuzz_translation_matches_evaluate(seed):
    """1000 random trees per seed, every node class, both targets, against evaluate."""
    rng, gen, variables, inputs = _fuzz_setup(seed)
    targets = ["numpy"] + (["casadi"] if _has_casadi() else [])
    classes, checked, bad = collections.Counter(), 0, []
    for i in range(1000):
        e = gen.any(int(rng.integers(1, 6)))
        try:
            layout = StateLayout(variables, inputs, collect_params([e], variables, inputs))
            fns = {t: lambdify(e, layout, t) for t in targets}
        except ValueError as err:  # a Variation (S² δq motif) is rejected by design
            assert "Variation" in str(err), err
            continue
        stack = [e]
        while stack:
            n = stack.pop()
            classes[type(n).__name__] += 1
            stack.extend(getattr(n, "nodes", None) or ())
        ev = Evaluator(env=sample_state(variables, inputs, rng), rng=np.random.default_rng(i))
        ref = ev(e)
        args = (
            layout.state_from_env(ev.env),
            layout.inputs_from_env(ev.env),
            layout.params_from_env(ev.env),
            layout.accel_from_env(ev.env),
        )
        for t, fn in fns.items():
            if not _close(fn(*args), ref):
                bad.append((t, str(e)))
        checked += 1
    assert not bad, bad[:3]
    assert checked > 900
    for cls in ("VVMul", "Transpose", "Vee", "Hat", "MVMul", "MMMul", "Cross", "SMMul"):
        assert classes[cls] > 50, (cls, classes)


@pytest.mark.parametrize("target", TARGETS)
@pytest.mark.parametrize(
    "build",
    [
        lambda a, b, M, N, s: VVMul(Transpose(Transpose(a)), Transpose(b)),
        lambda a, b, M, N, s: MVMul(MVMul(Transpose(a), M), N),
        lambda a, b, M, N, s: Transpose(MVMul(MVMul(Transpose(a), M), N)),
        lambda a, b, M, N, s: VVMul(MVMul(Transpose(a), M), b),
        lambda a, b, M, N, s: VVMul(Transpose(MVMul(Transpose(a), M)), Transpose(b)),
        lambda a, b, M, N, s: VVMul(SVMul(Transpose(a), s), b),
        lambda a, b, M, N, s: VVMul(VAdd(Transpose(a), MVMul(Transpose(b), M)), a),
        lambda a, b, M, N, s: Transpose(Transpose(M * N)),
    ],
    ids=[
        "TT-outer",
        "row-chain",
        "row-chain-T",
        "row-M-col",
        "outer-of-row",
        "s-row",
        "row-sum",
        "TT-M",
    ],
)
def test_row_and_transpose_shapes(build, target):
    a, b, M, N, s = Vector("a"), Vector("b"), Matrix("M"), Matrix("N"), Scalar("s")
    variables = SystemVariables(scalars=[s], vectors=[a])
    e = build(a, b, M, N, s)
    layout = StateLayout(variables, [], collect_params([e], variables))
    ev = Evaluator(env=sample_state(variables, [], np.random.default_rng(0)))
    ref = ev(e)
    got = lambdify(e, layout, target)(
        layout.state_from_env(ev.env), [], layout.params_from_env(ev.env)
    )
    assert _close(got, ref)


@pytest.mark.parametrize("target", TARGETS)
@pytest.mark.parametrize("value", [1e-20, -1e-20, 1e20, -2.5, -0.0, 5e-324, 1e300, 2**70, -3])
def test_numeric_literals_are_exact(value, target):
    s, x = Scalar("s"), Vector("x")
    layout = StateLayout(SystemVariables(scalars=[s], vectors=[x]))
    state = np.arange(1.0, 9.0)
    for e in (Number(value) * s, SVMul(x, Number(value))):
        if value != 0:  # -0.0 carries the Zero flag and is written 0.0
            assert repr(value) in to_source(e, layout)
        ref = evaluate(e, {"s": 1.0, "x": state[1:4].reshape(3, 1)})
        assert _close(lambdify(e, layout, target)(state, [], {}), ref, tol=1e-15)


def test_infinite_literal_generates_valid_code():
    s = Scalar("s")
    layout = StateLayout(SystemVariables(scalars=[s]))
    fn = lambdify(Number(float("inf")) * s, layout)  # or a ValueError here, naming it
    assert fn([2.0, 0.0], [], {}) == float("inf")


# ---------------------------------------------------------------------------
# 2. Identifiers and CSE
# ---------------------------------------------------------------------------

AWKWARD_NAMES = [
    "sum", "max", "np", "ca", "_t0", "t0", "_t1", "_hat", "hat", "_vec", "_arr",
    "_check_params", "_PARAMS", "LAYOUT", "copy", "lambda", "None", "out", "out_", "a",
    "u", "p", "state", "x^\\ast", "x_ast", "e_3", "α", "β", "ω̇", 'q"uote', "back\\slash\\",
    "__init__", "_", "float", "range", "mass", "s_", "1x", "s_1x", "m", "m ",
]  # fmt: skip


@pytest.mark.parametrize("target", TARGETS)
def test_awkward_symbol_names(target):
    syms = [Scalar(n) for n in AWKWARD_NAMES]
    x, q, R = Vector("x"), S2("q"), SO3("R")
    variables = SystemVariables(vectors=[x, q], matrices=[R])
    e = sum((Number(i + 1) * s for i, s in enumerate(syms[1:])), Number(1) * syms[0])
    e = e + Dot(x, q) * syms[0] + Dot(R.get_tangent_vector(), x) * syms[5]
    layout = StateLayout(variables, [], collect_params([e], variables))
    names = [s.name for s in layout]
    assert all(n.isidentifier() and not n.startswith("_") for n in names)
    assert len(set(names)) == len(names)
    ev = Evaluator(env=sample_state(variables, [], np.random.default_rng(1)))
    ref = ev(e)
    got = lambdify(e, layout, target)(
        layout.state_from_env(ev.env), [], layout.params_from_env(ev.env)
    )
    assert _close(got, ref, tol=1e-12)


def test_checked_in_modules_have_no_textual_common_subexpression():
    for path in sorted(GENERATED.glob("*.py")):
        tree = ast.parse(path.read_text())
        for fn in (n for n in tree.body if isinstance(n, ast.FunctionDef)):
            if fn.name.startswith("_"):
                continue
            seen = collections.Counter()
            for node in ast.walk(fn):
                if isinstance(node, (ast.Call, ast.BinOp)):
                    text = ast.unparse(node)
                    if text.count("(") > 1 and not text.startswith(("_", "np.zeros", "np.eye")):
                        seen[text] += 1
            dups = [t for t, n in seen.items() if n > 1]
            assert not dups, f"{path.name}:{fn.name} repeats {dups[:3]}"


# ---------------------------------------------------------------------------
# 3. StateLayout
# ---------------------------------------------------------------------------


def _big_system():
    s1, s2 = Scalar("s_1"), Scalar("s_2")
    x = Vector("x")
    q1, q2 = S2("q_1"), S2("q_2")
    R1, R2 = SO3("R_1"), SO3("R_2")
    variables = SystemVariables(scalars=[s1, s2], vectors=[q1, x, q2], matrices=[R1, R2])
    vel = [s1.t_diff(), s2.t_diff(), x.t_diff()]
    vel += [v.get_tangent_vector() for v in (q1, q2, R1, R2)]
    kinetic = [
        v * v if isinstance(v, Scalar) or v.type == ExprType.SCALAR else Dot(v, v) for v in vel
    ]
    L = HALF * sum(kinetic[1:], kinetic[0])
    return variables, L


def test_layout_order_and_roundtrip_with_several_manifolds():
    variables, _ = _big_system()
    layout = StateLayout(variables)
    keys = [s.key for s in layout.configuration]
    assert keys == ["s_1", "s_2", "q_1", "x", "q_2", "R_1", "R_2"]
    assert (layout.n_config, layout.n_state, layout.n_accel) == (29, 46, 17)
    env = sample_state(variables, [], np.random.default_rng(0))
    back = layout.env_from_arrays(layout.state_from_env(env), layout.accel_from_env(env))
    assert set(back) == set(env)
    for k in env:
        np.testing.assert_array_equal(np.asarray(back[k]), np.asarray(env[k]))


@pytest.mark.parametrize("target", TARGETS)
def test_kinematics_is_the_exact_flow_row_major(target):
    variables, L = _big_system()
    eom = compute_eom(L, Number(0), variables)
    layout = StateLayout(variables)
    mod = load_module(module_source(to_standard_form(eom, variables, []), layout, target))
    state = layout.state_from_env(sample_state(variables, [], np.random.default_rng(2)))
    h, fd = 1e-6, np.zeros(layout.n_config)
    for c, v in zip(layout.configuration, layout.velocity):
        w = state[v.slice]
        if c.manifold == "SO3":
            R = state[c.slice].reshape(3, 3)
            fd[c.slice] = ((R @ expm_so3(h * w) - R @ expm_so3(-h * w)) / (2 * h)).reshape(9)
        elif c.manifold == "S2":
            q = state[c.slice]
            fd[c.slice] = (expm_so3(h * w) @ q - expm_so3(-h * w) @ q) / (2 * h)
        else:
            fd[c.slice] = w
    assert np.max(np.abs(_flat(mod.kinematics(state)) - fd)) < 1e-8
    slots = layout.configuration + layout.velocity
    for d, s in zip(mod.LAYOUT["configuration"] + mod.LAYOUT["velocity"], slots, strict=True):
        assert tuple(d["slice"]) == (s.start, s.stop) and d["symbol"] == s.key


@pytest.mark.parametrize("target", TARGETS)
def test_scalar_vector_coupling_blocks(target):
    """1×3 and 3×1 mass blocks between R, R³, S² and SO(3) (CasADi transposes 1×3)."""
    s, x, q, R = Scalar("s"), Vector("x"), S2("q"), SO3("R")
    m, c = getScalars("m c", attr=["Constant"])
    e = Vector("e", attr=["Constant"])
    J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
    w, W = q.get_tangent_vector(), R.get_tangent_vector()
    sd, xd = s.t_diff(), x.t_diff()
    L = (
        HALF * m * sd * sd
        + HALF * Dot(xd, xd)
        + c * sd * (Dot(e, xd) + Dot(e, w) + Dot(e, W))
        + HALF * Dot(W, J * W)
        + HALF * Dot(w, w)
        + c * Dot(xd, R * W)
    )
    fs, F = Scalar("f_s"), Vector("F")
    variables = SystemVariables(scalars=[s], vectors=[x, q], matrices=[R])
    eom = compute_eom(L, s.delta() * fs + Dot(x.delta(), F), variables)
    sf = to_standard_form(eom, variables, [fs, F])
    exprs = [L] + [eq for _, eq in eom.values()]
    layout = StateLayout(variables, [fs, F], collect_params(exprs, variables, [fs, F]))
    mod = load_module(module_source(sf, layout, target))
    rng = np.random.default_rng(0)
    ev = Evaluator(env=sample_state(variables, [fs, F], rng), rng=rng)
    for sym in layout.params:
        ev(sym)
    state, p = layout.state_from_env(ev.env), layout.params_from_env(ev.env)
    Mgen = np.asarray(mod.mass(state, p), dtype=float)
    rows = [
        str(s.delta()),
        str(x.delta()),
        str(q.get_variation_vector()),
        str(R.get_variation_vector()),
    ]
    spans = [(0, 1), (1, 4), (4, 7), (7, 10)]
    for (r0, r1), key in zip(spans, rows):
        for acc in layout.acceleration:
            blk = sf[key].M.get(acc.key)
            got = Mgen[r0:r1, acc.start : acc.stop]
            ref = 0 if blk is None else np.asarray(ev(blk), float).reshape(got.shape)
            assert np.max(np.abs(got - ref)) < 1e-12
    assert np.max(np.abs(Mgen[0, 1:4])) > 0 and np.max(np.abs(Mgen[1:4, 0])) > 0


def test_omitted_acceleration_raises():
    s, f = Scalar("s"), Scalar("f")
    layout = StateLayout(SystemVariables(scalars=[s]), [f])
    fn = lambdify(s.t_diff().t_diff() + f, layout)
    with pytest.raises((TypeError, ValueError)):
        fn(np.array([1.0, 2.0]), [0.5], {})  # `a` defaults to None


def test_scalar_parameter_as_one_element_sequence():
    s, m = Scalar("s"), Scalar("m", attr=["Constant"])
    layout = StateLayout(SystemVariables(scalars=[s]), [], [m])
    fn = lambdify(m * Number(3) + s, layout)
    assert fn([1.0, 2.0], [], {"m": 1.5}) == 5.5
    # vectors accept any 3-entry sequence (`_vec`); a scalar given as [1.5] must either
    # be accepted as 1.5 or rejected, never become [1.5, 1.5, 1.5] + 1 = [2.5, 2.5, 2.5]
    try:
        got = fn([1.0, 2.0], [], {"m": [1.5]})
    except (TypeError, ValueError):
        return
    assert np.ndim(got) == 0 and float(got) == 5.5


# ---------------------------------------------------------------------------
# 4. Energy: hand-derived closed forms and the power balance
# ---------------------------------------------------------------------------


def _nonsymmetric_body():
    J = Matrix("J", attr=["Constant"])  # not symmetric: E = ½ Ωᵀ J Ω still
    R, M = SO3("R"), Vector("M")
    Om, eta = R.get_tangent_vector(), R.get_variation_vector()

    def energy(st, p):
        return 0.5 * st[9:12] @ p["J"] @ st[9:12]

    def power(st, u):
        return u[0:3] @ st[9:12]

    return Dot(Om, J * Om) * HALF, Dot(eta, M), SystemVariables(matrices=[R]), [M], energy, power


def _gyroscopic_particle():
    m, k, c = getScalars("m k c", attr=["Constant"])
    b, x, F = Vector("b", attr=["Constant"]), Vector("x"), Vector("F")
    v = x.t_diff()
    L = HALF * m * Dot(v, v) + c * Dot(v, Cross(b, x)) - HALF * k * Dot(x, x)

    def energy(st, p):  # the gyroscopic term is linear in ẋ and drops out
        return 0.5 * p["m"] * st[3:6] @ st[3:6] + 0.5 * p["k"] * st[0:3] @ st[0:3]

    def power(st, u):
        return u @ st[3:6]

    return L, Dot(x.delta(), F), SystemVariables(vectors=[x]), [F], energy, power


def _body_with_pendulum():
    """Rigid body R with a point mass on a rod q hung from the body point ρ (S² × SO(3))."""
    J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
    m, g, l = getScalars("m g l", attr=["Constant"])
    rho, e3 = Vector("\\rho", attr=["Constant"]), Vector("e_3", attr=["Constant"])
    R, q, M, f = SO3("R"), S2("q"), Vector("M"), Vector("f")
    Om, eta = R.get_tangent_vector(), R.get_variation_vector()
    xp = R * rho + q * l
    vp = xp.t_diff()
    L = HALF * Dot(Om, J * Om) + HALF * m * Dot(vp, vp) - m * g * Dot(xp, e3)

    def energy(st, p):
        qn, Rn, wn, Wn = st[0:3], st[3:12].reshape(3, 3), st[12:15], st[15:18]
        v = Rn @ np.cross(Wn, p["\\rho"]) + p["l"] * np.cross(wn, qn)
        h = (Rn @ p["\\rho"] + p["l"] * qn) @ p["e_3"]
        return 0.5 * Wn @ p["J"] @ Wn + 0.5 * p["m"] * v @ v + p["m"] * p["g"] * h

    def power(st, u):
        return u[0:3] @ st[15:18] + u[3:6] @ st[12:15]

    dW = Dot(eta, M) + Dot(q.get_variation_vector(), f)
    return L, dW, SystemVariables(vectors=[q], matrices=[R]), [M, f], energy, power


def _variable_inertia():
    m, k = getScalars("m k", attr=["Constant"])
    J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
    s, R, fs, M = Scalar("s"), SO3("R"), Scalar("f_s"), Vector("M")
    Om, eta = R.get_tangent_vector(), R.get_variation_vector()
    sd = s.t_diff()
    L = HALF * m * sd * sd + HALF * (One + s * s) * Dot(Om, J * Om) - HALF * k * s * s

    def energy(st, p):
        S, Sd, W = st[0], st[10], st[11:14]
        return 0.5 * p["m"] * Sd**2 + 0.5 * (1 + S * S) * W @ p["J"] @ W + 0.5 * p["k"] * S * S

    def power(st, u):
        return u[0] * st[10] + u[1:4] @ st[11:14]

    variables = SystemVariables(scalars=[s], matrices=[R])
    return L, s.delta() * fs + Dot(eta, M), variables, [fs, M], energy, power


def _two_spheres():
    m1, m2, k = getScalars("m_1 m_2 k", attr=["Constant"])
    b = Vector("b", attr=["Constant"])
    q1, q2, f1, f2 = S2("q_1"), S2("q_2"), Vector("f_1"), Vector("f_2")
    w1, w2 = q1.get_tangent_vector(), q2.get_tangent_vector()
    L = HALF * m1 * Dot(w1, w1) + HALF * m2 * Dot(w2, w2) + Dot(w1, b) - k * Dot(q1, q2)

    def energy(st, p):
        return (
            0.5 * p["m_1"] * st[6:9] @ st[6:9]
            + 0.5 * p["m_2"] * st[9:12] @ st[9:12]
            + p["k"] * st[0:3] @ st[3:6]
        )

    def power(st, u):
        return u[0:3] @ st[6:9] + u[3:6] @ st[9:12]

    dW = Dot(q1.get_variation_vector(), f1) + Dot(q2.get_variation_vector(), f2)
    return L, dW, SystemVariables(vectors=[q1, q2]), [f1, f2], energy, power


ENERGY_SYSTEMS = [
    _nonsymmetric_body,
    _gyroscopic_particle,
    _body_with_pendulum,
    _variable_inertia,
    _two_spheres,
]


def _state(layout, rng):
    st = np.zeros(layout.n_state)
    for c, v in zip(layout.configuration, layout.velocity):
        if c.manifold == "S2":
            q = rng.normal(size=3)
            q /= np.linalg.norm(q)
            w = rng.normal(size=3)
            st[c.slice], st[v.slice] = q, w - q * (q @ w)
        elif c.manifold == "SO3":
            st[c.slice], st[v.slice] = rand_rot(rng).reshape(9), rng.normal(size=3)
        else:
            st[c.slice], st[v.slice] = rng.normal(size=c.size), rng.normal(size=v.size)
    return st


def _params(layout, rng):
    p = {}
    for s in layout.param_slots:
        if s.shape == (3, 3):
            A = 0.3 * rng.normal(size=(3, 3))
            p[s.key] = np.diag([1.0, 2.0, 3.0]) + (A + A.T if s.symbol.is_symmetric else A)
        elif s.shape == (3,):
            p[s.key] = rng.normal(size=3)
        else:
            p[s.key] = 1.0 + abs(rng.normal())
    return p


@pytest.mark.parametrize("build", ENERGY_SYSTEMS, ids=lambda f: f.__name__.strip("_"))
def test_energy_closed_form_and_power_balance(build):
    L, dW, variables, inputs, closed, power = build()
    eom = compute_eom(L, dW, variables)
    exprs = [L] + [eq for _, eq in eom.values()]
    layout = StateLayout(variables, inputs, collect_params(exprs, variables, inputs))
    mod = load_module(
        module_source(to_standard_form(eom, variables, inputs), layout, lagrangian=L)
    )
    rng, h = np.random.default_rng(3), 1e-6
    for _ in range(5):
        st, p, u = _state(layout, rng), _params(layout, rng), rng.normal(size=layout.n_inputs)
        E = mod.energy(st, p)
        assert abs(E - closed(st, p)) < 1e-12 * max(1.0, abs(E))
        rate = np.concatenate([mod.kinematics(st), accelerations(mod, st, u, p)])
        dE = (mod.energy(st + h * rate, p) - mod.energy(st - h * rate, p)) / (2 * h)
        P = power(st, u)
        assert abs(dE - P) < 1e-7 * max(1.0, abs(P), float(np.max(np.abs(rate)))), (dE, P)


# ---------------------------------------------------------------------------
# 5. Generated modules
# ---------------------------------------------------------------------------


def test_generated_source_is_independent_of_hash_seed():
    code = (
        "import hashlib, importlib.util, sys;"
        f" spec = importlib.util.spec_from_file_location('r', {str(SCRIPT)!r});"
        " m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m);"
        " s = m.sources();"
        " print(hashlib.sha256(''.join(k + v for k, v in sorted(s.items())).encode()).hexdigest())"
    )
    digests = set()
    for seed in ("0", "1", "4242"):
        env = {**os.environ, "PYTHONHASHSEED": seed}
        proc = subprocess.run(
            [sys.executable, "-c", code], capture_output=True, text=True, cwd=ROOT, env=env
        )
        assert proc.returncode == 0, proc.stderr
        digests.add(proc.stdout.strip())
    assert len(digests) == 1


def test_check_mode_reports_stale_files(tmp_path):
    script = _load_script("regenerate_generated", ROOT / "scripts" / "regenerate_generated.py")
    for path in GENERATED.glob("*.py"):
        shutil.copy(path, tmp_path / path.name)
    assert script.main.__globals__["OUT"] == GENERATED
    script.OUT = tmp_path
    assert script.main(["--check"]) == 0  # the copy is up to date
    (tmp_path / "removed_system.py").write_text("# generated for a system no longer listed\n")
    assert script.main(["--check"]) == 1


# ---------------------------------------------------------------------------
# 6. Performance
# ---------------------------------------------------------------------------


def test_cse_keeps_shared_trees_linear():
    """A tree whose every level uses the previous one twice: 3·2^d nodes, d+O(1) lines."""
    x, s, a = Vector("x"), Scalar("s"), Vector("a")
    variables = SystemVariables(scalars=[s], vectors=[x])
    layout = StateLayout(variables, [], [a])
    e = x
    for _ in range(12):
        e = VAdd(Cross(e, a), SVMul(e, s))
    src = to_source(e, layout)
    assert src.count("\n") < 40
    env = sample_state(variables, [], np.random.default_rng(0))
    env["a"] = np.array([[0.1], [0.2], [0.3]])
    ref = evaluate(e, env)
    got = lambdify(e, layout)(layout.state_from_env(env), [], {"a": env["a"]})
    assert _close(got, ref, tol=1e-12)


def test_generated_numpy_code_avoids_np_cross():
    """``double_rigid_pendulum.bias`` makes 28 ``np.cross`` calls (283 µs); a component-wise
    ``_cross`` helper gives the same values in 47 µs."""
    src = (GENERATED / "double_rigid_pendulum.py").read_text()
    assert "np.cross(" not in src


def test_casadi_extra_is_not_an_exact_pin():
    text = (ROOT / "pyproject.toml").read_text()
    assert '"casadi==' not in text


def test_linearized_and_eom_modules_for_a_decoupled_system():
    """A variable absent from L: zero equations, an all-zero mass row, no inputs."""
    s, x, y = Scalar("s"), Vector("x"), Vector("y")
    m = Scalar("m", attr=["Constant"])
    variables = SystemVariables(scalars=[s], vectors=[x, y])
    L = HALF * m * Dot(x.t_diff(), x.t_diff())
    eom = compute_eom(L, Number(0), variables)
    layout = StateLayout(variables, [], [m])
    st = np.arange(1.0, layout.n_state + 1)
    sf = load_module(module_source(to_standard_form(eom, variables, []), layout))
    assert sf.input_matrix(st, {"m": 2.0}).shape == (7, 0)
    M = sf.mass(st, {"m": 2.0})
    np.testing.assert_array_equal(M, np.diag([0, -2, -2, -2, 0, 0, 0]))
    lin = load_module(module_source(linearize(eom, variables, []), layout))
    assert lin.B_lin(st, np.zeros(7), [], {"m": 2.0}).shape == (14, 0)
