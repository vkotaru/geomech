"""Arguments and helpers of the generated code (fixes F1–F5 of
``specs/reviews/codegen-review.md``).

- F1: ``state``, ``a`` and ``u`` are required exactly when a function's expressions use
  them; ``None`` then raises ``ValueError`` naming the argument. A function that does not
  use ``u`` (or ``a``) accepts ``None``. Generated functions never read ``None`` as zero;
  ``geomech.simulation.accelerations`` does (documented there).
- F2: a scalar parameter is converted by ``_scalar``: Python/NumPy scalars, 0-d and
  single-element arrays give a float; more entries raise ``ValueError`` naming it.
- F3: ``±inf`` and ``nan`` literals are written ``np.inf``/``np.nan`` (``ca.inf``), and a
  ``nan`` literal is inlined, not a parameter.
- F4: ``regenerate_generated.py --check`` reports stale and untracked files.
- F5: ``_cross`` equals ``np.cross`` bit for bit.
- Advisory: generation is linear in the number of node objects on shared trees.
"""

from __future__ import annotations

import ast
import importlib.util
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np
import pytest

from geomech import Cross, Number, Scalar, SVMul, SystemVariables, VAdd, Vector
from geomech.codegen import StateLayout, collect_params, lambdify, load_module, to_source
from geomech.codegen.source import namespace
from tests.codegen_systems import VBL_SYSTEMS, derived, modules, sample_points

ROOT = Path(__file__).resolve().parents[1]
GENERATED = ROOT / "examples" / "generated"


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


def _converted_arrays(module_text: str) -> dict[str, set[str]]:
    """``{function: arrays it converts with _arr}`` of a generated module."""
    out = {}
    for fn in ast.parse(module_text).body:
        if isinstance(fn, ast.FunctionDef) and not fn.name.startswith("_"):
            arrays = set()
            for node in ast.walk(fn):
                if (
                    isinstance(node, ast.Call)
                    and isinstance(node.func, ast.Name)
                    and node.func.id == "_arr"
                    and isinstance(node.args[0], ast.Name)
                ):
                    arrays.add(node.args[0].id)
            out[fn.name] = arrays
    return out


# ---------------------------------------------------------------------------
# F1: required arguments
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("target", TARGETS)
@pytest.mark.parametrize("build", VBL_SYSTEMS, ids=lambda b: b.__name__)
def test_none_raises_exactly_when_the_argument_is_used(build, target):
    """Every function of the eom and linearized modules: ``a=None``/``u=None`` raise
    ``ValueError`` naming the argument when the function uses it, and give the value of
    a zero array otherwise (the argument is then not read at all)."""
    from geomech.codegen import module_source

    d = derived(build)
    pt = sample_points(d, n=1, seed=4)[0]
    sources = [module_source(d.eom, d.layout, target)]
    if d.lin is not None:
        sources.append(module_source(d.lin, d.layout, target))
    checked = 0
    for text in sources:
        mod = load_module(text)
        for name, arrays in _converted_arrays(text).items():
            if name not in ("residual", "M_lin", "A_lin", "B_lin"):
                continue
            fn = getattr(mod, name)
            full = fn(pt.state, pt.a, pt.u, pt.p)
            for arg in ("a", "u"):
                args = {"state": pt.state, "a": pt.a, "u": pt.u, "p": pt.p, arg: None}
                if arg in arrays:
                    n = {"a": d.layout.n_accel, "u": d.layout.n_inputs}[arg]
                    with pytest.raises(
                        ValueError, match=rf"^{arg} must have {n} entries, got None"
                    ):
                        fn(**args)
                else:
                    np.testing.assert_array_equal(_flat(fn(**args)), _flat(full))
                checked += 1
    assert checked >= 2


def test_inputs_and_accelerations_are_used_where_expected():
    """The checks above are not vacuous: the residual of a forced system reads ``u`` and
    ``a``, its ``B_lin`` reads neither (the input enters linearly)."""
    d = derived(VBL_SYSTEMS[1])  # spring_mass_r: one scalar input
    from geomech.codegen import module_source

    eom = _converted_arrays(module_source(d.eom, d.layout))
    lin = _converted_arrays(module_source(d.lin, d.layout))
    assert {"state", "a", "u"} <= eom["residual"]
    assert "u" not in lin["B_lin"] and "a" not in lin["B_lin"]
    assert d.layout.n_inputs == 1  # the size-1 case that used to give NaN


@pytest.mark.parametrize("target", TARGETS)
def test_size_one_arguments_omitted(target):
    s, f = Scalar("s"), Scalar("f")
    layout = StateLayout(SystemVariables(scalars=[s]), [f])
    fn = lambdify(s.t_diff().t_diff() + f, layout, target)
    with pytest.raises(ValueError, match=r"^a must have 1 entries, got None \(a is required"):
        fn([1.0, 2.0], [0.5], {})
    with pytest.raises(ValueError, match=r"^u must have 1 entries, got None \(u is required"):
        fn([1.0, 2.0], None, {}, [3.0])
    assert float(_flat(fn([1.0, 2.0], [0.5], {}, [3.0]))[0]) == 3.5


@pytest.mark.parametrize("target", TARGETS)
def test_no_inputs_accepts_u_none(target):
    """A system without inputs, and an expression without inputs: ``u=None`` is fine."""
    s, f = Scalar("s"), Scalar("f")
    layout = StateLayout(SystemVariables(scalars=[s]))
    fn = lambdify(s * Number(2) + s.t_diff(), layout, target)
    assert float(_flat(fn([1.0, 3.0], None, None))[0]) == 5.0
    layout = StateLayout(SystemVariables(scalars=[s]), [f])
    fn = lambdify(s * Number(2), layout, target)
    assert float(_flat(fn([1.0, 3.0], None, {}))[0]) == 2.0


def test_standard_form_functions_take_no_u():
    """``mass``, ``bias``, ``input_matrix``, ``energy`` and ``kinematics`` have no ``u`` or
    ``a`` argument; ``geomech.simulation.accelerations`` is where ``u=None`` means zero."""
    from geomech.simulation import accelerations

    build = VBL_SYSTEMS[5]  # quadrotor: four inputs
    sf, _, _ = modules(build)
    pt = sample_points(derived(build), n=1, seed=1)[0]
    for name in ("mass", "bias", "input_matrix", "energy"):
        assert "u" not in getattr(sf, name).__code__.co_varnames[:2]
    np.testing.assert_array_equal(
        accelerations(sf, pt.state, None, pt.p),
        accelerations(sf, pt.state, np.zeros(derived(build).layout.n_inputs), pt.p),
    )


# ---------------------------------------------------------------------------
# F2: scalar parameters
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("target", TARGETS)
@pytest.mark.parametrize(
    "value",
    [1.5, np.float64(1.5), np.float32(1.5), np.array(1.5), [1.5], (1.5,), np.array([[1.5]])],
    ids=["float", "float64", "float32", "0-d", "list", "tuple", "1x1"],
)
def test_scalar_parameter_forms(value, target):
    s, m = Scalar("s"), Scalar("m", attr=["Constant"])
    layout = StateLayout(SystemVariables(scalars=[s]), [], [m])
    got = lambdify(m * Number(3) + s, layout, target)([1.0, 2.0], [], {"m": value})
    assert _flat(got).shape == (1,) and float(_flat(got)[0]) == 5.5
    if target == "numpy":
        assert isinstance(got, float)


@pytest.mark.parametrize("target", TARGETS)
@pytest.mark.parametrize("value", [[1.0, 2.0], np.ones(3), np.eye(3), []])
def test_scalar_parameter_with_several_entries_is_rejected(value, target):
    s, m = Scalar("s"), Scalar("m_1", attr=["Constant"])
    layout = StateLayout(SystemVariables(scalars=[s]), [], [m])
    fn = lambdify(m * Number(3) + s, layout, target)
    with pytest.raises(ValueError, match=r"parameter 'm_1' is a scalar, got \d+ entries"):
        fn([1.0, 2.0], [], {"m_1": value})


def test_scalar_parameter_in_a_module():
    build = VBL_SYSTEMS[6]  # double rigid pendulum: m1, m2, g scalars
    sf, _, _ = modules(build)
    pt = sample_points(derived(build), n=1, seed=2)[0]
    wrapped = {k: (np.array([v]) if np.ndim(v) == 0 else v) for k, v in pt.p.items()}
    for name in ("mass", "bias", "input_matrix", "energy"):
        got, ref = getattr(sf, name)(pt.state, wrapped), getattr(sf, name)(pt.state, pt.p)
        np.testing.assert_array_equal(got, ref)


@pytest.mark.skipif(not _has_casadi(), reason="casadi")
def test_casadi_scalar_parameter_may_be_symbolic():
    import casadi as ca

    s, m = Scalar("s"), Scalar("m", attr=["Constant"])
    layout = StateLayout(SystemVariables(scalars=[s]), [], [m])
    fn = lambdify(m * Number(3) + s, layout, "casadi")
    sym = ca.SX.sym("m")
    out = fn([1.0, 2.0], [], {"m": sym})
    assert float(ca.Function("f", [sym], [out])(1.5)) == 5.5
    assert float(fn([1.0, 2.0], [], {"m": ca.DM(1.5)})) == 5.5
    with pytest.raises(ValueError, match="parameter 'm' is a scalar, got 2 entries"):
        fn([1.0, 2.0], [], {"m": ca.SX.sym("m", 2)})


# ---------------------------------------------------------------------------
# F3: non-finite literals
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("target", TARGETS)
@pytest.mark.parametrize("value", [float("inf"), float("-inf"), float("nan")])
def test_non_finite_literals(value, target):
    s = Scalar("s")
    variables = SystemVariables(scalars=[s])
    e = Number(value) * s
    assert collect_params([e], variables) == []
    layout = StateLayout(variables)
    src = to_source(e, layout, target)
    if np.isnan(value):
        assert "np.nan" in src
    else:
        mod = "ca" if target == "casadi" else "np"
        assert ("-" if value < 0 else "") + f"{mod}.inf" in src
    got = float(_flat(lambdify(e, layout, target)([2.0, 0.0], [], {}))[0])
    np.testing.assert_equal(got, value * 2.0)


def test_nan_literal_is_a_number_literal():
    from geomech.core.base.expressions import is_number_literal

    assert is_number_literal(Number(float("nan")))
    assert is_number_literal(Number(float("inf")))
    assert not is_number_literal(Scalar("nan"))


# ---------------------------------------------------------------------------
# F4: --check reports stale and untracked files
# ---------------------------------------------------------------------------


def _script():
    spec = importlib.util.spec_from_file_location(
        "regenerate_generated_f4", ROOT / "scripts" / "regenerate_generated.py"
    )
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_check_reports_untracked_and_write_deletes_stale(tmp_path, capsys):
    if shutil.which("git") is None:
        pytest.skip("git")
    script = _script()
    out = tmp_path / "generated"
    out.mkdir()
    for path in GENERATED.glob("*.py"):
        shutil.copy(path, out / path.name)
    script.OUT = out
    git = ["git", "-C", str(out)]
    subprocess.run([*git, "init", "-q"], check=True)
    names = sorted(p.name for p in out.glob("*.py"))
    subprocess.run([*git, "add", *names[1:]], check=True)  # names[0] is not tracked
    assert script.main(["--check"]) == 1
    assert f"not tracked by git: {names[0]}" in capsys.readouterr().out
    subprocess.run([*git, "add", names[0]], check=True)
    assert script.main(["--check"]) == 0
    (out / "removed_system.py").write_text("# stale\n")
    assert script.main(["--check"]) == 1
    assert "stale (no longer generated): removed_system.py" in capsys.readouterr().out
    assert script.main([]) == 0
    assert not (out / "removed_system.py").exists()
    (out / names[0]).unlink()
    assert script.main(["--check"]) == 1
    assert f"missing: {names[0]}" in capsys.readouterr().out


# ---------------------------------------------------------------------------
# F5: the component-wise cross product
# ---------------------------------------------------------------------------


def test_cross_helper_equals_np_cross():
    cross = namespace("numpy")["_cross"]
    rng = np.random.default_rng(0)
    for scale in (1.0, 1e-200, 1e200):
        for _ in range(500):
            a, b = rng.normal(size=3) * scale, rng.normal(size=3)
            got = cross(a, b)
            assert got.shape == (3,)
            np.testing.assert_array_equal(got, np.cross(a, b))


def test_no_np_cross_in_generated_numpy_code():
    for path in sorted(GENERATED.glob("*.py")):
        assert "np.cross(" not in path.read_text(), path.name


# ---------------------------------------------------------------------------
# Advisory: generation cost on shared trees
# ---------------------------------------------------------------------------


def test_generation_is_linear_on_shared_trees():
    """Depth 20, 3·2²⁰ tree nodes in 61 objects: 8 s when every node was keyed by its
    ``str``, a few ms now."""
    x, s, a = Vector("x"), Scalar("s"), Vector("a")
    layout = StateLayout(SystemVariables(scalars=[s], vectors=[x]), [], [a])
    e = x
    for _ in range(20):
        e = VAdd(Cross(e, a), SVMul(e, s))
    t0 = time.perf_counter()
    src = to_source(e, layout)
    assert time.perf_counter() - t0 < 1.0
    assert src.count("_cross(") == 20
