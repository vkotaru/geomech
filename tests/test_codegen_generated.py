"""``examples/generated/`` (spec verification item 3) and the dependency rules.

- Regenerating the modules of the nine VBL systems reproduces the checked-in files byte
  for byte (CI also runs ``scripts/regenerate_generated.py`` and ``git diff --exit-code``).
- Each file imports and runs with ``geomech`` blocked: only ``numpy`` is needed.
- The files are ruff-clean (format and lint), when ruff is available.
- ``import geomech`` and ``import geomech.codegen`` import neither scipy nor casadi, and
  the CasADi target without casadi raises ``ImportError`` with the install hint.
"""

from __future__ import annotations

import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from tests.codegen_systems import derived, sample_points

ROOT = Path(__file__).resolve().parents[1]
GENERATED = ROOT / "examples" / "generated"


def _script():
    spec = importlib.util.spec_from_file_location(
        "regenerate_generated", ROOT / "scripts" / "regenerate_generated.py"
    )
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_regeneration_produces_no_diff():
    sources = _script().sources()
    assert len(sources) == 18
    on_disk = {p.name for p in GENERATED.glob("*.py")}
    assert on_disk == set(sources), f"stale or missing files: {on_disk ^ set(sources)}"
    drift = [name for name, src in sources.items() if (GENERATED / name).read_text() != src]
    assert not drift, (
        f"examples/generated is out of date for {drift}; run"
        " `python scripts/regenerate_generated.py`"
    )


_RUNNER = r"""
import importlib.abc, json, runpy, sys
import numpy as np

class _Block(importlib.abc.MetaPathFinder):
    def find_spec(self, name, path, target=None):
        if name == "geomech" or name.startswith("geomech."):
            raise ImportError("geomech is blocked in this check")

sys.meta_path.insert(0, _Block())
cases = json.loads(sys.argv[1])
out = {}
for path, (state, a, u, p) in cases.items():
    mod = runpy.run_path(path)
    assert not any(k.startswith("geomech") for k in sys.modules), "imported geomech"
    p = {k: np.array(v) if isinstance(v, list) else v for k, v in p.items()}
    names = [n for n in ("mass", "bias", "input_matrix", "energy") if n in mod]
    vals = {n: np.asarray(mod[n](state, p)).tolist() for n in names}
    for n in ("M_lin", "A_lin", "B_lin"):
        if n in mod:
            vals[n] = np.asarray(mod[n](state, a, u, p)).tolist()
    vals["kinematics"] = np.asarray(mod["kinematics"](state)).tolist()
    out[path] = vals
print(json.dumps(out))
"""


def _point(build):
    pt = sample_points(derived(build), n=1)[0]
    p = {k: (np.asarray(v).tolist() if np.ndim(v) else v) for k, v in pt.p.items()}
    return [pt.state.tolist(), pt.a.tolist(), pt.u.tolist(), p], pt


def test_generated_files_run_without_geomech():
    script = _script()
    from tests import vbl_systems
    from tests.codegen_systems import modules

    cases, expected = {}, {}
    for name in script.SYSTEMS:
        build = getattr(vbl_systems, name)
        args, pt = _point(build)
        sf_mod, _, lin_mod = modules(build)
        for fname, mod in ((f"{name}.py", sf_mod), (f"{name}_lin.py", lin_mod)):
            cases[str(GENERATED / fname)] = args
            expected[str(GENERATED / fname)] = (mod, pt)
    proc = subprocess.run(
        [sys.executable, "-I", "-c", _RUNNER, json.dumps(cases)],
        capture_output=True,
        text=True,
        cwd=str(GENERATED),
        check=False,
    )
    assert proc.returncode == 0, proc.stderr
    results = json.loads(proc.stdout)
    for path, vals in results.items():
        mod, pt = expected[path]
        for fn, val in vals.items():
            args = (pt.state,) if fn == "kinematics" else (pt.state, pt.p)
            if fn.endswith("_lin"):
                args = (pt.state, pt.a, pt.u, pt.p)
            np.testing.assert_allclose(val, getattr(mod, fn)(*args), rtol=1e-12, atol=1e-12)


def _ruff():
    exe = shutil.which("ruff") or str(Path(sys.executable).with_name("ruff"))
    return exe if Path(exe).exists() else None


@pytest.mark.parametrize("command", [["format", "--check"], ["check"]])
def test_generated_files_are_ruff_clean(command):
    ruff = _ruff()
    if ruff is None:
        pytest.skip("ruff is not installed (pip install -e '.[dev]')")
    proc = subprocess.run(
        [ruff, *command, str(GENERATED)], capture_output=True, text=True, cwd=ROOT, check=False
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr


def test_import_geomech_needs_no_optional_dependency():
    code = (
        "import sys, geomech, geomech.codegen;"
        " bad = [m for m in ('scipy', 'casadi') if m in sys.modules];"
        " assert not bad, bad"
    )
    proc = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, cwd=ROOT)
    assert proc.returncode == 0, proc.stderr


def test_missing_casadi_raises_with_install_hint(monkeypatch):
    from geomech.codegen import lambdify, load_module, module_source
    from tests.vbl_systems import point_mass_r3

    d = derived(point_mass_r3)
    monkeypatch.setitem(sys.modules, "casadi", None)
    with pytest.raises(ImportError, match=r"geomech\[casadi\]"):
        lambdify(d.system.L, d.layout, target="casadi")
    src = module_source(d.sf, d.layout, target="casadi")  # writing the source needs no casadi
    assert "import casadi as ca" in src
    with pytest.raises(ImportError, match=r"geomech\[casadi\]"):
        load_module(src)
