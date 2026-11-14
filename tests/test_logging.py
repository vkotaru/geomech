"""Phase 4.3: the ``geomech`` logger and the ``GEOMECH_DEBUG`` switch."""

from __future__ import annotations

import ast
import logging
import os
import subprocess
import sys
from pathlib import Path

import pytest

from geomech import compute_eom, linearize, to_standard_form
from tests.vbl_systems import rigid_pendulum_pivot, spherical_pendulum_s2

ROOT = Path(__file__).resolve().parent.parent


def _debug_messages(caplog, prefix):
    return [
        r.getMessage()
        for r in caplog.records
        if r.name == "geomech" and r.levelno == logging.DEBUG and r.getMessage().startswith(prefix)
    ]


def test_one_debug_line_per_compute_eom_stage(caplog):
    s = rigid_pendulum_pivot()
    with caplog.at_level(logging.DEBUG, logger="geomech"):
        compute_eom(s.L, s.dW, s.variables)
    lines = _debug_messages(caplog, "compute_eom: ")
    stages = [line.split(" ")[1] for line in lines]
    d_eta = "\\frac{d}{dt}(\\eta_{R})"
    assert stages == [
        "dL",
        "dW",
        "dS",
        "manifold_rules",
        "simplify",
        "ibp.prepare",
        f"ibp.collect[{d_eta}]",
        f"ibp.integrate[{d_eta}]",
        f"ibp.simplify[{d_eta}]",
        "ibp",
        "expand",
        "raw[\\eta_{R}]",
        "equation[\\eta_{R}]",
    ]
    for line in lines:
        head, _, text = line.partition(" nodes] ")
        assert head.split("[")[-1].isdigit(), line  # "[<count>"
        assert len(text) <= 200


def test_linearize_and_standard_form_log_per_row(caplog):
    s = spherical_pendulum_s2()
    eom = compute_eom(s.L, s.dW, s.variables)
    with caplog.at_level(logging.DEBUG, logger="geomech"):
        linearize(eom, s.variables, s.inputs)
        to_standard_form(eom, s.variables, s.inputs)
    lin = _debug_messages(caplog, "linearize: ")
    assert lin[0].startswith("linearize: dE[\\xi_{q}] [")
    assert any(line.startswith("linearize: raw[\\xi_{q}][") for line in lin)
    assert any(line.startswith("linearize: block[\\xi_{q}][") for line in lin)
    sf = _debug_messages(caplog, "to_standard_form: ")
    assert any(line.startswith("to_standard_form: f[\\xi_{q}] [") for line in sf)
    assert any(line.startswith("to_standard_form: M[\\xi_{q}][") for line in sf)


def test_silent_by_default(caplog):
    s = rigid_pendulum_pivot()
    with caplog.at_level(logging.INFO, logger="geomech"):
        compute_eom(s.L, s.dW, s.variables)
    assert [r for r in caplog.records if r.name.startswith("geomech")] == []


def test_geomech_debug_env_var_in_a_subprocess():
    code = (
        "from tests.vbl_systems import point_mass_r3\n"
        "from geomech import compute_eom\n"
        "s = point_mass_r3(); compute_eom(s.L, s.dW, s.variables)\n"
    )
    env = dict(os.environ, GEOMECH_DEBUG="1")
    env.pop("GEOMECH_TRACE", None)
    out = subprocess.run(
        [sys.executable, "-c", code], cwd=ROOT, env=env, capture_output=True, text=True
    )
    assert out.returncode == 0, out.stderr
    assert out.stdout == ""
    assert "geomech: compute_eom: dS [" in out.stderr
    assert "geomech: compute_eom: equation[\\delta{x}] [" in out.stderr
    # no rule trace unless GEOMECH_TRACE is set
    assert "geomech.rules" not in out.stderr


@pytest.mark.parametrize("value", ["", "0"])
def test_geomech_debug_off_values(value):
    code = (
        "import logging, geomech\n"
        "print(logging.getLogger('geomech').getEffectiveLevel() > logging.DEBUG)\n"
    )
    env = dict(os.environ, GEOMECH_DEBUG=value)
    out = subprocess.run(
        [sys.executable, "-c", code], cwd=ROOT, env=env, capture_output=True, text=True
    )
    assert out.stdout.strip() == "True", out.stderr


# Functions whose job is to print (display utilities); everything else logs.
_PRINTING_FUNCTIONS = {"print_eom", "print_tree", "summary", "show"}


def test_no_print_in_library_code():
    offenders = []
    for path in sorted((ROOT / "geomech").rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for func in ast.walk(tree):
            if not isinstance(func, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            if func.name in _PRINTING_FUNCTIONS:
                continue
            for node in ast.walk(func):
                if isinstance(node, ast.Call) and getattr(node.func, "id", None) == "print":
                    offenders.append(f"{path.relative_to(ROOT)}:{node.lineno} in {func.name}")
    assert offenders == []
