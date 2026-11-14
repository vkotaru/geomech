"""Phase 4.5 and 4.6: named rewrite rules, the rule trace, and un-shadowed submodules.

The rule names listed in each module docstring and the names passed to ``rule(...)`` in
the code must be the same set (both are parsed here), so the docstring stays the index
of what the code can do.
"""

from __future__ import annotations

import ast
import importlib
import logging
import os
import re
import subprocess
import sys
import types
from pathlib import Path

import pytest

import geomech
from geomech import S2, Cross, Dot, Scalar, Vector, full_simplify, set_trace
from geomech.utils.trace import rule, trace_enabled

ROOT = Path(__file__).resolve().parent.parent

RULE_MODULES = {
    "geomech.core.transformations.vector_identities": "vector",
    "geomech.core.transformations.simplification": "simplify",
    "geomech.core.transformations.pulling": "pull",
    "geomech.core.transformations.expansion": "expand",
    "geomech.core.transformations.manifold_rules": "manifold",
    "geomech.core.math.collecting": "collect",
    "geomech.core.math.ibp": "ibp",
}


def _source(module_name):
    return (ROOT / (module_name.replace(".", "/") + ".py")).read_text(encoding="utf-8")


def _docstring_names(module_name, prefix):
    doc = ast.get_docstring(ast.parse(_source(module_name))) or ""
    return set(re.findall(rf"^\s+({prefix}\.[A-Za-z0-9_]+)\b", doc, flags=re.MULTILINE))


def _code_names(module_name, prefix):
    """The literal names passed to ``rule(name, ...)`` (and forwarded through
    ``expand._distribute(..., name, expr)``), outside docstrings."""
    tree = ast.parse(_source(module_name))
    names = set()
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func = getattr(node.func, "id", None)
        args = node.args
        if func == "rule" and args and isinstance(args[0], ast.Constant):
            names.add(args[0].value)
        elif func == "_distribute" and len(args) > 6 and isinstance(args[6], ast.Constant):
            names.add(args[6].value)
    pattern = re.compile(rf"^{prefix}\.[A-Za-z0-9_]+$")
    bad = sorted(n for n in names if not pattern.match(str(n)))
    assert not bad, f"{module_name}: rule names without the '{prefix}.' prefix: {bad}"
    return names


@pytest.mark.parametrize("module_name, prefix", sorted(RULE_MODULES.items()))
def test_docstring_and_code_list_the_same_rules(module_name, prefix):
    documented = _docstring_names(module_name, prefix)
    used = _code_names(module_name, prefix)
    assert used, f"{module_name}: no rule(...) names found in the code"
    assert documented == used, (
        f"{module_name}: documented but not in code {sorted(documented - used)};"
        f" in code but not documented {sorted(used - documented)}"
    )


def test_every_rule_call_passes_a_literal_name():
    """``rule(name, ...)``: the name is a string literal (or the ``name`` parameter that
    ``expand._distribute`` forwards), so the parser above sees every name."""
    for module_name in RULE_MODULES:
        tree = ast.parse(_source(module_name))
        for node in ast.walk(tree):
            if isinstance(node, ast.Call) and getattr(node.func, "id", None) == "rule":
                first = node.args[0]
                assert isinstance(first, ast.Constant) or (
                    isinstance(first, ast.Name) and first.id == "name"
                ), f"{module_name}:{node.lineno}: rule() without a literal name"


def test_rule_names_are_unique_per_prefix():
    """No name is shared between modules (each prefix belongs to one module)."""
    prefixes = list(RULE_MODULES.values())
    assert len(prefixes) == len(set(prefixes))


# ---------------------------------------------------------------------------
# The trace switch
# ---------------------------------------------------------------------------


@pytest.fixture
def tracing():
    set_trace(True)
    yield
    set_trace(False)


def _s2_expression():
    """ω·q (vanishes: tangent ⊥ point) plus m·(q × (ω × q)) (the S² tangent recovery)."""
    q = S2("q")
    om = q.get_tangent_vector()
    m = Scalar("m")
    return q, om, Dot(Vector("a"), om) + Dot(om, q), Cross(q, Cross(om, q)) * m


def test_rule_returns_after_and_is_silent_when_off(caplog):
    assert not trace_enabled()
    a, b = Scalar("a"), Scalar("b")
    with caplog.at_level(logging.DEBUG, logger="geomech.rules"):
        assert rule("x.y", a, b) is b
    assert caplog.records == []


def test_set_trace_logs_each_firing(tracing, caplog):
    q, om, scalar, vector = _s2_expression()
    with caplog.at_level(logging.DEBUG, logger="geomech.rules"):
        assert str(full_simplify(scalar)) == "Dot(a,\\omega_{q})"
        assert str(full_simplify(vector)) == "\\omega_{q}*m"
    messages = [r.getMessage() for r in caplog.records if r.name == "geomech.rules"]
    assert "vector.dot_tangent_orthogonal: <v:\\omega_{q}, S2:q> ==> 0" in messages
    assert (
        "vector.s2_tangent_recovery: cross(S2:q, cross(v:\\omega_{q}, S2:q)) ==> v:\\omega_{q}"
        in messages
    )
    # a call that changes nothing is not logged
    assert all(" ==> " in m for m in messages)
    assert all(m.split(": ", 1)[0].count(".") == 1 for m in messages)


def test_trace_off_again(tracing, caplog):
    set_trace(False)
    _, _, scalar, vector = _s2_expression()
    with caplog.at_level(logging.DEBUG, logger="geomech.rules"):
        full_simplify(scalar)
        full_simplify(vector)
    assert not [r for r in caplog.records if r.name == "geomech.rules"]


def test_geomech_trace_env_var_in_a_subprocess():
    code = (
        "from geomech import S2, Dot, full_simplify\n"
        "q = S2('q'); full_simplify(Dot(q.get_tangent_vector(), q))\n"
    )
    env = dict(os.environ, GEOMECH_TRACE="1")
    env.pop("GEOMECH_DEBUG", None)
    out = subprocess.run(
        [sys.executable, "-c", code], cwd=ROOT, env=env, capture_output=True, text=True
    )
    assert out.returncode == 0, out.stderr
    assert (
        "geomech.rules: vector.dot_tangent_orthogonal: <v:\\omega_{q}, S2:q> ==> 0" in out.stderr
    )


# ---------------------------------------------------------------------------
# 4.6: the transformation submodules are not shadowed
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "module_name, func_name",
    [
        ("geomech.core.transformations.simplification", "simplify"),
        ("geomech.core.transformations.expansion", "expand"),
        ("geomech.core.transformations.pulling", "pull"),
        ("geomech.core.transformations.vector_identities", "vector_rules"),
        ("geomech.core.math.collecting", "collect"),
    ],
)
def test_submodule_is_the_module(module_name, func_name):
    module = importlib.import_module(module_name)
    ns = {}
    exec(f"import {module_name} as m", ns)
    assert ns["m"] is module
    assert isinstance(module, types.ModuleType)
    assert "rules" in module.__doc__.lower()
    # the old import path keeps working: the package exports the module's function under
    # its old name (the modules were renamed so that no function shadows its module)
    package = importlib.import_module(module_name.rsplit(".", 1)[0])
    assert getattr(package, func_name) is getattr(module, func_name)
    # the top-level name is the plain function
    assert getattr(geomech, func_name) is getattr(module, func_name)
    assert isinstance(getattr(geomech, func_name), types.FunctionType)


def test_callable_module_forwards():
    from geomech.core.transformations import (
        expand,
        pull,
        simplify,
        vector_rules,
    )
    from geomech.core.transformations import (
        full_simplify as fs,
    )

    q = S2("q")
    om = q.get_tangent_vector()
    expr = Dot(om, q) * Scalar("m")
    assert str(simplify(expr)) == str(geomech.simplify(expr))
    assert str(expand(expr)) == str(geomech.expand(expr))
    assert str(pull(expr)) == str(geomech.pull(expr))
    assert str(vector_rules(expr)) == str(geomech.vector_rules(expr))
    assert fs is geomech.full_simplify
