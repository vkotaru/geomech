"""Every rewrite goes through a named rule (phase 4 review, finding F6).

``tests/test_rule_names.py`` checks that the rule names in each module docstring and the
names passed to ``rule(...)`` are the same set, but it cannot see a rewrite that returns
a new node without calling ``rule()``. This guard runs the transformations on random
well-typed trees of the engine fuzzer with rule tracing on, and asserts that every output
that prints differently from its input logged at least one rule line.
"""

from __future__ import annotations

import logging
import sys
import warnings
from pathlib import Path

import pytest

from geomech import ExprType, collect, expand, full_simplify, pull, set_trace, simplify
from geomech.core.transformations.vector_identities import vector_rules

ROOT = Path(__file__).resolve().parent.parent
CASES = 1500


@pytest.fixture
def rules_log():
    lines = []

    class _H(logging.Handler):
        def emit(self, record):
            lines.append(record.getMessage())

    handler = _H()
    log = logging.getLogger("geomech.rules")
    log.addHandler(handler)
    set_trace(True)
    yield lines
    set_trace(False)
    log.removeHandler(handler)


def _fuzz_trees(seed):
    np = pytest.importorskip("numpy")
    sys.path.insert(0, str(ROOT / "specs/reviews/scripts"))
    try:
        import engine_fuzz as F
    finally:
        sys.path.pop(0)
    rng = np.random.default_rng(seed)
    gen = F.Gen(rng, F.Leaves())
    return [gen.any(int(rng.integers(1, 5))) for _ in range(CASES)]


def _first_vector_leaf(expr):
    stack = [expr]
    while stack:
        node = stack.pop()
        children = getattr(node, "nodes", None)
        if children:
            stack.extend(reversed(children))
        elif node.type == ExprType.VECTOR:
            return node
    return None


def test_every_rewrite_of_the_transformations_logs_a_rule(rules_log):
    unnamed, checked = [], 0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for expr in _fuzz_trees(seed=4):
            for name, fn in (
                ("expand", expand),
                ("pull", pull),
                ("vector_rules", vector_rules),
                ("simplify", simplify),
                ("full_simplify", full_simplify),
            ):
                rules_log.clear()
                try:
                    out = fn(expr)
                except Exception:  # noqa: BLE001 - a raising input rewrites nothing
                    continue
                checked += 1
                if str(out) != str(expr) and not rules_log:
                    unnamed.append(f"{name}: {expr!r} ==> {out!r}")
    assert checked > 5000
    assert not unnamed, unnamed[:5]


def test_every_rewrite_of_collect_logs_a_rule(rules_log):
    unnamed, checked = [], 0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for expr in _fuzz_trees(seed=5):
            if expr.type != ExprType.SCALAR:
                continue
            target = _first_vector_leaf(expr)
            if target is None:
                continue
            rules_log.clear()
            try:
                out = collect(expr, target)
            except Exception:  # noqa: BLE001
                continue
            checked += 1
            if str(out) != str(expr) and not rules_log:
                unnamed.append(f"collect[{target}]: {expr!r} ==> {out!r}")
    assert checked > 200
    assert not unnamed, unnamed[:5]
