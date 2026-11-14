"""Numeric snapshot of every pipeline output, for changes that rewrite pinned equation strings.

A simplifier change may print an equation differently; it must not change its value. This
script records, at one commit, every ``compute_eom`` / ``to_standard_form`` / ``linearize``
/ ``integrate_by_parts`` call of the test suite plus a chain over every builder system, and
later replays the same calls with the working tree and compares.

  python specs/reviews/scripts/pinned_snapshot.py snapshot --out base.json
  python specs/reviews/scripts/pinned_snapshot.py compare --base base.json [--out cmp.json]

``snapshot``
  1. Runs the test suite (``tests/``, quiet) with the three pipeline functions wrapped; each
     top-level call pickles its arguments together with the id of the test that made it.
     Identical calls (same function, same pickled arguments) are kept once with all ids.
  2. Adds a chain ``compute_eom → to_standard_form → linearize`` for every builder system
     of ``tests/vbl_systems.py``, ``tests/phase1_systems.py`` and ``tests/test_eom_oracle.py``.
  3. Replays every call outside pytest and stores, per call: ``str`` of every output
     (equation, standard-form block, linearization block), its value at five fixed states
     (seeds 1001..1005), and for ``compute_eom`` the oracle ``check_eom`` (seed 0, 5 trials),
     or the exception class and message.

``compare``
  Replays the stored calls with the current code and reports, per output, whether the
  string changed, the maximum relative difference of the values over the five states
  (``max |new − old| / max(1, max |old|)``), and ``check_eom`` before/after. Any value
  difference above ``1e-12`` or a ``check_eom`` above ``max(old, 1e-6)`` is a failure,
  with one exception: a ``linearize`` block of an S² system may differ as a matrix when
  every row ``−M δẋ + A δx + B δu`` is identical along admissible directions (see
  ``linearization_rhs``); those are listed separately.

Fixed states. The configuration variables get a point of the manifold
(``geomech.numeric.sample_state``: |q| = 1, ω ⊥ q, ω̇ ⊥ q, R ∈ SO(3)); every other leaf, time
derivative chain and variation symbol gets a value drawn from a generator seeded by its
printed name, so the value does not depend on the order in which a tree is traversed and
two differently written trees of the same expression see the same numbers. Leaves respect
their flags (S² / unit norm, SO(3) / orthogonal, symmetric, skew-symmetric, zero, identity;
a ``Scalar`` with ``value`` uses it).
"""

from __future__ import annotations

import argparse
import base64
import functools
import hashlib
import json
import os
import pickle
import sys
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402
import pytest  # noqa: E402

import geomech  # noqa: E402
import geomech.core.math  # noqa: E402
import geomech.core.math.ibp  # noqa: E402
import geomech.dynamics  # noqa: E402
import geomech.dynamics.eom  # noqa: E402
import geomech.dynamics.standard_form  # noqa: E402
import geomech.dynamics.vbl  # noqa: E402
from geomech import S2, SO3, Matrix, Scalar  # noqa: E402
from geomech.core.base.expressions import IdentityMatrix  # noqa: E402
from geomech.core.operations.calculus import TimeDerivative, Variation  # noqa: E402
from geomech.numeric import (  # noqa: E402
    Evaluator,
    _dstate_values,
    check_eom,
    linearized_rhs,
    rand_rot,
    sample_direction,
    sample_state,
)

FUNCS = {
    "compute_eom": geomech.dynamics.eom.compute_eom,
    "to_standard_form": geomech.dynamics.standard_form.to_standard_form,
    "linearize": geomech.dynamics.vbl.linearize,
    "integrate_by_parts": geomech.core.math.ibp.integrate_by_parts,
}
STATE_SEEDS = (1001, 1002, 1003, 1004, 1005)
VALUE_RTOL = 1e-12
PINNED_FILES = (
    "tests/test_manifold_pipelines.py",
    "tests/test_standard_form.py",
    "tests/test_s2_variation_fix.py",
    "tests/test_collect_fix.py",
    "tests/test_scalar_adversarial.py",
    "tests/test_eom_oracle.py",
)


# ---------------------------------------------------------------------------
# Recording (inside pytest)
# ---------------------------------------------------------------------------


class Recorder:
    def __init__(self):
        self.current = "collect"
        self.depth = 0
        self.calls = {}  # (func, sha) -> {"func", "args", "ids"}

    def wrap(self, name, func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            top = self.depth == 0
            if top:
                try:
                    blob = pickle.dumps((args, kwargs))
                except Exception:  # noqa: BLE001 - unpicklable call (a lambda in kwargs)
                    blob = None
                if blob is not None:
                    key = (name, hashlib.sha256(blob).hexdigest())
                    entry = self.calls.setdefault(key, {"func": name, "args": blob, "ids": []})
                    if self.current not in entry["ids"]:
                        entry["ids"].append(self.current)
            self.depth += 1
            try:
                return func(*args, **kwargs)
            finally:
                self.depth -= 1

        return wrapper

    # pytest hooks (tryfirst: the id is set before setup_method/fixtures run)
    @pytest.hookimpl(tryfirst=True)
    def pytest_runtest_setup(self, item):
        self.current = item.nodeid

    def pytest_runtest_teardown(self, item):
        self.current = "between-tests"


def _install(recorder):
    wrapped = {n: recorder.wrap(n, f) for n, f in FUNCS.items()}
    mods = [
        geomech,
        geomech.core.math,
        geomech.core.math.ibp,
        geomech.dynamics,
        geomech.dynamics.eom,
        geomech.dynamics.standard_form,
        geomech.dynamics.vbl,
    ]
    for mod in mods:
        for n, f in FUNCS.items():
            if getattr(mod, n, None) is f:
                setattr(mod, n, wrapped[n])


def _uninstall():
    mods = [
        geomech,
        geomech.core.math,
        geomech.core.math.ibp,
        geomech.dynamics,
        geomech.dynamics.eom,
        geomech.dynamics.standard_form,
        geomech.dynamics.vbl,
    ]
    for mod in mods:
        for n, f in FUNCS.items():
            if hasattr(mod, n):
                setattr(mod, n, f)


def record_suite():
    rec = Recorder()
    _install(rec)
    try:
        code = pytest.main(
            [os.path.join(ROOT, "tests"), "-q", "-p", "no:cacheprovider", "-x", "--no-header"],
            plugins=[rec],
        )
    finally:
        _uninstall()
    print(f"[snapshot] pytest exit code {code}; {len(rec.calls)} distinct pipeline calls")
    return list(rec.calls.values())


def builder_systems():
    from tests import phase1_systems, test_eom_oracle, vbl_systems

    seen, out = set(), []

    def add(mod, f):
        key = f"{f.__module__}.{f.__name__}"
        if key not in seen:
            seen.add(key)
            out.append((key, f))

    for name in dir(vbl_systems):
        f = getattr(vbl_systems, name)
        if callable(f) and getattr(f, "__module__", "") == vbl_systems.__name__:
            if name not in ("System",) and not name.startswith("_"):
                add(vbl_systems, f)
    for f in test_eom_oracle.SYSTEMS:
        add(test_eom_oracle, f)
    for lst in ("S2_VARIATION_SYSTEMS", "COLLECT_SYSTEMS", "ADVERSARIAL_SYSTEMS"):
        for f in getattr(phase1_systems, lst):
            add(phase1_systems, f)
    for name in (
        "stiff_point_mass",
        "so3_outer_product_on_velocity",
        "s2_outer_product_on_velocity",
    ):
        add(phase1_systems, getattr(phase1_systems, name))
    return out


# ---------------------------------------------------------------------------
# Fixed numeric states
# ---------------------------------------------------------------------------


def _children(node):
    nodes = getattr(node, "nodes", None)
    if nodes:
        return list(nodes)
    return []


def _name_rng(seed, key, kind):
    h = hashlib.sha256(f"{seed}|{kind}|{key}".encode()).hexdigest()
    return np.random.default_rng(int(h[:16], 16))


def _base_of(node):
    while isinstance(node, TimeDerivative):
        node = node.expr
    if isinstance(node, Variation):
        node = node.expr
    return node


def _leaf_value(node, rng):
    if isinstance(node, Scalar):
        return float(rng.normal())
    if isinstance(node, Matrix):
        if isinstance(node, SO3) or node.is_orthogonal:
            return rand_rot(rng)
        A = rng.normal(size=(3, 3))
        if node.is_symmetric:
            return A + A.T
        if node.is_skew_symmetric:
            return A - A.T
        return A
    v = rng.normal(size=(3, 1))
    if isinstance(node, S2) or node.is_unit_norm:
        v = v / np.linalg.norm(v)
    return v


def _plain_value(node, rng):
    """Value of a derived symbol (time derivative / variation): no structure."""
    if isinstance(node, Scalar) or getattr(node, "type", None) == geomech.ExprType.SCALAR:
        return float(rng.normal())
    if isinstance(node, Matrix) or getattr(node, "type", None) == geomech.ExprType.MATRIX:
        return rng.normal(size=(3, 3))
    return rng.normal(size=(3, 1))


def fill_env(env, exprs, seed):
    stack = [e for e in exprs if e is not None and hasattr(e, "type")]
    while stack:
        n = stack.pop()
        if isinstance(n, (TimeDerivative, Variation)):
            key = str(n)
            if key not in env:
                env[key] = _plain_value(_base_of(n), _name_rng(seed, key, "derived"))
            # Evaluator never looks inside a derived symbol
            continue
        kids = _children(n)
        if kids:
            stack.extend(kids)
            continue
        if n.is_zero or (isinstance(n, Matrix) and n.flags.is_identity):
            continue
        if isinstance(n, Scalar) and n.value is not None:
            continue
        key = str(n)
        if key not in env:
            env[key] = _leaf_value(n, _name_rng(seed, key, type(n).__name__))
    return env


def state_env(variables, seed):
    if variables is None:
        return {}
    try:
        return sample_state(variables, [], np.random.default_rng(seed))
    except Exception:  # noqa: BLE001 - odd variables object in an error-path test
        return {}


def value_of(expr, env):
    try:
        v = Evaluator(env=dict(env), rng=np.random.default_rng(0))(expr)
    except Exception as e:  # noqa: BLE001
        return {"error": f"{type(e).__name__}: {e}"[:300]}
    return np.asarray(v, dtype=float).reshape(-1).tolist()


# ---------------------------------------------------------------------------
# Replay
# ---------------------------------------------------------------------------


def _refresh(obj, seen=None):
    """Give unpickled identity leaves the current IdentityMatrix flags."""
    if seen is None:
        seen = set()
    if id(obj) in seen:
        return
    seen.add(id(obj))
    if isinstance(obj, Matrix) and obj.flags is not None and obj.flags.is_identity:
        obj.attr = list(IdentityMatrix.attr)
        obj.flags = IdentityMatrix.flags
        return
    if isinstance(obj, dict):
        for v in obj.values():
            _refresh(v, seen)
        return
    if isinstance(obj, (list, tuple)):
        for v in obj:
            _refresh(v, seen)
        return
    for attr in ("nodes", "expr", "left", "right"):
        ch = getattr(obj, attr, None) if hasattr(obj, "type") else None
        if ch is not None:
            _refresh(ch, seen)
    for attr in ("scalars", "vectors", "matrices"):
        ch = getattr(obj, attr, None)
        if ch is not None:
            _refresh(ch, seen)
    if hasattr(obj, "__dict__") and not hasattr(obj, "type"):
        for v in vars(obj).values():
            _refresh(v, seen)


def _variables_of(func, args, kwargs):
    if func == "integrate_by_parts":
        return None
    if func == "compute_eom":
        return args[2] if len(args) > 2 else kwargs.get("variables")
    return args[1] if len(args) > 1 else kwargs.get("variables")


def outputs_of(func, out):
    """``{label: expr}`` of every output tree."""
    res = {}
    if func == "compute_eom":
        for k, (_, eq) in out.items():
            res[f"E[{k}]"] = eq
    elif func == "to_standard_form":
        for k, sf in out.items():
            for a, b in sf.M.items():
                res[f"M[{k}][{a}]"] = b
            res[f"f[{k}]"] = sf.f
            for a, b in sf.G.items():
                res[f"G[{k}][{a}]"] = b
    elif func == "integrate_by_parts":
        res["ibp"] = out
    else:
        for name, blocks in (("M", out.M), ("A", out.A), ("B", out.B)):
            for row, cols in blocks.items():
                for col, b in cols.items():
                    res[f"{name}[{row}][{col}]"] = b
    return res


def evaluate_outputs(outs, variables, extra=()):
    vals = {k: [] for k in outs}
    for seed in STATE_SEEDS:
        env = state_env(variables, seed)
        fill_env(env, list(outs.values()) + list(extra), seed)
        for k, e in outs.items():
            vals[k].append(value_of(e, env))
    return vals


def linearization_rhs(lin, variables, inputs, outs):
    """``rhs[row]`` = −M δẋ + A δx + B δu at the fixed states, for an admissible direction.

    The direction (``geomech.numeric.sample_direction``: ξ ⊥ q, δω = ξ × ω + w with w ⊥ q,
    ...) satisfies the linearized constraints. On S² the blocks themselves are not unique:
    two equations that agree on the sphere (``E + c (q·ω)``) have variations that differ
    by ``c (δq·ω + q·δω)``, which is zero only along admissible directions. The value of
    the whole row along such a direction is what ``finite_difference_check`` tests and
    what must not change.
    """
    vals = {}
    for i, seed in enumerate(STATE_SEEDS):
        env = state_env(variables, seed)
        fill_env(env, list(outs.values()), seed)
        try:
            ev = Evaluator(env=dict(env), rng=np.random.default_rng(0))
            raw, pert = sample_direction(
                variables, inputs, ev.env, np.random.default_rng(seed + 7)
            )
            dvals = _dstate_values(lin, variables, ev.env, pert)
            for row in lin.rows:
                v = linearized_rhs(lin, row, ev, pert, dvals)
                vals.setdefault(f"rhs[{row}]", []).append(np.asarray(v, float).tolist())
        except Exception as e:  # noqa: BLE001
            vals.setdefault("rhs:error", []).append(f"{type(e).__name__}: {e}"[:200])
    return vals


def run_call(func, args, kwargs):
    t0 = time.perf_counter()
    try:
        out = FUNCS[func](*args, **kwargs)
    except Exception as e:  # noqa: BLE001
        return None, {"error": f"{type(e).__name__}: {str(e)[:300]}"}, time.perf_counter() - t0
    return out, None, time.perf_counter() - t0


def replay_entry(entry):
    args, kwargs = pickle.loads(entry["args"])
    _refresh((args, kwargs))
    func = entry["func"]
    variables = _variables_of(func, args, kwargs)
    out, err, _ = run_call(func, args, kwargs)
    rec = {"func": func, "ids": entry["ids"]}
    if err is not None:
        rec.update(err)
        return rec
    outs = outputs_of(func, out)
    rec["strings"] = {k: str(v) for k, v in outs.items()}
    rec["values"] = evaluate_outputs(outs, variables)
    if func == "linearize":
        inputs = args[2] if len(args) > 2 else kwargs.get("inputs") or []
        rec["rhs"] = linearization_rhs(out, variables, list(inputs or []), outs)
    if func == "compute_eom":
        L, dW = args[0], args[1] if len(args) > 1 else kwargs.get("inf_work")
        params = kwargs.get("params", args[3] if len(args) > 3 else ())
        try:
            rec["check_eom"] = float(check_eom(L, dW, variables, out, params=params))
        except Exception as e:  # noqa: BLE001
            rec["check_eom"] = f"{type(e).__name__}: {str(e)[:200]}"
    return rec


def replay_chain(key, build):
    s = build()
    rec = {"func": "chain", "ids": [f"system:{key}"]}
    eom, err, dt = run_call("compute_eom", (s.L, s.dW, s.variables), {})
    if err is not None:
        rec.update(err)
        return rec
    outs = outputs_of("compute_eom", eom)
    try:
        rec["check_eom"] = float(check_eom(s.L, s.dW, s.variables, eom, inputs=s.inputs))
    except Exception as e:  # noqa: BLE001
        rec["check_eom"] = f"{type(e).__name__}: {str(e)[:200]}"
    sf, err2, _ = run_call("to_standard_form", (eom, s.variables, list(s.inputs)), {})
    if err2 is None:
        outs.update({"sf:" + k: v for k, v in outputs_of("to_standard_form", sf).items()})
    else:
        rec["sf_error"] = err2["error"]
    lin, err3, _ = run_call("linearize", (eom, s.variables, list(s.inputs)), {})
    if err3 is None:
        outs.update({"lin:" + k: v for k, v in outputs_of("linearize", lin).items()})
    else:
        rec["lin_error"] = err3["error"]
    rec["strings"] = {k: str(v) for k, v in outs.items()}
    rec["values"] = evaluate_outputs(outs, s.variables)
    if err3 is None:
        lin_outs = {k: v for k, v in outs.items() if k.startswith("lin:")}
        rec["rhs"] = linearization_rhs(lin, s.variables, list(s.inputs), lin_outs)
    rec["seconds"] = dt
    return rec


def snapshot(out_path):
    calls = record_suite()
    records = []
    for entry in calls:
        rec = replay_entry(entry)
        rec["args"] = base64.b64encode(entry["args"]).decode()
        records.append(rec)
    for key, build in builder_systems():
        records.append(replay_chain(key, build))
    head = os.popen(f"git -C {ROOT} rev-parse --short HEAD").read().strip()
    with open(out_path, "w") as fh:
        json.dump({"commit": head, "records": records}, fh)
    print(f"[snapshot] {len(records)} records written to {out_path}")


# ---------------------------------------------------------------------------
# Compare
# ---------------------------------------------------------------------------


def _max_rel(old, new):
    if isinstance(old, dict) or isinstance(new, dict):
        if isinstance(old, dict) and isinstance(new, dict):
            return 0.0
        return float("inf")
    a, b = np.asarray(old, dtype=float), np.asarray(new, dtype=float)
    if a.shape != b.shape:
        return float("inf")
    return float(
        np.max(np.abs(a - b), initial=0.0) / max(1.0, float(np.max(np.abs(a), initial=0.0)))
    )


def compare(base_path, out_path=None):
    from tests import phase1_systems, test_eom_oracle, vbl_systems  # noqa: F401

    base = json.load(open(base_path))
    builders = dict(builder_systems())
    changed, failures, same, errors_kept = [], [], 0, 0
    blocks_equivalent = []
    eom_checks = []
    for old in base["records"]:
        if old["func"] == "chain":
            key = old["ids"][0][len("system:") :]
            new = replay_chain(key, builders[key])
        else:
            entry = {"func": old["func"], "ids": old["ids"], "args": base64.b64decode(old["args"])}
            new = replay_entry(entry)
        tag = old["ids"][0]
        if "error" in old:
            if "error" in new and new["error"].split(":")[0] == old["error"].split(":")[0]:
                errors_kept += 1
            else:
                failures.append(
                    {
                        "ids": old["ids"],
                        "why": "error changed",
                        "old": old["error"],
                        "new": new.get("error", "now succeeds"),
                    }
                )
            continue
        if "error" in new:
            failures.append({"ids": old["ids"], "why": "now raises", "new": new["error"]})
            continue
        if set(old["strings"]) != set(new["strings"]):
            failures.append(
                {
                    "ids": old["ids"],
                    "why": "output keys differ",
                    "old": sorted(old["strings"]),
                    "new": sorted(new["strings"]),
                }
            )
            continue
        rhs_ok = True
        for k in old.get("rhs", {}):
            if k == "rhs:error":
                continue
            d = max(_max_rel(o, n) for o, n in zip(old["rhs"][k], new.get("rhs", {}).get(k, [])))
            if d > VALUE_RTOL or k not in new.get("rhs", {}):
                rhs_ok = False
                failures.append(
                    {"ids": old["ids"], "why": f"linearization row differs {d:.3e}", "key": k}
                )
        for k in old["strings"]:
            diffs = [_max_rel(o, n) for o, n in zip(old["values"][k], new["values"][k])]
            d = max(diffs, default=0.0)
            is_block = old["func"] == "linearize" or k.startswith("lin:")
            if d > VALUE_RTOL and is_block and rhs_ok and old.get("rhs"):
                blocks_equivalent.append({"ids": old["ids"], "key": k, "max_rel_diff": d})
            elif d > VALUE_RTOL:
                failures.append(
                    {
                        "ids": old["ids"],
                        "why": f"value differs {d:.3e}",
                        "key": k,
                        "old": old["strings"][k],
                        "new": new["strings"][k],
                    }
                )
            if old["strings"][k] != new["strings"][k]:
                changed.append(
                    {
                        "ids": old["ids"],
                        "key": k,
                        "old": old["strings"][k],
                        "new": new["strings"][k],
                        "max_rel_diff": d,
                    }
                )
            else:
                same += 1
        if "check_eom" in old:
            o, n = old["check_eom"], new.get("check_eom")
            eom_checks.append((tag, o, n))
            if isinstance(o, float) and (not isinstance(n, float) or n > max(o, 1e-6)):
                failures.append({"ids": old["ids"], "why": f"check_eom {o} -> {n}"})
    print(f"[compare] base commit {base['commit']}: {len(base['records'])} records")
    print(
        f"[compare] outputs unchanged: {same}; strings changed: {len(changed)};"
        f" expected errors kept: {errors_kept}; failures: {len(failures)}"
    )
    worst = max((c["max_rel_diff"] for c in changed), default=0.0)
    print(f"[compare] worst value difference among changed strings: {worst:.3e}")
    print(
        f"[compare] linearization blocks that differ as matrices but give identical rows"
        f" along admissible directions: {len(blocks_equivalent)}"
        f" ({len({tuple(b['ids']) for b in blocks_equivalent})} calls)"
    )
    better = sum(
        1 for _, o, n in eom_checks if isinstance(o, float) and isinstance(n, float) and n < o
    )
    print(f"[compare] check_eom: {len(eom_checks)} compared, {better} lower than before")
    for f in failures[:40]:
        print("  FAIL", json.dumps(f)[:600])
    if out_path:
        with open(out_path, "w") as fh:
            json.dump(
                {
                    "changed": changed,
                    "failures": failures,
                    "check_eom": eom_checks,
                    "blocks_equivalent": blocks_equivalent,
                },
                fh,
                indent=1,
            )
    return 1 if failures else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("snapshot")
    s.add_argument("--out", required=True)
    c = sub.add_parser("compare")
    c.add_argument("--base", required=True)
    c.add_argument("--out")
    a = ap.parse_args()
    if a.cmd == "snapshot":
        snapshot(a.out)
        return 0
    return compare(a.base, a.out)


if __name__ == "__main__":
    sys.exit(main())
