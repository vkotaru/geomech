"""Regenerate ``tests/golden/``: the reference output of every system of ``geomech.systems``.

    python scripts/regenerate_golden.py            # rewrite the files, delete stale ones
    python scripts/regenerate_golden.py --check    # exit 1 if a stored equation changed

One JSON file per system (``<name>.json`` at the builder's defaults, and
``<name>-2.json`` for a parametric system at size 2) holds, for every equation of
``compute_eom`` and every block of ``to_standard_form`` (``M`` per acceleration, ``f``,
``G`` per input), the printed string (``str``) and its value at ``N_STATES`` fixed states,
and the system's default ``params`` and ``x0``.
The states are stored in the file (every leaf the expressions use: configuration,
velocity, acceleration, inputs and random parameter values), drawn once from the seeds
``SEEDS`` with ``geomech.numeric.sample_state``, so a check does not depend on the sampler.

A check re-derives every system and evaluates the new expressions at the stored states:

- a value that differs by more than ``RTOL · max(1, |stored value|)`` is a change of the
  equations and fails, and so does a change of a default parameter or of ``x0``;
- a printed string that differs while every value agrees is a change of form only; it
  fails too, with the instruction to check the numbers and regenerate the files.

``RTOL = 1e-9``: two printings of one equation evaluate to the same value up to the
rounding of a few dozen floating-point operations (about 1e-14 relative), and a changed
term moves a value by its own size (of order one at the random states).
"""

from __future__ import annotations

import argparse
import inspect
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "tests" / "golden"

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from geomech import compute_eom, to_standard_form  # noqa: E402
from geomech.numeric import Evaluator, sample_state  # noqa: E402
from geomech.systems import get, list_systems  # noqa: E402

SEEDS = (11, 12, 13)
N_STATES = len(SEEDS)
RTOL = 1e-9
REGENERATE = "python scripts/regenerate_golden.py"


def entries() -> list[tuple[str, str, dict]]:
    """``(file stem, system name, builder kwargs)`` of every golden file."""
    out = []
    for info in list_systems():
        out.append((info.name, info.name, {}))
        for param in info.parameters:
            out.append((f"{info.name}-2", info.name, {param: 2}))
    return out


def blocks(eom: dict, sf: dict) -> dict[str, Any]:
    """``{label: expression}`` of every equation and standard-form block, in a fixed order."""
    out: dict[str, Any] = {}
    for key, (_, eq) in eom.items():
        out[f"eom[{key}]"] = eq
    for key, s in sf.items():
        for acc, m in s.M.items():
            out[f"M[{key}][{acc}]"] = m
        out[f"f[{key}]"] = s.f
        for u, gk in s.G.items():
            out[f"G[{key}][{u}]"] = gk
    return out


def derive(name: str, kwargs: dict):
    s = get(name, **kwargs)
    eom = compute_eom(s.L, s.dW, s.variables)
    sf = to_standard_form(eom, s.variables, list(s.inputs))
    return s, blocks(eom, sf)


def _value(expr, env) -> list[float]:
    val = Evaluator(env=dict(env), rng=np.random.default_rng(0))(expr)
    return np.atleast_1d(np.asarray(val, dtype=float)).reshape(-1).tolist()


def _to_json(v):
    return float(v) if np.isscalar(v) else np.asarray(v, dtype=float).tolist()


def _from_json(v):
    return float(v) if isinstance(v, (int, float)) else np.asarray(v, dtype=float)


def _numeric(env: dict) -> dict:
    return {k: _from_json(v) for k, v in env.items()}


def states(s, exprs) -> list[dict]:
    """The stored states: every leaf of *exprs* at each seed of ``SEEDS``."""
    out = []
    for seed in SEEDS:
        rng = np.random.default_rng(seed)
        ev = Evaluator(env=sample_state(s.variables, list(s.inputs), rng), rng=rng)
        for e in exprs:
            ev(e)
        out.append({k: _to_json(v) for k, v in sorted(ev.env.items())})
    return out


def golden(name: str, kwargs: dict) -> dict:
    """The content of one golden file."""
    s, exprs = derive(name, kwargs)
    envs = states(s, [s.L, *exprs.values()])
    num = [_numeric(env) for env in envs]
    return {
        "system": name,
        "kwargs": kwargs,
        "rtol": RTOL,
        "params": {k: _to_json(v) for k, v in sorted(s.params.items())},
        "x0": np.asarray(s.x0, dtype=float).tolist(),
        "states": envs,
        "blocks": {
            label: {"string": str(e), "values": [_value(e, env) for env in num]}
            for label, e in exprs.items()
        },
    }


def _close(got, ref) -> bool:
    got, ref = np.asarray(got, dtype=float), np.asarray(ref, dtype=float)
    if got.shape != ref.shape:
        return False
    return not ref.size or float(np.max(np.abs(got - ref))) <= RTOL * max(
        1.0, float(np.max(np.abs(ref)))
    )


def compare_defaults(stored: dict, s) -> list[str]:
    """Changes of the default ``params`` and ``x0`` of the system *s* against *stored*."""
    out = []
    old = stored.get("params")
    if old is None or set(old) != set(s.params):
        out.append(f"params: stored {sorted(old or {})}, now {sorted(s.params)}")
    else:
        out += [
            f"params[{k!r}]: stored {old[k]}, now {_to_json(s.params[k])}"
            for k in sorted(old)
            if not _close(_to_json(s.params[k]), old[k])
        ]
    if "x0" not in stored or not _close(s.x0, stored["x0"]):
        out.append(f"x0: stored {stored.get('x0')}, now {np.asarray(s.x0).tolist()}")
    return out


def compare(stored: dict, exprs: dict) -> tuple[list[str], list[str]]:
    """``(value changes, string-only changes)`` of re-derived *exprs* against *stored*.

    The tolerance is always ``RTOL``; a file that carries another ``rtol`` is reported as a
    change (a hand-edited tolerance cannot loosen the check).
    """
    envs = [_numeric(env) for env in stored["states"]]
    rtol = RTOL
    values, strings = [], []
    if stored.get("rtol") != RTOL:
        values.append(f"rtol: the file says {stored.get('rtol')!r}, the check uses {RTOL!r}")
    old, new = stored["blocks"], exprs
    for label in sorted(set(old) ^ set(new)):
        where = "no longer derived" if label in old else "new"
        values.append(f"{label}: {where}")
    for label in [k for k in old if k in new]:
        e = new[label]
        for i, (env, ref) in enumerate(zip(envs, old[label]["values"])):
            try:
                got = np.asarray(_value(e, env))
            except Exception as exc:  # noqa: BLE001  (reported as a change)
                values.append(f"{label}: cannot evaluate at state {i}: {exc}")
                break
            ref = np.asarray(ref)
            if got.shape != ref.shape:
                values.append(f"{label}: shape {got.shape} != stored {ref.shape}")
                break
            err = float(np.max(np.abs(got - ref))) if ref.size else 0.0
            if err > rtol * max(1.0, float(np.max(np.abs(ref)))):
                values.append(
                    f"{label}: value changed at state {i} (|Δ| = {err:.3e}):"
                    f" stored {old[label]['string']!r}, now {str(e)!r}"
                )
                break
        else:
            if str(e) != old[label]["string"]:
                strings.append(f"{label}: {old[label]['string']!r} -> {str(e)!r}")
    return values, strings


def check_one(stem: str, name: str, kwargs: dict) -> tuple[list[str], list[str]]:
    path = OUT / f"{stem}.json"
    if not path.exists():
        return [f"{stem}: missing golden file {path.name}"], []
    stored = json.loads(path.read_text())
    s, exprs = derive(name, kwargs)
    values, strings = compare(stored, exprs)
    values = compare_defaults(stored, s) + values
    return [f"{stem}: {m}" for m in values], [f"{stem}: {m}" for m in strings]


def _dump(data: dict) -> str:
    return json.dumps(data, indent=1, ensure_ascii=False) + "\n"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=inspect.cleandoc(__doc__).splitlines()[0])
    parser.add_argument("--check", action="store_true", help="report changes, write nothing")
    args = parser.parse_args(argv)
    todo = entries()
    stems = {stem for stem, _, _ in todo}
    stale = sorted(p.name for p in OUT.glob("*.json") if p.stem not in stems)
    if args.check:
        values, strings = [], []
        for stem, name, kwargs in todo:
            v, s = check_one(stem, name, kwargs)
            values += v
            strings += s
        for m in values:
            print("value changed:", m)
        for m in strings:
            print("printed form changed (values agree):", m)
        if stale:
            print("stale (no longer generated): " + ", ".join(stale))
        if values or strings or stale:
            if strings and not values:
                print(f"check the numbers above, then regenerate with `{REGENERATE}`")
            return 1
        print(f"checked {len(todo)} golden files in {OUT}")
        return 0
    OUT.mkdir(parents=True, exist_ok=True)
    for stem, name, kwargs in todo:
        (OUT / f"{stem}.json").write_text(_dump(golden(name, kwargs)))
    for fname in stale:
        (OUT / fname).unlink()
    if stale:
        print("deleted stale: " + ", ".join(stale))
    print(f"wrote {len(todo)} golden files in {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
