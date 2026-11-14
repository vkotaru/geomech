"""Regenerate ``examples/generated/``: NumPy modules for the nine systems of
``tests/vbl_systems.py``, one from ``to_standard_form`` (``<name>.py``, with ``energy``) and
one from ``linearize`` (``<name>_lin.py``).

    python scripts/regenerate_generated.py            # rewrite the files, delete stale ones
    python scripts/regenerate_generated.py --check    # exit 1 if anything would change

``--check`` reports a module that would change or is missing, a stale ``*.py`` (on disk
but no longer generated) and, inside a git work tree, a generated file that git does not
track (output that was never committed). The CI lint job runs ``--check``; the test suite
compares the file set through ``tests/test_codegen_generated.py``.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "examples" / "generated"

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from geomech import compute_eom, linearize, to_standard_form  # noqa: E402
from geomech.codegen import StateLayout, collect_params, module_source  # noqa: E402
from tests import vbl_systems  # noqa: E402

SYSTEMS = [
    "point_mass_r3",
    "spring_mass_r",
    "free_rigid_body",
    "rigid_pendulum_pivot",
    "spherical_pendulum_s2",
    "quadrotor",
    "double_rigid_pendulum",
    "elastic_spherical_pendulum",
    "quadrotor_with_load",
]


def _description(build) -> str:
    """The first paragraph of the builder's docstring, on one line per sentence chunk."""
    doc = (build.__doc__ or build.__name__).strip().split("\n\n")[0]
    return " ".join(doc.split()).replace("``", "")


def derive(name: str):
    """``(system, eom, layout)`` of the builder *name* in ``tests/vbl_systems.py``."""
    s = getattr(vbl_systems, name)()
    eom = compute_eom(s.L, s.dW, s.variables)
    exprs = [s.L] + [eq for _, eq in eom.values()]
    layout = StateLayout(s.variables, s.inputs, collect_params(exprs, s.variables, s.inputs))
    return s, eom, layout


def sources() -> dict[str, str]:
    """``{file name: source}`` of every generated module."""
    out = {}
    for name in SYSTEMS:
        s, eom, layout = derive(name)
        desc = _wrap(_description(getattr(vbl_systems, name)))
        sf = to_standard_form(eom, s.variables, s.inputs)
        lin = linearize(eom, s.variables, s.inputs)
        out[f"{name}.py"] = module_source(sf, layout, lagrangian=s.L, description=desc)
        out[f"{name}_lin.py"] = module_source(lin, layout, description=desc)
    return out


def _wrap(text: str, width: int = 88) -> str:
    lines, line = [], ""
    for word in text.split():
        if line and len(line) + 1 + len(word) > width:
            lines.append(line)
            line = word
        else:
            line = f"{line} {word}" if line else word
    return "\n".join(lines + [line])


def untracked(out: Path, names) -> list[str]:
    """The files of *names* in *out* that git does not track (``[]`` outside a git work
    tree or without git)."""
    try:
        proc = subprocess.run(
            ["git", "-C", str(out), "ls-files", "--", *names],
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError:
        return []
    if proc.returncode != 0:
        return []
    tracked = {Path(line).name for line in proc.stdout.splitlines()}
    return [n for n in names if n not in tracked]


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true", help="report drift, write nothing")
    args = parser.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)
    srcs = sources()
    drift, missing = [], []
    for fname, src in srcs.items():
        path = OUT / fname
        if path.exists() and path.read_text() == src:
            continue
        (drift if path.exists() else missing).append(fname)
        if not args.check:
            path.write_text(src)
    stale = sorted(p.name for p in OUT.glob("*.py") if p.name not in srcs)
    if not args.check:
        for name in stale:
            (OUT / name).unlink()
    new = untracked(OUT, sorted(srcs))
    problems = [
        ("out of date", drift),
        ("missing", missing),
        ("stale (no longer generated)", stale),
        ("not tracked by git", new),
    ]
    if args.check:
        failed = [(what, names) for what, names in problems if names]
        for what, names in failed:
            print(f"{what}: " + ", ".join(names))
        if failed:
            return 1
        print(f"checked {len(srcs)} modules in {OUT}")
        return 0
    if stale:
        print("deleted stale: " + ", ".join(stale))
    if new:
        print("not tracked by git (add them): " + ", ".join(new))
    print(f"wrote {len(srcs)} modules in {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
