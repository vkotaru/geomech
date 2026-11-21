"""Derivation time, size and numerical check of every catalog entry, for the appendix.

    python sim/catalog_numbers.py

Writes ``sim/verification_rows.tex``: one row per entry of ``geomech.systems`` (in the
order of the appendix) and one per model of ``geomech.indexed.systems``. For a library
entry: the median wall time of ``compute_eom`` over three runs, the number of additive
terms of its equations, and the ``check_eom`` error. For an indexed model: the time of
``compute_indexed_eom`` (derived once for every n) and the largest ``check_eom``
error of the equations instantiated for n = 1, ..., 6 against the Lagrangian unrolled for
that n. The entry names are read from
``07_appendix.tex``.
"""

from __future__ import annotations

import re
import statistics
import time
from pathlib import Path

import geomech.indexed.systems as I
import geomech.systems as S
from geomech import check_eom, compute_eom
from geomech.core.transformations.utils import addends

HERE = Path(__file__).resolve().parent
APPENDIX = HERE.parent / "07_appendix.tex"
SIZES = range(1, 7)

INDEXED = [
    ("spherical_chain", "Spherical chain", r"$(\Sph)^n$", "$2n$"),
    ("rigid_chain", "Rigid chain", r"$(\SOthree)^n$", "$3n$"),
    (
        "quadrotor_flexible_cable",
        "Quadrotor, flexible cable",
        r"$\RealR^3 \times \SOthree \times (\Sph)^n$",
        "$6 + 2n$",
    ),
    (
        "multi_quadrotor_point_mass",
        "$n$ quadrotors, point mass",
        r"$\RealR^3 \times (\SOthree \times \Sph)^n$",
        "$3 + 5n$",
    ),
    (
        "multi_quadrotor_rigid_load",
        "$n$ quadrotors, rigid body",
        r"$\RealR^3 \times \SOthree \times (\SOthree \times \Sph)^n$",
        "$6 + 5n$",
    ),
]


def entries():
    """(title, library name) in the order of the appendix; without the paper source
    (``07_appendix.tex`` next to ``sim/``), every library system under its own name."""
    if not APPENDIX.exists():
        return [(name, name) for name in S.names()]
    text = APPENDIX.read_text()
    return re.findall(r"\\subsubsection\{(.+?)\}\n\\builder\{(\w+)\}", text)


def timed(f, runs=3):
    out, times = None, []
    for _ in range(runs):
        t = time.perf_counter()
        out = f()
        times.append(time.perf_counter() - t)
    return out, statistics.median(times)


def sci(x):
    m, e = f"{x:.1e}".split("e")
    return f"${m} \\times 10^{{{int(e)}}}$"


def secs(t):
    return f"{t:.3f}" if t < 10 else f"{t:.1f}"


_TEX = {
    "×": r" \times ",
    "SO(3)": r"\SOthree",
    "S²": r"\Sph",
    "S¹": r"\Sone",
    "R³": r"\RealR^3",
    "R": r"\RealR",
    "ⁿ": "^n",
    "²": "^2",
    "³": "^3",
}
_TOKEN = re.compile("|".join(re.escape(k) for k in sorted(_TEX, key=len, reverse=True)))


def manifold_tex(text):
    """The library's configuration string (``R³ × SO(3)``) in the paper's macros."""
    out = _TOKEN.sub(lambda m: _TEX[m.group(0)], text)
    return "$" + re.sub(r" +", " ", out).strip() + "$"


def main():
    rows = []
    for title, name in entries():
        s = S.get(name)
        eom, t = timed(lambda: compute_eom(s.L, s.dW, s.variables))
        terms = sum(len(addends(eq)) for _, eq in eom.values())
        err = check_eom(s.L, s.dW, s.variables, eom, inputs=s.inputs)
        rows.append(
            f"{title} & {manifold_tex(s.configuration)} & {s.dof} & {secs(t)} & "
            f"{terms} & {sci(err)} \\\\"
        )
        print(f"{name:<36} {t:7.3f} s {terms:5d} terms  check_eom {err:.1e}")
    rows.append("\\midrule")
    for name, title, config, dof in INDEXED:
        model = I.MODELS[name]()
        eom, t = timed(model.eom)
        worst = 0.0
        for n in SIZES:
            L, dW, X, inputs = model.unroll(n)
            worst = max(worst, check_eom(L, dW, X, eom.instantiate(n), inputs=inputs))
        rows.append(f"{title} & {config} & {dof} & {secs(t)} & -- & {sci(worst)} \\\\")
        print(f"{name:<36} {t:7.3f} s  check_eom ≤ {worst:.1e} (n = 1..6)")
    (HERE / "verification_rows.tex").write_text(
        "% written by sim/catalog_numbers.py\n" + "\n".join(rows) + "\n"
    )


if __name__ == "__main__":
    main()
