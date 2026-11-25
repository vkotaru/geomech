"""Derivation time and size against n for three parametric entries (the scaling table).

    python sim/scaling_numbers.py

Prints one LaTeX row per n: the median wall time of compute_eom over three runs and the
number of additive terms of the derived equations, for the n-link spherical chain, the
n-link rigid chain and n quadrotors carrying a point mass.
"""

import statistics
import time

import geomech.systems as S
from geomech import compute_eom
from geomech.core.transformations.utils import addends

CASES = [
    ("n_link_spherical_chain", "n_links"),
    ("n_link_rigid_chain", "n_links"),
    ("multi_quadrotor_point_mass", "n_quads"),
]


def measure(name, key, n):
    s = S.get(name, **{key: n})
    times = []
    for _ in range(3):
        t = time.perf_counter()
        eom = compute_eom(s.L, s.dW, s.variables)
        times.append(time.perf_counter() - t)
    return statistics.median(times), sum(len(addends(eq)) for _, eq in eom.values())


def fmt(t):
    return f"{t:.3f}" if t < 1 else f"{t:.2f}"


for n in (1, 2, 3, 4, 5, 10):
    cells = []
    for name, key in CASES:
        t, terms = measure(name, key, n)
        cells += [fmt(t), str(terms)]
    print(f"{n:<2} & " + " & ".join(cells) + " \\\\")
