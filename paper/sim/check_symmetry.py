"""Check that the standard-form mass matrix of every library system is symmetric.

    python sim/check_symmetry.py

For each system, the generated model is simulated unforced for one second from its initial
state, and the relative asymmetry max|M - M^T| / max(1, max|M|) is evaluated at about 20 of
the stored states (at x0 alone if the simulation fails). The script fails above 1e-14.
"""

from __future__ import annotations

import warnings

import numpy as np

from geomech import systems, to_standard_form
from geomech.codegen.module import generate_module, load_module
from geomech.simulation import simulate

TOL = 1e-14


def asymmetry(name: str) -> float:
    system = systems.get(name)
    sf = to_standard_form(system.derive(), system.variables, list(system.inputs))
    src = generate_module(sf, system.layout(), lagrangian=system.L)
    mod = load_module(src, name=f"sym_{name}")
    try:
        states = simulate(mod, system.x0, 1.0, system.params).state
    except Exception as err:  # report and fall back to x0
        print(f"{name}: simulation failed ({err}); checking x0 only")
        states = np.atleast_2d(system.x0)
    worst = 0.0
    for x in states[:: max(1, len(states) // 20)]:
        M = mod.mass(x, system.params)
        worst = max(worst, float(np.abs(M - M.T).max() / max(1.0, np.abs(M).max())))
    return worst


def main() -> None:
    warnings.filterwarnings("ignore")
    worst = 0.0
    for name in systems.names():
        w = asymmetry(name)
        worst = max(worst, w)
        print(f"{name:45s} {w:.2e}")
    print(f"worst {worst:.2e}")
    if worst > TOL:
        raise SystemExit(f"a mass matrix is asymmetric by {worst:.1e}")


if __name__ == "__main__":
    main()
