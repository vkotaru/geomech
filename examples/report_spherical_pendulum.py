"""Spherical pendulum on S² (Lee et al. Eq. 5.14): derive, verify, simulate, report.

An S² point has no attitude, so the pendulum is drawn as a rod from the pivot to the bob
plus a point mass, with the path of the bob traced.

    python examples/report_spherical_pendulum.py [out.html]

Writes ``spherical_pendulum.html`` in the current directory by default.
"""

import inspect
import sys
import time

import numpy as np

from geomech import (
    S2,
    Dot,
    Scalar,
    SystemVariables,
    Vector,
    compute_eom,
    getScalars,
    to_standard_form,
)
from geomech.codegen import StateLayout, collect_params, load_module, module_source
from geomech.core.base.expressions import ZeroVector
from geomech.simulation import simulate
from geomech.viz import PointMass, Report, Rod, Scene, Trace


def model():
    """L = ½ m l² ω·ω − m g l q·e₃, δW = ξ·f."""
    half = Scalar("0.5", value=0.5, attr=["Constant"])
    m, g, l = getScalars("m g l", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    q, f = S2("q"), Vector("f")
    omega, xi = q.get_tangent_vector(), q.get_variation_vector()
    L = half * m * l * l * Dot(omega, omega) - m * g * l * Dot(q, e3)
    dW = Dot(xi, f)
    return L, dW, SystemVariables(vectors=[q]), [f], dict(q=q, l=l)


def main(out="spherical_pendulum.html"):
    L, dW, variables, inputs, s = model()
    t0 = time.perf_counter()
    eom = compute_eom(L, dW, variables)
    derive_time = time.perf_counter() - t0
    sf = to_standard_form(eom, variables, inputs)
    layout = StateLayout(variables, inputs, collect_params([L, dW], variables, inputs))
    plant = load_module(module_source(sf, layout, lagrangian=L), "spherical_pendulum")

    p = {"e3": np.array([0.0, 0.0, 1.0]), "g": 9.81, "l": 1.0, "m": 1.0}
    q0 = np.array([np.sin(1.0), 0.0, -np.cos(1.0)])  # 1 rad from hanging
    w0 = np.array([0.3, 2.0, -0.2])
    w0 -= (w0 @ q0) * q0  # tangent at q0
    traj = simulate(plant, layout.pack(q=q0, omega_q=w0), 10.0, p)

    bob = s["q"] * s["l"]
    scene = Scene(layout)  # drawn with traj.p, the simulated parameters
    scene.add(Rod(ZeroVector, bob, radius=0.012))
    scene.add(PointMass(bob, radius=0.08))
    scene.add(Trace(bob))

    report = Report(
        "Spherical pendulum",
        system=(L, dW, variables, inputs),
        eom=eom,
        standard_form=sf,
        scene=scene,
        trajectory=traj,
        derive_time=derive_time,
        units={"m": "kg", "g": "m/s²", "l": "m"},
        source=inspect.getsource(model),
        description="A point mass on a massless rod, on S² (Lee et al. Eq. 5.14).",
    )
    path = report.write(out)
    print(f"check_eom error {report.check_error:.2e}, energy drift {report.energy_drift():.2e}")
    print(f"wrote {path} ({path.stat().st_size / 1e3:.0f} kB)")
    return report


if __name__ == "__main__":
    main(*sys.argv[1:2])
