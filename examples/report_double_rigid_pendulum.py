"""Double rigid pendulum on SO(3) × SO(3): derive, verify, simulate and write a report.

Each link is a slender box with a non-square cross-section and an anisotropic inertia
(the inertia of that box about its centre of mass). Both links start with an angular
velocity that has a component along the link's own axis (body e₃), so they spin about
their long axes while they swing: the tinted +e₁ face and the body axes in the report
make that spin visible.

    python examples/report_double_rigid_pendulum.py [out.html]

Writes ``double_rigid_pendulum.html`` in the current directory by default.
"""

import inspect
import sys
import time

import numpy as np

from geomech import (
    SO3,
    Dot,
    Matrix,
    Scalar,
    SystemVariables,
    Vector,
    compute_eom,
    getScalars,
    to_standard_form,
)
from geomech.codegen import StateLayout, collect_params, load_module, module_source
from geomech.core.base.expressions import ZeroVector
from geomech.numeric import expm_so3
from geomech.simulation import simulate
from geomech.viz import Box, Report, RigidBody, Scene, Trace


def model():
    """Two links hinged in series, inertias about their centres of mass."""
    half = Scalar("0.5", value=0.5, attr=["Constant"])
    J1 = Matrix("J1", attr=["Constant", "SymmetricMatrix"])
    J2 = Matrix("J2", attr=["Constant", "SymmetricMatrix"])
    rho1 = Vector("\\rho_1", attr=["Constant"])  # pivot -> centre of mass of link 1
    rho2 = Vector("\\rho_2", attr=["Constant"])  # joint -> centre of mass of link 2
    l1 = Vector("l_1", attr=["Constant"])  # pivot -> joint, in link 1's frame
    m1, m2, g = getScalars("m1 m2 g", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    M1, M2 = Vector("M1"), Vector("M2")
    R1, R2 = SO3("R1"), SO3("R2")
    Om1, Om2 = R1.get_tangent_vector(), R2.get_tangent_vector()
    eta1, eta2 = R1.get_variation_vector(), R2.get_variation_vector()
    x1 = R1 * rho1
    x2 = R1 * l1 + R2 * rho2
    v1, v2 = x1.t_diff(), x2.t_diff()
    T = (
        half * Dot(Om1, J1 * Om1)
        + half * Dot(Om2, J2 * Om2)
        + half * m1 * Dot(v1, v1)
        + half * m2 * Dot(v2, v2)
    )
    V = m1 * g * Dot(x1, e3) + m2 * g * Dot(x2, e3)
    dW = Dot(eta1, M1) + Dot(eta2, M2)
    variables = SystemVariables(matrices=[R1, R2])
    return T - V, dW, variables, [M1, M2], dict(R1=R1, R2=R2, l1=l1, rho2=rho2)


def box_inertia(m, lx, ly, lz):
    return m / 12 * np.diag([ly**2 + lz**2, lx**2 + lz**2, lx**2 + ly**2])


def main(out="double_rigid_pendulum.html"):
    L, dW, variables, inputs, s = model()
    t0 = time.perf_counter()
    eom = compute_eom(L, dW, variables)
    derive_time = time.perf_counter() - t0
    sf = to_standard_form(eom, variables, inputs)
    layout = StateLayout(variables, inputs, collect_params([L, dW], variables, inputs))
    plant = load_module(module_source(sf, layout, lagrangian=L), "double_rigid_pendulum")

    link1, link2 = (0.12, 0.05, 1.0), (0.12, 0.05, 0.8)
    p = {
        "J1": box_inertia(1.0, *link1),
        "J2": box_inertia(0.8, *link2),
        "\\rho_1": np.array([0.0, 0.0, -0.5]),
        "\\rho_2": np.array([0.0, 0.0, -0.4]),
        "l_1": np.array([0.0, 0.0, -1.0]),
        "m1": 1.0,
        "m2": 0.8,
        "g": 9.81,
        "e3": np.array([0.0, 0.0, 1.0]),
    }
    # each Ω has a component along its link (body e3): a spin about the link's own axis
    x0 = layout.pack(
        R1=expm_so3([0.7, 0.3, 0.0]),
        R2=expm_so3([1.1, -0.4, 0.3]),
        Omega_R1=[0.4, -0.6, 4.0],
        Omega_R2=[-0.8, 0.5, -5.0],
    )
    traj = simulate(plant, x0, 10.0, p)

    scene = Scene(layout)  # drawn with traj.p, the simulated parameters
    scene.add(RigidBody(s["R1"], origin=ZeroVector, shape=Box(*link1), offset=[0, 0, -0.5]))
    scene.add(RigidBody(s["R2"], origin=s["R1"] * s["l1"], shape=Box(*link2), offset=[0, 0, -0.4]))
    scene.add(Trace(s["R1"] * s["l1"] + s["R2"] * (s["rho2"] * 2)))  # the tip of link 2

    report = Report(
        "Double rigid pendulum",
        system=(L, dW, variables, inputs),
        eom=eom,
        standard_form=sf,
        scene=scene,
        trajectory=traj,
        derive_time=derive_time,
        units={
            "m1": "kg",
            "m2": "kg",
            "g": "m/s²",
            "J1": "kg m²",
            "J2": "kg m²",
            "l_1": "m",
            "\\rho_1": "m",
            "\\rho_2": "m",
        },  # fmt: skip
        source=inspect.getsource(model),
        description="Two rigid links on SO(3) × SO(3), each spinning about its own axis.",
    )
    path = report.write(out)
    print(f"check_eom error {report.check_error:.2e}, energy drift {report.energy_drift():.2e}")
    print(f"wrote {path} ({path.stat().st_size / 1e3:.0f} kB)")
    return report


if __name__ == "__main__":
    main(*sys.argv[1:2])
