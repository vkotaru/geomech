"""Quadrotor with a cable-suspended load on R³ × S² × SO(3): derive, verify, stabilise
with LQR about hover, simulate and write a report.

The quadrotor is a box with two arms and four rotor cylinders at the arm ends (one rigid
body, drawn as several shapes with the same attitude), the cable a rod and the load a
sphere. Needs scipy for ``lqr`` (``pip install "geomech[sim]"``).

    python examples/report_quadrotor_load.py [out.html]

Writes ``quadrotor_load.html`` in the current directory by default.
"""

import inspect
import sys
import time

import numpy as np

from geomech import (
    S2,
    SO3,
    Cross,
    Dot,
    Matrix,
    Scalar,
    SystemVariables,
    Vector,
    compute_eom,
    getScalars,
    linearize,
    to_standard_form,
)
from geomech.codegen import StateLayout, collect_params, load_module, module_source
from geomech.numeric import expm_so3
from geomech.simulation import lqr, simulate
from geomech.viz import Box, Cylinder, PointMass, Report, RigidBody, Rod, Scene, Trace


def model():
    """Load at x_L, cable direction q (x_Q = x_L − l q), quadrotor attitude R.

    L = ½ m_Q ẋ_Q·ẋ_Q + ½ m_L ẋ_L·ẋ_L + ½ Ω·JΩ − m_Q g e₃·x_Q − m_L g e₃·x_L,
    δW = f (R e₃)·(δx_L − l ξ × q) + M·η.
    """
    half = Scalar("0.5", value=0.5, attr=["Constant"])
    mQ, mL, g, l = getScalars("m_Q m_L g l", attr=["Constant"])
    J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
    e3 = Vector("e3", attr=["Constant"])
    xL, q, R = Vector("x_L"), S2("q"), SO3("R")
    xi = q.get_variation_vector()
    Om, eta = R.get_tangent_vector(), R.get_variation_vector()
    f, M = Scalar("f"), Vector("M")
    xQ = xL - q * l
    vQ, vL = xQ.t_diff(), xL.t_diff()
    L = (
        half * mQ * Dot(vQ, vQ)
        + half * mL * Dot(vL, vL)
        + half * Dot(Om, J * Om)
        - mQ * g * Dot(e3, xQ)
        - mL * g * Dot(e3, xL)
    )
    dW = Dot(xL.delta() - Cross(xi, q) * l, (R * e3) * f) + Dot(eta, M)
    variables = SystemVariables(vectors=[xL, q], matrices=[R])
    return L, dW, variables, [f, M], dict(xL=xL, xQ=xQ, R=R)


def main(out="quadrotor_load.html"):
    L, dW, variables, inputs, s = model()
    t0 = time.perf_counter()
    eom = compute_eom(L, dW, variables)
    derive_time = time.perf_counter() - t0
    sf = to_standard_form(eom, variables, inputs)
    layout = StateLayout(variables, inputs, collect_params([L, dW], variables, inputs))
    plant = load_module(module_source(sf, layout, lagrangian=L), "quadrotor_load")
    plant_lin = load_module(module_source(linearize(eom, variables, inputs), layout), "ql_lin")

    p = {
        "J": np.diag([0.02, 0.025, 0.04]),
        "e3": np.array([0.0, 0.0, 1.0]),
        "g": 9.81,
        "l": 1.0,
        "m_L": 0.3,
        "m_Q": 1.5,
    }
    # hover: load at the origin, cable hanging (q from the quadrotor to the load is −e3)
    x_eq = layout.pack(x_L=np.zeros(3), q=[0, 0, -1], R=np.eye(3))
    u_eq = np.r_[(p["m_Q"] + p["m_L"]) * p["g"], 0, 0, 0]
    ctrl = lqr(plant_lin, x_eq, u_eq, p, Q=1.0, R=1.0, module=plant)

    q0 = np.array([np.sin(0.5), 0.0, -np.cos(0.5)])  # load swung 0.5 rad
    w0 = np.array([0.0, 1.0, 0.0])
    w0 -= (w0 @ q0) * q0
    x0 = layout.pack(
        x_L=[0.6, -0.4, 0.3], q=q0, R=expm_so3([0.2, -0.3, 0.8]), omega_q=w0, Omega_R=[0, 0, 1.5]
    )
    traj = simulate(plant, x0, 10.0, p, u=ctrl)

    R, xQ, xL = s["R"], s["xQ"], s["xL"]
    arm = 0.25
    scene = Scene(layout, params=p)
    scene.add(RigidBody(R, origin=xQ, shape=Box(0.16, 0.12, 0.06), color="link1", name="body"))
    scene.add(RigidBody(R, xQ, Box(2 * arm, 0.025, 0.02), axes=False, color="cable", name="arm 1"))
    scene.add(RigidBody(R, xQ, Box(0.025, 2 * arm, 0.02), axes=False, color="cable", name="arm 2"))
    for k, off in enumerate(
        [(arm, 0, 0.025), (0, arm, 0.025), (-arm, 0, 0.025), (0, -arm, 0.025)]
    ):
        color = "link2" if k == 0 else "rotor"  # the rotor on +e1 is marked
        scene.add(
            RigidBody(
                R,
                xQ,
                Cylinder(0.09, 0.012, "z"),
                offset=off,
                axes=False,
                color=color,
                name=f"rotor {k + 1}",
            )  # fmt: skip
        )
    scene.add(Rod(xQ, xL, radius=0.006))
    scene.add(PointMass(xL, radius=0.07))
    scene.add(Trace(xL))
    scene.add(Trace(xQ, color="link1"))

    report = Report(
        "Quadrotor with a suspended load",
        system=(L, dW, variables, inputs),
        eom=eom,
        standard_form=sf,
        scene=scene,
        trajectory=traj,
        params=p,
        derive_time=derive_time,
        units={"m_Q": "kg", "m_L": "kg", "g": "m/s²", "l": "m", "J": "kg m²"},
        source=inspect.getsource(model),
        description="LQR about hover from a swinging, yawing start (forced: the energy"
        " changes by design).",
    )
    path = report.write(out)
    print(f"check_eom error {report.check_error:.2e}")
    print(f"wrote {path} ({path.stat().st_size / 1e3:.0f} kB)")
    return report


if __name__ == "__main__":
    main(*sys.argv[1:2])
