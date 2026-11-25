"""Closed-loop simulation for the figure in Section IV of the paper.

Quadrotor with a cable-suspended load (``geomech.systems.quadrotor_suspended_load``,
configuration space R^3 x SO(3) x S^2). The script derives the equations of motion,
linearizes them about hover with the variation-based linearization, designs an LQR gain
from that linearization, and simulates the nonlinear closed loop with the Lie-group
integrator (RKMK4) of ``geomech.simulation``.

Run from this directory with geomech installed (``pip install -e "<geomech>[sim]"``):

    python quadrotor_load_lqr.py          # writes quadrotor_load_lqr.dat

The data file has one row per 20th stored state (every 20 ms): time, |x_L - x_L,eq|,
the cable angle from hanging (deg), the attitude error |log(R_eq^T R)| (deg), the thrust
f and |M|. The script prints the largest constraint drift over the run.
"""

import numpy as np

from geomech import linearize, to_standard_form
from geomech.codegen import load_module, module_source
from geomech.numeric import expm_so3, log_so3
from geomech.simulation import lqr, simulate
from geomech.systems import get

sys = get("quadrotor_suspended_load")
eom = sys.derive()
sf = to_standard_form(eom, sys.variables, sys.inputs)
lin = linearize(eom, sys.variables, sys.inputs)
layout = sys.layout()
plant = load_module(module_source(sf, layout, lagrangian=sys.L), "plant")
plant_lin = load_module(module_source(lin, layout), "plant_lin")
p = dict(sys.params)

# hover: load below the quadrotor, cable along -e3, level attitude
mQ, mL, g = p["m_Q"], p["m_L"], p["g"]
x_eq = layout.pack(x_L=np.zeros(3), q=np.array([0.0, 0.0, -1.0]), R=np.eye(3))
u_eq = np.r_[(mQ + mL) * g, np.zeros(3)]
n_red = 3 + 2 + 3  # reduced configuration perturbation: x_L, xi (2), eta
Q = np.r_[[10.0] * 3, [10.0] * 2, [10.0] * 3, [1.0] * n_red]
ctrl = lqr(plant_lin, x_eq, u_eq, p, Q=Q, R=np.r_[0.1, [10.0] * 3], module=plant)

# initial condition: load displaced by 0.5 m, cable 30 deg from hanging, quadrotor tilted 20 deg
q0 = np.array([np.sin(np.deg2rad(30.0)), 0.0, -np.cos(np.deg2rad(30.0))])
R0 = expm_so3(np.deg2rad(20.0) * np.array([0.0, 1.0, 1.0]) / np.sqrt(2.0))
x0 = layout.pack(x_L=np.array([0.5, -0.3, 0.2]), q=q0, R=R0)

T = 8.0
traj = simulate(plant, x0, T, p, u=ctrl, method="rkmk4", dt=1e-3)

rows = []
stride = 20
for k in range(0, len(traj.t), stride):
    s = layout.unpack(traj.state[k])
    xL = np.asarray(s["x_L"], dtype=float).ravel()
    q = np.asarray(s["q"], dtype=float).ravel()
    R = np.asarray(s["R"], dtype=float).reshape(3, 3)
    ang_q = np.rad2deg(np.arccos(np.clip(-q[2], -1.0, 1.0)))
    ang_R = np.rad2deg(np.linalg.norm(log_so3(R)))
    u = traj.u[k]
    rows.append([traj.t[k], np.linalg.norm(xL), ang_q, ang_R, u[0], np.linalg.norm(u[1:4])])
rows = np.array(rows)
np.savetxt(
    "quadrotor_load_lqr.dat",
    rows,
    header="t xL_err q_angle_deg R_angle_deg f M_norm",
    comments="",
    fmt="%.6e",
)
print(
    "closed-loop eigenvalues (max real part):",
    np.max(ctrl.closed_loop_eigenvalues.real),
)
print("final |x_L|, cable angle, attitude angle:", rows[-1, 1:4])
print("max constraint drift over the run:", traj.max_drift())
