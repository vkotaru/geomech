"""Stabilise an inverted rigid pendulum on SO(3) with LQR designed from the linearization.

The whole chain, from the Lagrangian to a closed-loop simulation:

1. derive the equations with ``compute_eom`` (Lee et al. Eq. 6.8, inertia about the pivot),
2. split them into ``M a + f + G u = 0`` (``to_standard_form``) and linearize them
   (``linearize``, ``M δẋ = A δx + B δu``),
3. generate NumPy code for both (``geomech.codegen``),
4. design an LQR gain about the inverted equilibrium (``geomech.simulation.lqr``), whose
   attitude error is ``η = vee(log(R_eqᵀ R))``,
5. simulate the nonlinear closed loop with the Lie-group integrator (RKMK4).

The pendulum starts 20° from upright about a random axis and is at rest. Run with
``python examples/inverted_rigid_pendulum_lqr.py`` (needs scipy: ``pip install
'geomech[sim]'``).
"""

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
    linearize,
    to_standard_form,
)
from geomech.codegen import StateLayout, collect_params, load_module, module_source
from geomech.numeric import expm_so3
from geomech.simulation import lqr, simulate

# --- 1. The model ----------------------------------------------------------------------
HALF = Scalar("0.5", value=0.5, attr=["Constant"])
J = Matrix("J", attr=["Constant", "SymmetricMatrix"])  # inertia about the pivot
rho = Vector("\\rho", attr=["Constant"])  # pivot to centre of mass, body frame
m, g = getScalars("m g", attr=["Constant"])
e3 = Vector("e3", attr=["Constant"])
M = Vector("M")  # control torque, body frame
R = SO3("R")
Om, eta = R.get_tangent_vector(), R.get_variation_vector()

L = HALF * Dot(Om, J * Om) - m * g * Dot(R * rho, e3)
dW = Dot(eta, M)
variables = SystemVariables(matrices=[R])

# --- 2. Equations, standard form, linearization ----------------------------------------
eom = compute_eom(L, dW, variables)
sf = to_standard_form(eom, variables, [M])
lin = linearize(eom, variables, [M])

# --- 3. Generated code -----------------------------------------------------------------
# one layout, built with the variables and inputs of to_standard_form and linearize, serves
# both modules; collect_params lists the constants (here of L and dW)
params = collect_params([L, dW], variables, [M])
layout = StateLayout(variables, [M], params)
plant = load_module(module_source(sf, layout, lagrangian=L), "pendulum")
plant_lin = load_module(module_source(lin, layout), "pendulum_lin")

p = {
    "J": np.diag([0.30, 0.30, 0.02]),
    "\\rho": np.array([0.0, 0.0, -0.5]),  # hanging down at R = I
    "e3": np.array([0.0, 0.0, 1.0]),
    "g": 9.81,
    "m": 1.0,
}

# --- 4. LQR about the inverted equilibrium ---------------------------------------------
R_eq = np.diag([1.0, -1.0, -1.0])  # a half turn about e1: R_eq ρ points up
x_eq = layout.pack(R=R_eq)  # Omega_R defaults to zero
# module=plant checks that (x_eq, u_eq) is an equilibrium of the nonlinear equations
ctrl = lqr(plant_lin, x_eq, np.zeros(3), p, Q=np.r_[[100.0] * 3, [10.0] * 3], R=1.0, module=plant)
print("reduced state :", ctrl.labels)
print("gain K        :\n", np.array2string(ctrl.K, precision=3, suppress_small=True))
eig = ctrl.closed_loop_eigenvalues
print("closed loop   :", np.array2string(np.sort(eig.real), precision=3))

# --- 5. Closed-loop simulation ---------------------------------------------------------
rng = np.random.default_rng(0)
axis = rng.normal(size=3)
axis /= np.linalg.norm(axis)
R0 = R_eq @ expm_so3(np.deg2rad(20.0) * axis)
x0 = layout.pack(R=R0, Omega_R=np.zeros(3))

open_loop = simulate(plant, x0, 1.0, p)
traj = simulate(plant, x0, 5.0, p, u=ctrl)

eta_norm = np.array([np.linalg.norm(ctrl.error(x)[:3]) for x in traj.state])
print(f"\nopen loop: |eta| goes from {np.rad2deg(eta_norm[0]):.1f} deg to", end=" ")
print(f"{np.rad2deg(np.linalg.norm(ctrl.error(open_loop.state[-1])[:3])):.1f} deg in 1 s")
print("closed loop:")
for t in (0.0, 1.0, 2.0, 3.0, 4.0, 5.0):
    k = int(np.argmin(np.abs(traj.t - t)))
    print(f"  t = {traj.t[k]:.1f} s   |eta| = {eta_norm[k]:.3e} rad")
print(f"max |RᵀR − I| over the run: {traj.max_drift():.1e}")
assert eta_norm[-1] < 1e-3, eta_norm[-1]
