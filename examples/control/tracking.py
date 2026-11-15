"""Stretch examples: time-varying LQR tracking along a trajectory, designed from the
variation-based linearization evaluated along the reference.

    python examples/control/tracking.py quadrotor_figure_eight
    python examples/control/tracking.py suspended_load_circle
    python examples/control/tracking.py --all

**Reduced time-varying linearization.** The linearization module (``linearize`` →
``generate_module``) is evaluated at every reference point ``(x_r(t), a_r(t), u_r(t))``
with ``geomech.simulation.reduced_linearization`` (``basis``: the projection of e₁ on the
tangent plane, completed by ``q × b₁``, smooth while q stays away from ±e₁). Its
equilibrium reduction is exact only for a constant basis and ``ω = 0``; along a moving
reference three terms are added for each S² point q with reference rate ω:

- the column map ``T_c``: the velocity perturbation is not tangent,
  ``q·δω = −(q × ω)·ξ`` (from ``q·ω = 0``), so ``δω = B z_ω − q (q × ω)ᵀ B z_ξ``;
- the time-varying basis: ``δx = T_c(t) z`` gives ``δẋ = T_c ż + Ṫ_c z``;
- the row projection: the S² equation is ``Bᵀ E = 0`` and the reference satisfies it with
  a normal residual ``E_r = q (q·E_r)``; with ``δB ᵀq = −Bᵀδq`` and ``δq = −q̂ ξ`` the
  linearized row gains ``(q·E_r) Bᵀ q̂ B z_ξ``.

With ``Tᵀ`` the row projection (``T`` the equilibrium reduction) and ``M``, ``A``, ``B`` the
generated matrices (``δE = −M δv̇ + A δx + B δu``):
``Tᵀ M T_c ż = (Tᵀ A T_c − Tᵀ M Ṫ_c + C) z + Tᵀ B δu``. ``validate_reduction`` compares the
result with central differences of the nonlinear error dynamics along the reference.

The Riccati differential equation ``−Ṗ = AᵀP + PA − P B R⁻¹ BᵀP + Q``, ``P(t_f) = Q``, is
integrated backward with ``scipy.integrate.solve_ivp`` and ``u = u_r(t) − R⁻¹B(t)ᵀP(t) z``
is simulated with the nonlinear generated model (RKMK4, dt = 1 ms). This lives in the
example, not in the library.

**References.**

- ``quadrotor_figure_eight``: the quadrotor follows a 3-D figure-eight
  ``x_d = (1.5 sin wt, 0.75 sin 2wt, 0.5 sin² wt)``, ``w = 2π/8``. The reference is built
  by differential flatness with analytic derivatives (``flat_reference``): thrust vector
  ``F = m (ẍ_d + g e₃)``, ``b₃ = F/|F|``, ``b₂ = (b₃ × e₁)/|b₃ × e₁|``, ``b₁ = b₂ × b₃``
  (zero yaw), ``R_d = [b₁ b₂ b₃]``; ``Ṙ_d`` and ``R̈_d`` from the jerk and snap through the
  derivatives of the normalisations, ``Ω_d = vee(R_dᵀṘ_d)``,
  ``Ω̇_d = vee(skew(R_dᵀR̈_d))``, thrust ``f = F·b₃`` and moment ``M = JΩ̇_d + Ω_d × JΩ_d``.
  It is checked independently: ``‖Ṙ_d − R_dΩ̂_d‖`` with ``Ṙ_d`` from finite differences of
  ``R_d`` and ``∫Ω̇_d dt`` against ``ΔΩ_d``.
- ``suspended_load_circle``: the quadrotor with a suspended load. The reference is a
  trajectory of the nonlinear model itself (not a flatness construction): the hover LQR
  of ``run.py`` driven with a set point moving on a 1 m circle (period 5 s), recorded with
  its inputs, so it is feasible by construction.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from control_helper import RESULTS, build, equilibrium_input, rotate  # noqa: E402

from geomech.numeric import expm_so3, hat, log_so3  # noqa: E402
from geomech.simulation import (  # noqa: E402
    accelerations,
    lqr,
    model_of,
    reduced_linearization,
    simulate,
)

DT = 1e-3
E1, E3 = np.eye(3)[0], np.eye(3)[2]


def vee(A) -> np.ndarray:
    return np.array([A[2, 1], A[0, 2], A[1, 0]])


def basis(q) -> np.ndarray:
    """A smooth tangent basis at *q* (valid away from q = ±e₁)."""
    q = np.asarray(q, dtype=float).reshape(3)
    b1 = E1 - (E1 @ q) * q
    b1 /= np.linalg.norm(b1)
    return np.column_stack([b1, np.cross(q, b1)])


class Reference:
    """A reference sampled every ``dt / 2`` (so that every RKMK4 stage time is a sample)."""

    def __init__(self, t, X, U, A_acc):
        self.t, self.X, self.U, self.A_acc = t, X, U, A_acc
        self.h = t[1] - t[0]

    def index(self, t) -> int:
        return int(np.clip(round((t - self.t[0]) / self.h), 0, len(self.t) - 1))


# ---------------------------------------------------------------------------
# Time-varying reduced linearization
# ---------------------------------------------------------------------------


def _s2_blocks(model, lin):
    """Per S² variable: (var, full ξ slice, full δω slice)."""
    states, nv = lin["states"], len(model.variables)
    out = []
    for i, var in enumerate(model.variables):
        if var.manifold == "S2":
            out.append((var, slice(*states[i]["slice"]), slice(*states[nv + i]["slice"])))
    return out


def column_map(model, lin, x, T, coupling=True):
    """``T_c``: the equilibrium reduction ``T`` plus, for each S² point, the normal part of
    the velocity perturbation ``−q (q × ω)ᵀ B`` in the ξ columns of the δω rows."""
    Tc = T.copy()
    if coupling:
        for var, xi, dom in _s2_blocks(model, lin):
            q, om = x[var.cfg], x[var.vel]
            cols = np.flatnonzero(np.abs(T[xi]).sum(axis=0))
            Tc[dom, cols] += -np.outer(q, np.cross(q, om)) @ T[xi][:, cols]
    return Tc


def curvature(model, lin, plant, x, a, u, p, T):
    """``C``: the S² rows' term ``(q·E_r) Bᵀ q̂ B`` in the ξ columns (``E_r`` the equation
    residual at the reference, normal to the sphere)."""
    n = T.shape[1]
    C = np.zeros((n, n))
    r = plant.mass(x, p) @ a + plant.bias(x, p) + plant.input_matrix(x, p) @ u
    for var, xi, dom in _s2_blocks(model, lin):
        q = x[var.cfg]
        nq = float(q @ r[var.acc])
        cx = np.flatnonzero(np.abs(T[xi]).sum(axis=0))
        cw = np.flatnonzero(np.abs(T[dom]).sum(axis=0))
        B = T[xi][:, cx]
        C[np.ix_(cw, cx)] += nq * B.T @ hat(q) @ B
    return C


def reduced_along(plant, plant_lin, ref, p, k_list=None, terms=("Tdot", "coupling", "curv")):
    """``(A_k, B_k)`` of the reduced error dynamics at the reference samples *k_list*."""
    model = model_of(plant_lin)
    lin = plant_lin.LINEARIZATION
    ks = range(len(ref.t)) if k_list is None else k_list

    def parts(k):
        red = reduced_linearization(plant_lin, ref.X[k], ref.U[k], p, a=ref.A_acc[k], basis=basis)
        Tc = column_map(model, lin, ref.X[k], red.T, coupling="coupling" in terms)
        return red, Tc

    out_A, out_B = [], []
    for k in ks:
        red, Tc = parts(k)
        km, kp = max(k - 1, 0), min(k + 1, len(ref.t) - 1)
        if "Tdot" in terms:
            Tm, Tp = parts(km)[1], parts(kp)[1]
            Td = (Tp - Tm) / ((kp - km) * ref.h)
        else:
            Td = np.zeros_like(Tc)
        T = red.T
        Mr = T.T @ red.M_full @ Tc
        rhs = T.T @ red.A_full @ Tc - T.T @ red.M_full @ Td
        if "curv" in terms:
            rhs = rhs + curvature(model, lin, plant, ref.X[k], ref.A_acc[k], ref.U[k], p, T)
        out_A.append(np.linalg.solve(Mr, rhs))
        out_B.append(np.linalg.solve(Mr, T.T @ red.B_full))
    return np.array(out_A), np.array(out_B)


def error(model, x, xr) -> np.ndarray:
    """Reduced error of *x* about the reference state *xr* (the coordinates of ``lqr``,
    with the basis ``basis(q_ref)``)."""
    cfg, vel = [], []
    for var in model.variables:
        if var.manifold == "S2":
            B = basis(xr[var.cfg])
            cfg.append(B.T @ np.cross(xr[var.cfg], x[var.cfg]))
            vel.append(B.T @ (x[var.vel] - xr[var.vel]))
        elif var.manifold == "SO3":
            cfg.append(log_so3(xr[var.cfg].reshape(3, 3).T @ x[var.cfg].reshape(3, 3)))
            vel.append(x[var.vel] - xr[var.vel])
        else:
            cfg.append(x[var.cfg] - xr[var.cfg])
            vel.append(x[var.vel] - xr[var.vel])
    return np.concatenate(cfg + vel)


def validate_reduction(plant, plant_lin, ref, p, samples):
    """Max relative difference between the reduced ``A(t)`` (with the listed correction
    terms) and central differences of the nonlinear error dynamics
    ``ż = d/dt error(x(t), x_r(t))`` along the reference, at the reference samples."""
    model = model_of(plant)

    def step(x, u, h):
        a = accelerations(plant, x, u, p)
        y = x.copy()
        for v in model.variables:
            if v.manifold == "S2":
                y[v.cfg] = expm_so3(h * x[v.vel]) @ x[v.cfg]
            elif v.manifold == "SO3":
                y[v.cfg] = (x[v.cfg].reshape(3, 3) @ expm_so3(h * x[v.vel])).ravel()
            else:
                y[v.cfg] = x[v.cfg] + h * x[v.vel]
            y[v.vel] = x[v.vel] + h * a[v.acc]
        return y

    def to_state(z, xr):
        x, k, n = xr.copy(), 0, len(z) // 2
        for v in model.variables:
            d = (
                2
                if v.manifold == "S2"
                else (v.cfg.stop - v.cfg.start if v.manifold != "SO3" else 3)
            )
            c, w = z[k : k + d], z[n + k : n + k + d]
            k += d
            if v.manifold == "S2":
                B = basis(xr[v.cfg])
                q = expm_so3(B @ c) @ xr[v.cfg]
                om = xr[v.vel] + B @ w
                om = om - xr[v.cfg] * (om @ q) / (xr[v.cfg] @ q)  # ω·q = 0, Bᵀ(ω − ω_r) = w
                x[v.cfg], x[v.vel] = q, om
            elif v.manifold == "SO3":
                x[v.cfg] = (xr[v.cfg].reshape(3, 3) @ expm_so3(c)).ravel()
                x[v.vel] = xr[v.vel] + w
            else:
                x[v.cfg], x[v.vel] = xr[v.cfg] + c, xr[v.vel] + w
        return x

    def zdot(z, k, h=1e-6):
        xr, ur = ref.X[k], ref.U[k]
        x = to_state(z, xr)
        ep = error(model, step(x, ur, h), step(xr, ur, h))
        em = error(model, step(x, ur, -h), step(xr, ur, -h))
        return (ep - em) / (2 * h)

    variants = {
        "equilibrium reduction only": (),
        "+ basis rate (−TᵀMṪ)": ("Tdot",),
        "+ normal velocity coupling": ("Tdot", "coupling"),
        "+ curvature term (all)": ("Tdot", "coupling", "curv"),
    }
    out = {name: 0.0 for name in variants}
    for k in samples:
        n = len(error(model, ref.X[k], ref.X[k]))
        eps = 1e-5
        A_fd = np.column_stack(
            [(zdot(eps * e, k) - zdot(-eps * e, k)) / (2 * eps) for e in np.eye(n)]
        )
        scale = np.abs(A_fd).max()
        for name, terms in variants.items():
            A = reduced_along(plant, plant_lin, ref, p, [k], terms)[0][0]
            out[name] = max(out[name], float(np.abs(A - A_fd).max() / scale))
    return out


def tvlqr(plant, plant_lin, ref: Reference, p, Q, R):
    """``K(t_k)`` on the reference grid from the backward Riccati differential equation."""
    from scipy.integrate import solve_ivp

    As, Bs = reduced_along(plant, plant_lin, ref, p)
    n = As.shape[1]
    Qm, Rm = np.diag(Q), np.diag(R)
    Rinv = np.linalg.inv(Rm)

    def interp(arr, t):
        s = (t - ref.t[0]) / ref.h
        k = int(np.clip(np.floor(s), 0, len(ref.t) - 2))
        w = s - k
        return (1 - w) * arr[k] + w * arr[k + 1]

    def rhs(t, y):
        P = y.reshape(n, n)
        A, B = interp(As, t), interp(Bs, t)
        dP = -(A.T @ P + P @ A - P @ B @ Rinv @ B.T @ P + Qm)
        return dP.reshape(-1)

    sol = solve_ivp(
        rhs,
        (ref.t[-1], ref.t[0]),
        Qm.reshape(-1),
        t_eval=ref.t[::-1],
        rtol=1e-8,
        atol=1e-10,
        method="LSODA",
    )
    Ps = sol.y.T[::-1].reshape(-1, n, n)
    return np.array([Rinv @ B.T @ P for B, P in zip(Bs, Ps)])


def track(name, sys_name, ref, x0, Q, R, horizon, units, extra=None, dt=DT, save=True):
    sys_, layout, plant, plant_lin, _ = build(sys_name, {})
    p = dict(sys_.params)
    model = model_of(plant)
    Ks = tvlqr(plant, plant_lin, ref, p, Q, R)

    def ctrl(t, x):
        k = ref.index(t)
        return ref.U[k] - Ks[k] @ error(model, x, ref.X[k])

    traj = simulate(plant, x0, horizon, p, u=ctrl, dt=dt, output_every=max(1, round(1e-2 / dt)))
    pos = layout.configuration[0].slice  # x (quadrotor) or x_L (load)
    att = layout.slot("R").slice
    ks = [ref.index(t) for t in traj.t]
    pos_err = np.linalg.norm(traj.state[:, pos] - ref.X[ks][:, pos], axis=1)
    att_err = np.array(
        [
            np.linalg.norm(log_so3(ref.X[k][att].reshape(3, 3).T @ x[att].reshape(3, 3)))
            for k, x in zip(ks, traj.state)
        ]
    )
    half = traj.t >= horizon / 2
    ref_pos = ref.X[:, pos]
    peaks = {}
    for slot in layout.input_slots:
        dev = np.abs(traj.u[:, slot.slice] - ref.U[ks][:, slot.slice])
        peaks[slot.key] = {"peak_deviation": float(dev.max()), "unit": units[slot.key]}
    metrics = {
        "name": name,
        "system": sys_name,
        "horizon": horizon,
        "dt": dt,
        "reference_radius": float(np.max(np.linalg.norm(ref_pos - ref_pos.mean(axis=0), axis=1))),
        "reference_speed_max": float(
            np.max(np.linalg.norm(np.gradient(ref_pos, ref.h, axis=0), axis=1))
        ),
        "initial_position_error": float(pos_err[0]),
        "initial_attitude_error": float(att_err[0]),
        "rms_position_error_second_half": float(np.sqrt(np.mean(pos_err[half] ** 2))),
        "final_position_error": float(pos_err[-1]),
        "final_attitude_error": float(att_err[-1]),
        "peak_inputs": peaks,
        "constraint_drift": traj.max_drift(),
        **(extra or {}),
    }
    if save:
        RESULTS.mkdir(exist_ok=True)
        (RESULTS / f"tracking_{name}.json").write_text(json.dumps(metrics, indent=1) + "\n")
        np.savez_compressed(
            RESULTS / f"tracking_{name}.npz",
            t=traj.t,
            state=traj.state,
            u=traj.u,
            position_error=pos_err,
            attitude_error=att_err,
        )
    return metrics


# ---------------------------------------------------------------------------
# Quadrotor: figure-eight by differential flatness (analytic derivatives)
# ---------------------------------------------------------------------------

PERIOD_FIG8 = 8.0


def figure_eight(t, k=0, period=PERIOD_FIG8):
    """The k-th time derivative of ``(1.5 sin wt, 0.75 sin 2wt, 0.5 sin² wt)``
    (``0.5 sin² wt = 0.25 − 0.25 cos 2wt``)."""
    w = 2 * np.pi / period
    s = k * np.pi / 2
    z = -0.25 * (2 * w) ** k * np.cos(2 * w * t + s) + (0.25 if k == 0 else 0.0)
    return np.array(
        [1.5 * w**k * np.sin(w * t + s), 0.75 * (2 * w) ** k * np.sin(2 * w * t + s), z]
    )


def _unit(v, vd, vdd):
    """``u = v/|v|`` and its first two time derivatives."""
    n = np.linalg.norm(v)
    u = v / n
    P = np.eye(3) - np.outer(u, u)
    ud = P @ vd / n
    nd = u @ vd
    udd = P @ vdd / n - (np.outer(ud, u) + np.outer(u, ud)) @ vd / n - ud * nd / n
    return u, ud, udd


def flat_reference(t, m, g, J, period=PERIOD_FIG8):
    """``(x, v, a, R, Ω, Ω̇, f, M)`` of the quadrotor on the figure-eight at time *t*."""
    x, v, a, j, s = (figure_eight(t, k, period) for k in range(5))
    F, Fd, Fdd = m * (a + g * E3), m * j, m * s
    b3, b3d, b3dd = _unit(F, Fd, Fdd)
    c, cd, cdd = np.cross(b3, E1), np.cross(b3d, E1), np.cross(b3dd, E1)
    b2, b2d, b2dd = _unit(c, cd, cdd)
    b1 = np.cross(b2, b3)
    b1d = np.cross(b2d, b3) + np.cross(b2, b3d)
    b1dd = np.cross(b2dd, b3) + 2 * np.cross(b2d, b3d) + np.cross(b2, b3dd)
    R = np.column_stack([b1, b2, b3])
    Rd = np.column_stack([b1d, b2d, b3d])
    Rdd = np.column_stack([b1dd, b2dd, b3dd])
    Om = vee(R.T @ Rd)
    X = R.T @ Rdd
    Omd = vee(0.5 * (X - X.T))
    f = F @ b3
    M = J @ Omd + np.cross(Om, J @ Om)
    return x, v, a, R, Om, Omd, f, M


def flatness_checks(m, g, J, period=PERIOD_FIG8, horizon=16.0):
    """Independent kinematic checks of the flat reference: ``max ‖Ṙ_d − R_dΩ̂_d‖`` with
    ``Ṙ_d`` by 5-point differences of ``R_d`` (h = 1e-3), ``max |∫Ω̇_d dt − ΔΩ_d|`` by the
    trapezoidal rule (h = 1e-4), and ``max |ẍ_d − FD(ẋ_d)|``; with the magnitudes."""
    R = lambda t: flat_reference(t, m, g, J, period)[3]  # noqa: E731
    h = 1e-3
    ts = np.linspace(0.01, horizon - 0.01, 400)
    kin = max(
        np.abs(
            (-R(t + 2 * h) + 8 * R(t + h) - 8 * R(t - h) + R(t - 2 * h)) / (12 * h)
            - R(t) @ hat(flat_reference(t, m, g, J, period)[4])
        ).max()
        for t in ts
    )
    tt = np.arange(0.0, horizon, 1e-4)
    data = [flat_reference(t, m, g, J, period) for t in tt]
    Om = np.array([d[4] for d in data])
    Omd = np.array([d[5] for d in data])
    integ = np.concatenate(
        [np.zeros((1, 3)), np.cumsum(0.5 * (Omd[1:] + Omd[:-1]) * 1e-4, axis=0)]
    )
    integral_err = float(np.max(np.abs(integ - (Om - Om[0]))))
    return {
        "kinematic_error_Rdot": float(kin),
        "integral_error_Omega": integral_err,
        "max_Omega": float(np.abs(Om).max()),
        "max_Omega_dot": float(np.abs(Omd).max()),
    }


def quadrotor_figure_eight(horizon=16.0, dt=DT, save=True):
    sys_, layout, plant, plant_lin, _ = build("quadrotor", {})
    p = dict(sys_.params)
    m, g, J = p["m"], p["g"], np.asarray(p["J"])
    tt = np.arange(0.0, horizon + dt, dt / 2)
    X, U, Acc = [], [], []
    for t in tt:
        x, v, a, R, Om, Omd, f, M = flat_reference(t, m, g, J)
        X.append(layout.pack(x=x, R=R, x_dot=v, Omega_R=Om))
        U.append(np.r_[f, M])
        Acc.append(np.r_[a, Omd])
    ref = Reference(tt, np.array(X), np.array(U), np.array(Acc))
    checks = flatness_checks(m, g, J, horizon=horizon)
    x0v = layout.unpack(ref.X[0])
    x0 = layout.pack(
        x=x0v["x"] + np.array([0.3, -0.3, -0.3]),
        R=rotate(x0v["R"], 20, E3),
        x_dot=x0v["x_dot"],
        Omega_R=x0v["Omega_R"],
    )
    Q = [10.0] * 6 + [1.0] * 6
    return track(
        "quadrotor_figure_eight",
        "quadrotor",
        ref,
        x0,
        Q,
        [1.0] * 4,
        horizon,
        {"f": "N", "M": "N·m"},
        extra={"flatness_checks": checks},
        dt=dt,
        save=save,
    )


# ---------------------------------------------------------------------------
# Quadrotor with a suspended load: a reference simulated with the nonlinear model
# ---------------------------------------------------------------------------

PERIOD_CIRCLE = 5.0


def load_reference(horizon, dt=DT, period=PERIOD_CIRCLE):
    sys_, layout, plant, plant_lin, _ = build("quadrotor_suspended_load", {})
    p = dict(sys_.params)
    x_eq = layout.pack(x_L=np.zeros(3), q=-E3, R=np.eye(3))
    u_eq, _ = equilibrium_input(plant, x_eq, p)
    hover = lqr(plant_lin, x_eq, u_eq, p, Q=1.0, R=1.0, module=plant)

    def setpoint(t):
        w = 2 * np.pi / period
        s = min(t / 2.0, 1.0)  # ramp the radius in over 2 s
        return s * np.array([np.cos(w * t) - 1.0, np.sin(w * t), 0.0])

    def nominal(t, x):
        xe = x_eq.copy()
        xe[layout.slot("x_L").slice] = setpoint(t)
        return u_eq - hover.K @ hover.error(x) + hover.K @ hover.error(xe)

    # the nominal run, sampled every dt/2: simulate at dt/2 and keep every state
    nom = simulate(plant, x_eq, horizon, p, u=nominal, dt=dt / 2)
    U = np.array([nominal(t, x) for t, x in zip(nom.t, nom.state)])
    Acc = np.array([accelerations(plant, x, u, p) for x, u in zip(nom.state, U)])
    return Reference(nom.t, nom.state, U, Acc), layout, plant, plant_lin, p


def suspended_load_circle(horizon=12.0, dt=DT, save=True, validate=True):
    ref, layout, plant, plant_lin, p = load_reference(horizon, dt)
    x0v = layout.unpack(ref.X[0])
    q0 = expm_so3(np.deg2rad(15) * np.array([1.0, 0.0, 0.0])) @ x0v["q"]
    x0 = layout.pack(x_L=x0v["x_L"] + np.array([0.3, 0.2, -0.2]), q=q0, R=rotate(x0v["R"], 20, E3))
    model = model_of(plant)
    q_vel = next(v for v in model.variables if v.manifold == "S2").vel
    extra = {
        "reference": "simulated (hover LQR with a set point on a 1 m circle, period 5 s)",
        "max_reference_omega_q": float(np.max(np.linalg.norm(ref.X[:, q_vel], axis=1))),
    }
    if validate:
        n = len(ref.t)
        samples = [int(n * f) for f in (0.3, 0.5, 0.7, 0.9)]
        extra["reduction_vs_finite_differences"] = validate_reduction(
            plant, plant_lin, ref, p, samples
        )
    Q = [10.0] * 3 + [10.0] * 2 + [1.0] * 3 + [1.0] * 8
    return track(
        "suspended_load_circle",
        "quadrotor_suspended_load",
        ref,
        x0,
        Q,
        [1.0] * 4,
        horizon,
        {"f": "N", "M": "N·m"},
        extra=extra,
        dt=dt,
        save=save,
    )


TRACKING = {
    "quadrotor_figure_eight": quadrotor_figure_eight,
    "suspended_load_circle": suspended_load_circle,
}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("names", nargs="*")
    parser.add_argument("--all", action="store_true")
    args = parser.parse_args(argv)
    for name in list(TRACKING) if args.all else args.names:
        m = TRACKING[name]()
        print(json.dumps(m, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
