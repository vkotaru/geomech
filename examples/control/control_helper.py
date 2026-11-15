"""Shared machinery of the LQR examples in ``examples/control/``.

One example is an ``Example`` record: a system of ``geomech.systems`` (used as it is),
an equilibrium configuration, the weights, an initial perturbation and a horizon. ``run``
takes it through the whole chain:

1. ``compute_eom`` → ``to_standard_form`` and ``linearize``;
2. NumPy code generation of the plant (standard form, with ``energy``) and of the
   linearization, from one ``StateLayout``;
3. the equilibrium input ``u_eq``: the least-squares solution of ``bias + G u = 0`` at the
   equilibrium configuration (S² rows on the tangent plane), with the residual reported;
4. the reduced linearization ``δẋ = A δx + B δu`` from
   ``geomech.simulation.reduced_linearization`` (the coordinates of ``lqr``: S²
   perturbations in a 2-D tangent basis, SO(3) as ``vee(log(R_eqᵀR))``), optionally
   without cyclic coordinates (``drop``: a configuration perturbation whose columns of
   ``A`` and ``B`` are zero, such as a wheel angle);
5. controllability: the controllable subspace (an orthonormal Krylov basis of
   ``[B, AB, A²B, …]``), the eigenvalues of ``A`` on its complement (the uncontrollable
   modes) and a PBH cross-check. A fully controllable system is designed with
   ``geomech.simulation.lqr`` (which also checks the equilibrium); otherwise LQR is
   designed on the controllable subspace;
6. the closed loop ``u = u_eq − K δx(state)`` simulated with the nonlinear generated model
   (RKMK4 on the manifold, dt = 1 ms, no input saturation);
7. metrics: the initial error per group (position m, attitude rad, velocity), the
   normalized error ``ν(t) = √(V(z(t)) / V(z(0)))`` with ``V = zᵀPz`` the LQR cost-to-go,
   the settling time ``t_s`` (last exit of ν from the 2 % band), convergence (ν ≤ 2 %
   over the last 10 % of the horizon), the peak input deviation per input with its unit,
   the constraint drift. Written to ``results/<name>.json``; the time histories to
   ``results/<name>.npz``.

Run the examples with ``python examples/control/run.py <name>`` or ``--all``.
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

import numpy as np

from geomech import linearize, to_standard_form
from geomech.codegen import load_module, module_source
from geomech.numeric import expm_so3
from geomech.simulation import (
    LQRController,
    lqr,
    model_of,
    reduced_linearization,
    simulate,
    tangent_basis,
)
from geomech.systems import get

RESULTS = Path(__file__).resolve().parent / "results"
DT = 1e-3
OUTPUT_EVERY = 10  # store every 10th step (100 Hz) for the metrics and the saved histories
SETTLE_FRACTION = 0.02


@dataclass
class Example:
    """One closed-loop example.

    Attributes:
        name: registry name.
        system: name in ``geomech.systems``; ``kwargs`` its builder arguments.
        equilibrium: ``f(sys, layout) -> dict`` of slot values of the equilibrium
            configuration (``layout.pack`` keys; velocities are zero).
        perturbation: ``f(sys, layout, x_eq) -> x0``, the initial state.
        Q: the state weight, as a function ``f(labels) -> diagonal`` of the reduced
            coordinates (after ``drop``).
        R: the input weight (scalar or diagonal).
        horizon: simulated time (s).
        drop: labels of cyclic configuration coordinates removed from the design.
        units: the unit of each input (``"N"`` or ``"N·m"``), for the peak-input report.
        description: what the example shows; notes: caveats.
    """

    name: str
    system: str
    equilibrium: Callable
    perturbation: Callable
    Q: Callable
    R: Any
    horizon: float
    kwargs: dict = field(default_factory=dict)
    drop: tuple = ()
    units: dict = field(default_factory=dict)  # input name → "N" or "N·m"
    description: str = ""
    notes: str = ""


# ---------------------------------------------------------------------------
# Building blocks
# ---------------------------------------------------------------------------

_BUILT: dict = {}


def build(system: str, kwargs: dict):
    """``(sys, layout, plant, plant_lin, derive_seconds)``, cached per system."""
    key = (system, tuple(sorted(kwargs.items())))
    if key not in _BUILT:
        sys_ = get(system, **kwargs)
        t0 = time.perf_counter()
        eom = sys_.derive()
        t_derive = time.perf_counter() - t0
        inputs = list(sys_.inputs)
        sf = to_standard_form(eom, sys_.variables, inputs)
        lin = linearize(eom, sys_.variables, inputs)
        layout = sys_.layout()
        tag = f"{system}_{abs(hash(key)) % 10**8}"
        plant = load_module(module_source(sf, layout, lagrangian=sys_.L), f"plant_{tag}")
        plant_lin = load_module(module_source(lin, layout), f"lin_{tag}")
        _BUILT[key] = (sys_, layout, plant, plant_lin, t_derive)
    return _BUILT[key]


def equilibrium_input(plant, x_eq, p) -> tuple[np.ndarray, float]:
    """``(u_eq, relative residual)``: least squares of ``bias + G u = 0`` at ``x_eq`` with
    the S² rows projected on their tangent planes."""
    model = model_of(plant)
    bias = np.asarray(plant.bias(x_eq, p), dtype=float)
    G = np.asarray(plant.input_matrix(x_eq, p), dtype=float).reshape(len(bias), -1)
    P = np.eye(len(bias))
    for var in model.s2:
        q = x_eq[var.cfg]
        P[var.acc, var.acc] = np.eye(3) - np.outer(q, q)
    u, *_ = np.linalg.lstsq(P @ G, -P @ bias, rcond=None)
    r = P @ (bias + G @ u)
    scale = max(1.0, float(np.max(np.abs(bias))), float(np.max(np.abs(G @ u))))
    return u, float(np.max(np.abs(r))) / scale


def controllable_subspace(A, B, tol=1e-9) -> np.ndarray:
    """An orthonormal basis of the controllable subspace (block Krylov iteration with
    re-orthogonalisation; ``tol`` relative to the largest singular value of ``B``)."""
    n = A.shape[0]
    scale = max(np.linalg.norm(B, 2), 1e-300)
    U, s, _ = np.linalg.svd(B, full_matrices=False)
    V = U[:, s > tol * scale]
    W = V
    while V.shape[1] < n and W.shape[1]:
        W = A @ W
        W = W - V @ (V.T @ W)
        W = W - V @ (V.T @ W)
        U, s, _ = np.linalg.svd(W, full_matrices=False)
        W = U[:, s > tol * max(np.linalg.norm(A, 2), 1.0) * scale] if s.size else W[:, :0]
        V = np.hstack([V, W])
    return V


def pbh_uncontrollable(A, B, tol=1e-7) -> list:
    """The eigenvalues λ of ``A`` with ``rank [A − λI, B] < n`` (PBH test)."""
    n = A.shape[0]
    out = []
    for lam in np.linalg.eigvals(A):
        s = np.linalg.svd(np.hstack([A - lam * np.eye(n), B]), compute_uv=False)
        if s[-1] < tol * max(1.0, s[0]):
            out.append(complex(lam))
    return out


def uncontrollable_eigenvalues(A, V) -> list:
    """The eigenvalues of ``A`` on the quotient by the controllable subspace (an
    ``A``-invariant subspace with orthonormal basis *V*): the uncontrollable modes."""
    n = A.shape[0]
    if V.shape[1] == n:
        return []
    U, _, _ = np.linalg.svd(V, full_matrices=True)
    W = U[:, V.shape[1] :]
    return [complex(z) for z in np.linalg.eigvals(W.T @ A @ W)]


def label_group(label: str) -> str:
    """``"velocity"``, ``"attitude"`` (S¹, S² and SO(3) configuration errors, rad) or
    ``"position"`` (R³ and R configuration errors, m) for a reduced coordinate."""
    if is_velocity(label):
        return "velocity"
    if label.startswith(("gamma_", "xi_", "eta_")):
        return "attitude"
    return "position"


def settling(t, nu, fraction=SETTLE_FRACTION, tail=0.1):
    """``(t_s, converged)`` of the normalized error *nu* (``nu[0] = 1``): ``t_s`` is the
    last exit time from the band ``nu ≤ fraction`` (the first sample after the last one
    above it; ``None`` if the last sample is above it), and *converged* requires
    ``nu ≤ fraction`` over the whole final *tail* fraction of the horizon."""
    above = np.flatnonzero(nu > fraction)
    if not above.size:
        ts = float(t[0])
    else:
        k = above[-1] + 1
        ts = float(t[k]) if k < len(t) else None
    window = t >= t[-1] - tail * (t[-1] - t[0])
    return ts, bool(ts is not None and np.all(nu[window] <= fraction))


# ---------------------------------------------------------------------------
# Design and simulation
# ---------------------------------------------------------------------------


@dataclass
class Design:
    ctrl: LQRController
    keep: np.ndarray  # indices of the reduced coordinates used by the controller
    labels: list
    n_reduced: int
    rank: int
    uncontrollable: list
    cl_max_real: float
    cl_min_damping: float  # smallest damping ratio −Re λ / |λ| of the controlled modes
    cl_fastest: float  # the largest |λ| of the closed loop (RKMK4 needs |λ| dt below ~2.8)
    used_library_lqr: bool
    V: np.ndarray  # orthonormal basis of the controllable subspace (kept coordinates)
    P: np.ndarray  # Riccati solution in the coordinates Vᵀ z

    def lyapunov(self, z_keep) -> float:
        """``zᵀPz`` of the reduced error (on the controllable subspace): the LQR cost-to-go,
        a weighted norm of the error that the closed loop decreases."""
        zc = self.V.T @ z_keep
        return float(zc @ self.P @ zc)


def design(ex: Example, sys_, layout, plant, plant_lin, x_eq, u_eq) -> Design:
    p = dict(sys_.params)
    red = reduced_linearization(plant_lin, x_eq, u_eq, p)
    A, B, T, bases, labels = red.A, red.B, red.T, red.bases, red.labels
    keep = np.array([i for i, lab in enumerate(labels) if lab not in ex.drop])
    dropped = np.setdiff1d(np.arange(len(labels)), keep)
    if dropped.size:
        assert np.allclose(A[:, dropped], 0) and np.allclose(B[dropped], 0), "not cyclic"
    Ak, Bk = A[np.ix_(keep, keep)], B[keep]
    kept_labels = [labels[i] for i in keep]
    Q = np.diag(np.asarray(ex.Q(kept_labels), dtype=float))
    n_u = Bk.shape[1]
    Rm = np.eye(n_u) * ex.R if np.ndim(ex.R) == 0 else np.diag(ex.R)
    V = controllable_subspace(Ak, Bk)
    rank = V.shape[1]
    unc = uncontrollable_eigenvalues(Ak, V)
    assert len(pbh_uncontrollable(Ak, Bk)) >= (1 if unc else 0)  # PBH agrees
    if rank == len(keep) and not dropped.size:
        ctrl = lqr(plant_lin, x_eq, u_eq, p, Q=Q, R=Rm, module=plant)
        assert np.allclose(ctrl.A, A) and np.allclose(ctrl.B, B)
        library = True
        K_full = ctrl.K
        V = np.eye(len(keep))
        P = ctrl.P
    else:
        from scipy.linalg import solve_continuous_are

        Ac, Bc, Qc = V.T @ Ak @ V, V.T @ Bk, V.T @ Q @ V
        P = solve_continuous_are(Ac, Bc, Qc, Rm)
        Kc = np.linalg.solve(Rm, Bc.T @ P) @ V.T
        K_full = np.zeros((n_u, len(labels)))
        K_full[:, keep] = Kc
        model = model_of(plant_lin)
        axes = {v.key: model.axis(v, p, "lqr") for v in model.s1}  # S¹ errors: angles
        ctrl = LQRController(K_full, P, A, B, T, x_eq, u_eq, labels, model, bases, axes)
        library = False
    eig = np.linalg.eigvals(Ak - Bk @ K_full[:, keep])
    eig_c = np.linalg.eigvals(V.T @ (Ak - Bk @ K_full[:, keep]) @ V)  # controlled modes
    cl = float(np.max(eig_c.real))
    damping = float(np.min(-eig_c.real / np.abs(eig_c)))
    fastest = float(np.max(np.abs(eig)))
    return Design(
        ctrl, keep, kept_labels, len(keep), rank, unc, cl, damping, fastest, library, V, P
    )


def _common(ex, d, residual, u_eq, T, dt, t_derive):
    return {
        "name": ex.name,
        "system": ex.system,
        "kwargs": ex.kwargs,
        "description": ex.description,
        "notes": ex.notes,
        "stabilizable": all(z.real < -1e-9 for z in d.uncontrollable),
        "equilibrium_residual": residual,
        "u_eq": u_eq.tolist(),
        "n_reduced": d.n_reduced,
        "dropped": [lab for lab in d.ctrl.labels if lab not in d.labels],
        "controllability_rank": d.rank,
        "uncontrollable_eigenvalues": [[z.real, z.imag] for z in d.uncontrollable],
        "closed_loop_max_real": d.cl_max_real,
        "closed_loop_min_damping": d.cl_min_damping,
        "closed_loop_fastest": d.cl_fastest,
        "library_lqr": d.used_library_lqr,
        "horizon": T,
        "dt": dt,
        "derive_seconds": t_derive,
    }


def initial_errors(d: Design, z0) -> dict:
    """Initial error per group (Euclidean norm of the coordinates of the group)."""
    out = {}
    for group in ("position", "attitude", "velocity"):
        idx = [i for i, lab in enumerate(d.labels) if label_group(lab) == group]
        out[group] = float(np.linalg.norm(z0[idx])) if idx else None
    return out


def run(ex: Example, horizon: float | None = None, dt: float = DT, save: bool = True) -> dict:
    """Run one example; return (and, with *save*, write) its metrics.

    Errors: the reduced error ``z`` of ``LQRController.error`` (without dropped
    coordinates). The initial error is reported per group (position m, attitude rad,
    velocity in the units of each coordinate); the settling metrics use the normalized
    error ``ν(t) = √(V(z(t)) / V(z(0)))`` with ``V = zᵀPz`` the LQR cost-to-go on the
    controllable subspace.
    """
    sys_, layout, plant, plant_lin, t_derive = build(ex.system, ex.kwargs)
    p = dict(sys_.params)
    x_eq = layout.pack(**ex.equilibrium(sys_, layout))
    u_eq, residual = equilibrium_input(plant, x_eq, p)
    d = design(ex, sys_, layout, plant, plant_lin, x_eq, u_eq)
    x0 = ex.perturbation(sys_, layout, x_eq)
    T = ex.horizon if horizon is None else horizon
    metrics = _common(ex, d, residual, u_eq, T, dt, t_derive)
    z0 = d.ctrl.error(x0)[d.keep]
    metrics["initial_error"] = initial_errors(d, z0)
    t0 = time.perf_counter()
    failure = None
    try:
        with np.errstate(all="ignore"):
            traj = simulate(
                plant, x0, T, p, u=d.ctrl, dt=dt, output_every=max(1, round(1e-2 / dt))
            )
    except (ValueError, np.linalg.LinAlgError, FloatingPointError) as exc:
        failure = str(exc).splitlines()[0][:200]
    metrics["simulate_seconds"] = time.perf_counter() - t0
    if failure is not None:
        metrics.update(converged=False, failure=failure, settling_time=None)
        if save:
            RESULTS.mkdir(exist_ok=True)
            (RESULTS / f"{ex.name}.json").write_text(json.dumps(metrics, indent=1) + "\n")
        return metrics
    Z = np.array([d.ctrl.error(x)[d.keep] for x in traj.state])
    v0 = d.lyapunov(Z[0])
    nu = np.sqrt(np.maximum([d.lyapunov(z) for z in Z], 0.0) / v0)
    ts, converged = settling(traj.t, nu)
    peaks = {}
    for slot in layout.input_slots:
        dev = np.abs(traj.u[:, slot.slice] - u_eq[slot.slice])
        peaks[slot.key] = {
            "peak_deviation": float(dev.max()),
            "peak_abs": float(np.abs(traj.u[:, slot.slice]).max()),
            "unit": ex.units.get(slot.key, ""),
        }
    E = traj.energy()
    metrics.update(
        final_normalized_error=float(nu[-1]),
        max_normalized_error=float(nu.max()),
        settling_time=ts,
        converged=converged,
        peak_inputs=peaks,
        constraint_drift=traj.max_drift(),
        energy_change=float(E[-1] - E[0]),
    )
    if save:
        RESULTS.mkdir(exist_ok=True)
        (RESULTS / f"{ex.name}.json").write_text(json.dumps(metrics, indent=1) + "\n")
        np.savez_compressed(
            RESULTS / f"{ex.name}.npz",
            t=traj.t,
            state=traj.state,
            u=traj.u,
            z=Z,
            normalized_error=nu,
            x_eq=x_eq,
            u_eq=u_eq,
            labels=np.array(d.labels),
        )
    return metrics


# ---------------------------------------------------------------------------
# Helpers for the example definitions
# ---------------------------------------------------------------------------

UP = np.array([0.0, 0.0, 1.0])
DOWN = -UP
FLIP = np.diag([1.0, -1.0, -1.0])  # a half turn about e1: a link hanging along −e3 points up


def tilt(q, angle_deg, axis) -> np.ndarray:
    """*q* rotated by *angle_deg* about *axis*."""
    axis = np.asarray(axis, float) / np.linalg.norm(axis)
    return expm_so3(np.deg2rad(angle_deg) * axis) @ np.asarray(q, float)


def rotate(R, angle_deg, axis) -> np.ndarray:
    """``R expm(angle · axis^)`` (a body-frame rotation of *R*)."""
    axis = np.asarray(axis, float) / np.linalg.norm(axis)
    return np.asarray(R, float) @ expm_so3(np.deg2rad(angle_deg) * axis)


def tangent(q, v) -> np.ndarray:
    q = np.asarray(q, float)
    v = np.asarray(v, float)
    return v - (v @ q) * q


def is_velocity(label: str) -> bool:
    """Whether a reduced coordinate is a velocity perturbation (``δẋ``, ``δω``, ``δΩ``)."""
    return "_dot" in label or label.startswith(("delta_omega", "delta_Omega"))


def weights(config=10.0, velocity=1.0, **by_prefix):
    """``Q(labels)``: *config* for configuration perturbations, *velocity* for velocity
    perturbations, overridden by the first matching label prefix in *by_prefix*."""

    def Q(labels):
        out = []
        for lab in labels:
            for prefix, w in by_prefix.items():
                if lab.startswith(prefix):
                    out.append(w)
                    break
            else:
                out.append(velocity if is_velocity(lab) else config)
        return out

    return Q


__all__ = [
    "Example",
    "build",
    "design",
    "run",
    "equilibrium_input",
    "controllable_subspace",
    "pbh_uncontrollable",
    "uncontrollable_eigenvalues",
    "settling",
    "label_group",
    "tilt",
    "rotate",
    "tangent",
    "weights",
    "tangent_basis",
]
