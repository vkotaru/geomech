"""Integrators on the state manifold of a generated model: ``simulate`` and ``Trajectory``.

Two methods:

- ``"rkmk4"`` (default): fixed-step Runge–Kutta–Munthe-Kaas of order 4. S² and SO(3)
  points are moved by rotations, so ``|q| = 1`` and ``RᵀR = I`` hold to round-off with no
  projection; every other slot (R, R³ and all velocities) gets classical RK4.
- ``"solve_ivp"``: ``scipy.integrate.solve_ivp`` (DOP853 by default) on the embedded
  state; the outputs are projected onto the manifold and the drift before projection is
  reported.

scipy is imported only by ``method="solve_ivp"``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable

import numpy as np
from numpy.typing import ArrayLike

from geomech.simulation.model import (
    STATE_TOL,
    Model,
    check_state,
    model_of,
    solve_accelerations,
)

__all__ = ["simulate", "Trajectory", "constraint_drift", "project_state"]

_METHODS = ("rkmk4", "solve_ivp")


# ---------------------------------------------------------------------------
# Trajectory
# ---------------------------------------------------------------------------


@dataclass
class Trajectory:
    """The result of ``simulate``.

    Attributes:
        t: the output times, shape ``(N,)``.
        state: the states, shape ``(N, n_state)`` (``LAYOUT`` order). With
            ``method="solve_ivp"`` they are projected onto the manifold.
        constraint_drift: ``{name: (N,) array}`` of the constraint violations of the
            integrated (unprojected) states: ``"|q|-1"`` and ``"q.omega_q"`` (``|q·ω|``)
            per S² point, ``"R^T R-I"`` (Frobenius norm) per SO(3) point; empty for a
            system on R and R³ only.
        method: ``"rkmk4"`` or ``"solve_ivp"``.
        p: the parameters used.
        module: the generated module.
        info: integrator details (``n_steps``, ``dt`` or ``nfev``, ``rtol``, …).
    """

    t: np.ndarray
    state: np.ndarray
    constraint_drift: dict
    method: str
    p: dict
    module: Any = field(repr=False)
    info: dict = field(default_factory=dict)

    def energy(self, fn: Callable | None = None) -> np.ndarray:
        """``fn(state, p)`` at every output state, ``(N,)``.

        *fn* defaults to the module's ``energy`` (``generate_module(...,
        lagrangian=L)``); ``geomech.codegen.energy_function(L, layout)`` gives one for any
        Lagrangian.
        """
        if fn is None:
            fn = getattr(self.module, "energy", None)
            if fn is None:
                raise TypeError(
                    "Trajectory.energy: the module has no energy(); generate it with"
                    " generate_module(..., lagrangian=L) or pass energy_function(L, layout)"
                )
        # The generated functions take one state at a time; at about 20 µs per state
        # (double rigid pendulum) this loop is a few percent of the simulation that
        # produced the states, so it is not vectorised.
        return np.fromiter((float(fn(x, self.p)) for x in self.state), float, len(self.state))

    def max_drift(self) -> float:
        """The largest constraint violation over the trajectory (0 without S²/SO(3))."""
        return max((float(np.max(v)) for v in self.constraint_drift.values()), default=0.0)

    def __getitem__(self, key: str) -> np.ndarray:
        """The values of one state slot over time, by printed symbol (``"q"``,
        ``"\\\\omega_{q}"``) or generated name (``"omega_q"``): ``(N,)``, ``(N, 3)`` or
        ``(N, 3, 3)`` for an SO(3) point."""
        lay = self.module.LAYOUT
        for slot in lay["configuration"] + lay["velocity"]:
            if key in (slot["symbol"], slot["name"]):
                vals = self.state[:, slice(*slot["slice"])]
                shape = tuple(slot["shape"])
                return vals[:, 0] if shape == () else vals.reshape((len(self.t), *shape))
        names = [s["name"] for s in lay["configuration"] + lay["velocity"]]
        raise KeyError(f"Trajectory: no state slot {key!r} (slots: {names})")


def constraint_drift(model: Model, states: np.ndarray) -> dict:
    """``{name: (N,) array}`` of the manifold constraint violations of *states*."""
    states = np.atleast_2d(states)
    out = {}
    for var in model.s2:
        q, w = states[:, var.cfg], states[:, var.vel]
        name = _name(model, var.key)
        out[f"|{name}|-1"] = np.abs(np.linalg.norm(q, axis=1) - 1.0)
        out[f"{name}.{_name(model, var.vel_key)}"] = np.abs(np.einsum("ij,ij->i", q, w))
    for var in model.so3:
        R = states[:, var.cfg].reshape(-1, 3, 3)
        name = _name(model, var.key)
        RtR = np.einsum("nki,nkj->nij", R, R)
        out[f"{name}^T {name}-I"] = np.linalg.norm(RtR - np.eye(3), axis=(1, 2))
    return out


def _name(model: Model, key: str) -> str:
    for slot in model.layout["configuration"] + model.layout["velocity"]:
        if slot["symbol"] == key:
            return slot["name"]
    return key


def project_state(model: Model, state: np.ndarray) -> np.ndarray:
    """*state* projected onto the state manifold: each S² point normalised and its
    velocity made tangent (``ω ← ω − (q·ω) q``), each SO(3) point replaced by its polar
    factor (the nearest rotation, ``U Vᵀ`` of the SVD). A projected state passes
    ``check_state``, so it can start another ``simulate``."""
    x = np.array(state, dtype=float)
    for var in model.s2:
        q = x[var.cfg] / np.linalg.norm(x[var.cfg])
        x[var.cfg] = q
        x[var.vel] -= (q @ x[var.vel]) * q
    for var in model.so3:
        U, _, Vt = np.linalg.svd(x[var.cfg].reshape(3, 3))
        R = U @ Vt
        if np.linalg.det(R) < 0:  # not reachable from a rotation by small drift
            U[:, -1] = -U[:, -1]
            R = U @ Vt
        x[var.cfg] = R.reshape(9)
    return x


# ---------------------------------------------------------------------------
# simulate
# ---------------------------------------------------------------------------


def _input_function(model: Model, u) -> Callable[[float, np.ndarray], np.ndarray]:
    if u is None:
        zero = np.zeros(model.n_inputs)
        return lambda t, x: zero
    if callable(u):

        def ufun(t, x):
            val = u(t, x)
            if val is None:
                raise TypeError(
                    f"simulate: u(t, state) returned None at t = {t:.6g} (a controller"
                    " missing its return statement?); return the inputs, or pass u=None"
                    " for zero inputs"
                )
            try:
                return model.as_input(val, "u(t, state)")
            except ValueError as exc:
                raise ValueError(f"{exc} (at t = {t:.6g})") from None

        return ufun
    const = model.as_input(u)
    return lambda t, x: const


def simulate(
    module: Any,
    x0: ArrayLike,
    t_span: float | tuple[float, float],
    p: dict,
    u: ArrayLike | Callable | None = None,
    method: str = "rkmk4",
    dt: float = 1e-3,
    *,
    rtol: float = 1e-10,
    atol: float = 1e-12,
    ivp_method: str = "DOP853",
    tol: float = STATE_TOL,
    output_every: int = 1,
) -> Trajectory:
    """Integrate a generated standard-form model from *x0* over *t_span*.

    Args:
        module: a module from ``geomech.codegen.generate_module`` of a
            ``to_standard_form`` result (``mass``, ``bias``, ``input_matrix``,
            ``kinematics``).
        x0: the initial state, flat in ``LAYOUT`` order. It must lie on the manifold
            (``|q| = 1``, ``q·ω = 0``, ``R`` orthogonal with ``det R = +1``, to *tol*),
            else ``ValueError`` naming the slot.
        t_span: ``(t0, tf)``, or ``tf`` for ``(0, tf)``.
        p: the parameters, keyed by printed symbol.
        u: ``None`` (zero inputs), a constant 1-D array, or a callable ``u(t, state)``
            returning the inputs (e.g. an ``LQRController``). The inputs must be finite;
            a callable that returns ``None`` raises ``TypeError``. A callable is
            evaluated at the RKMK4 stage states and, with ``"solve_ivp"``, at trial
            steps that may be rejected and at times that are not monotone: it should be
            a pure function of ``(t, state)``.
        method: ``"rkmk4"`` (fixed step, default) or ``"solve_ivp"``.
        dt: the step of ``"rkmk4"`` (shortened so that it divides ``tf − t0``) and the
            output spacing of ``"solve_ivp"``.
        rtol: passed to ``scipy.integrate.solve_ivp``.
        atol: passed to ``scipy.integrate.solve_ivp``.
        ivp_method: the ``solve_ivp`` method (``"DOP853"``, ``"RK45"``, …).
        tol: the tolerance of the *x0* check.
        output_every: store every ``output_every``-th step of the ``dt`` grid (and
            always the last), to thin long runs; the integration itself is unchanged.

    The accelerations come from ``geomech.simulation.accelerations`` (the S² rows solved
    on the tangent plane). The stored states are on the manifold (``project_state``), so
    ``traj.state[-1]`` is a valid *x0* for the next segment. See ``docs/simulation.md``.
    """
    if method not in _METHODS:
        raise ValueError(f"simulate: method must be one of {_METHODS}, got {method!r}")
    model = model_of(module).require("mass", "bias", "input_matrix", "kinematics")
    x0 = check_state(module, x0, tol=tol, what="x0")
    t0, tf = _span(t_span)
    if not (dt > 0 and np.isfinite(dt)):
        raise ValueError(f"simulate: dt must be positive, got {dt}")
    if isinstance(output_every, bool) or not isinstance(output_every, (int, np.integer)):
        raise ValueError(f"simulate: output_every must be a positive int, got {output_every!r}")
    if output_every < 1:
        raise ValueError(f"simulate: output_every must be a positive int, got {output_every}")
    n_steps = max(1, int(np.ceil((tf - t0) / dt - 1e-9)))
    keep = np.unique(np.r_[np.arange(0, n_steps + 1, int(output_every)), n_steps])
    ufun = _input_function(model, u)
    p = dict(p or {})
    solve_accelerations(model, x0, ufun(t0, x0), p, check=True)  # singular system: raise now
    if method == "rkmk4":
        t, states = _rkmk4(model, x0, t0, tf, n_steps, p, ufun, keep)
        drift = constraint_drift(model, states)
        info = {"n_steps": n_steps, "dt": (tf - t0) / n_steps, "nfev": 4 * n_steps}
    else:
        t, raw, info = _solve_ivp(
            model, x0, t0, tf, n_steps, p, ufun, rtol, atol, ivp_method, keep
        )
        drift = constraint_drift(model, raw)
        states = np.array([project_state(model, x) for x in raw])
    info["output_every"] = int(output_every)
    return Trajectory(t, states, drift, method, p, module, info)


def _span(t_span) -> tuple[float, float]:
    if np.ndim(t_span) == 0:
        t0, tf = 0.0, float(t_span)
    else:
        if len(t_span) != 2:
            raise ValueError(f"simulate: t_span must be (t0, tf) or tf, got {t_span!r}")
        t0, tf = float(t_span[0]), float(t_span[1])
    if not (np.isfinite(t0) and np.isfinite(tf)):
        raise ValueError(f"simulate: t_span must be finite, got ({t0}, {tf})")
    if not tf > t0:
        raise ValueError(f"simulate: t_span must have tf > t0, got ({t0}, {tf})")
    return t0, tf


# ---------------------------------------------------------------------------
# RKMK4
# ---------------------------------------------------------------------------


class _LieRK:
    """The product-group form of the state used by RKMK4.

    An algebra element is ``(e, s, r)``: ``e`` the increments of the vector-space slots
    (R/R³ configuration entries, then every velocity entry), ``s`` one 3-vector per S²
    point (the rotation acting on ``q`` from the left) and ``r`` one 3-vector per SO(3)
    point (acting on ``R`` from the right).
    """

    def __init__(self, model: Model):
        self.model = model
        idx = [np.arange(v.cfg.start, v.cfg.stop) for v in model.euclid]
        vel = [np.arange(v.vel.start, v.vel.stop) for v in model.euclid]
        self.eu_cfg = np.concatenate(idx) if idx else np.zeros(0, dtype=int)
        self.eu_vel = np.concatenate(vel) if vel else np.zeros(0, dtype=int)
        self.m = self.eu_cfg.size
        self.nc = model.n_config

    def generator(self, t, x, p, ufun):
        """The algebra element ``A(x)`` with ``ẋ = A(x)·x``: rates of the vector slots,
        ``ω`` of each S² point (``q̇ = ω × q``) and ``Ω`` of each SO(3) point
        (``Ṙ = R Ω̂``)."""
        a = solve_accelerations(self.model, x, ufun(t, x), p)
        e = np.concatenate([x[self.eu_vel], a])
        s = np.array([x[v.vel] for v in self.model.s2]).reshape(-1, 3)
        r = np.array([x[v.vel] for v in self.model.so3]).reshape(-1, 3)
        return e, s, r

    def act(self, U, x):
        """``exp(U) · x``: translate the vector slots, rotate each ``q`` to
        ``expm(s) q`` and each ``R`` to ``R expm(r)``."""
        e, s, r = U
        y = x.copy()
        y[self.eu_cfg] += e[: self.m]
        y[self.nc :] += e[self.m :]
        for j, v in enumerate(self.model.s2):
            y[v.cfg] = _expm(s[j]) @ x[v.cfg]
        for j, v in enumerate(self.model.so3):
            y[v.cfg] = (x[v.cfg].reshape(3, 3) @ _expm(r[j])).reshape(9)
        return y

    @staticmethod
    def dexpinv(U, G, h):
        """``h · dexp⁻¹_U(G)`` truncated after the ``ad²`` term (Bernoulli numbers
        1, −½, 1/12; the ``ad³`` coefficient is zero). Left action (S²):
        ``G − ½[U, G] + 1/12 [U, [U, G]]``; right action (SO(3)): ``dexp⁻¹_{−U}``, the
        sign of the odd term flipped. ``[a, b] = a × b`` on so(3) ≅ R³; the vector part
        commutes."""
        e, s, r = U
        ge, gs, gr = G
        if s.size:
            us = _cross(s, gs)
            gs = gs - 0.5 * us + _cross(s, us) / 12.0
        if r.size:
            ur = _cross(r, gr)
            gr = gr + 0.5 * ur + _cross(r, ur) / 12.0
        return h * ge, h * gs, h * gr


def _cross(a, b):
    """Row-wise ``a × b`` of two ``(n, 3)`` arrays (``np.cross`` is slow on tiny arrays)."""
    a0, a1, a2 = a[:, 0], a[:, 1], a[:, 2]
    b0, b1, b2 = b[:, 0], b[:, 1], b[:, 2]
    return np.stack([a1 * b2 - a2 * b1, a2 * b0 - a0 * b2, a0 * b1 - a1 * b0], axis=1)


def _expm(v):
    """Rodrigues' formula for one 3-vector (``geomech.numeric.expm_so3``, inlined)."""
    x, y, z = v
    th2 = x * x + y * y + z * z
    if th2 < 1e-16:
        a, b = 1.0 - th2 / 6.0, 0.5 - th2 / 24.0
    else:
        th = np.sqrt(th2)
        a, b = np.sin(th) / th, (1.0 - np.cos(th)) / th2
    K = np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])
    return np.eye(3) + a * K + b * (K @ K)


def _lin(*terms):
    """Σ cᵢ Uᵢ for (c, U) pairs, componentwise."""
    return tuple(sum(c * U[k] for c, U in terms) for k in range(3))


def _rkmk4(model, x0, t0, tf, n_steps, p, ufun, keep):
    """Runge–Kutta–Munthe-Kaas with the classical RK4 tableau (Munthe-Kaas 1999).

    With ``y = exp(U)·y_n``, the curve ``U(t)`` in the Lie algebra solves
    ``U̇ = dexp⁻¹_U(A(exp(U)·y_n))``, ``U(t_n) = 0``. One step applies classical RK4 to
    that ODE and maps back with the exponential:

        K₁ = h F(0),         K₂ = h F(K₁/2),     K₃ = h F(K₂/2),     K₄ = h F(K₃),
        U  = (K₁ + 2K₂ + 2K₃ + K₄)/6,            y_{n+1} = exp(U)·y_n,

    with ``F(V) = dexp⁻¹_V(A(exp(V)·y_n))`` at the stage times ``t, t+h/2, t+h/2, t+h``.
    Since ``U = O(h)``, truncating ``dexp⁻¹`` after the ``ad_U²`` term changes ``F`` by
    ``O(h⁴)`` (``B₃ = 0``), so each ``Kᵢ`` by ``O(h⁵)``: the local error stays
    ``O(h⁵)`` and the method has the order 4 of the underlying RK scheme (Munthe-Kaas,
    "High order Runge–Kutta methods on manifolds", Appl. Numer. Math. 29 (1999); Hairer,
    Lubich, Wanner, *Geometric Numerical Integration*, §IV.8.2). On the vector-space
    slots every bracket vanishes and the step is classical RK4. ``q`` and ``R`` are only
    ever multiplied by rotation matrices, so their constraints hold to round-off.

    ``ω`` is a vector slot, translated while ``q`` is rotated, so ``q·ω`` would drift at
    the truncation error. After each step ``ω ← ω − (q·ω) q``: a change of the size of
    the local error, so the order is still 4, and ``q·ω`` stays at round-off.

    Only the steps in *keep* (sorted indices into ``0..n_steps``) are stored.
    """
    lie = _LieRK(model)
    h = (tf - t0) / n_steps
    t = t0 + h * np.arange(n_steps + 1)
    t[-1] = tf
    out = np.empty((keep.size, model.n_state))
    out[0] = x = x0
    store = np.zeros(n_steps + 1, dtype=bool)
    store[keep] = True
    k = 1
    s2 = [(v.cfg, v.vel) for v in model.s2]
    zero = (np.zeros(lie.m + model.n_accel), np.zeros((len(model.s2), 3)))
    zero += (np.zeros((len(model.so3), 3)),)
    for n in range(n_steps):
        tn = t[n]
        K1 = lie.dexpinv(zero, lie.generator(tn, x, p, ufun), h)
        U2 = _lin((0.5, K1))
        K2 = lie.dexpinv(U2, lie.generator(tn + 0.5 * h, lie.act(U2, x), p, ufun), h)
        U3 = _lin((0.5, K2))
        K3 = lie.dexpinv(U3, lie.generator(tn + 0.5 * h, lie.act(U3, x), p, ufun), h)
        K4 = lie.dexpinv(K3, lie.generator(tn + h, lie.act(K3, x), p, ufun), h)
        U = _lin((1 / 6, K1), (1 / 3, K2), (1 / 3, K3), (1 / 6, K4))
        x = lie.act(U, x)
        for cfg, vel in s2:  # keep ω tangent (see above)
            q = x[cfg]
            x[vel] -= (q @ x[vel]) * q
        if store[n + 1]:
            out[k] = x
            k += 1
    return t[keep], out


# ---------------------------------------------------------------------------
# solve_ivp
# ---------------------------------------------------------------------------


def _solve_ivp(model, x0, t0, tf, n_steps, p, ufun, rtol, atol, ivp_method, keep):
    try:
        from scipy.integrate import solve_ivp
    except ImportError as exc:  # pragma: no cover - depends on the environment
        raise ImportError(
            "simulate(method='solve_ivp') needs scipy: pip install 'geomech[sim]'"
        ) from exc
    kin = model.module.kinematics

    def rhs(t, y):
        return np.concatenate([kin(y), solve_accelerations(model, y, ufun(t, y), p)])

    t_eval = np.linspace(t0, tf, n_steps + 1)[keep]
    sol = solve_ivp(rhs, (t0, tf), x0, method=ivp_method, t_eval=t_eval, rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(f"simulate: solve_ivp failed: {sol.message}")
    info = {"nfev": int(sol.nfev), "rtol": rtol, "atol": atol, "ivp_method": ivp_method}
    return sol.t, sol.y.T, info
