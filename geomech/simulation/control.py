"""LQR about an equilibrium from the variation-based linearization (``linearize``).

``lqr(lin_module, x_eq, u_eq, p, Q, R)`` evaluates the generated ``M_lin``, ``A_lin`` and
``B_lin`` at the equilibrium, reduces the S² perturbations to a 2-D tangent basis, solves
the continuous algebraic Riccati equation and returns an ``LQRController``: the gain ``K``
and the state error ``δx(state)`` on the manifold, callable as ``u(t, state)`` for
``simulate``.

The perturbation coordinates are those of ``linearize`` (``docs/vbl.md``):

==========  ==========================  =========================================
Variable    Configuration error         Velocity error
==========  ==========================  =========================================
R, R³       ``x − x_eq``                ``ẋ − ẋ_eq``
S²          ``Bᵀ (q_eq × q)`` (2-D)     ``Bᵀ (ω − ω_eq)`` (2-D)
SO(3)       ``vee(log(R_eqᵀ R))``       ``Ω − Ω_eq``
==========  ==========================  =========================================

``δq = ξ × q`` with ``ξ ⊥ q``, so ``ξ = q_eq × q`` to first order, and ``δR = R η̂`` so
``R = R_eq exp(η̂)``. ``B = [b₁ b₂]`` is an orthonormal basis of the tangent plane at
``q_eq`` (``tangent_basis``).

scipy (``scipy.linalg.solve_continuous_are``) is imported only when ``lqr`` is called.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from numpy.typing import ArrayLike

from geomech.numeric import log_so3
from geomech.simulation.model import STATE_TOL, Model, check_state, model_of

__all__ = ["lqr", "LQRController", "tangent_basis"]


def tangent_basis(q) -> np.ndarray:
    """An orthonormal basis ``[b₁ b₂]`` (3×2) of the plane orthogonal to the unit vector
    *q*, with ``b₁ × b₂ = q``. ``b₁`` is the coordinate axis least aligned with *q*,
    projected and normalised, so the basis is a deterministic function of *q*."""
    q = np.asarray(q, dtype=float).reshape(3)
    q = q / np.linalg.norm(q)
    e = np.eye(3)[int(np.argmin(np.abs(q)))]
    b1 = e - (e @ q) * q
    b1 /= np.linalg.norm(b1)
    b2 = np.cross(q, b1)
    return np.column_stack([b1, b2])


@dataclass
class LQRController:
    """The LQR gain about an equilibrium and the state error it acts on.

    ``u = u_eq − K δx(state)``; call the controller as ``u(t, state)`` (the signature
    ``simulate`` expects).

    Attributes:
        K: the gain, ``(n_inputs, n_reduced)``.
        P: the solution of the Riccati equation.
        A, B: the reduced linearization ``δẋ = A δx + B δu`` (``A = M⁻¹ A_lin`` and
            ``B = M⁻¹ B_lin`` after the S² reduction).
        T: ``(n_perturbation, n_reduced)``, block diagonal: ``I`` for the R, R³ and SO(3)
            perturbations, the tangent basis ``B_q`` at ``q_eq`` for each S² ``ξ`` and ``δω``.
        x_eq, u_eq: the equilibrium.
        labels: the name of each reduced coordinate (``"xi_q[0]"``, ``"eta_R[2]"``, …).
    """

    K: np.ndarray
    P: np.ndarray
    A: np.ndarray
    B: np.ndarray
    T: np.ndarray
    x_eq: np.ndarray
    u_eq: np.ndarray
    labels: list
    model: Model = field(repr=False)
    bases: dict = field(repr=False, default_factory=dict)

    @property
    def closed_loop_eigenvalues(self) -> np.ndarray:
        """The eigenvalues of ``A − B K``."""
        return np.linalg.eigvals(self.A - self.B @ self.K)

    def error(self, state) -> np.ndarray:
        """The reduced state error ``δx`` of *state* about ``x_eq`` (see the module
        docstring)."""
        x = np.asarray(state, dtype=float).reshape(-1)
        xe = self.x_eq
        cfg, vel = [], []
        for var in self.model.variables:
            if var.manifold == "S2":
                Bq = self.bases[var.key]
                cfg.append(Bq.T @ np.cross(xe[var.cfg], x[var.cfg]))
                vel.append(Bq.T @ (x[var.vel] - xe[var.vel]))
            elif var.manifold == "SO3":
                Re, R = xe[var.cfg].reshape(3, 3), x[var.cfg].reshape(3, 3)
                cfg.append(log_so3(Re.T @ R))
                vel.append(x[var.vel] - xe[var.vel])
            else:
                cfg.append(x[var.cfg] - xe[var.cfg])
                vel.append(x[var.vel] - xe[var.vel])
        return np.concatenate(cfg + vel)

    def __call__(self, t: float, state) -> np.ndarray:
        return self.u_eq - self.K @ self.error(state)


def _weight(W, n, name) -> np.ndarray:
    W = np.asarray(W, dtype=float)
    if W.ndim == 0:
        return float(W) * np.eye(n)
    if W.ndim == 1:
        W = np.diag(W)
    if W.shape != (n, n):
        raise ValueError(f"lqr: {name} must be a scalar, ({n},) or ({n}, {n}), got {W.shape}")
    return W


def lqr(
    lin_module: Any,
    x_eq: ArrayLike,
    u_eq: ArrayLike,
    p: dict,
    Q: ArrayLike,
    R: ArrayLike,
    *,
    a_eq: ArrayLike | None = None,
    tol: float = STATE_TOL,
) -> LQRController:
    """The continuous-time LQR gain of a generated linearization about an equilibrium.

    Evaluates ``M_lin``, ``A_lin``, ``B_lin`` of *lin_module* (``generate_module`` of a
    ``LinearizedSystem``) at ``(x_eq, a_eq, u_eq)``, reduces them to the tangent
    coordinates ``δx = Tᵀ δx_full`` (``Mr = Tᵀ M T``, ``Ar = Tᵀ A T``, ``Br = Tᵀ B``),
    forms ``Ā = Mr⁻¹ Ar``, ``B̄ = Mr⁻¹ Br``, solves ``scipy.linalg.solve_continuous_are``
    and returns ``K = R⁻¹ B̄ᵀ P`` in an ``LQRController``.

    The S² reduction: the 3-D perturbations ``ξ`` and ``δω`` of an S² point satisfy
    ``ξ · q = δω · q = 0`` at an equilibrium (where ``ω = 0``), so each has 2 degrees of
    freedom. In 3-D the normal direction is an uncontrollable mode with a zero
    eigenvalue, so the Riccati equation has no stabilising solution. ``T`` maps the 2-D
    tangent coordinates in the basis ``tangent_basis(q_eq)`` to the 3-D ones, and the S²
    rows (whose equations hold only on the tangent plane) are projected with ``Tᵀ``.
    A singular S² mass block (``Hat(q) Hat(q) m l²``) becomes invertible in this basis.

    Args:
        lin_module: the generated module of ``linearize`` (``M_lin``, ``A_lin``,
            ``B_lin``, ``LINEARIZATION``).
        x_eq: the equilibrium state (checked like ``simulate``'s ``x0``).
        u_eq: the equilibrium input.
        p: the parameters.
        Q: the weight of the reduced state: a scalar (times the identity), a vector (the
            diagonal) or a matrix.
        R: the weight of the input, in the same forms.
        a_eq: the accelerations at ``x_eq`` (zeros by default).
        tol: the tolerance of the ``x_eq`` check.
    """
    try:
        from scipy.linalg import solve_continuous_are
    except ImportError as exc:  # pragma: no cover - depends on the environment
        raise ImportError("lqr needs scipy: pip install 'geomech[sim]'") from exc
    model = model_of(lin_module).require("M_lin", "A_lin", "B_lin")
    lin = getattr(lin_module, "LINEARIZATION", None)
    if not isinstance(lin, dict):
        raise TypeError("lqr: the module has no LINEARIZATION; generate it from linearize")
    x_eq = check_state(lin_module, x_eq, tol=tol, what="x_eq")
    u_eq = model.as_input(u_eq, "u_eq")
    a_eq = np.zeros(model.n_accel) if a_eq is None else np.asarray(a_eq, float).reshape(-1)

    T, bases, labels = _reduction(model, lin, x_eq)
    M = lin_module.M_lin(x_eq, a_eq, u_eq, p)
    A = lin_module.A_lin(x_eq, a_eq, u_eq, p)
    B = lin_module.B_lin(x_eq, a_eq, u_eq, p)
    Mr = T.T @ M @ T
    try:
        Abar = np.linalg.solve(Mr, T.T @ A @ T)
        Bbar = np.linalg.solve(Mr, T.T @ B)
    except np.linalg.LinAlgError as exc:
        raise np.linalg.LinAlgError(
            f"lqr: the reduced M_lin is singular at x_eq (condition number"
            f" {np.linalg.cond(Mr):.3g})"
        ) from exc
    n = T.shape[1]
    Qm = _weight(Q, n, "Q")
    Rm = _weight(R, model.n_inputs, "R")
    P = solve_continuous_are(Abar, Bbar, Qm, Rm)
    K = np.linalg.solve(Rm, Bbar.T @ P)
    return LQRController(K, P, Abar, Bbar, T, x_eq, u_eq, labels, model, bases)


def _reduction(model: Model, lin: dict, x_eq: np.ndarray):
    """``T``, the S² bases and the reduced labels; checks the LINEARIZATION order."""
    nv = len(model.variables)
    states = lin["states"]
    if len(states) != 2 * nv or [r["slice"] for r in lin["rows"]] != [s["slice"] for s in states]:
        raise ValueError(
            "lqr: LINEARIZATION does not have one configuration and one velocity"
            " perturbation per variable with matching row blocks"
        )
    n_full = states[-1]["slice"][1]
    bases = {v.key: tangent_basis(x_eq[v.cfg]) for v in model.s2}
    cols, labels = [], []
    for k, entry in enumerate(states):
        var = model.variables[k % nv]
        a, b = entry["slice"]
        name = _ident(entry["symbol"])
        if var.manifold == "S2":
            block = np.zeros((n_full, 2))
            block[a:b] = bases[var.key]
            labels += [f"{name}[{i}]" for i in range(2)]
        else:
            block = np.zeros((n_full, b - a))
            block[a:b] = np.eye(b - a)
            labels += [name] if b - a == 1 else [f"{name}[{i}]" for i in range(b - a)]
        cols.append(block)
    return np.hstack(cols), bases, labels


def _ident(symbol: str) -> str:
    """A readable identifier for a perturbation symbol (``\\frac{d}{dt}(\\delta{x})`` →
    ``delta_x_dot``)."""
    while (m := re.fullmatch(r"\\frac\{d\}\{dt\}\((.*)\)", symbol)) is not None:
        symbol = m.group(1) + "_dot"
    name = re.sub(r"[^0-9A-Za-z_]+", "_", symbol.replace("\\", "")).strip("_")
    return re.sub(r"_+", "_", name)
