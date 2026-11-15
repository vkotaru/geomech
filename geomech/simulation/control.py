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
S¹          ``atan2(e·(q_eq × q),       ``ω − ω_eq`` (1-D)
            q_eq·q)`` (1-D)
S²          ``Bᵀ (q_eq × q)`` (2-D)     ``Bᵀ (ω − ω_eq)`` (2-D)
SO(3)       ``vee(log(R_eqᵀ R))``       ``Ω − Ω_eq``
==========  ==========================  =========================================

``δq = ξ × q`` with ``ξ ⊥ q``, so ``ξ = q_eq × q`` to first order, and ``δR = R η̂`` so
``R = R_eq exp(η̂)``. On S¹ ``δq = γ e × q`` already has one coordinate, the angle ``γ``
from ``q_eq`` to ``q`` about the axis ``e`` (its tangent basis is ``e × q``, of size 1),
so ``T`` is the identity there. ``B = [b₁ b₂]`` is an orthonormal basis of the tangent plane at
``q_eq`` (``tangent_basis``).

scipy (``scipy.linalg.solve_continuous_are``) is imported only when ``lqr`` is called.
"""

from __future__ import annotations

import re
import warnings
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from numpy.typing import ArrayLike

from geomech.numeric import log_so3
from geomech.simulation.model import (
    STATE_TOL,
    Model,
    check_state,
    model_of,
    solve_accelerations,
)

__all__ = [
    "lqr",
    "LQRController",
    "tangent_basis",
    "reduced_linearization",
    "ReducedLinearization",
]


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
        axes: ``{printed S¹ symbol: unit axis}`` from the ``p`` of ``lqr``.
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
    axes: dict = field(repr=False, default_factory=dict)

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
            elif var.manifold == "S1":  # the angle from q_eq to q about e
                qe, q = xe[var.cfg], x[var.cfg]
                e = self.axes[var.key]
                cfg.append(np.array([np.arctan2(e @ np.cross(qe, q), qe @ q)]))
                vel.append(x[var.vel] - xe[var.vel])
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
    module: Any = None,
    eq_tol: float = 1e-8,
) -> LQRController:
    """The continuous-time LQR gain of a generated linearization about an equilibrium.

    Evaluates ``M_lin``, ``A_lin``, ``B_lin`` of *lin_module* (``generate_module`` of a
    ``LinearizedSystem``) at ``(x_eq, a_eq, u_eq)``, reduces them to the tangent
    coordinates ``δx = Tᵀ δx_full`` (``Mr = Tᵀ M T``, ``Ar = Tᵀ A T``, ``Br = Tᵀ B``),
    forms ``Ā = Mr⁻¹ Ar``, ``B̄ = Mr⁻¹ Br``, solves ``scipy.linalg.solve_continuous_are``
    and returns ``K = R⁻¹ B̄ᵀ P`` in an ``LQRController``.

    An S¹ point needs no reduction: its perturbations ``γ`` and ``δω`` are scalars.
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
        module: the generated standard-form module of the same system
            (``generate_module`` of ``to_standard_form``, with the same layout: the same
            variables, velocities and inputs in the same order, checked). When
            given, ``(x_eq, u_eq)`` is checked to be an equilibrium, see below. When
            omitted, it is not checked and a ``UserWarning`` says so.
        eq_tol: the relative tolerance of the equilibrium check.

    The gain is only meaningful about an equilibrium: the linearization is evaluated at
    ``(x_eq, a_eq, u_eq)`` and the S² reduction assumes ``ω_eq = 0``. At a point that is
    not an equilibrium ``lqr`` still returns a gain with stable closed-loop eigenvalues,
    but the nonlinear loop settles somewhere else. With *module*, ``lqr`` therefore
    requires every velocity of ``x_eq`` to be zero (to *tol*) and the equation residual
    ``r = mass·a_eq + bias + G u_eq`` (S² rows projected on the tangent plane) to satisfy
    ``‖r‖∞ ≤ eq_tol · scale``, with ``scale`` the largest of ``‖bias‖∞``,
    ``‖G u_eq‖∞`` and ``‖A_lin T‖∞`` (the force per unit perturbation, which is nonzero
    at an equilibrium held by a potential even when the bias itself is zero there).
    Otherwise it raises ``ValueError`` with the residual and the accelerations
    ``accelerations(module, x_eq, u_eq, p)``.

    Raises:
        ValueError: bad arguments; a ``LINEARIZATION`` whose variables or inputs are not
            those of its ``LAYOUT`` in the same order; *module* with another layout than
            *lin_module* (the message names both variable and input orders); or
            ``(x_eq, u_eq)`` is not an equilibrium of *module*.
        np.linalg.LinAlgError: the reduced ``M_lin`` is singular.
    """
    try:
        from scipy.linalg import solve_continuous_are
    except ImportError as exc:  # pragma: no cover - depends on the environment
        raise ImportError("lqr needs scipy: pip install 'geomech[sim]'") from exc
    model = model_of(lin_module).require("M_lin", "A_lin", "B_lin")
    lin = getattr(lin_module, "LINEARIZATION", None)
    if not isinstance(lin, dict):
        raise TypeError("lqr: the module has no LINEARIZATION; generate it from linearize")
    x_eq = check_state(lin_module, x_eq, tol=tol, what="x_eq", p=p)
    u_eq = model.as_input(u_eq, "u_eq")
    a_eq = np.zeros(model.n_accel) if a_eq is None else np.asarray(a_eq, float).reshape(-1)
    if module is not None:
        _check_same_layout(model, model_of(module).require("mass", "bias", "input_matrix"))

    red = _reduce(lin_module, model, lin, x_eq, a_eq, u_eq, p, "lqr")
    T, bases, labels, A = red.T, red.bases, red.labels, red.A_full
    Abar, Bbar = red.A, red.B
    n = T.shape[1]
    Qm = _weight(Q, n, "Q")
    Rm = _weight(R, model.n_inputs, "R")
    if module is None:
        warnings.warn(
            "lqr: (x_eq, u_eq) is not checked to be an equilibrium; pass module=<the"
            " standard-form module> to check it (a gain about a non-equilibrium point"
            " stabilises the loop somewhere else)",
            UserWarning,
            stacklevel=2,
        )
    else:
        _check_equilibrium(model, module, x_eq, u_eq, a_eq, p, tol, eq_tol, A @ T)
    P = solve_continuous_are(Abar, Bbar, Qm, Rm)
    K = np.linalg.solve(Rm, Bbar.T @ P)
    axes = {v.key: model.axis(v, p, "lqr") for v in model.s1}
    return LQRController(K, P, Abar, Bbar, T, x_eq, u_eq, labels, model, bases, axes)


@dataclass
class ReducedLinearization:
    """The linearization ``δẋ = A δx + B δu`` in the reduced coordinates of ``lqr``.

    Attributes:
        A, B: ``Mr⁻¹ Tᵀ A_lin T`` and ``Mr⁻¹ Tᵀ B_lin`` with ``Mr = Tᵀ M_lin T``.
        T: ``(n_perturbation, n_reduced)``, block diagonal: ``I`` for the R, R³ and SO(3)
            perturbations, the tangent basis of each S² point for its ``ξ`` and ``δω``.
        labels: the name of each reduced coordinate (``"xi_q[0]"``, ``"eta_R[2]"``, …),
            configuration coordinates first, then velocities.
        bases: ``{printed S² symbol: 3×2 tangent basis}`` at the state.
        M_full, A_full, B_full: ``M_lin``, ``A_lin``, ``B_lin`` evaluated at the state.
    """

    A: np.ndarray
    B: np.ndarray
    T: np.ndarray
    labels: list
    bases: dict
    M_full: np.ndarray = field(repr=False)
    A_full: np.ndarray = field(repr=False)
    B_full: np.ndarray = field(repr=False)


def reduced_linearization(
    lin_module: Any,
    x: ArrayLike,
    u: ArrayLike,
    p: dict,
    *,
    a: ArrayLike | None = None,
    basis=tangent_basis,
    tol: float = STATE_TOL,
) -> ReducedLinearization:
    """The generated linearization at ``(x, a, u)`` in the reduced coordinates of ``lqr``.

    This is the matrix pair ``lqr`` designs with: ``M_lin``, ``A_lin``, ``B_lin`` of
    *lin_module* evaluated at the state *x* (checked like ``simulate``'s ``x0``), the
    accelerations *a* (zeros by default) and the input *u*, reduced with ``T`` (the S²
    perturbations ``ξ`` and ``δω`` in the 2-D basis ``basis(q)``, ``tangent_basis`` by
    default) to ``A = Mr⁻¹ Tᵀ A_lin T`` and ``B = Mr⁻¹ Tᵀ B_lin``, ``Mr = Tᵀ M_lin T``.
    Use it to check controllability or to design another controller with the same
    coordinates as ``LQRController.error``.

    The reduction assumes ``δω ⊥ q``, a constant basis and an S² equation residual of zero,
    which hold at an equilibrium (``ω = 0``). Along a moving reference the reduced
    matrices miss three terms: the time-varying basis (``−Tᵀ M Ṫ``), the normal component
    of the velocity perturbation (``q·δω = −(q × ω)·ξ``) and the row projection of the
    equation's normal residual ``E_r = q (q·E_r)`` (``(q·E_r) Bᵀ q̂ B`` in the ξ columns of
    the S² rows); see ``examples/control/tracking.py`` for the corrected time-varying
    version.

    Raises ``np.linalg.LinAlgError`` when the reduced ``M_lin`` is singular and the
    errors of ``lqr`` for a module without ``LINEARIZATION`` or with a layout mismatch.
    """
    model = model_of(lin_module).require("M_lin", "A_lin", "B_lin")
    lin = getattr(lin_module, "LINEARIZATION", None)
    if not isinstance(lin, dict):
        raise TypeError(
            "reduced_linearization: the module has no LINEARIZATION; generate it from linearize"
        )
    x = check_state(lin_module, x, tol=tol, what="x", p=p)
    u = model.as_input(u, "u")
    a = np.zeros(model.n_accel) if a is None else np.asarray(a, float).reshape(-1)
    return _reduce(lin_module, model, lin, x, a, u, p, "reduced_linearization", basis)


def _reduce(lin_module, model, lin, x, a, u, p, caller, basis=tangent_basis):
    T, bases, labels = _reduction(model, lin, x, basis)
    M = lin_module.M_lin(x, a, u, p)
    A = lin_module.A_lin(x, a, u, p)
    B = lin_module.B_lin(x, a, u, p)
    Mr = T.T @ M @ T
    try:
        Abar = np.linalg.solve(Mr, T.T @ A @ T)
        Bbar = np.linalg.solve(Mr, T.T @ B)
    except np.linalg.LinAlgError as exc:
        raise np.linalg.LinAlgError(
            f"{caller}: the reduced M_lin is singular at x_eq (condition number"
            f" {np.linalg.cond(Mr):.3g})"
        ) from exc
    return ReducedLinearization(Abar, Bbar, T, labels, bases, M, A, B)


def _names(symbols) -> str:
    return "[" + ", ".join(f"'{s}'" for s in symbols) + "]"


def _layout_order(layout: dict) -> tuple[list[str], list[str]]:
    """The configuration variables and the inputs of a ``LAYOUT``, in order."""
    return (
        [c["symbol"] for c in layout["configuration"]],
        [u["symbol"] for u in layout["inputs"]],
    )


def _check_same_layout(lin_model: Model, model: Model) -> None:
    """Raise ``ValueError`` unless the two modules have the same layout (variables,
    velocities, inputs and slices; the parameters may differ)."""
    keys = (
        "n_state",
        "n_config",
        "n_accel",
        "n_inputs",
        "configuration",
        "velocity",
        "acceleration",
        "inputs",
    )
    a, b = lin_model.layout, model.layout
    if all(a.get(k) == b.get(k) for k in keys):
        return
    (lv, lu), (mv, mu) = _layout_order(a), _layout_order(b)
    raise ValueError(
        "lqr: module does not have the layout of lin_module: lin_module has the variables"
        f" {_names(lv)} and inputs {_names(lu)}, module has the variables {_names(mv)} and"
        f" inputs {_names(mu)}; build one"
        " StateLayout with the SystemVariables and the inputs (same order) passed to"
        " to_standard_form and linearize, and generate both modules from it"
    )


def _check_equilibrium(lin_model, module, x_eq, u_eq, a_eq, p, tol, eq_tol, AT) -> None:
    """Raise ``ValueError`` unless ``(x_eq, u_eq)`` is an equilibrium of *module* (see
    ``lqr``)."""
    model = model_of(module).require("mass", "bias", "input_matrix")
    vel = x_eq[model.n_config :]
    if np.max(np.abs(vel), initial=0.0) > tol:
        raise ValueError(
            "lqr: x_eq is not an equilibrium: its velocities must be zero, got"
            f" {np.array2string(vel, precision=6)}"
        )
    mod = module
    bias = mod.bias(x_eq, p)
    Gu = mod.input_matrix(x_eq, p) @ u_eq
    r = mod.mass(x_eq, p) @ a_eq + bias + Gu
    for var in model.s2:
        q = x_eq[var.cfg]
        r[var.acc] -= (q @ r[var.acc]) * q
    scale = max(np.max(np.abs(bias), initial=0.0), np.max(np.abs(Gu), initial=0.0))
    scale = max(scale, np.max(np.abs(AT), initial=0.0))
    res = float(np.max(np.abs(r), initial=0.0))
    if not res <= eq_tol * scale:
        a = solve_accelerations(model, x_eq, u_eq, p, check=True)
        raise ValueError(
            f"lqr: (x_eq, u_eq) is not an equilibrium of module: the equation residual is"
            f" {res:.3g} (relative {res / scale if scale > 0 else np.inf:.3g} > eq_tol ="
            f" {eq_tol:g}) and accelerations(module, x_eq, u_eq, p) ="
            f" {np.array2string(a, precision=6)}, not a_eq"
        )


def _reduction(model: Model, lin: dict, x_eq: np.ndarray, basis=tangent_basis):
    """``T``, the S² bases and the reduced labels; checks the LINEARIZATION order.

    The rows, states and inputs of ``LINEARIZATION`` are matched to the layout by symbol
    (the kinematic row ``kin:<variable>`` of each variable, ``δu`` of each input), not
    only by position."""
    nv = len(model.variables)
    states = lin["states"]
    if len(states) != 2 * nv or [r["slice"] for r in lin["rows"]] != [s["slice"] for s in states]:
        raise ValueError(
            "lqr: LINEARIZATION does not have one configuration and one velocity"
            " perturbation per variable with matching row blocks"
        )
    lay_vars, lay_inputs = _layout_order(model.layout)
    lin_vars = [r["label"][len("kin:") :] for r in lin["rows"][:nv]]
    lin_inputs = [u["symbol"] for u in lin["inputs"]]
    if lin_vars != lay_vars or lin_inputs != ["\\delta{" + u + "}" for u in lay_inputs]:
        raise ValueError(
            f"lqr: the LINEARIZATION of lin_module has the variables {_names(lin_vars)} and"
            f" inputs {_names(lin_inputs)}, but its LAYOUT has the variables"
            f" {_names(lay_vars)} and inputs {_names(lay_inputs)}; regenerate it with a"
            " StateLayout built from the SystemVariables and the inputs (same order) passed"
            " to linearize"
        )
    n_full = states[-1]["slice"][1]
    bases = {v.key: np.asarray(basis(x_eq[v.cfg]), dtype=float) for v in model.s2}
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
