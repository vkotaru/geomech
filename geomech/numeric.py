"""Numeric evaluation of expression trees and finite-difference checks of derived equations.

Three layers, each built on the previous one:

1. ``hat``, ``vee``, ``rand_rot``, ``expm_so3``, ``log_so3``: the SO(3) primitives on
   numpy arrays.
2. ``Evaluator`` / ``evaluate``: evaluate any expression tree given values for its leaves.
   Leaves are keyed by their printed name, ``str(leaf)``; velocities and accelerations by
   ``str(TimeDerivative(leaf))`` and so on. Unlisted leaves are drawn at random.
3. Two finite-difference oracles that need no closed-form reference:

   - ``check_eom`` / ``assert_eom`` compare the output of ``compute_eom`` with the
     Euler-Lagrange residual of the *numeric* Lagrangian at random trajectory points.
   - ``finite_difference_check`` compares the output of ``linearize`` with a central
     difference of the equations along a constraint-respecting perturbation.

All vectors are column arrays of shape ``(3, 1)``, matrices are ``(3, 3)`` and scalars are
Python floats; the conventions follow Lee et al. (``S2.t_diff()`` is ``ω × q`` and
``SO3.t_diff()`` is ``R Ω̂``).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Sequence

import numpy as np

from geomech.core.base import S2, SO3, Expr, Matrix, Scalar, Vector
from geomech.core.operations import (
    Add,
    Cross,
    Dot,
    Hat,
    MAdd,
    MMMul,
    Mul,
    MVMul,
    SMMul,
    SVMul,
    TimeDerivative,
    Transpose,
    VAdd,
    Variation,
    Vee,
    VVMul,
)

if TYPE_CHECKING:
    from geomech.dynamics import LinearizedSystem, SystemVariables

__all__ = [
    "hat",
    "vee",
    "rand_rot",
    "expm_so3",
    "log_so3",
    "Evaluator",
    "evaluate",
    "sample_state",
    "sample_direction",
    "perturbed_state",
    "linearized_rhs",
    "finite_difference_check",
    "check_eom",
    "assert_eom",
    "TermResolution",
    "DEFAULT_H",
]

# ---------------------------------------------------------------------------
# SO(3) primitives
# ---------------------------------------------------------------------------


def hat(v: np.ndarray) -> np.ndarray:
    """Skew-symmetric matrix of ``v`` so that ``hat(v) @ w == v × w``."""
    v = np.asarray(v, dtype=float).reshape(3)
    return np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]], dtype=float)


def vee(M: np.ndarray) -> np.ndarray:
    """Inverse of ``hat``: the ``(3, 1)`` vector of a skew-symmetric matrix."""
    return np.array([[M[2, 1]], [M[0, 2]], [M[1, 0]]])


def rand_rot(rng: np.random.Generator) -> np.ndarray:
    """Random rotation matrix (QR of a Gaussian matrix, determinant forced to +1)."""
    A = rng.normal(size=(3, 3))
    Q, _ = np.linalg.qr(A)
    if np.linalg.det(Q) < 0:
        Q[:, 0] = -Q[:, 0]
    return Q


def expm_so3(v):
    """Rodrigues formula: ``exp(hat(v))`` for ``v`` in R³."""
    v = np.asarray(v, dtype=float).reshape(3)
    theta = np.linalg.norm(v)
    K = hat(v)
    if theta < 1e-12:
        return np.eye(3) + K + 0.5 * K @ K
    return np.eye(3) + np.sin(theta) / theta * K + (1 - np.cos(theta)) / theta**2 * (K @ K)


def log_so3(R) -> np.ndarray:
    """Principal logarithm of a rotation: the ``(3,)`` vector ``v`` with ``|v| ≤ π`` and
    ``expm_so3(v) == R``.

    The angle is ``atan2(|vee(R − Rᵀ)|/2, (tr R − 1)/2)``, accurate over the whole range.
    Near the identity the axis comes from the skew part with the series of ``θ/sin θ``;
    near ``θ = π`` (where the skew part vanishes) it comes from the symmetric part
    ``(R + Rᵀ)/2 = cos θ I + (1 − cos θ) n nᵀ``, with the sign of ``n`` taken from the skew
    part. At exactly ``θ = π`` both ``±π n`` are logarithms; ``+π n`` with the largest
    component of ``n`` positive is returned. ``R`` must be a rotation (not checked).
    """
    R = np.asarray(R, dtype=float).reshape(3, 3)
    w = 0.5 * np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]])
    s = np.linalg.norm(w)  # sin θ
    c = 0.5 * (np.trace(R) - 1.0)  # cos θ
    theta = np.arctan2(s, c)
    if theta < 1e-4:
        return w * (1.0 + theta**2 / 6.0 + 7.0 * theta**4 / 360.0)
    if c > -0.9:
        return w * (theta / s)
    S = 0.5 * (R + R.T) - c * np.eye(3)  # (1 − cos θ) n nᵀ
    k = int(np.argmax(np.diag(S)))
    n = S[:, k] / np.sqrt(S[k, k] * (1.0 - c))
    if n @ w < 0:
        n = -n
    return theta * n


def _cross(a, b):
    return np.cross(np.asarray(a).reshape(3), np.asarray(b).reshape(3)).reshape(3, 1)


def _perp(v, q, rng):
    """A random vector orthogonal to the unit vector *q* (or *v* projected if given)."""
    if v is None:
        v = rng.normal(size=(3, 1))
    q = q.reshape(3, 1)
    return v - q * float((q.T @ v)[0, 0])


# ---------------------------------------------------------------------------
# Evaluator
# ---------------------------------------------------------------------------


class Evaluator:
    """Evaluate an expression tree numerically.

    Leaves are keyed by ``str(node)``: a value found in ``env`` is used, otherwise a random
    one is drawn (and stored in ``env`` so that repeated evaluations are consistent).
    ``TimeDerivative`` chains over a leaf are independent values keyed by their printed
    name; a ``Variation`` must have a value in ``env`` (a perturbation is never random).

    Random draws respect the leaf's declared structure, because the engine does: a
    ``Scalar`` with a ``value`` takes that value; ``SO3`` and ``OrthogonalMatrix`` leaves are
    rotations; ``S2`` and ``UnitNorm`` vectors have unit norm; ``SymmetricMatrix`` and
    ``SkewSymmetricMatrix`` leaves are (skew-)symmetric; the zero and identity constants
    are resolved by type.
    """

    def __init__(self, env=None, rng=None):
        self.env = {} if env is None else env
        self.rng = np.random.default_rng(7) if rng is None else rng

    def leaf(self, node):
        # Zero / identity constants are resolved by type, not by name: the scalar,
        # vector and matrix zeros all print as "0".
        if node.is_zero:
            if isinstance(node, Matrix):
                return np.zeros((3, 3))
            if isinstance(node, Vector):
                return np.zeros((3, 1))
            return 0.0
        if isinstance(node, Matrix) and node.flags.is_identity:
            return np.eye(3)
        key = str(node)
        if key in self.env:
            return self.env[key]
        rng = self.rng
        if isinstance(node, Scalar):
            val = float(node.value) if node.value is not None else float(rng.normal())
        elif isinstance(node, SO3) or (isinstance(node, Matrix) and node.flags.is_orthogonal):
            val = rand_rot(rng)
        elif isinstance(node, Matrix):
            A = rng.normal(size=(3, 3))
            if node.is_symmetric:
                A = A + A.T
            elif node.flags.is_skew_symmetric:
                A = A - A.T
            val = A
        elif isinstance(node, Vector):  # S2 leaves carry the UnitNorm flag
            v = rng.normal(size=(3, 1))
            val = v / np.linalg.norm(v) if node.is_unit_norm else v
        else:
            raise TypeError(f"unknown leaf {type(node)} {key}")
        self.env[key] = val
        return val

    def _derived_leaf(self, node):
        """TimeDerivative chains over a leaf (or over a Variation of a leaf)."""
        key = str(node)
        if key in self.env:
            return self.env[key]
        inner = node
        while isinstance(inner, TimeDerivative):
            inner = inner.expr
        if isinstance(inner, Variation):
            raise KeyError(f"no value for the perturbation symbol '{key}' in env")
        if isinstance(inner, Scalar):
            self.env[key] = float(self.rng.normal())
        elif isinstance(inner, Vector):
            self.env[key] = self.rng.normal(size=(3, 1))
        elif isinstance(inner, Matrix):
            self.env[key] = self.rng.normal(size=(3, 3))
        else:
            raise ValueError(f"TimeDerivative of non-leaf survived: {node}")
        return self.env[key]

    def __call__(self, node):
        if isinstance(node, (Add, VAdd, MAdd)):
            return sum(self(n) for n in node.nodes)
        if isinstance(node, (Mul, SVMul, SMMul)):
            return self(node.left) * self(node.right)
        if isinstance(node, (MVMul, MMMul, VVMul)):
            out = self(node.left) @ self(node.right)
            if out.shape == (1, 1):
                return float(out[0, 0])
            return out
        if isinstance(node, Dot):
            return float((self(node.left).T @ self(node.right))[0, 0])
        if isinstance(node, Cross):
            return _cross(self(node.left), self(node.right))
        if isinstance(node, Hat):
            return hat(self(node.expr))
        if isinstance(node, Vee):
            return vee(self(node.expr))
        if isinstance(node, Transpose):
            return self(node.expr).T
        if isinstance(node, TimeDerivative):
            return self._derived_leaf(node)
        if isinstance(node, Variation):
            key = str(node)
            if key not in self.env:
                raise KeyError(f"no value for the perturbation symbol '{key}' in env")
            return self.env[key]
        return self.leaf(node)


def evaluate(expr: Expr, env: dict | None = None, *, seed: int = 0) -> Any:
    """Evaluate ``expr`` with the leaf values in ``env``; draw the rest from ``seed``.

    ``env`` maps the *printed name* of a leaf, ``str(leaf)``, to its value: a float for a
    ``Scalar``, a ``(3, 1)`` array for a ``Vector``/``S2`` and a ``(3, 3)`` array for a
    ``Matrix``/``SO3``. Velocities and accelerations are leaves of their own, keyed by
    ``str(TimeDerivative(x))``, ``str(TimeDerivative(TimeDerivative(x)))``, and the
    tangent vectors of the manifold leaves by ``str(q.get_tangent_vector())``. Every leaf
    that is not in ``env`` is drawn once from ``numpy.random.default_rng(seed)``, so two
    calls with the same ``env`` and ``seed`` return the same value; ``env`` itself is not
    modified. A ``Variation`` symbol has no natural random value and must be in ``env``.

    Example:
        >>> x = Vector("x"); e3 = Vector("e3", attr=["Constant"])
        >>> v = x.t_diff()
        >>> env = {str(x): np.array([[1.0], [0.0], [0.0]]),
        ...        str(v): np.array([[0.0], [2.0], [0.0]]),
        ...        "e3": np.array([[0.0], [0.0], [1.0]])}
        >>> evaluate(Dot(Cross(x, v), e3), env)
        2.0

    Returns a float for a scalar expression, a ``(3, 1)`` array for a vector expression
    and a ``(3, 3)`` array for a matrix expression.
    """
    ev = Evaluator(env=dict(env or {}), rng=np.random.default_rng(seed))
    return ev(expr)


# ---------------------------------------------------------------------------
# Trajectory points and tangent perturbations (linearization check)
# ---------------------------------------------------------------------------


def _kinematic_symbols(variables):
    """(config var, [position key, velocity key, acceleration key]) per configuration variable."""
    out = []
    for s in variables.scalars:
        out.append((s, [str(s), str(s.t_diff()), str(s.t_diff().t_diff())]))
    for v in variables.vectors:
        if isinstance(v, S2):
            w = v.get_tangent_vector()
            out.append((v, [str(v), str(w), str(w.t_diff())]))
        else:
            out.append((v, [str(v), str(v.t_diff()), str(v.t_diff().t_diff())]))
    for m in variables.matrices:
        W = m.get_tangent_vector()
        out.append((m, [str(m), str(W), str(W.t_diff())]))
    return out


def sample_state(
    variables: SystemVariables, inputs: Sequence[Expr], rng: np.random.Generator
) -> dict[str, Any]:
    """Random trajectory point on the manifold: |q| = 1, q·ω = 0, q·ω̇ = 0, R ∈ SO(3).

    Returns an ``env`` for ``Evaluator`` holding the configuration, velocity and
    acceleration of every configuration variable and a value for every input.
    """
    env: dict[str, Any] = {}
    for var, (kq, kv, ka) in _kinematic_symbols(variables):
        if isinstance(var, Scalar):
            env[kq], env[kv], env[ka] = (float(rng.normal()) for _ in range(3))
        elif isinstance(var, S2):
            q = rng.normal(size=(3, 1))
            q /= np.linalg.norm(q)
            env[kq] = q
            env[kv] = _perp(None, q, rng)
            env[ka] = _perp(None, q, rng)
        elif isinstance(var, SO3):
            env[kq] = rand_rot(rng)
            env[kv] = rng.normal(size=(3, 1))
            env[ka] = rng.normal(size=(3, 1))
        else:
            env[kq], env[kv], env[ka] = (rng.normal(size=(3, 1)) for _ in range(3))
    for u in inputs:
        env[str(u)] = float(rng.normal()) if isinstance(u, Scalar) else rng.normal(size=(3, 1))
    return env


def sample_direction(
    variables: SystemVariables, inputs: Sequence[Expr], env: dict, rng: np.random.Generator
) -> tuple[dict, dict]:
    """Random perturbation direction tangent to the constraints.

    Returns ``(raw, pert)``. ``raw`` holds the generators used by ``perturbed_state``
    (ξ, w, a for S²; η, δΩ, δΩ̇ for SO(3); δx, δẋ, δẍ for R³; δs, δṡ, δs̈ for R; δu).
    ``pert`` maps the printed perturbation symbols of ``linearize`` (states, velocity
    dstates and input variations) to their numeric values; for S² the induced
    velocity direction is δω = ξ × ω + w and δω̇ = ξ × ω̇ + a.
    """
    raw: dict[str, Any] = {}
    pert: dict[str, Any] = {}
    for var, (kq, kv, ka) in _kinematic_symbols(variables):
        if isinstance(var, Scalar):
            ds, dv, da = (float(rng.normal()) for _ in range(3))
            raw[kq] = (ds, dv, da)
            d = var.delta()
            pert[str(d)], pert[str(d.t_diff())], pert[str(d.t_diff().t_diff())] = ds, dv, da
        elif isinstance(var, S2):
            q = env[kq]
            xi, w, a = (_perp(None, q, rng) for _ in range(3))
            raw[kq] = (xi, w, a)
            omega = var.get_tangent_vector()
            pert[str(var.get_variation_vector())] = xi
            pert[str(Variation(omega))] = _cross(xi, env[kv]) + w
            pert[str(Variation(omega).t_diff())] = _cross(xi, env[ka]) + a
        elif isinstance(var, SO3):
            eta, dW, ddW = (rng.normal(size=(3, 1)) for _ in range(3))
            raw[kq] = (eta, dW, ddW)
            Om = var.get_tangent_vector()
            pert[str(var.get_variation_vector())] = eta
            pert[str(Variation(Om))] = dW
            pert[str(Variation(Om).t_diff())] = ddW
        else:
            dx, ddx, dddx = (rng.normal(size=(3, 1)) for _ in range(3))
            raw[kq] = (dx, ddx, dddx)
            d = var.delta()
            pert[str(d)], pert[str(d.t_diff())], pert[str(d.t_diff().t_diff())] = dx, ddx, dddx
    for u in inputs:
        du = float(rng.normal()) if isinstance(u, Scalar) else rng.normal(size=(3, 1))
        raw[str(u)] = du
        pert[str(u.delta())] = du
    return raw, pert


def perturbed_state(
    variables: SystemVariables, inputs: Sequence[Expr], env: dict, raw: dict, eps: float
) -> dict:
    """Copy of *env* moved by *eps* along *raw* while staying on the manifold.

    R, R³: linear; SO(3): R exp(ε η̂), Ω + ε δΩ, Ω̇ + ε δΩ̇;
    S²: q(ε) = exp(ε ξ̂) q, ω(ε) = exp(ε ξ̂)(ω + ε w), ω̇(ε) = exp(ε ξ̂)(ω̇ + ε a).
    """
    out = dict(env)
    for var, (kq, kv, ka) in _kinematic_symbols(variables):
        g0, g1, g2 = raw[kq]
        if isinstance(var, S2):
            E = expm_so3(eps * g0)
            out[kq] = E @ env[kq]
            out[kv] = E @ (env[kv] + eps * g1)
            out[ka] = E @ (env[ka] + eps * g2)
        elif isinstance(var, SO3):
            out[kq] = env[kq] @ expm_so3(eps * g0)
            out[kv] = env[kv] + eps * g1
            out[ka] = env[ka] + eps * g2
        else:
            out[kq] = env[kq] + eps * g0
            out[kv] = env[kv] + eps * g1
            out[ka] = env[ka] + eps * g2
    for u in inputs:
        out[str(u)] = env[str(u)] + eps * raw[str(u)]
    return out


def _dstate_values(lin, variables, env, pert):
    """Numeric δẋ: velocity dstates from *pert*, configuration dstates from the kinematics."""
    vals = {}
    n_cfg = len(lin.states) // 2
    for i, (cfg, dcfg) in enumerate(zip(lin.states[:n_cfg], lin.dstates[:n_cfg])):
        vel = lin.states[n_cfg + i]
        var = (variables.scalars + variables.vectors + variables.matrices)[i]
        if isinstance(var, S2):
            q, w = env[str(var)], env[str(var.get_tangent_vector())]
            P = q @ q.T
            vals[str(dcfg)] = P @ hat(w) @ pert[str(cfg)] + (np.eye(3) - P) @ pert[str(vel)]
        elif isinstance(var, SO3):
            W = env[str(var.get_tangent_vector())]
            vals[str(dcfg)] = -hat(W) @ pert[str(cfg)] + pert[str(vel)]
        else:
            vals[str(dcfg)] = pert[str(vel)]
    for d in lin.dstates[n_cfg:]:
        vals[str(d)] = pert[str(d)]
    return vals


def _as_array(x):
    return np.atleast_2d(np.asarray(x, dtype=float))


def linearized_rhs(
    lin: LinearizedSystem, row: str, ev: Evaluator, pert: dict, dvals: dict
) -> np.ndarray:
    """Numeric −M δv̇ + A δx + B δu for one row, blocks evaluated with *ev*."""
    total = 0.0
    for col, block in lin.M[row].items():
        total = total - _as_array(ev(block)) @ _as_array(dvals[col]).reshape(-1, 1)
    for col, block in lin.A[row].items():
        total = total + _as_array(ev(block)) @ _as_array(pert[col]).reshape(-1, 1)
    for col, block in lin.B[row].items():
        total = total + _as_array(ev(block)) @ _as_array(pert[col]).reshape(-1, 1)
    return np.asarray(total, dtype=float).reshape(-1)


def finite_difference_check(
    lin, eom, variables, inputs, trials=5, eps=1e-5, rtol=1e-6, seed=0
) -> float:
    """Compare (E(ε) − E(−ε)) / 2ε with the linearization for every dynamic row.

    The error of each row is measured relative to ``max(1, |rhs|)`` and asserted to be
    below *rtol* at *trials* random base states and directions. Returns the worst error.
    """
    rng = np.random.default_rng(seed)
    worst = 0.0
    for _ in range(trials):
        base = sample_state(variables, inputs, rng)
        ev = Evaluator(env=base, rng=rng)
        for _, (_, eq) in eom.items():
            ev(eq)  # populate every constant leaf so all evaluations share them
        raw, pert = sample_direction(variables, inputs, ev.env, rng)
        env_p = perturbed_state(variables, inputs, ev.env, raw, eps)
        env_m = perturbed_state(variables, inputs, ev.env, raw, -eps)
        ev_p, ev_m = Evaluator(env=env_p, rng=rng), Evaluator(env=env_m, rng=rng)
        dvals = _dstate_values(lin, variables, ev.env, pert)
        for row, (_, eq) in eom.items():
            fd = (_as_array(ev_p(eq)) - _as_array(ev_m(eq))).reshape(-1) / (2 * eps)
            rhs = linearized_rhs(lin, row, ev, pert, dvals)
            scale = max(1.0, float(np.max(np.abs(rhs))))
            err = float(np.max(np.abs(fd - rhs))) / scale
            worst = max(worst, err)
            assert err < rtol, f"row {row}: finite difference {fd} vs linearization {rhs}"
    return worst


# ---------------------------------------------------------------------------
# Euler-Lagrange oracle for compute_eom
# ---------------------------------------------------------------------------


def _dkey(sym):
    return str(TimeDerivative(sym))


def _is_expr(obj):
    return hasattr(obj, "type")


def _config_vars(variables):
    return list(variables.scalars) + list(variables.vectors) + list(variables.matrices)


def _variation_key(var):
    """Printed name of the variation vector that keys the equation of *var* in ``eom``."""
    if isinstance(var, (S2, SO3)):
        return str(var.get_variation_vector())
    return str(var.delta())


def _zero_variations(variables):
    return {
        _variation_key(v): 0.0 if isinstance(v, Scalar) else np.zeros((3, 1)) for v in variables
    }


def _random_like(leaf, rng):
    """A random value with the shape of *leaf* (symmetric if the leaf is declared so)."""
    if isinstance(leaf, Scalar):
        return float(rng.normal())
    if isinstance(leaf, Matrix):
        A = rng.normal(size=(3, 3))
        return A + A.T if leaf.is_symmetric else A
    return rng.normal(size=(3, 1))


class _Coordinates:
    """Coordinates of the numeric Lagrangian at a random trajectory point.

    ``state`` maps ``("q", var)`` and ``("v", var)`` to the coordinates the finite
    differences act on: ``(s, ṡ)`` on R, ``(x, ẋ)`` on R³, ``(q, q̇)`` on S² with ``|q| = 1``
    and ``q · q̇ = 0``, ``(R, Ω)`` on SO(3), and ``(p, ṗ)`` for a prescribed parameter.
    ``rate`` maps the same keys to their time derivative at the point, so that
    ``d/dt f = (f(state + h·rate) − f(state − h·rate)) / 2h``. ``accel`` holds the
    acceleration leaves the library equation refers to (``ω̇``, ``Ω̇``, ``s̈``, ``ẍ``, ``p̈``).
    """

    def __init__(self, variables, params, rng):
        self.vars = _config_vars(variables)
        self.params = list(params)
        self.state, self.rate, self.accel = {}, {}, {}
        for var in self.vars + self.params:
            self._sample(var, rng)

    def _sample(self, var, rng):
        kq, kv = ("q", var), ("v", var)
        if isinstance(var, S2):
            q = rng.normal(size=(3, 1))
            q /= np.linalg.norm(q)
            w, wd = _perp(None, q, rng), _perp(None, q, rng)
            self.state[kq], self.state[kv] = q, _cross(w, q)
            self.rate[kq], self.rate[kv] = _cross(w, q), _cross(wd, q) + _cross(w, _cross(w, q))
            self.accel[_dkey(var.get_tangent_vector())] = wd
        elif isinstance(var, SO3):
            R, Om, Omd = rand_rot(rng), rng.normal(size=(3, 1)), rng.normal(size=(3, 1))
            self.state[kq], self.state[kv] = R, Om
            self.rate[kq], self.rate[kv] = R @ hat(Om), Omd
            self.accel[_dkey(var.get_tangent_vector())] = Omd
        else:
            x, xd, xdd = (_random_like(var, rng) for _ in range(3))
            self.state[kq], self.state[kv] = x, xd
            self.rate[kq], self.rate[kv] = xd, xdd
            self.accel[_dkey(var.t_diff())] = xdd

    def env(self, state):
        """``Evaluator`` env (printed-name keys) of the configuration and velocity leaves."""
        env = {}
        for var in self.vars + self.params:
            q, v = state[("q", var)], state[("v", var)]
            if isinstance(var, S2):
                env[str(var)], env[str(var.get_tangent_vector())] = q, _cross(q, v)
            elif isinstance(var, SO3):
                env[str(var)], env[str(var.get_tangent_vector())] = q, v
            else:
                env[str(var)], env[_dkey(var)] = q, v
        return env

    def full_env(self):
        """Env at the base point including the accelerations (for evaluating an equation)."""
        return {**self.env(self.state), **self.accel}

    def along(self, h):
        """The state moved a time ``h`` along the trajectory (first order)."""
        return {k: v + h * self.rate[k] for k, v in self.state.items()}


def _shift(state, key, idx, h):
    """Copy of *state* with entry *idx* of ``state[key]`` moved by *h*."""
    out = dict(state)
    val = state[key]
    if np.isscalar(val):
        out[key] = val + h
    else:
        val = np.array(val, dtype=float, copy=True)
        val.flat[idx] += h
        out[key] = val
    return out


def _grad(fun, state, key, h):
    """Central-difference gradient of the scalar *fun* with respect to ``state[key]``."""
    val = state[key]
    if np.isscalar(val):
        return (fun(_shift(state, key, 0, h)) - fun(_shift(state, key, 0, -h))) / (2 * h)
    g = np.zeros_like(np.asarray(val, dtype=float))
    for i in range(g.size):
        g.flat[i] = (fun(_shift(state, key, i, h)) - fun(_shift(state, key, i, -h))) / (2 * h)
    return g


def _ddt(fun, co, h):
    """Central-difference time derivative of *fun* along the trajectory through *co*."""
    return (fun(co.along(h)) - fun(co.along(-h))) / (2 * h)


# Default finite-difference step of ``check_eom`` / ``assert_eom``.  ``d/dt ∂L/∂v`` nests
# two central differences with the same step, so its roundoff is about ε|L|/h² with
# ε ≈ 1e-16 (not the ε/h of a single difference) while its truncation is about h².  At
# h = 1e-5 the roundoff term is ≈ 1e-6 for |L| ≈ 1, within a factor 10 of the 1e-5 gate
# and large enough for a correct equation to fail it at some trajectory points; at
# h = 1e-4 roundoff ≈ 1e-8 and truncation ≈ 1e-8, the minimum of the two for this scheme.
DEFAULT_H = 1e-4


def _constant_env(lagrangian, inf_work, inputs, eom, co, rng):
    """Every leaf value shared by all evaluations at this point.

    The trajectory point, the inputs and every constant parameter are drawn once so that
    the numeric Lagrangian and the library equation see the same numbers.
    """
    zeros = _zero_variations(co.vars)
    ev = Evaluator(env={**co.full_env(), **zeros}, rng=rng)
    for expr in [lagrangian, *inputs, *(eq for _, eq in eom.values())]:
        ev(expr)
    if _is_expr(inf_work):
        ev(inf_work)
    return {k: v for k, v in ev.env.items() if k not in zeros}


def _forcing(inf_work, consts, co, var, rng):
    """Generalized force on *var*: ``inf_work`` with its variation a basis vector, others 0."""
    if not _is_expr(inf_work):
        return 0.0
    key = _variation_key(var)
    basis = [1.0] if isinstance(var, Scalar) else [np.eye(3)[:, [k]] for k in range(3)]
    f = []
    for b in basis:
        env = {**consts, **_zero_variations(co.vars), key: b}
        f.append(float(Evaluator(env=env, rng=rng)(inf_work)))
    return f[0] if isinstance(var, Scalar) else np.array(f).reshape(3, 1)


def _euler_lagrange(lnum, co, var, h):
    """``(-g, P)``: minus the unforced Euler-Lagrange residual and the comparison projector.

    R, R³: ``g = d/dt ∂L/∂v − ∂L/∂q``;
    S²:    ``g`` of the constrained Lagrangian ``L̃(q, q̇) = L(q, q × q̇)``, mapped to the
           coefficient of ``ξ`` by ``−q × g`` and compared on the tangent plane;
    SO(3): ``g = d/dt ∂L/∂Ω − (∂L/∂Ω) × Ω − vee(A − Aᵀ)`` with ``A = Rᵀ ∂L/∂R``.
    """
    kq, kv = ("q", var), ("v", var)
    dL_dv = lambda state: _grad(lnum, state, kv, h)
    if isinstance(var, SO3):
        R, Om = co.state[kq], co.state[kv]
        A = R.T @ _grad(lnum, co.state, kq, h)
        g = _ddt(dL_dv, co, h) - _cross(dL_dv(co.state), Om) - vee(A - A.T)
        return -g, np.eye(3)
    g = _ddt(dL_dv, co, h) - _grad(lnum, co.state, kq, h)
    if isinstance(var, S2):
        q = co.state[kq]
        return -_cross(q, g), np.eye(3) - q @ q.T
    return -g, (1.0 if isinstance(var, Scalar) else np.eye(3))


def _scale(e_true, proj):
    """``max(1, |P E_true|)``, the normalisation of the relative error."""
    if np.isscalar(proj):
        return max(1.0, abs(float(e_true)))
    return max(1.0, float(np.max(np.abs(proj @ e_true))))


def _tangent_norm(val, proj):
    """``(|P t|, |t|)`` for a term value *val* (floats for a scalar equation)."""
    if np.isscalar(proj):
        return abs(float(val)), abs(float(val))
    t = np.asarray(val, dtype=float).reshape(3, 1)
    return float(np.max(np.abs(proj @ t))), float(np.max(np.abs(t)))


def _relative_error(e_lib, e_true, proj):
    """``max |P (E_lib − E_true)| / max(1, |P E_true|)``."""
    if np.isscalar(proj):
        return abs(float(e_lib) - float(e_true)) / _scale(e_true, proj)
    d = proj @ (np.asarray(e_lib, dtype=float).reshape(3, 1) - e_true)
    return float(np.max(np.abs(d))) / _scale(e_true, proj)


# A term whose tangent part is below this fraction of its own size or of the equation's
# scale is parallel to q, or zero on the constraint set (``q × (ω × (e₃ × q))`` with
# ω ⊥ q), to within roundoff: it has no content on the tangent plane and nothing to
# resolve. Such a term is also far below the finite-difference noise (≈ 1e-8 · scale).
_ZERO_TERM_RTOL = 1e-10


def _terms(expr):
    """The additive terms of an equation (the nodes of its top-level sum)."""
    return list(expr.nodes) if isinstance(expr, (Add, VAdd)) else [expr]


class TermResolution:
    """How well one additive term of an equation is resolved by ``check_eom``.

    ``ratio`` is the best value over the trials of ``2 |P t| / max(1, |P E_true|)``: the
    relative error the check would report if the term had the wrong sign. The term is
    *resolved* at tolerance ``rtol`` when ``ratio >= rtol``; a term with ``ratio < rtol``
    could be wrong (sign, missing factor) without the check noticing. ``magnitude``
    (``|P t|``) and ``scale`` (``max(1, |P E_true|)``, so that the resolution at ``rtol`` is
    ``rtol · scale``) are the values at that best trial. ``tangent`` is ``False`` when the
    term is parallel to ``q`` or zero (to within roundoff, ``1e-10`` of the equation's
    scale) at every trial, so that there is nothing to resolve.
    """

    __slots__ = ("term", "ratio", "magnitude", "scale", "tangent")

    def __init__(self, term):
        self.term = term
        self.ratio = 0.0
        self.magnitude = 0.0
        self.scale = 1.0
        self.tangent = False

    def update(self, tangent_norm, full_norm, scale):
        if tangent_norm <= _ZERO_TERM_RTOL * max(full_norm, scale):
            return
        self.tangent = True
        ratio = 2.0 * tangent_norm / scale
        if ratio > self.ratio:
            self.ratio, self.magnitude, self.scale = ratio, tangent_norm, scale

    def unresolved(self, rtol):
        return self.tangent and self.ratio < rtol

    def __repr__(self):
        return (
            f"TermResolution({self.term!s}, ratio={self.ratio:.2e}, "
            f"magnitude={self.magnitude:.2e}, scale={self.scale:.2e}, tangent={self.tangent})"
        )


def _eom_errors(lagrangian, inf_work, variables, eom, inputs, params, trials, seed, h):
    """``(worst, terms)`` over *trials* random trajectory points.

    ``worst[key]`` is the largest relative error of the equation keyed by *key*;
    ``terms[key]`` lists a ``TermResolution`` for each additive term of that equation.
    """
    rng = np.random.default_rng(seed)
    worst, terms = {}, {}
    for _ in range(trials):
        co = _Coordinates(variables, params, rng)
        consts = _constant_env(lagrangian, inf_work, inputs, eom, co, rng)

        def lnum(state):
            return float(Evaluator(env={**consts, **co.env(state)}, rng=rng)(lagrangian))

        for var in co.vars:
            key = _variation_key(var)
            if key not in eom:
                raise ValueError(f"eom has no equation for the variation '{key}' of {var}")
            e_true, proj = _euler_lagrange(lnum, co, var, h)
            e_true = e_true + _forcing(inf_work, consts, co, var, rng)
            ev = Evaluator(env=dict(consts), rng=rng)
            e_lib = ev(eom[key][1])
            err = _relative_error(e_lib, e_true, proj)
            worst[key] = max(worst.get(key, 0.0), err)
            scale = _scale(e_true, proj)
            if key not in terms:
                terms[key] = [TermResolution(t) for t in _terms(eom[key][1])]
            for res, t in zip(terms[key], _terms(eom[key][1])):
                res.update(*_tangent_norm(ev(t), proj), scale)
    return worst, terms


def check_eom(
    lagrangian: Expr,
    inf_work: Expr | None,
    variables: SystemVariables,
    eom: dict[str, tuple[Expr, Expr]],
    *,
    inputs: Sequence[Expr] = (),
    params: Sequence[Expr] = (),
    trials: int = 5,
    seed: int = 0,
    h: float = DEFAULT_H,
) -> float:
    """Worst relative error of ``eom`` against a finite-difference Euler-Lagrange residual.

    ``compute_eom`` returns, for each configuration variable, the coefficient ``E`` of its
    variation vector in ``δ∫L dt + ∫δW dt`` after integration by parts, so that the
    equation of motion is ``E = 0``. By this convention ``E = −(Euler-Lagrange) + forcing``:
    for a point mass ``E = −m ẍ + f``. This function evaluates the *numeric* Lagrangian with
    ``Evaluator`` and forms the Euler-Lagrange residual ``g`` by central differences, with no
    use of the library's variational rules:

    - R, R³: ``g = d/dt ∂L/∂v − ∂L/∂q``, ``E_true = −g + f``;
    - S²: constrained Lagrangian ``L̃(q, q̇) = L(q, q × q̇)`` on ``|q| = 1`` (``q̇ = ω × q``),
      ``E_true = −q × g + f``, compared on the tangent plane only (Lee et al., Ch. 5);
    - SO(3): body-frame Euler-Poincaré form ``g = d/dt ∂L/∂Ω − (∂L/∂Ω) × Ω − vee(A − Aᵀ)``
      with ``A = Rᵀ ∂L/∂R``, ``E_true = −g + f`` (Lee et al., Ch. 6).

    The forcing ``f`` is ``inf_work`` evaluated with the variable's variation set to each
    basis vector and the other variations zero. ``d/dt`` is a directional derivative along
    the sampled velocity and acceleration, so no integration is performed: the library's
    ``E(q, v, v̇, u)`` and ``E_true`` must agree at random accelerations.

    Args:
        lagrangian: the Lagrangian passed to ``compute_eom``.
        inf_work: the virtual work passed to ``compute_eom`` (``None`` or a non-expression
            means no forcing).
        variables: the ``SystemVariables`` passed to ``compute_eom``.
        eom: the result of ``compute_eom``, ``{str(variation): (variation, E)}``.
        inputs: input leaves (``Scalar``/``Vector``); drawn once per trial like any
            non-constant leaf, listed for symmetry with ``linearize``.
        params: prescribed time-varying symbols that are neither constant nor
            configuration variables (for example a cable length ``l(t)``). They are given a
            random value, rate and acceleration; a non-constant parameter *not* listed here
            is treated as frozen, which makes a correct equation look wrong.
        trials: number of random trajectory points.
        seed: seed of the random draws.
        h: finite-difference step. ``d/dt ∂L/∂v`` nests two central differences, so the
            noise is about ``h² + 1e-16 |L| / h²``: ≈ 1e-8 at the default 1e-4 for a
            Lagrangian of order one, ≈ 1e-6 at 1e-5 (see ``DEFAULT_H``).

    Returns:
        The largest error over trials and variables, each measured as
        ``max |P (E_lib − E_true)| / max(1, |P E_true|)`` with ``P`` the tangent projector
        on S² and the identity elsewhere. A correct equation of order one gives ~1e-7;
        a wrong sign or a dropped term of order one gives 0.1 or more.

    Resolution. The error is relative to the size of the whole equation, so a term that
    is small next to a large one (a stiff spring ``k = 1e7`` next to gravity, physical
    constants put in ``Scalar(..., value=...)``) moves the error by only
    ``2 |P t| / max(1, |P E_true|)`` when it is wrong, which can be below any gate. The
    finite-difference noise itself also scales with ``|L|``. ``check_eom`` reports the
    relative error only; ``assert_eom`` additionally checks that every term of the
    equation is resolved at ``rtol`` and names the ones that are not. There is always a
    grey zone: a term whose tangent part is below the noise of the residual cannot be
    verified by any normalisation.
    """
    errors, _ = _eom_errors(lagrangian, inf_work, variables, eom, inputs, params, trials, seed, h)
    return max(errors.values(), default=0.0)


def assert_eom(
    lagrangian: Expr,
    inf_work: Expr | None,
    variables: SystemVariables,
    eom: dict[str, tuple[Expr, Expr]],
    *,
    inputs: Sequence[Expr] = (),
    params: Sequence[Expr] = (),
    trials: int = 5,
    seed: int = 0,
    h: float = DEFAULT_H,
    rtol: float = 1e-5,
    require_resolved: bool = True,
) -> float:
    """``check_eom`` that raises ``AssertionError`` naming the worst variable above *rtol*.

    With ``require_resolved`` (the default) it also raises when some additive term of an
    equation is *unresolved*: its tangent part is so small next to the rest of the
    equation that a wrong sign in it would move the relative error by less than *rtol*
    (``2 |P t| / max(1, |P E_true|) < rtol`` at every trial), so the check could not have
    noticed it. The message names the term, its size and the scale of the equation. Terms
    parallel to ``q`` on S² or identically zero have nothing to resolve and are ignored.
    Pass ``require_resolved=False`` to accept such an equation on the relative error alone,
    or scale the parameters so that the terms are of comparable size.

    Returns the worst relative error when every equation passes.
    """
    errors, terms = _eom_errors(
        lagrangian, inf_work, variables, eom, inputs, params, trials, seed, h
    )
    if not errors:
        return 0.0
    key, err = max(errors.items(), key=lambda kv: kv[1])
    assert err < rtol, (
        f"equation for '{key}' disagrees with the Euler-Lagrange residual: "
        f"relative error {err:.3e} >= {rtol:.1e}"
    )
    if require_resolved:
        for key, resolutions in terms.items():
            bad = [r for r in resolutions if r.unresolved(rtol)]
            assert not bad, (
                f"equation for '{key}' passes at relative error {errors[key]:.3e} but "
                + ", ".join(
                    f"the term '{r.term}' (tangent size {r.magnitude:.2e}) is below the "
                    f"resolution {rtol * r.scale:.2e} = rtol * max(1, |E_true|)"
                    for r in bad
                )
                + "; a wrong sign or a missing factor there would not have been detected."
                " Scale the parameters so the terms are of comparable size, or pass"
                " require_resolved=False"
            )
    return err
