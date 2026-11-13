"""Numeric evaluation of geomech expression trees for test cross-checks.

``Evaluator`` assigns a random value to every leaf (keyed by its printed name) unless a
value is injected through ``env``; ``hat``, ``vee``, ``rand_rot`` and ``expm_so3`` are the
SO(3) primitives. ``finite_difference_check`` perturbs a sampled trajectory point along a
constraint-respecting direction and compares the central difference of each equation of
motion with the symbolic linearization produced by ``linearize``.
"""

from __future__ import annotations

import numpy as np

from geomech import (
    S2,
    SO3,
    Add,
    Cross,
    Dot,
    Hat,
    MAdd,
    Matrix,
    MMMul,
    Mul,
    MVMul,
    Scalar,
    SMMul,
    SVMul,
    TimeDerivative,
    Transpose,
    VAdd,
    Variation,
    Vector,
    Vee,
    VVMul,
)

# ---------------------------------------------------------------------------
# SO(3) primitives
# ---------------------------------------------------------------------------


def hat(v):
    v = np.asarray(v, dtype=float).reshape(3)
    return np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]], dtype=float)


def vee(M):
    return np.array([[M[2, 1]], [M[0, 2]], [M[1, 0]]])


def rand_rot(rng):
    A = rng.normal(size=(3, 3))
    Q, _ = np.linalg.qr(A)
    if np.linalg.det(Q) < 0:
        Q[:, 0] = -Q[:, 0]
    return Q


def expm_so3(v):
    """Rodrigues formula: exp(Hat(v)) for v in R³."""
    v = np.asarray(v, dtype=float).reshape(3)
    theta = np.linalg.norm(v)
    K = hat(v)
    if theta < 1e-12:
        return np.eye(3) + K + 0.5 * K @ K
    return np.eye(3) + np.sin(theta) / theta * K + (1 - np.cos(theta)) / theta**2 * (K @ K)


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
        elif isinstance(node, SO3):
            val = rand_rot(rng)
        elif isinstance(node, S2):
            v = rng.normal(size=(3, 1))
            val = v / np.linalg.norm(v)
        elif isinstance(node, Matrix):
            A = rng.normal(size=(3, 3))
            val = A + A.T if node.is_symmetric else A
        elif isinstance(node, Vector):
            val = rng.normal(size=(3, 1))
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


# ---------------------------------------------------------------------------
# Trajectory points and tangent perturbations
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


def sample_state(variables, inputs, rng):
    """Random trajectory point on the manifold: |q| = 1, q·ω = 0, q·ω̇ = 0, R ∈ SO(3)."""
    env = {}
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


def sample_direction(variables, inputs, env, rng):
    """Random perturbation direction tangent to the constraints.

    Returns ``(raw, pert)``. ``raw`` holds the generators used by ``perturbed_state``
    (ξ, w, a for S²; η, δΩ, δΩ̇ for SO(3); δx, δẋ, δẍ for R³; δs, δṡ, δs̈ for R; δu).
    ``pert`` maps the printed perturbation symbols of ``linearize`` (states, velocity
    dstates and input variations) to their numeric values; for S² the induced
    velocity direction is δω = ξ × ω + w and δω̇ = ξ × ω̇ + a.
    """
    raw, pert = {}, {}
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
            dx, dv, da = (rng.normal(size=(3, 1)) for _ in range(3))
            raw[kq] = (dx, dv, da)
            d = var.delta()
            pert[str(d)], pert[str(d.t_diff())], pert[str(d.t_diff().t_diff())] = dx, dv, da
    for u in inputs:
        du = float(rng.normal()) if isinstance(u, Scalar) else rng.normal(size=(3, 1))
        raw[str(u)] = du
        pert[str(u.delta())] = du
    return raw, pert


def perturbed_state(variables, inputs, env, raw, eps):
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


def linearized_rhs(lin, row, ev, pert, dvals):
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
