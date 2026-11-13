"""Independent Euler-Lagrange check for compute_eom output.

For each configuration variable the library equation E (coefficient of the variation
after IBP) must equal  -g + forcing, where g is the Euler-Lagrange residual computed
here by central finite differences of the *numeric* Lagrangian:

  R / R3 : g = d/dt dL/dv - dL/dq
  S2     : q-curve with |q|=1, w = q x qdot (so qdot = w x q); constrained EL in R^3,
           g = d/dt dLt/dqdot - dLt/dq with Lt(q, qdot) = L(q, q x qdot); E_true = -q x g + f
           (only the tangent components of E are meaningful).
  SO3    : dR = R hat(eta), dOmega = etadot + Omega x eta  =>
           g = d/dt dL/dOmega - (dL/dOmega) x Omega - vee(A - A^T),  A = R^T dL/dR.

The d/dt is a directional derivative along the sampled velocity/acceleration, so no
integration is needed: E and -g + f must agree for *random* accelerations.
"""

import os
import sys

# Run from the repository root (or anywhere): make `geomech` and `tests` importable.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
import numpy as np

from geomech import S2, SO3, Matrix, Scalar, TimeDerivative
from tests.numeric_helpers import Evaluator, hat, rand_rot, vee


def cross(a, b):
    return np.cross(np.asarray(a).reshape(3), np.asarray(b).reshape(3)).reshape(3, 1)


def _dkey(sym):
    return str(TimeDerivative(sym))


class Coords:
    """Coordinate state for the numeric Lagrangian and its time direction."""

    def __init__(self, variables, rng, extra_params=()):
        self.vars = list(variables.scalars) + list(variables.vectors) + list(variables.matrices)
        self.rng = rng
        # 'state' holds the numeric coordinates; 'vel' their d/dt; 'acc' d/dt of vel (where needed)
        self.state, self.vel, self.acc = {}, {}, {}
        for v in self.vars:
            if isinstance(v, S2):
                q = rng.normal(size=(3, 1))
                q /= np.linalg.norm(q)
                P = np.eye(3) - q @ q.T
                w = P @ rng.normal(size=(3, 1))
                wd = P @ rng.normal(size=(3, 1))
                qd = cross(w, q)
                qdd = cross(wd, q) + cross(w, cross(w, q))
                self.state[("q", v)] = q
                self.state[("qd", v)] = qd
                self.vel[("q", v)] = qd
                self.vel[("qd", v)] = qdd
                self.extra = getattr(self, "extra", {})
                self.extra[v] = (w, wd)
            elif isinstance(v, SO3):
                R = rand_rot(rng)
                Om = rng.normal(size=(3, 1))
                Omd = rng.normal(size=(3, 1))
                self.state[("R", v)] = R
                self.state[("Om", v)] = Om
                self.vel[("R", v)] = R @ hat(Om)
                self.vel[("Om", v)] = Omd
            elif isinstance(v, Scalar):
                s, sd, sdd = (float(rng.normal()) for _ in range(3))
                self.state[("s", v)] = s
                self.state[("sd", v)] = sd
                self.vel[("s", v)] = sd
                self.vel[("sd", v)] = sdd
            else:
                x, xd, xdd = (rng.normal(size=(3, 1)) for _ in range(3))
                self.state[("x", v)] = x
                self.state[("xd", v)] = xd
                self.vel[("x", v)] = xd
                self.vel[("xd", v)] = xdd
        # non-constant parameters p(t): value, pdot, pddot
        for p in extra_params:
            if isinstance(p, Scalar):
                a, b, c = (float(rng.normal()) for _ in range(3))
            elif isinstance(p, Matrix):
                a, b, c = (rng.normal(size=(3, 3)) for _ in range(3))
                if p.is_symmetric:
                    a, b, c = (m + m.T for m in (a, b, c))
            else:
                a, b, c = (rng.normal(size=(3, 1)) for _ in range(3))
            self.state[("p", p)] = a
            self.state[("pd", p)] = b
            self.vel[("p", p)] = b
            self.vel[("pd", p)] = c

    def env(self, state):
        """Library env (printed-name keys) from a coordinate state."""
        env = {}
        for v in self.vars:
            if isinstance(v, S2):
                q, qd = state[("q", v)], state[("qd", v)]
                env[str(v)] = q
                env[str(v.get_tangent_vector())] = cross(q, qd)
            elif isinstance(v, SO3):
                env[str(v)] = state[("R", v)]
                env[str(v.get_tangent_vector())] = state[("Om", v)]
            elif isinstance(v, Scalar):
                env[str(v)] = state[("s", v)]
                env[_dkey(v)] = state[("sd", v)]
            else:
                env[str(v)] = state[("x", v)]
                env[_dkey(v)] = state[("xd", v)]
        for (kind, p), val in state.items():
            if kind == "p":
                env[str(p)] = val
            if kind == "pd":
                env[_dkey(p)] = val
        return env

    def full_env(self):
        """Env at the base point including accelerations (for evaluating E)."""
        env = self.env(self.state)
        for v in self.vars:
            if isinstance(v, S2):
                w, wd = self.extra[v]
                env[_dkey(v.get_tangent_vector())] = wd
            elif isinstance(v, SO3):
                env[_dkey(v.get_tangent_vector())] = self.vel[("Om", v)]
            elif isinstance(v, Scalar):
                env[_dkey(v.t_diff())] = self.vel[("sd", v)]
            else:
                env[_dkey(v.t_diff())] = self.vel[("xd", v)]
        for (kind, p), val in self.vel.items():
            if kind == "pd":
                env[_dkey(p.t_diff())] = val
        return env


def _shift(state, key, idx, h):
    s = dict(state)
    val = (
        np.array(state[key], dtype=float, copy=True) if not np.isscalar(state[key]) else state[key]
    )
    if np.isscalar(val):
        s[key] = val + h
    else:
        val.flat[idx] += h
        s[key] = val
    return s


def _grad(fun, state, key, h=1e-5):
    val = state[key]
    if np.isscalar(val):
        return (fun(_shift(state, key, 0, h)) - fun(_shift(state, key, 0, -h))) / (2 * h)
    g = np.zeros_like(np.asarray(val, dtype=float))
    for i in range(g.size):
        g.flat[i] = (fun(_shift(state, key, i, h)) - fun(_shift(state, key, i, -h))) / (2 * h)
    return g


def _along(state, vel, h):
    s = {}
    for k, v in state.items():
        s[k] = v + h * vel[k] if k in vel else v
    return s


def _ddt(fun, state, vel, h=1e-5):
    return (fun(_along(state, vel, h)) - fun(_along(state, vel, -h))) / (2 * h)


def check_system(L, dW, variables, inputs, eom, trials=3, seed=0, extra_params=(), verbose=True):
    """Return the worst |E_lib - (-g + f)| over trials (tangent part on S2)."""
    rng = np.random.default_rng(seed)
    worst = 0.0
    for t in range(trials):
        co = Coords(variables, rng, extra_params)
        base_env = co.full_env()
        ev = Evaluator(env=base_env, rng=rng)
        for u in inputs:
            ev(u)
        for _, (_, eq) in eom.items():
            ev(eq)  # draw every constant leaf once
        ev(L)
        if hasattr(dW, "type"):
            zero_env = dict(ev.env)
            for w in co.vars:
                wk = (
                    str(w.delta())
                    if isinstance(w, Scalar) or not isinstance(w, (S2, SO3))
                    else str(w.get_variation_vector())
                )
                zero_env[wk] = 0.0 if isinstance(w, Scalar) else np.zeros((3, 1))
            Evaluator(env=zero_env, rng=rng)(dW)
            for k, val in zero_env.items():
                if k not in ev.env and not k.startswith("\\delta"):
                    ev.env[k] = val
        consts = dict(ev.env)

        def Lnum(state):
            e = dict(consts)
            e.update(co.env(state))
            return float(Evaluator(env=e, rng=rng)(L))

        results = {}
        for v in co.vars:
            if isinstance(v, S2):
                dLdqd = lambda st: _grad(Lnum, st, ("qd", v))
                g = _ddt(dLdqd, co.state, co.vel) - _grad(Lnum, co.state, ("q", v))
                q = co.state[("q", v)]
                Etrue = -cross(q, g)
                key = str(v.get_variation_vector())
                proj = np.eye(3) - q @ q.T
            elif isinstance(v, SO3):
                dLdOm = lambda st: _grad(Lnum, st, ("Om", v))
                R = co.state[("R", v)]
                Om = co.state[("Om", v)]
                G = _grad(Lnum, co.state, ("R", v))
                A = R.T @ G
                g = _ddt(dLdOm, co.state, co.vel) - cross(dLdOm(co.state), Om) - vee(A - A.T)
                Etrue = -g
                key = str(v.get_variation_vector())
                proj = np.eye(3)
            elif isinstance(v, Scalar):
                dLdsd = lambda st: _grad(Lnum, st, ("sd", v))
                g = _ddt(dLdsd, co.state, co.vel) - _grad(Lnum, co.state, ("s", v))
                Etrue = -g
                key = str(v.delta())
                proj = 1.0
            else:
                dLdxd = lambda st: _grad(Lnum, st, ("xd", v))
                g = _ddt(dLdxd, co.state, co.vel) - _grad(Lnum, co.state, ("x", v))
                Etrue = -g
                key = str(v.delta())
                proj = np.eye(3)
            # forcing: evaluate dW with this variation = basis vector, others zero
            if hasattr(dW, "type"):
                f = []
                basis = [1.0] if isinstance(v, Scalar) else [np.eye(3)[:, [k]] for k in range(3)]
                for b in basis:
                    e = dict(consts)
                    for w in co.vars:
                        wk = (
                            str(w.delta())
                            if isinstance(w, (Scalar,)) or not isinstance(w, (S2, SO3))
                            else str(w.get_variation_vector())
                        )
                        e[wk] = 0.0 if isinstance(w, Scalar) else np.zeros((3, 1))
                    e[key] = b
                    f.append(float(Evaluator(env=e, rng=rng)(dW)))
                f = f[0] if isinstance(v, Scalar) else np.array(f).reshape(3, 1)
            else:
                f = 0.0
            Etrue = Etrue + f
            Elib = Evaluator(env=dict(consts), rng=rng)(eom[key][1])
            Elib = np.asarray(Elib, dtype=float).reshape(-1, 1) if not np.isscalar(Elib) else Elib
            if np.isscalar(proj):
                err = abs(float(Elib) - float(Etrue))
                scale = max(1.0, abs(float(Etrue)))
            else:
                d = proj @ (np.asarray(Elib).reshape(3, 1) - Etrue)
                err = float(np.max(np.abs(d)))
                scale = max(1.0, float(np.max(np.abs(proj @ Etrue))))
            results[key] = (err / scale, Elib, Etrue)
            worst = max(worst, err / scale)
        if verbose:
            for k, (e, a, b) in results.items():
                print(f"  trial {t} {k}: rel err {e:.2e}")
                if e > 1e-4:
                    print("     lib :", np.asarray(a).reshape(-1))
                    print("     true:", np.asarray(b).reshape(-1))
    return worst
