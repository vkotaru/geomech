"""Check the quadrotor and flexible-cable equations of geomech against the source papers.

    python sim/check_references.py

Uses the harness of ``check_appendix.py`` (random admissible states, ``compute_eom``
output, tangent projection of S² equations). The reference equations are coded from the
papers, not from the appendix:

* Lee, Leok, McClamroch, "Geometric tracking control of a quadrotor UAV on SE(3)", CDC
  2010 (arXiv:1003.2005), Eqs. (2)-(5): ẋ = v, m v̇ = m g e3 - f R e3, Ṙ = R Ω̂,
  J Ω̇ + Ω × J Ω = M, with e3 pointing down. A state of geomech (e3 up) is mapped to that
  frame by the rotation D = diag(1, -1, -1) about e1, which reverses e3:
  x' = D x, R' = D R D, Ω' = D Ω, J' = D J D, M' = D M, and f unchanged.

* Kotaru, Wu, Sreenath, "Differential-flatness and control of quadrotor(s) with a payload
  suspended through flexible cable(s)", ICC 2018 (arXiv:1711.04895), Eqs. (1)-(5):
  x_i = x_{i-1} + l_i q_i with x_0 the quadrotor and x_n the load, tension magnitudes T_i,
      m_Q (ẍ_0 + g e3) = f R e3 + T_1 q_1,
      m_j (ẍ_j + g e3) = -T_j q_j + T_{j+1} q_{j+1},   j = 1, ..., n-1,
      m_n (ẍ_n + g e3) = -T_n q_n,
      J Ω̇ + Ω̂ J Ω = M.
  The tensions are eliminated in two ways: (a) by d'Alembert's principle, projecting the
  residuals r_i onto the virtual displacements δx_i = δx_0 + Σ_{j<=i} l_j ξ_j × q_j that
  (1) admits, which gives Σ_i r_i for δx_0 and l_k q_k × Σ_{i>=k} r_i for ξ_k, for
  arbitrary random T_i; (b) by solving (2)-(5) with q_j·ω̇_j = 0 for the accelerations and
  the tensions at a random state and input, and evaluating the equations of geomech there.
"""

from __future__ import annotations

import numpy as np
from check_appendix import (
    SEEDS,
    TOL,
    Evaluator,
    State,
    cross,
    derived,
    dot,
    project,
    qddot_s2,
    random_params,
    relative,
)

import geomech.systems as S
from geomech.numeric import _kinematic_symbols, hat

D = np.diag([1.0, -1.0, -1.0])  # rotation by π about e1: reverses e3


def states(system):
    """Random states at the default and at random parameters, seed by seed."""
    for randomize in (False, True):
        for seed in SEEDS:
            rng = np.random.default_rng(seed)
            params = random_params(system, rng) if randomize else None
            yield State(system, rng, params), rng


# ---------------------------------------------------------------- Lee et al. 2010


def lee2010():
    system = S.get("quadrotor")
    eqs = derived(system)
    worst = {"x (Eq. 3)": 0.0, "R (Eq. 5)": 0.0, "kinematics (Eqs. 2, 4)": 0.0}
    for x, _ in states(system):
        m, g, e3, J, f, M = x["m"], x["g"], x["e3"], x["J"], x["f"], x["M"]
        X, V, A = x.q("x")
        R, Om, Omd = x.q("R")
        # the state in the frame of Lee et al.
        Xp, Vp, Ap = D @ X, D @ V, D @ A
        Rp, Op, Odp, Jp, Mp = D @ R @ D, D @ Om, D @ Omd, D @ J @ D, D @ M
        # residuals of Eqs. (3) and (5), as printed
        r3 = m * Ap - (m * g * e3 - f * Rp @ e3)
        r5 = Jp @ Odp + cross(Op, Jp @ Op) - Mp
        # Eqs. (2) and (4) in that frame: ẋ' = D ẋ, Ṙ' = R' hat(Ω')
        kin = np.max(np.abs(D @ (R @ hat(Om)) @ D - Rp @ hat(Op)))
        ev = Evaluator(env=dict(x.env))
        # back to the geomech frame; compute_eom returns forces minus inertia (E = 0)
        worst["x (Eq. 3)"] = max(worst["x (Eq. 3)"], relative(ev(eqs["x"]), -D.T @ r3))
        worst["R (Eq. 5)"] = max(worst["R (Eq. 5)"], relative(ev(eqs["R"]), -D.T @ r5))
        worst["kinematics (Eqs. 2, 4)"] = max(worst["kinematics (Eqs. 2, 4)"], kin)
        assert np.allclose(Xp[2], -X[2]) and np.allclose(Vp[2], -V[2])
    for k, v in worst.items():
        print(f"Lee 2010 quadrotor, {k:<24} worst relative difference {v:.1e}")
    return max(worst.values())


# ---------------------------------------------------------- Kotaru et al. 2018


def kotaru_residuals(x, n, T):
    """Residuals r_0, ..., r_n of Eqs. (2)-(4) and r_R of Eq. (5) at the state of x."""
    mQ, g, e3, J, f, M = x["m_Q"], x["g"], x["e3"], x["J"], x["f"], x["M"]
    _, _, a0 = x.q("x_Q")
    R, Om, Omd = x.q("R")
    q = {j: x.q(f"q{j}")[0] for j in range(1, n + 1)}
    acc = [a0]
    for j in range(1, n + 1):
        acc.append(acc[-1] + x[f"l{j}"] * qddot_s2(*x.q(f"q{j}")))  # Eq. (1), twice
    r = [mQ * (acc[0] + g * e3) - (f * R @ e3 + T[1] * q[1])]
    for j in range(1, n):
        r.append(x[f"m{j}"] * (acc[j] + g * e3) - (-T[j] * q[j] + T[j + 1] * q[j + 1]))
    r.append(x[f"m{n}"] * (acc[n] + g * e3) - (-T[n] * q[n]))
    rR = J @ Omd + cross(Om, J @ Om) - M
    return r, rR


def eliminated(x, n, T):
    """The tension-free equations, by d'Alembert's principle, as forces minus inertia."""
    r, rR = kotaru_residuals(x, n, T)
    out = {"x_Q": -sum(r), "R": -rR}
    for k in range(1, n + 1):
        out[f"q{k}"] = -x[f"l{k}"] * cross(x.q(f"q{k}")[0], sum(r[k:]))
    return out


def solve_kotaru(x, system, n):
    """Accelerations and tensions solving Eqs. (2)-(5) with q_j · ω̇_j = 0 at the state of x;
    sets the accelerations of x and returns the tensions."""
    keys = {str(v): k for v, k in _kinematic_symbols(system.variables)}
    acc_keys = [keys["x_Q"][2]] + [keys[f"q{j}"][2] for j in range(1, n + 1)] + [keys["R"][2]]

    def refresh():
        x.kin = {str(v): tuple(x.env[k] for k in ks) for v, ks in keys_list}

    keys_list = _kinematic_symbols(system.variables)

    def residual(z):
        for i, k in enumerate(acc_keys):
            x.env[k] = z[3 * i : 3 * i + 3].reshape(3, 1)
        refresh()
        T = {j: z[3 * len(acc_keys) + j - 1] for j in range(1, n + 1)}
        r, rR = kotaru_residuals(x, n, T)
        tangent = [dot(x.q(f"q{j}")[0], x.q(f"q{j}")[2]) for j in range(1, n + 1)]
        return np.concatenate([np.ravel(v) for v in r] + [np.ravel(rR), tangent])

    size = 3 * len(acc_keys) + n
    b = residual(np.zeros(size))
    A = np.column_stack([residual(np.eye(size)[i]) - b for i in range(size)])
    z = np.linalg.solve(A, -b)
    left = residual(z)  # also leaves the accelerations of x at the solution
    assert np.max(np.abs(left)) < 1e-9, left
    return z[3 * len(acc_keys) :]


def kotaru2018(n):
    system = S.get("quadrotor_flexible_cable", n_links=n)
    eqs = derived(system)
    worst_a = worst_t = worst_b = 0.0
    for x, rng in states(system):
        ev = Evaluator(env=dict(x.env))
        T1 = {j: rng.normal() for j in range(1, n + 1)}
        T2 = {j: 10 * rng.normal() for j in range(1, n + 1)}
        want, other = eliminated(x, n, T1), eliminated(x, n, T2)
        for v, eq in eqs.items():
            got, w, o = project(system, x, v, ev(eq), want[v], other[v])
            worst_a = max(worst_a, relative(got, w))
            worst_t = max(worst_t, relative(o, w))  # independent of the tensions
        # (b) forward dynamics of the Newton-Euler model, then the equations of geomech
        scale = max(np.max(np.abs(Evaluator(env=dict(x.env))(eq))) for eq in eqs.values())
        solve_kotaru(x, system, n)
        ev = Evaluator(env=dict(x.env))
        for v, eq in eqs.items():
            (res,) = project(system, x, v, ev(eq))
            worst_b = max(worst_b, float(np.max(np.abs(res))) / max(1.0, scale))
    print(
        f"Kotaru 2018 flexible cable (n = {n}): projection {worst_a:.1e}, "
        f"tension independence {worst_t:.1e}, forward dynamics {worst_b:.1e}"
    )
    return max(worst_a, worst_t, worst_b)


def main():
    worst = max([lee2010()] + [kotaru2018(n) for n in (1, 3, 5)])
    print(f"all references: {worst:.1e}")
    if worst > TOL:
        raise SystemExit(f"a reference model differs from geomech by {worst:.1e}")


if __name__ == "__main__":
    main()
