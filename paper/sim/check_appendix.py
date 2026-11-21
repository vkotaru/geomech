"""Check the hand-arranged equations of the appendix against the output of geomech.

    python sim/check_appendix.py

For each entry below, the equations as displayed in the appendix are written in NumPy
(``hand``) and compared, at random admissible states, with the equations that
``compute_eom`` derives for the library system: the relative difference of every
equation is printed and must be at round-off level. An equation of a point of S² is
compared after projection onto its tangent plane. Each displayed equation ``lhs = rhs``
is coded as ``rhs - lhs``, the sign convention of ``compute_eom`` (``E = 0``).

Every displayed equation of the appendix (ENTRIES and NEW_ENTRIES) is checked at the
default parameters and again at random admissible parameters, at the sizes the appendix
states. An equation displayed as
a nonzero multiple of the coefficient that compute_eom returns is scaled by the factor
declared in SCALES. ``check_references.py`` compares the quadrotor and the flexible cable
with the models of the cited papers.
"""

from __future__ import annotations

import re

import numpy as np

import geomech.systems as S
from geomech import S2, compute_eom
from geomech.numeric import Evaluator, _kinematic_symbols, sample_state

SEEDS = range(5)
TOL = 1e-12


def cross(a, b):
    return np.cross(a.ravel(), b.ravel()).reshape(3, 1)


def dot(a, b):
    return float(a.ravel() @ b.ravel())


def qddot_s2(q, w, wd):
    return cross(wd, q) - dot(w, w) * q


class State:
    """Values of one random state: configuration, velocity and acceleration of each
    variable by name, parameters and inputs by printed name."""

    def __init__(self, system, rng, params=None):
        params = params or {
            str(k): (np.asarray(v, float).reshape(3, -1) if np.ndim(v) else float(v))
            for k, v in system.params.items()
        }
        self.env = sample_state(system.variables, system.inputs, rng, env=params)
        self.kin = {
            str(var): tuple(self.env[k] for k in keys)
            for var, keys in _kinematic_symbols(system.variables)
        }

    def __getitem__(self, name):
        return self.env[name]

    def q(self, name):
        return self.kin[name]


# --------------------------------------------------------------------- entries


def planar_pendulum(x):
    m, g, l, tau = x["m"], x["g"], x["l"], x["\\tau"]
    e2, e3 = x["e2"], x["e3"]
    q, w, wd = x.q("q")
    return {"q": tau - (m * l * l * wd + m * g * l * dot(e3, cross(e2, q)))}


def double_planar_pendulum(x):
    m1, m2, l1, l2, g = x["m1"], x["m2"], x["l1"], x["l2"], x["g"]
    e2, e3 = x["e2"], x["e3"]
    q1, w1, wd1 = x.q("q1")
    q2, w2, wd2 = x.q("q2")
    c, s = dot(q1, q2), dot(q1, cross(e2, q2))
    E1 = (m1 + m2) * l1 * l1 * wd1 + m2 * l1 * l2 * c * wd2 + m2 * l1 * l2 * s * w2 * w2
    E1 += (m1 + m2) * g * l1 * dot(e3, cross(e2, q1))
    E2 = m2 * l1 * l2 * c * wd1 + m2 * l2 * l2 * wd2 - m2 * l1 * l2 * s * w1 * w1
    E2 += m2 * g * l2 * dot(e3, cross(e2, q2))
    return {"q1": x["\\tau_1"] - E1, "q2": x["\\tau_2"] - E2}


def cart_planar_pendulum(x):
    M, m, g, l, f = x["M"], x["m"], x["g"], x["l"], x["f"]
    e1, e2, e3 = x["e1"], x["e2"], x["e3"]
    s, sd, sdd = x.q("s")
    q, w, wd = x.q("q")
    c = dot(e1, cross(e2, q))
    Es = (M + m) * sdd + m * l * c * wd - m * l * w * w * dot(e1, q) + (M + m) * g * dot(e1, e3)
    Eq = m * l * c * sdd + m * l * l * wd + m * g * l * dot(e3, cross(e2, q))
    return {"s": f - Es, "q": -Eq}


def _spacecraft(x, panels):
    """Body equation ``Ḣ + Ω × H = 0`` and the wheel and panel equations."""
    Jb = x["J_b"]
    R, Om, Omd = x.q("R")
    H = Jb @ Om
    Hd = Jb @ Omd
    out = {}
    for i in (1, 2, 3):
        I, b = x[f"I_{i}"], x[f"b{i}"]
        _, w, wd = x.q(f"q_{i}")
        H = H + I * (w + dot(b, Om)) * b
        Hd = Hd + I * (wd + dot(b, Omd)) * b
        out[f"q_{i}"] = x[f"\\tau_{i}"] - I * (wd + dot(b, Omd))
    for k in range(1, panels + 1):
        a, mp, d = x[f"a{k}"], x[f"m_p{k}"], x[f"d_{k}"]
        Ia, In, Iw = x[f"Ia_{k}"], x[f"In_{k}"], x[f"Iw_{k}"]
        n, w, wd = x.q(f"n_{k}")
        v = cross(a, n)  # w_k in the text: the third panel axis
        H = H + mp * d * d * cross(a, cross(Om, a)) + Ia * (w + dot(a, Om)) * a
        H = H + In * dot(n, Om) * n + Iw * dot(v, Om) * v
        nd, vd = w * v, -w * n  # body-frame rates of n_k and a_k × n_k
        Hd = Hd + mp * d * d * cross(a, cross(Omd, a)) + Ia * (wd + dot(a, Omd)) * a
        Hd = Hd + In * ((dot(nd, Om) + dot(n, Omd)) * n + dot(n, Om) * nd)
        Hd = Hd + Iw * ((dot(vd, Om) + dot(v, Omd)) * v + dot(v, Om) * vd)
        out[f"n_{k}"] = x[f"\\mu_{k}"] - (
            Ia * (wd + dot(a, Omd)) - (In - Iw) * dot(n, Om) * dot(v, Om)
        )
    out["R"] = -(Hd + cross(Om, H))
    return out


def spacecraft_reaction_wheels_s1(x):
    return _spacecraft(x, 0)


def spacecraft_solar_arrays(x):
    return _spacecraft(x, 2)


def hose(n):
    def f(x):
        g, e3 = x["g"], x["e3"]
        x0, v0, a0 = x.q("x_0")
        qs = [x.q(f"q{j}") for j in range(1, n + 1)]
        acc = [a0]
        for j, (q, w, wd) in enumerate(qs, start=1):
            acc.append(acc[-1] + x[f"l{j}"] * qddot_s2(q, w, wd))
        Mk = [x[f"m{k}"] for k in range(n + 1)]
        Mk[0] += x["m_Q0"]
        Mk[n] += x[f"m_Q{n}"]
        F = {k: x[f"f{k}"] * x.q(f"R{k}")[0] @ e3 for k in (0, n)}
        out = {"x_0": F[0] + F[n] - sum(Mk[k] * (acc[k] + g * e3) for k in range(n + 1))}
        for j in range(1, n + 1):
            q = qs[j - 1][0]
            tail = sum(Mk[k] * (acc[k] + g * e3) for k in range(j, n + 1))
            out[f"q{j}"] = x[f"l{j}"] * cross(q, F[n] - tail)
        for k in (0, n):
            J = x[f"J{k}"]
            _, Om, Omd = x.q(f"R{k}")
            out[f"R{k}"] = x[f"M{k}"] - (J @ Omd + cross(Om, J @ Om))
        return out

    return f


# The entries below transcribe the remaining displayed equations of 07_appendix.tex (and
# the main-text equations it refers to: sp_result, rp_result, ql_xL--ql_R), each coded as
# ``rhs - lhs``. Where a displayed equation is a nonzero multiple of the coefficient that
# compute_eom returns (a ``= 0`` equation, or one divided by a length), the factor is
# declared in SCALES below and the displayed form is multiplied by it before comparison.


def body_acc(R, Om, Omd, r):
    """d²/dt²(R r) of eq:body_acc."""
    return R @ (cross(Omd, r) + cross(Om, cross(Om, r)))


def euler(J, Om, Omd):
    return J @ Omd + cross(Om, J @ Om)


def point_mass(x):
    _, _, xdd = x.q("x")
    return {"x": x["F"] - x["m"] * xdd}


def spring_mass(x):
    s, _, sdd = x.q("x")
    return {"x": x["f"] - (x["m"] * sdd + x["k"] * s)}


def spherical_pendulum(x):  # eq:sp_result
    m, g, l, e3 = x["m"], x["g"], x["l"], x["e3"]
    q, w, wd = x.q("q")
    return {"q": x["M"] - (m * l * l * wd + m * g * l * cross(q, e3))}


def particle_on_rotating_sphere(x):  # eq:app_rot_sphere
    m, r, W0, e3, C = x["m"], x["L"], x["\\Omega"], x["e3"], x["C"]
    q, w, wd = x.q("q")
    rhs = 2 * m * r * r * W0 * dot(q, e3) * cross(w, q) + m * r * r * W0 * W0 * cross(q, C @ q)
    return {"q": x["M"] + rhs - m * r * r * wd}


def double_spherical_pendulum(x):  # eq:app_dsp1, eq:app_dsp2
    m1, m2, l1, l2, g, e3 = x["m1"], x["m2"], x["l1"], x["l2"], x["g"], x["e3"]
    q1, w1, wd1 = x.q("q1")
    q2, w2, wd2 = x.q("q2")
    E1 = (m1 + m2) * l1**2 * wd1 - m2 * l1 * l2 * cross(q1, cross(q2, wd2))
    E1 = E1 - m2 * l1 * l2 * dot(w2, w2) * cross(q1, q2) + (m1 + m2) * g * l1 * cross(q1, e3)
    E2 = m2 * l2**2 * wd2 - m2 * l1 * l2 * cross(q2, cross(q1, wd1))
    E2 = E2 - m2 * l1 * l2 * dot(w1, w1) * cross(q2, q1) + m2 * g * l2 * cross(q2, e3)
    return {"q1": x["M1"] - E1, "q2": x["M2"] - E2}


def n_link_spherical_chain(n):  # eq:app_nlink
    def f(x):
        g, e3 = x["g"], x["e3"]
        m = {i: x[f"m{i}"] for i in range(1, n + 1)}
        l = {i: x[f"l{i}"] for i in range(1, n + 1)}
        kin = {i: x.q(f"q{i}") for i in range(1, n + 1)}
        qdd = {j: qddot_s2(*kin[j]) for j in kin}
        out = {}
        for k in range(1, n + 1):
            qk = kin[k][0]
            E = sum(
                l[k] * l[j] * sum(m[i] for i in range(max(k, j), n + 1)) * cross(qk, qdd[j])
                for j in range(1, n + 1)
            )
            E = E + sum(m[i] for i in range(k, n + 1)) * g * l[k] * cross(qk, e3)
            out[f"q{k}"] = x[f"M{k}"] - E
        return out

    return f


def elastic_spherical_pendulum(x):
    m, g, k, l0, e3 = x["m"], x["g"], x["k"], x["l0"], x["e3"]
    l, ld, ldd = x.q("l")
    q, w, wd = x.q("q")
    El = m * (ldd - l * dot(w, w)) + m * g * dot(q, e3) + k * (l - l0)
    Eq = m * l * l * wd + 2 * m * l * ld * w + m * g * l * cross(q, e3)
    return {"l": x["f_l"] - El, "q": x["M_q"] - Eq}


def cart_spherical_pendulum(x):
    # the library names the cart mass M and the rail axis e_c = e1
    mc, m, g, l, ec, e3 = x["M"], x["m"], x["g"], x["l"], x["e1"], x["e3"]
    s, sd, sdd = x.q("s")
    q, w, wd = x.q("q")
    Es = (mc + m) * sdd + m * l * dot(ec, qddot_s2(q, w, wd)) + (mc + m) * g * dot(ec, e3)
    Eq = m * l * l * wd + m * l * sdd * cross(q, ec) + m * g * l * cross(q, e3)
    return {"s": x["f"] - Es, "q": -Eq}


def free_rigid_body(x):
    _, Om, Omd = x.q("R")
    return {"R": x["M"] - euler(x["J"], Om, Omd)}


def heavy_top(x):
    R, Om, Omd = x.q("R")
    m, g, l, e3 = x["m"], x["g"], x["l"], x["e3"]
    return {"R": x["M"] - (euler(x["J"], Om, Omd) + m * g * l * cross(e3, R.T @ e3))}


def rigid_pendulum(x):  # eq:rp_result
    R, Om, Omd = x.q("R")
    m, g, rho, e3 = x["m"], x["g"], x["\\rho"], x["e3"]
    return {"R": x["M"] - (euler(x["J"], Om, Omd) + m * g * cross(rho, R.T @ e3))}


def double_rigid_pendulum(x):  # eq:app_drp1, eq:app_drp2
    m1, m2, g, e3 = x["m1"], x["m2"], x["g"], x["e3"]
    rho1, rho2, d1 = x["\\rho_1"], x["\\rho_2"], x["l_1"]
    R1, O1, Od1 = x.q("R1")
    R2, O2, Od2 = x.q("R2")
    x1dd = body_acc(R1, O1, Od1, rho1)
    x2dd = body_acc(R1, O1, Od1, d1) + body_acc(R2, O2, Od2, rho2)
    E1 = euler(x["J1"], O1, Od1) + m1 * cross(rho1, R1.T @ (x1dd + g * e3))
    E1 = E1 + m2 * cross(d1, R1.T @ (x2dd + g * e3))
    E2 = euler(x["J2"], O2, Od2) + m2 * cross(rho2, R2.T @ (x2dd + g * e3))
    return {"R1": x["M1"] - E1, "R2": x["M2"] - E2}


def n_link_rigid_chain(n):  # eq:app_rigid_chain; the library names d_k l_k
    def f(x):
        g, e3 = x["g"], x["e3"]
        kin = {k: x.q(f"R{k}") for k in range(1, n + 1)}
        d = {k: x.env.get(f"l_{k}") for k in range(1, n + 1)}  # d_n does not enter L
        xdd = {}
        for k in range(1, n + 1):
            xdd[k] = sum(body_acc(*kin[i], d[i]) for i in range(1, k)) + body_acc(
                *kin[k], x[f"\\rho_{k}"]
            )
        tau = {k: x[f"\\tau_{k}"] for k in range(1, n + 1)}
        tau[n + 1] = np.zeros((3, 1))
        out = {}
        for k in range(1, n + 1):
            R, Om, Omd = kin[k]
            mk, rho = x[f"m{k}"], x[f"\\rho_{k}"]
            E = euler(x[f"J{k}"], Om, Omd) + mk * cross(rho, R.T @ (xdd[k] + g * e3))
            if k < n:
                tail = sum(x[f"m{j}"] * (xdd[j] + g * e3) for j in range(k + 1, n + 1))
                E = E + cross(d[k], R.T @ tail)
            out[f"R{k}"] = R.T @ (tau[k] - tau[k + 1]) - E
        return out

    return f


def two_rigid_bodies_spherical_joint(x):
    m1, m2, g, e3 = x["m1"], x["m2"], x["g"], x["e3"]
    rho1, rho2 = x["\\rho_1"], x["\\rho_2"]
    _, _, xdd = x.q("x")
    R1, O1, Od1 = x.q("R1")
    R2, O2, Od2 = x.q("R2")
    x2dd = xdd + body_acc(R1, O1, Od1, rho1) + body_acc(R2, O2, Od2, rho2)
    return {
        "x": x["F"] - (m1 * (xdd + g * e3) + m2 * (x2dd + g * e3)),
        "R1": x["M1"] - (euler(x["J1"], O1, Od1) + m2 * cross(rho1, R1.T @ (x2dd + g * e3))),
        "R2": x["M2"] - (euler(x["J2"], O2, Od2) + m2 * cross(rho2, R2.T @ (x2dd + g * e3))),
    }


def rigid_body(x):
    m, g, e3 = x["m"], x["g"], x["e3"]
    _, _, xdd = x.q("x")
    _, Om, Omd = x.q("R")
    return {"x": x["F"] - m * (xdd + g * e3), "R": x["M"] - euler(x["J"], Om, Omd)}


def spacecraft_reaction_wheels(x):  # eq:app_sc_R, eq:app_sc_wheel
    _, Om, Omd = x.q("R")
    I = {i: x[f"I_{i}"] for i in (1, 2, 3)}
    b = {i: x[f"b{i}"] for i in (1, 2, 3)}
    th = {i: x.q(f"\\theta_{i}") for i in (1, 2, 3)}
    Jp = x["J_b"] + sum(I[i] * b[i] @ b[i].T for i in (1, 2, 3))
    ER = Jp @ Omd + sum(I[i] * th[i][2] * b[i] for i in (1, 2, 3))
    ER = ER + cross(Om, Jp @ Om + sum(I[i] * th[i][1] * b[i] for i in (1, 2, 3)))
    out = {"R": -ER}
    for i in (1, 2, 3):
        out[f"\\theta_{i}"] = x[f"\\tau_{i}"] - I[i] * (dot(b[i], Omd) + th[i][2])
    return out


def biped_3d(x):
    mH, m, g, l, a, J, e3 = x["m_H"], x["m"], x["g"], x["l"], x["a"], x["J"], x["e3"]
    R1, O1, Od1 = x.q("R1")
    R2, O2, Od2 = x.q("R2")
    xHdd = l * body_acc(R1, O1, Od1, e3)
    x1dd = a * body_acc(R1, O1, Od1, e3)
    x2dd = xHdd - (l - a) * body_acc(R2, O2, Od2, e3)
    G1 = l * mH * (xHdd + g * e3) + a * m * (x1dd + g * e3) + l * m * (x2dd + g * e3)
    E1 = euler(J, O1, Od1) + cross(e3, R1.T @ G1)
    E2 = euler(J, O2, Od2) - (l - a) * m * cross(e3, R2.T @ (x2dd + g * e3))
    return {"R1": R1.T @ (x["M_A"] - x["M_H"]) - E1, "R2": R2.T @ x["M_H"] - E2}


def quadrotor(x):
    m, g, e3 = x["m"], x["g"], x["e3"]
    _, _, xdd = x.q("x")
    R, Om, Omd = x.q("R")
    return {
        "x": x["f"] * R @ e3 - m * (xdd + g * e3),
        "R": x["M"] - euler(x["J"], Om, Omd),
    }


def quadrotor_offset_com(x):  # eq:app_offset_x, eq:app_offset_R
    m, g, e3, rho = x["m"], x["g"], x["e3"], x["\\rho"]
    _, _, xdd = x.q("x")
    R, Om, Omd = x.q("R")
    xcdd = xdd + body_acc(R, Om, Omd, rho)
    return {
        "x": x["f"] * R @ e3 - m * (xcdd + g * e3),
        "R": x["M"] - (euler(x["J"], Om, Omd) + m * cross(rho, R.T @ (xcdd + g * e3))),
    }


def quadrotor_suspended_load(x):  # eq:ql_xL, eq:ql_q, eq:ql_R
    mQ, mL, l, g, e3 = x["m_Q"], x["m_L"], x["l"], x["g"], x["e3"]
    _, _, xLdd = x.q("x_L")
    q, w, wd = x.q("q")
    R, Om, Omd = x.q("R")
    fRe3 = x["f"] * R @ e3
    Ex = (mQ + mL) * (xLdd + g * e3) - mQ * l * (cross(wd, q) - dot(w, w) * q)
    Eq = mQ * l * wd - cross(q, mQ * (xLdd + g * e3) - fRe3)
    return {"x_L": fRe3 - Ex, "q": -Eq, "R": x["M"] - euler(x["J"], Om, Omd)}


def quadrotor_elastic_load(x):
    mQ, mL, k, l0, g, e3 = x["m_Q"], x["m_L"], x["k"], x["l0"], x["g"], x["e3"]
    l, ld, ldd = x.q("l")
    _, _, xLdd = x.q("x_L")
    q, w, wd = x.q("q")
    R, Om, Omd = x.q("R")
    qd = cross(w, q)
    xQdd = xLdd - ldd * q - 2 * ld * qd - l * qddot_s2(q, w, wd)
    FQ = mQ * (xQdd + g * e3) - x["f"] * R @ e3
    return {
        "x_L": -(mL * (xLdd + g * e3) + FQ),
        "l": k * (l - l0) - dot(q, FQ),
        "q": -cross(q, FQ),
        "R": x["M"] - euler(x["J"], Om, Omd),
    }


def quadrotor_flexible_cable(n):
    def f(x):
        mQ, g, e3 = x["m_Q"], x["g"], x["e3"]
        _, _, xQdd = x.q("x_Q")
        R, Om, Omd = x.q("R")
        kin = {j: x.q(f"q{j}") for j in range(1, n + 1)}
        m = {i: x[f"m{i}"] for i in range(1, n + 1)}
        xdd = {
            i: xQdd + sum(x[f"l{j}"] * qddot_s2(*kin[j]) for j in range(1, i + 1))
            for i in range(1, n + 1)
        }
        Fi = {i: m[i] * (xdd[i] + g * e3) for i in xdd}
        out = {"x_Q": x["f"] * R @ e3 - (mQ * (xQdd + g * e3) + sum(Fi.values()))}
        for k in range(1, n + 1):
            out[f"q{k}"] = -cross(kin[k][0], sum(Fi[i] for i in range(k, n + 1)))
        out["R"] = x["M"] - euler(x["J"], Om, Omd)
        return out

    return f


def multi_quadrotor_point_mass(n):
    def f(x):
        mL, g, e3 = x["m_L"], x["g"], x["e3"]
        _, _, xLdd = x.q("x_L")
        out, F = {}, {}
        for i in range(1, n + 1):
            q, w, wd = x.q(f"q{i}")
            R, Om, Omd = x.q(f"R{i}")
            xQdd = xLdd - x[f"l{i}"] * qddot_s2(q, w, wd)
            F[i] = x[f"m{i}"] * (xQdd + g * e3) - x[f"f{i}"] * R @ e3
            out[f"q{i}"] = -cross(q, F[i])
            out[f"R{i}"] = x[f"M{i}"] - euler(x[f"J{i}"], Om, Omd)
        out["x_L"] = -(mL * (xLdd + g * e3) + sum(F.values()))
        return out

    return f


def multi_quadrotor_rigid_load(n):
    def f(x):
        mL, g, e3 = x["m_L"], x["g"], x["e3"]
        _, _, xLdd = x.q("x_L")
        RL, OL, OdL = x.q("R_L")
        out, F = {}, {}
        for i in range(1, n + 1):
            q, w, wd = x.q(f"q{i}")
            R, Om, Omd = x.q(f"R{i}")
            rho = x[f"\\rho_{i}"]
            xQdd = xLdd + body_acc(RL, OL, OdL, rho) - x[f"l{i}"] * qddot_s2(q, w, wd)
            F[i] = x[f"m{i}"] * (xQdd + g * e3) - x[f"f{i}"] * R @ e3
            out[f"q{i}"] = -cross(q, F[i])
            out[f"R{i}"] = x[f"M{i}"] - euler(x[f"J{i}"], Om, Omd)
        out["x_L"] = -(mL * (xLdd + g * e3) + sum(F.values()))
        torque = sum(cross(x[f"\\rho_{i}"], RL.T @ F[i]) for i in range(1, n + 1))
        out["R_L"] = -(euler(x["J_L"], OL, OdL) + torque)
        return out

    return f


ENTRIES = [
    ("planar_pendulum", {}, planar_pendulum),
    ("double_planar_pendulum", {}, double_planar_pendulum),
    ("cart_planar_pendulum", {}, cart_planar_pendulum),
    ("spacecraft_reaction_wheels_s1", {}, spacecraft_reaction_wheels_s1),
    ("spacecraft_solar_arrays", {}, spacecraft_solar_arrays),
    *[("multi_quadrotor_hose", {"n_links": n}, hose(n)) for n in (1, 2, 3, 4)],
]


NEW_ENTRIES = [
    ("point_mass", {}, point_mass),
    ("spring_mass", {}, spring_mass),
    ("spherical_pendulum", {}, spherical_pendulum),
    ("particle_on_rotating_sphere", {}, particle_on_rotating_sphere),
    ("double_spherical_pendulum", {}, double_spherical_pendulum),
    *[("n_link_spherical_chain", {"n_links": n}, n_link_spherical_chain(n)) for n in range(1, 6)],
    ("elastic_spherical_pendulum", {}, elastic_spherical_pendulum),
    ("cart_spherical_pendulum", {}, cart_spherical_pendulum),
    ("free_rigid_body", {}, free_rigid_body),
    ("heavy_top", {}, heavy_top),
    ("rigid_pendulum", {}, rigid_pendulum),
    ("double_rigid_pendulum", {}, double_rigid_pendulum),
    *[("n_link_rigid_chain", {"n_links": n}, n_link_rigid_chain(n)) for n in range(1, 6)],
    ("two_rigid_bodies_spherical_joint", {}, two_rigid_bodies_spherical_joint),
    ("rigid_body", {}, rigid_body),
    ("spacecraft_reaction_wheels", {}, spacecraft_reaction_wheels),
    ("biped_3d", {}, biped_3d),
    ("quadrotor", {}, quadrotor),
    ("quadrotor_offset_com", {}, quadrotor_offset_com),
    ("quadrotor_suspended_load", {}, quadrotor_suspended_load),
    ("quadrotor_elastic_load", {}, quadrotor_elastic_load),
    *[
        ("quadrotor_flexible_cable", {"n_links": n}, quadrotor_flexible_cable(n))
        for n in range(1, 6)
    ],
    *[
        ("multi_quadrotor_point_mass", {"n_quads": n}, multi_quadrotor_point_mass(n))
        for n in (1, 2, 3)
    ],
    *[
        ("multi_quadrotor_rigid_load", {"n_quads": n}, multi_quadrotor_rigid_load(n))
        for n in (1, 2, 3)
    ],
]

# Displayed equations that are a nonzero multiple of the coefficient compute_eom returns,
# as (entry, variable regex) -> factor(x, variable); the factor multiplies the displayed form.
# Each factor is the constant by which the variation of the configuration enters: for a
# cable q with x_Q = x_L - l q, the coefficient of ξ is -l q × (·); for a link with
# x_i = x_Q + l q + ..., it is l q × (·); for an elastic cable, δl enters as -q δl.
SCALES: dict[tuple[str, str], object] = {
    ("quadrotor_suspended_load", "q"): lambda x, v: x["l"],
    ("quadrotor_elastic_load", "q"): lambda x, v: -x["l"],
    ("quadrotor_elastic_load", "l"): lambda x, v: -1.0,
    ("quadrotor_flexible_cable", r"q\d+"): lambda x, v: x[f"l{v[1:]}"],
    ("multi_quadrotor_point_mass", r"q\d+"): lambda x, v: -x[f"l{v[1:]}"],
    ("multi_quadrotor_rigid_load", r"q\d+"): lambda x, v: -x[f"l{v[1:]}"],
}

NOT_DISPLAYED: list[str] = []  # every library entry is displayed and checked

KEEP = {"e2", "e3"}  # the S¹ axis e2 and gravity along e3 stay fixed; e_c (e1) may tilt


def random_params(system, rng, keep=KEEP):
    """Random admissible parameters: positive scalars, unit vectors where the default is a
    unit vector, other vectors at random, symmetric positive definite matrices where the
    default is symmetric."""
    out = {}
    for k, v in system.params.items():
        name, v = str(k), np.asarray(v, float)
        if name in keep:
            out[name] = v.reshape(3, -1) if v.ndim else float(v)
        elif v.ndim == 0:
            out[name] = float(rng.uniform(0.5, 2.0))
        elif v.ndim == 1:
            u = rng.normal(size=(3, 1))
            out[name] = u / np.linalg.norm(u) if np.isclose(np.linalg.norm(v), 1) else u
        else:
            A = rng.normal(size=(3, 3))
            out[name] = A @ A.T + 0.5 * np.eye(3) if np.allclose(v, v.T) else A
    return out


def relative(a, b):
    a, b = np.atleast_1d(np.asarray(a, float)), np.atleast_1d(np.asarray(b, float))
    return float(np.max(np.abs(a - b)) / max(1.0, np.max(np.abs(b))))


def derived(system):
    """{variable name: equation} of compute_eom, keyed by the variable of each variation."""
    eom = compute_eom(system.L, system.dW, system.variables)
    # each key of eom names its variation, which ends in the variable: \\xi_{q1}, \\delta{s}
    return {re.search(r"\{([^{}]+)\}$", key).group(1): eq for key, (_, eq) in eom.items()}


def project(system, x, vname, *vals):
    """Tangent projection of the equation of an S² variable, identity otherwise."""
    if vname not in {str(v) for v in system.variables.vectors if isinstance(v, S2)}:
        return vals
    q = x.q(vname)[0]
    P = np.eye(3) - q @ q.T
    return tuple(P @ v for v in vals)


def check(name, kwargs, hand, randomize=False, per_eq=None):
    system = S.get(name, **kwargs)
    eqs = derived(system)
    worst = 0.0
    for seed in SEEDS:
        rng = np.random.default_rng(seed)
        params = random_params(system, rng) if randomize else None
        x = State(system, rng, params)
        want = hand(x)
        assert set(want) == set(eqs), (name, sorted(want), sorted(eqs))
        for vname, eq in eqs.items():
            w = want[vname]
            for (ename, pattern), factor in SCALES.items():
                if ename == name and re.fullmatch(pattern, vname):
                    w = factor(x, vname) * w
            got, w = project(system, x, vname, Evaluator(env=dict(x.env))(eq), w)
            r = relative(got, w)
            worst = max(worst, r)
            if per_eq is not None:
                per_eq[vname] = max(per_eq.get(vname, 0.0), r)
    size = "".join(f" (n = {v})" for v in kwargs.values())
    mode = "random params" if randomize else "default params"
    print(f"{name + size:<44} {mode:<15} worst relative difference {worst:.1e}")
    return worst


LISTING_OUTPUT = (
    "(Cross(q,Cross(\\frac{d}{dt}(\\omega_{q}),q))*m*l*l*(-1)+Cross(q,e3)*l*m*g*(-1)+tau)"
)


def check_listing():
    """The spherical-pendulum listing of the Examples section prints LISTING_OUTPUT."""
    from geomech import Dot, SystemVariables, Vector, getScalars

    m, g, l = getScalars("m g l", attr=["Constant"])
    e3 = Vector("e3", attr=["Constant"])
    q, tau = S2("q"), Vector("tau")
    xi = q.get_variation_vector()
    x = l * q
    v = x.t_diff()
    L = 0.5 * m * Dot(v, v) - m * g * Dot(x, e3)
    eom = compute_eom(L, Dot(xi, tau), SystemVariables(vectors=[q]))
    out = str(eom["\\xi_{q}"][1])
    print(f"spherical-pendulum listing output matches: {out == LISTING_OUTPUT}")
    if out != LISTING_OUTPUT:
        raise SystemExit(f"the listing prints {out}")


def main():
    check_listing()
    worst = 0.0
    for name, kwargs, hand in ENTRIES + NEW_ENTRIES:
        per_eq: dict[str, float] = {}
        worst = max(worst, check(name, kwargs, hand, False, per_eq))
        worst = max(worst, check(name, kwargs, hand, True, per_eq))
        bad = {k: f"{v:.1e}" for k, v in per_eq.items() if v > TOL}
        if bad:
            print(f"    equations above {TOL:.0e}: {bad}")
    if NOT_DISPLAYED:
        print(f"not displayed, not checked: {', '.join(NOT_DISPLAYED)}")
    print(f"all entries: {worst:.1e}")
    if worst > TOL:
        raise SystemExit(f"an appendix equation differs from geomech by {worst:.1e}")


if __name__ == "__main__":
    main()
