"""The systems of ``geomech.systems`` whose equations Lee et al. give in closed form agree
with those equations numerically (spec ``specs/systems-library.md``, verification item 2).

Each test evaluates the library's equation ``E`` (``compute_eom``: ``E = −(Euler–Lagrange)
+ forcing``, ``docs/conventions.md``) and the book's equation, rewritten in the same
convention, at ``N_STATES`` random admissible states (``sample_state``: ``|q| = 1``,
``ω ⊥ q``, ``ω̇ ⊥ q``, ``R ∈ SO(3)``) with the system's default parameters, and compares
them on the tangent plane for S² (``(I − qqᵀ)E``) and in full elsewhere. The hand
derivation of each rewriting is in the test's docstring.
"""

from __future__ import annotations

import functools

import numpy as np
import pytest

from geomech import compute_eom
from geomech.numeric import Evaluator, hat, sample_state
from geomech.systems import get

N_STATES = 5
TOL = 1e-10
#: the worst relative error of each comparison, reported by ``test_report_book_errors``
WORST: dict[str, float] = {}


@functools.cache
def derived(name, **kwargs):
    s = get(name, **kwargs)
    return s, compute_eom(s.L, s.dW, s.variables)


def states(name, seed=0, **kwargs):
    """``(env, value)`` pairs: ``value(key_or_symbol)`` returns a flat numpy value."""
    s, eom = derived(name, **kwargs)
    rng = np.random.default_rng(seed)
    for _ in range(N_STATES):
        env = sample_state(s.variables, list(s.inputs), rng)
        for k, v in s.params.items():
            env[k] = (
                float(v)
                if np.isscalar(v)
                else np.asarray(v, float).reshape((3, 1) if np.ndim(v) == 1 else (3, 3))
            )
        ev = Evaluator(env=env, rng=rng)

        def value(x, ev=ev):
            key = x if isinstance(x, str) else str(x)
            val = ev.env[key] if key in ev.env else ev(x)
            arr = np.asarray(val, dtype=float)
            return arr.reshape(3) if arr.shape == (3, 1) else arr

        yield eom, value


def lib_equation(eom, var):
    return eom[
        str(var.get_variation_vector() if hasattr(var, "get_variation_vector") else var.delta())
    ][1]


def projector(q):
    return np.eye(3) - np.outer(q, q)


def assert_same(what, got, ref, P=None):
    P = np.eye(len(np.atleast_1d(ref))) if P is None else P
    err = float(np.max(np.abs(P @ (got - ref)))) / max(1.0, float(np.max(np.abs(P @ ref))))
    WORST[what] = max(WORST.get(what, 0.0), err)
    assert err < TOL, f"{what}: library {got} vs book {ref} (relative error {err:.2e})"


def w_s2(q):
    """``(q, ω, ω̇)`` symbols of an S² point."""
    w = q.get_tangent_vector()
    return q, w, w.t_diff()


# ---------------------------------------------------------------------------
# S²
# ---------------------------------------------------------------------------


def test_spherical_pendulum_eq_5_38():
    """Lee et al. eq. 5.38 (from eq. 5.14): ``m l² ω̇ + m g l S(q) e₃ = 0``.

    With ``δW = ξ·f`` the library's equation is ``E = −(m l² ω̇ + m g l q × e₃) + f``
    on the tangent plane of q.
    """
    for eom, val in states("spherical_pendulum"):
        s = derived("spherical_pendulum")[0]
        q, w, wd = w_s2(s.symbols["q"])
        qv, wdv = val(q), val(wd)
        m, g, l = val("m"), val("g"), val("l")
        book = m * l * l * wdv + m * g * l * np.cross(qv, val("e3"))
        ref = -book + val(s.symbols["f"])
        assert_same("spherical_pendulum", val(lib_equation(eom, q)), ref, projector(qv))


def test_particle_on_rotating_sphere_eq_5_43():
    """Lee et al. eq. 5.43, in (q, q̇) form: ``EL = m L² q̈ + m L² |q̇|² q
    + 2 m L² (I − qqᵀ) S(Ωe₃) q̇ − m L² Ω² (I − qqᵀ) C q = 0``.

    The library varies ``q`` by ``δq = ξ × q`` and returns the coefficient of ξ. The
    (q, q̇) principle gives ``δS = ∫ δq·(−EL + F) dt`` for the force ``F`` with
    ``ξ·f = δq·F``, and ``δq·v = (ξ × q)·v = ξ·(q × v)``, so
    ``E = −q × EL + f`` on the tangent plane, with ``q̇ = ω × q`` and
    ``q̈ = ω̇ × q + ω × (ω × q)``. ``C = diag(1, 1, 0)`` and ``e₃ = (0, 0, 1)`` as in the book.
    """
    for eom, val in states("particle_on_rotating_sphere"):
        s = derived("particle_on_rotating_sphere")[0]
        q, w, wd = w_s2(s.symbols["q"])
        qv, wv, wdv = val(q), val(w), val(wd)
        m, L, Om = val("m"), val("L"), val("\\Omega")
        C, e3 = val("C"), val("e3")
        qd = np.cross(wv, qv)
        qdd = np.cross(wdv, qv) + np.cross(wv, qd)
        P = projector(qv)
        EL = (
            m * L**2 * qdd
            + m * L**2 * (qd @ qd) * qv
            + 2 * m * L**2 * P @ hat(Om * e3) @ qd
            - m * L**2 * Om**2 * P @ C @ qv
        )
        ref = -np.cross(qv, EL) + val(s.symbols["f"])
        assert_same("particle_on_rotating_sphere", val(lib_equation(eom, q)), ref, P)


def _eq_5_71(val, s):
    """Rows of Lee et al. eq. 5.71 (double spherical pendulum, (q, ω) form):

    ``(m₁+m₂)L₁² ω̇₁ + m₂L₁L₂ Sᵀ(q₁)S(q₂) ω̇₂ + m₂L₁L₂ |ω₂|² S(q₂) q₁ + (m₁+m₂) g L₁ S(q₁) e₃``,
    ``m₂L₁L₂ Sᵀ(q₂)S(q₁) ω̇₁ + m₂L₂² ω̇₂ + m₂L₁L₂ |ω₁|² S(q₁) q₂ + m₂ g L₂ S(q₂) e₃``.
    """
    q1, w1, wd1 = (val(x) for x in w_s2(s.symbols["qs"][0]))
    q2, w2, wd2 = (val(x) for x in w_s2(s.symbols["qs"][1]))
    m1, m2, L1, L2, g, e3 = (val(k) for k in ("m1", "m2", "l1", "l2", "g", "e3"))
    S1, S2_ = hat(q1), hat(q2)
    row1 = (
        (m1 + m2) * L1**2 * wd1
        + m2 * L1 * L2 * S1.T @ S2_ @ wd2
        + m2 * L1 * L2 * (w2 @ w2) * S2_ @ q1
        + (m1 + m2) * g * L1 * S1 @ e3
    )
    row2 = (
        m2 * L1 * L2 * S2_.T @ S1 @ wd1
        + m2 * L2**2 * wd2
        + m2 * L1 * L2 * (w1 @ w1) * S1 @ q2
        + m2 * g * L2 * S2_ @ e3
    )
    return (q1, row1), (q2, row2)


@pytest.mark.parametrize(
    "name, kwargs",
    [("double_spherical_pendulum", {}), ("n_link_spherical_chain", {"n_links": 2})],
)
def test_double_spherical_pendulum_eq_5_71(name, kwargs):
    """Lee et al. eq. 5.71, rows ``Rᵢ = 0``; the library's ``Eᵢ = −Rᵢ + fᵢ`` on the tangent
    plane of qᵢ (``δW = Σ ξᵢ·fᵢ``)."""
    for eom, val in states(name, **kwargs):
        s = derived(name, **kwargs)[0]
        for (qv, row), q, f in zip(_eq_5_71(val, s), s.symbols["qs"], s.symbols["inputs"]):
            assert_same(name, val(lib_equation(eom, q)), -row + val(f), projector(qv))


# ---------------------------------------------------------------------------
# SO(3) and SE(3)
# ---------------------------------------------------------------------------


def _euler(val, J, Om, Omd):
    return J @ Omd + np.cross(Om, J @ Om)


def test_free_rigid_body_eq_6_16():
    """Lee et al. eq. 6.16: ``J Ω̇ + Ω × JΩ = 0``; with ``δW = η·M``,
    ``E = −(J Ω̇ + Ω × JΩ) + M``."""
    for eom, val in states("free_rigid_body"):
        s = derived("free_rigid_body")[0]
        R = s.symbols["R"]
        Om = R.get_tangent_vector()
        ref = -_euler(val, val("J"), val(Om), val(Om.t_diff())) + val(s.symbols["M"])
        assert_same("free_rigid_body", val(lib_equation(eom, R)), ref)


def test_rigid_pendulum_eq_6_20():
    """Lee et al. eq. 6.20 (eq. 6.8 for the potential ``U = −m g e₃ᵀRρ₀``):
    ``J Ω̇ + Ω × JΩ − m g ρ₀ × Rᵀe₃ = 0``, ρ₀ from the centre of mass to the pivot. The
    library's ρ goes from the pivot to the centre of mass, ``ρ₀ = −ρ``:
    ``E = −(J Ω̇ + Ω × JΩ + m g ρ × Rᵀe₃) + M``."""
    for eom, val in states("rigid_pendulum"):
        s = derived("rigid_pendulum")[0]
        R = s.symbols["R"]
        Om = R.get_tangent_vector()
        rho0 = -val("\\rho")
        book = _euler(val, val("J"), val(Om), val(Om.t_diff()))
        book = book - val("m") * val("g") * np.cross(rho0, val(R).T @ val("e3"))
        assert_same("rigid_pendulum", val(lib_equation(eom, R)), -book + val(s.symbols["M"]))


def test_heavy_top_eq_6_20():
    """Lee et al. eq. 6.20 with the centre of mass at ``l e₃`` from the pivot in the body
    frame, ``ρ₀ = −l e₃``: ``E = −(J Ω̇ + Ω × JΩ + m g l e₃ × Rᵀe₃) + M``."""
    for eom, val in states("heavy_top"):
        s = derived("heavy_top")[0]
        R = s.symbols["R"]
        Om = R.get_tangent_vector()
        rho0 = -val("l") * val("e3")
        book = _euler(val, val("J"), val(Om), val(Om.t_diff()))
        book = book - val("m") * val("g") * np.cross(rho0, val(R).T @ val("e3"))
        assert_same("heavy_top", val(lib_equation(eom, R)), -book + val(s.symbols["M"]))


def test_rigid_body_eqs_7_39_7_40():
    """Lee et al. eqs. 7.39–7.40: ``J Ω̇ + Ω × JΩ = 0``, ``m ẍ + m g e₃ = 0``; with
    ``δW = δx·F + η·M``: ``E_R = −(J Ω̇ + Ω × JΩ) + M``, ``E_x = −(m ẍ + m g e₃) + F``."""
    for eom, val in states("rigid_body"):
        s = derived("rigid_body")[0]
        R, x = s.symbols["R"], s.symbols["x"]
        Om = R.get_tangent_vector()
        ref_R = -_euler(val, val("J"), val(Om), val(Om.t_diff())) + val(s.symbols["M"])
        xdd = val(x.t_diff().t_diff())
        ref_x = -(val("m") * xdd + val("m") * val("g") * val("e3")) + val(s.symbols["F"])
        assert_same("rigid_body", val(lib_equation(eom, R)), ref_R)
        assert_same("rigid_body", val(lib_equation(eom, x)), ref_x)


def test_spacecraft_reaction_wheels_eqs_9_61_9_62():
    """Lee et al. eqs. 9.61–9.62 (symmetric wheels, spin axes through their centres of
    mass, no translation): ``J' Ω̇ + Σ αᵢ sᵢ ω̇ᵢ + Ω × (J'Ω + Σ αᵢ sᵢ ωᵢ) = M`` and
    ``αᵢ sᵢ·Ω̇ + αᵢ ω̇ᵢ = τᵢ``, with ``J'`` the inertia of the whole spacecraft (wheels
    locked), αᵢ the axial wheel inertia and ωᵢ the spin rate relative to the bus.

    Here ``sᵢ = bᵢ``, ``αᵢ = Iᵢ``, ``ωᵢ = θ̇ᵢ`` and ``J' = J_b + Σ Iᵢ bᵢbᵢᵀ`` (``J_b`` holds
    everything but the wheels' axial inertias); the bus gets no external moment (``M = 0``:
    the motor torques are internal and appear only in the wheel equations). So
    ``E_R = −(J' Ω̇ + Σ Iᵢ bᵢ θ̈ᵢ + Ω × (J'Ω + Σ Iᵢ bᵢ θ̇ᵢ))`` and
    ``E_θᵢ = −(Iᵢ bᵢ·Ω̇ + Iᵢ θ̈ᵢ) + τᵢ``."""
    name = "spacecraft_reaction_wheels"
    for eom, val in states(name):
        s = derived(name)[0]
        R = s.symbols["R"]
        Om = R.get_tangent_vector()
        Omv, Omd = val(Om), val(Om.t_diff())
        bs = [val(b) for b in s.symbols["b"]]
        Is = [val(f"I_{i}") for i in (1, 2, 3)]
        ths = s.symbols["theta"]
        wd = [float(val(th.t_diff())) for th in ths]
        wdd = [float(val(th.t_diff().t_diff())) for th in ths]
        Jp = val("J_b") + sum(I * np.outer(b, b) for I, b in zip(Is, bs))
        h = Jp @ Omv + sum(I * w * b for I, w, b in zip(Is, wd, bs))
        ref_R = -(Jp @ Omd + sum(I * a * b for I, a, b in zip(Is, wdd, bs)) + np.cross(Omv, h))
        assert_same(name, val(lib_equation(eom, R)), ref_R)
        for th, tau, I, b, a in zip(ths, s.symbols["tau"], Is, bs, wdd):
            ref = -(I * b @ Omd + I * a) + float(val(tau))
            assert_same(name, np.atleast_1d(val(lib_equation(eom, th))), np.atleast_1d(ref))


#: (system, builder kwargs) → the equation its test above compares with; every system
#: whose ``reference.equation`` is not empty must be listed here with that equation.
COVERED = {
    ("spherical_pendulum", ()): "eq. 5.38",
    ("particle_on_rotating_sphere", ()): "eq. 5.43",
    ("double_spherical_pendulum", ()): "eq. 5.71",
    ("n_link_spherical_chain", (("n_links", 2),)): "eq. 5.71",
    ("free_rigid_body", ()): "eq. 6.16",
    ("rigid_pendulum", ()): "eq. 6.20",
    ("heavy_top", ()): "eq. 6.20",
    ("rigid_body", ()): "eqs. 7.39–7.40",
    ("spacecraft_reaction_wheels", ()): "eqs. 9.61–9.62",
}


def test_every_named_equation_has_a_book_test():
    """``Reference.equation`` names exactly the equations compared in this file (the docs
    table shows it as "Compared with")."""
    from geomech.systems import list_systems

    named = {}
    for info in list_systems():
        sizes = [()] + [((p, n),) for p in info.parameters for n in (1, 2, 3)]
        for kw in sizes:
            eq = get(info.name, **dict(kw)).reference.equation
            if eq:
                named[(info.name, kw)] = eq
    assert named == COVERED


def test_report_book_errors():
    """Runs last in this file: every comparison above ran and stayed below ``TOL``."""
    assert WORST, "the book comparisons did not run"
    for name, err in sorted(WORST.items()):
        print(f"{name:32s} worst relative error {err:.1e}")
        assert err < TOL
