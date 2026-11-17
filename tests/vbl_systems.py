"""Systems used by the variation-based linearization tests.

Each function returns ``System(L, dW, variables, inputs, symbols)`` where ``symbols`` is a
dict of the leaves the tests need (configuration variables, parameters, inputs). The
definitions mirror those pinned in ``tests/test_manifold_pipelines.py`` and
``tests/test_scalar_adversarial.py`` so that the linearizations can be derived by hand from
the pinned equations of motion.

The nine systems are those of ``geomech.systems`` (same symbols, same Lagrangian); each
builder here returns the library's system with the symbols its tests use.
"""

from __future__ import annotations

from typing import NamedTuple

from geomech import Scalar, SystemVariables, systems

HALF = Scalar("0.5", value=0.5, attr=["Constant"])


class System(NamedTuple):
    L: object
    dW: object
    variables: SystemVariables
    inputs: list
    symbols: dict


def _from_library(s, keys):
    """The library system *s* (``geomech.systems``) as a ``System`` with the symbols *keys*."""
    return System(s.L, s.dW, s.variables, list(s.inputs), {k: s.symbols[k] for k in keys.split()})


def point_mass_r3():
    """Point mass in R³: m ẍ = F.  L = ½ m ẋ·ẋ, δW = δx·F."""
    return _from_library(systems.point_mass(), "m x F")


def spring_mass_r():
    """Spring-mass on R: m s̈ + k s = f (``TestSpringMassPipeline``)."""
    return _from_library(systems.spring_mass(), "m k x f")


def free_rigid_body():
    """Free rigid body on SO(3), Lee et al. Eq. 6.16 (``TestFreeRigidBodyPipeline``)."""
    return _from_library(systems.free_rigid_body(), "J M R")


def rigid_pendulum_pivot():
    """Rigid pendulum on SO(3), inertia about the pivot, Lee et al. Eq. 6.8
    (``TestRigidPendulumPivotPipeline``)."""
    return _from_library(systems.rigid_pendulum(), "J rho m g e3 M R")


def spherical_pendulum_s2():
    """Spherical pendulum on S², ω-form, Lee et al. Eq. 5.14
    (``TestSphericalPendulumPipeline``)."""
    return _from_library(systems.spherical_pendulum(), "m g l e3 q M")


def quadrotor():
    """Quadrotor on R³ × SO(3) (Lee, Leok, McClamroch CDC 2010; Wu & Sreenath 2015).

    L = ½ m ẋ·ẋ + ½ Ω·JΩ − m g e₃·x,  δW = f (R e₃)·δx + M·η,  inputs [f (Scalar), M].
    """
    return _from_library(systems.quadrotor(), "m g J e3 x R f M")


def double_rigid_pendulum():
    """Double rigid pendulum on SO(3) × SO(3) (``TestDoubleRigidPendulumMassMatrixExact``)."""
    return _from_library(systems.double_rigid_pendulum(), "J1 J2 rho1 rho2 l1 m1 m2 g R1 R2")


def elastic_spherical_pendulum():
    """Elastic spherical pendulum on R × S² (``tests/test_scalar_adversarial.py``).

    KE = ½ m (l̇² + l² ω·ω), PE = m g l (q·e₃) + ½ k (l − l₀)², δW = f_l δl + ξ·M_q.
    """
    return _from_library(systems.elastic_spherical_pendulum(), "m g k l0 e3 l q f_l M_q")


def quadrotor_with_load():
    """Quadrotor with a cable-suspended point-mass load on R³ × S² × SO(3)
    (Sreenath, Lee, Kumar CDC 2013).

    Load position x_L, cable direction q from load to quadrotor, cable length l,
    x_Q = x_L − l q.  L = ½ m_Q ẋ_Q·ẋ_Q + ½ m_L ẋ_L·ẋ_L + ½ Ω·JΩ − m_Q g e₃·x_Q − m_L g e₃·x_L,
    δW = f (R e₃)·(δx_L − l ξ × q) + M·η.
    """
    return _from_library(systems.quadrotor_suspended_load(), "mQ mL g l J e3 xL q R f M")
