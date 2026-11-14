"""Physical parameters and initial states for the simulation tests (``tests/test_simulation*.py``).

``SETUPS[name]`` gives, for each of the nine systems of ``tests/vbl_systems.py``, a
parameter dict with positive masses, positive-definite inertias and ``e3 = (0, 0, 1)``,
and an initial state on the manifold (``|q| = 1``, ``q·ω = 0``, ``R ∈ SO(3)``) with
moderate velocities, so that every system moves visibly over 10 s and its energy is far
from zero. ``standard_module(build)``/``linearized_module(build)`` are the generated
modules, cached by ``tests/codegen_systems.py``.
"""

from __future__ import annotations

import numpy as np

from geomech.numeric import expm_so3
from tests import vbl_systems
from tests.codegen_systems import modules

E3 = np.array([0.0, 0.0, 1.0])
J_A = np.array([[1.0, 0.1, 0.0], [0.1, 2.0, 0.2], [0.0, 0.2, 3.0]])
J_B = np.array([[0.5, 0.0, 0.05], [0.0, 0.7, 0.0], [0.05, 0.0, 0.4]])


def unit(v) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v)


def tangent(q, v) -> np.ndarray:
    """*v* projected onto the tangent plane at the unit vector *q*."""
    v = np.asarray(v, dtype=float)
    return v - (v @ q) * q


def rot(v) -> np.ndarray:
    """``expm_so3(v)`` flattened row-major."""
    return expm_so3(np.asarray(v, dtype=float)).reshape(9)


Q0 = unit([np.sin(1.0), 0.0, -np.cos(1.0)])  # 1 rad from hanging (q = −e₃ is down)
W0 = tangent(Q0, [0.3, 1.5, -0.2])

SETUPS = {
    "point_mass_r3": (
        {"m": 1.5},
        np.r_[[0.1, -0.2, 0.3], [1.0, 0.5, -0.4]],
    ),
    "spring_mass_r": (
        {"m": 1.5, "k": 4.0},
        np.r_[[0.4], [0.3]],
    ),
    "free_rigid_body": (
        {"J": J_A},
        np.r_[rot([0.1, 0.2, 0.3]), [0.3, -1.2, 0.8]],
    ),
    "rigid_pendulum_pivot": (
        {"J": J_A, "\\rho": np.array([0.0, 0.1, -0.5]), "e3": E3, "g": 9.81, "m": 1.0},
        np.r_[rot([0.5, 0.3, 0.0]), [0.2, 0.5, -0.3]],
    ),
    "spherical_pendulum_s2": (
        {"e3": E3, "g": 9.81, "l": 1.0, "m": 1.0},
        np.r_[Q0, W0],
    ),
    "quadrotor": (
        {"J": J_A, "e3": E3, "g": 9.81, "m": 2.0},
        np.r_[[0.0, 0.5, 1.0], rot([0.2, -0.1, 0.4]), [0.3, -0.2, 0.5], [0.4, 0.3, -0.6]],
    ),
    "double_rigid_pendulum": (
        {
            "J1": J_A,
            "J2": J_B,
            "\\rho_1": np.array([0.0, 0.0, -0.5]),
            "\\rho_2": np.array([0.05, 0.0, -0.4]),
            "l_1": np.array([0.0, 0.0, -1.0]),
            "m1": 1.0,
            "m2": 0.8,
            "g": 9.81,
            "e3": E3,
        },
        np.r_[rot([0.4, 0.1, 0.0]), rot([-0.3, 0.5, 0.2]), [0.2, 0.3, -0.1], [-0.4, 0.2, 0.3]],
    ),
    "elastic_spherical_pendulum": (
        {"e3": E3, "g": 9.81, "k": 40.0, "l0": 1.0, "m": 1.0},
        np.r_[[1.1], Q0, [0.2], W0],
    ),
    "quadrotor_with_load": (
        {"J": J_A, "e3": E3, "g": 9.81, "l": 1.0, "m_L": 0.3, "m_Q": 1.5},
        np.r_[
            [0.0, 0.2, 0.5],
            Q0 * -1.0,
            rot([0.1, 0.2, -0.3]),
            [0.2, -0.1, 0.3],
            tangent(-Q0, [0.5, -0.8, 0.4]),
            [0.3, -0.2, 0.4],
        ],
    ),
}

BUILDERS = {name: getattr(vbl_systems, name) for name in SETUPS}


def standard_module(name: str):
    """The generated standard-form module (with ``energy``) of the system *name*."""
    return modules(BUILDERS[name])[0]


def linearized_module(name: str):
    """The generated linearization module of the system *name*."""
    return modules(BUILDERS[name])[2]
