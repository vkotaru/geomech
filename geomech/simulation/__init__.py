"""Simulation of generated models on their state manifold, and LQR from ``linearize``.

``accelerations`` solves the standard form for the accelerations (S² rows on the tangent
plane); ``simulate`` integrates with a Runge–Kutta–Munthe-Kaas method (``"rkmk4"``) or
``scipy.integrate.solve_ivp`` and returns a ``Trajectory``; ``lqr`` designs a gain about an
equilibrium from a generated linearization. Simulation consumes only the modules of
``geomech.codegen``. scipy (extra ``geomech[sim]``) is imported only by
``method="solve_ivp"`` and ``lqr``. See ``docs/simulation.md``.
"""

from geomech.simulation.control import (
    LQRController,
    ReducedLinearization,
    lqr,
    reduced_linearization,
    tangent_basis,
)
from geomech.simulation.integrators import Trajectory, constraint_drift, simulate
from geomech.simulation.model import accelerations, check_state, model_of

__all__ = [
    "accelerations",
    "check_state",
    "simulate",
    "Trajectory",
    "constraint_drift",
    "lqr",
    "LQRController",
    "tangent_basis",
    "reduced_linearization",
    "ReducedLinearization",
    "model_of",
]
