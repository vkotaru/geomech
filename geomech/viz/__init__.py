"""3-D visualizer and HTML report for systems derived with geomech.

``Scene`` declares rigid bodies, point masses, rods and traces with geomech expressions of
the configuration variables and parameters, evaluated through ``geomech.codegen.lambdify``
on the system's ``StateLayout``. ``Report`` puts the animation of a ``Trajectory``, the
Lagrangian, the equations of motion, the standard form, the verification numbers and the
parameters into one self-contained HTML page (three.js and MathJax from pinned CDN URLs).
``latex_names`` gives readable LaTeX for the leaf names (``J1`` → ``J_1``). See
``docs/visualization.md``.

``import geomech`` does not import this package, and this package does not import scipy.
"""

from geomech.viz.latex import latex_names, pretty_name
from geomech.viz.report import Report
from geomech.viz.scene import (
    Box,
    Cylinder,
    PointMass,
    RigidBody,
    Rod,
    Scene,
    SceneFrames,
    Sphere,
    Trace,
)

__all__ = [
    "Scene",
    "SceneFrames",
    "RigidBody",
    "PointMass",
    "Rod",
    "Trace",
    "Box",
    "Cylinder",
    "Sphere",
    "Report",
    "latex_names",
    "pretty_name",
]
