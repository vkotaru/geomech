"""The 3-D scene of a system: rigid bodies, point masses, rods and traces.

Geometry is declared with geomech expressions of the configuration variables and the
parameters, and evaluated with ``geomech.codegen.lambdify`` on the system's
``StateLayout``: the visualizer has no kinematics of its own.

.. code-block:: python

    scene = Scene(layout, params=p)
    scene.add(RigidBody(R1, origin=ZeroVector, shape=Box(0.06, 0.03, 1.0), offset=[0, 0, -0.5]))
    scene.add(PointMass(x_L, radius=0.05))
    scene.add(Rod(x_Q, x_L))
    scene.add(Trace(R1 * l1 + R2 * l2))
    frames = scene.frames(traj)      # positions and rotations per stored state
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Callable, Sequence

import numpy as np

if TYPE_CHECKING:
    from geomech.codegen import StateLayout

__all__ = [
    "Box",
    "Cylinder",
    "Sphere",
    "RigidBody",
    "PointMass",
    "Rod",
    "Trace",
    "Scene",
    "SceneFrames",
]

#: Named colours, resolved by the report from its theme tokens (``--viz-<name>``).
PALETTE = ("link1", "link2", "link3", "link4", "load", "cable", "rotor", "trace")
_DEFAULT_BODY = ("link1", "link2", "link3", "link4")

_ORTHO_TOL = 1e-6


# ---------------------------------------------------------------------------
# Shapes
# ---------------------------------------------------------------------------


def _positive(value, what):
    try:
        v = float(value)
    except (TypeError, ValueError):
        raise TypeError(f"{what} must be a number, got {type(value).__name__}") from None
    if not (np.isfinite(v) and v > 0):
        raise ValueError(f"{what} must be positive and finite, got {value!r}")
    return v


@dataclass(frozen=True)
class Box:
    """A box with edges ``lx``, ``ly``, ``lz`` along the body axes e₁, e₂, e₃.

    Give it a non-square cross-section (``lx != ly`` for a link along e₃) so that a spin
    about the link's own axis is visible.
    """

    lx: float
    ly: float
    lz: float

    def __post_init__(self):
        for name in ("lx", "ly", "lz"):
            object.__setattr__(self, name, _positive(getattr(self, name), f"Box: {name}"))

    @property
    def extent(self) -> float:
        return max(self.lx, self.ly, self.lz)

    def to_dict(self) -> dict:
        return {"kind": "box", "size": [self.lx, self.ly, self.lz]}


@dataclass(frozen=True)
class Cylinder:
    """A cylinder of ``radius`` and ``length`` along the body axis ``"x"``, ``"y"`` or ``"z"``."""

    radius: float
    length: float
    axis: str = "z"

    def __post_init__(self):
        object.__setattr__(self, "radius", _positive(self.radius, "Cylinder: radius"))
        object.__setattr__(self, "length", _positive(self.length, "Cylinder: length"))
        if self.axis not in ("x", "y", "z"):
            raise ValueError(f"Cylinder: axis must be 'x', 'y' or 'z', got {self.axis!r}")

    @property
    def extent(self) -> float:
        return max(2 * self.radius, self.length)

    def to_dict(self) -> dict:
        return {"kind": "cylinder", "radius": self.radius, "length": self.length,
                "axis": self.axis}  # fmt: skip


@dataclass(frozen=True)
class Sphere:
    """A sphere of ``radius``."""

    radius: float

    def __post_init__(self):
        object.__setattr__(self, "radius", _positive(self.radius, "Sphere: radius"))

    @property
    def extent(self) -> float:
        return 2 * self.radius

    def to_dict(self) -> dict:
        return {"kind": "sphere", "radius": self.radius}


_SHAPES = (Box, Cylinder, Sphere)


# ---------------------------------------------------------------------------
# Scene elements
# ---------------------------------------------------------------------------


def _check_color(color):
    if color is None:
        return None
    if not isinstance(color, str):
        raise TypeError(f"color must be a palette name or '#rrggbb', got {type(color).__name__}")
    if color in PALETTE:
        return color
    if len(color) == 7 and color[0] == "#":
        try:
            int(color[1:], 16)
            return color.lower()
        except ValueError:
            pass
    raise ValueError(f"color must be one of {list(PALETTE)} or '#rrggbb', got {color!r}")


@dataclass
class RigidBody:
    """A rigid body drawn as a solid shape with its body frame.

    Args:
        R: the attitude, an ``SO3`` configuration variable (or a matrix expression of the
            configuration variables and parameters that evaluates to a rotation).
        origin: the body's reference point in the inertial frame, a vector expression
            (``ZeroVector`` for a body hinged at the origin).
        shape: a ``Box``, ``Cylinder`` or ``Sphere``, rotated by ``R``.
        offset: where the shape's centre sits in the body frame: the shape is centred at
            ``origin + R offset``. A 3-sequence of numbers or a constant vector expression
            (a parameter such as ``\\rho``); ``None`` is zero.
        axes: draw the body frame (e₁ red, e₂ green, e₃ blue) at the shape's centre.
        color: a palette name (``"link1"`` … ``"link4"``, ``"load"``, ``"rotor"``,
            ``"cable"``, ``"trace"``) or ``"#rrggbb"``; ``None`` picks the next link colour.
        name: a label (defaults to ``str(R)``).
    """

    R: Any
    origin: Any
    shape: Any
    offset: Any = None
    axes: bool = True
    color: str | None = None
    name: str | None = None

    def __post_init__(self):
        if not isinstance(self.shape, _SHAPES):
            raise TypeError(
                f"RigidBody: shape must be a Box, Cylinder or Sphere, got"
                f" {type(self.shape).__name__}"
            )
        self.color = _check_color(self.color)
        if self.name is None:
            self.name = str(self.R)


@dataclass
class PointMass:
    """A point mass drawn as a sphere at ``position`` (a vector expression)."""

    position: Any
    radius: float = 0.05
    color: str | None = None
    name: str | None = None

    def __post_init__(self):
        self.radius = _positive(self.radius, "PointMass: radius")
        self.color = _check_color(self.color)
        if self.name is None:
            self.name = str(self.position)


@dataclass
class Rod:
    """A rod (cable, massless link) between two vector expressions ``a`` and ``b``."""

    a: Any
    b: Any
    radius: float = 0.008
    color: str | None = None
    name: str | None = None

    def __post_init__(self):
        self.radius = _positive(self.radius, "Rod: radius")
        self.color = _check_color(self.color)
        if self.name is None:
            self.name = f"{self.a} - {self.b}"


@dataclass
class Trace:
    """The path of a point (a vector expression), drawn up to the current time."""

    point: Any
    color: str | None = None
    name: str | None = None

    def __post_init__(self):
        self.color = _check_color(self.color)
        if self.name is None:
            self.name = str(self.point)


# ---------------------------------------------------------------------------
# Frames
# ---------------------------------------------------------------------------


@dataclass
class SceneFrames:
    """The evaluated scene, one entry per stored state.

    Attributes:
        t: the times, ``(N,)``.
        bodies: per ``RigidBody``, a dict with ``"position"`` (``(N, 3)``, the body's
            ``origin``), ``"rotation"`` (``(N, 3, 3)``) and ``"center"`` (``(N, 3)``, the
            shape's centre ``origin + R offset``).
        points: per ``PointMass``, its position ``(N, 3)``.
        rods: per ``Rod``, its end points ``(N, 2, 3)``.
        traces: per ``Trace``, the traced point ``(N, 3)``.
    """

    t: np.ndarray
    bodies: list = field(default_factory=list)
    points: list = field(default_factory=list)
    rods: list = field(default_factory=list)
    traces: list = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.t)

    def bounds(self) -> tuple[np.ndarray, np.ndarray]:
        """The smallest and largest inertial coordinates of every evaluated point."""
        pts = [b["position"] for b in self.bodies] + [b["center"] for b in self.bodies]
        pts += list(self.points) + [r.reshape(-1, 3) for r in self.rods] + list(self.traces)
        if not pts:
            return np.zeros(3), np.zeros(3)
        allp = np.concatenate([np.asarray(p).reshape(-1, 3) for p in pts])
        return allp.min(axis=0), allp.max(axis=0)


class _Evaluated:
    """One vector or matrix expression compiled against a layout.

    A bare configuration slot is read from the state directly (vectorised); anything
    else is ``lambdify``'d and called per state.
    """

    def __init__(self, expr, layout, kind: str, owner: str):
        from geomech.codegen import lambdify

        self.expr = expr
        self.shape = (3,) if kind == "vector" else (3, 3)
        self.slice = None
        self.fn: Callable | None = None
        slot = layout.slot(expr) if not getattr(expr, "nodes", None) else None
        if slot is not None and slot.kind == "configuration":
            self.slice = slot.slice
        else:
            self.fn = lambdify(expr, layout)

    def __call__(self, states: np.ndarray, p: dict) -> np.ndarray:
        n = len(states)
        if self.slice is not None:
            return states[:, self.slice].reshape((n, *self.shape)).copy()
        assert self.fn is not None
        out = np.empty((n, *self.shape))
        for i, x in enumerate(states):
            out[i] = np.asarray(self.fn(x, None, p), dtype=float).reshape(self.shape)
        return out


class _Constant:
    """A fixed inertial point (``origin=None`` or three numbers)."""

    def __init__(self, value: np.ndarray):
        self.value = value

    def __call__(self, states: np.ndarray, p: dict) -> np.ndarray:
        return np.tile(self.value, (len(states), 1))


def _numeric3(value, owner, what) -> np.ndarray:
    if value is None:
        return np.zeros(3)
    try:
        arr = np.asarray(value, dtype=float).reshape(-1)
    except (TypeError, ValueError):
        arr = np.zeros(0)
    if arr.size != 3 or not np.all(np.isfinite(arr)):
        raise TypeError(
            f"{owner}: the {what} must be a vector expression, None or 3 finite numbers,"
            f" got {type(value).__name__} {value!r}"
        )
    return arr


def _leaf_kinds(expr, layout, owner: str, what: str):
    """Raise ``ValueError`` naming the first leaf that is not a configuration variable,
    a parameter of *layout* or a literal."""
    from geomech.codegen.layout import base_symbol, is_literal_leaf
    from geomech.core.operations.calculus import TimeDerivative

    def walk(node):
        if isinstance(node, TimeDerivative) and base_symbol(node) is not None:
            leaf = node
        elif getattr(node, "nodes", None):
            for child in node.nodes:
                walk(child)
            return
        else:
            leaf = node
        if is_literal_leaf(leaf):
            return
        slot = layout.slot(leaf)
        if slot is None:
            raise ValueError(
                f"{owner}: the {what} uses '{leaf}', which is neither a configuration"
                f" variable nor a parameter of the layout (params: {layout.param_names});"
                " add it to the StateLayout params or replace it by a number"
            )
        if slot.kind not in ("configuration", "param"):
            raise ValueError(
                f"{owner}: the {what} uses the {slot.kind} '{leaf}'; scene geometry may use"
                " configuration variables and parameters only"
            )

    walk(expr)


def _as_expr(value, owner, what):
    from geomech.core.base.expressions import Expr

    if isinstance(value, Expr):
        return value
    raise TypeError(f"{owner}: {what} must be a geomech expression, got {type(value).__name__}")


class Scene:
    """A set of bodies, point masses, rods and traces drawn from a trajectory.

    Args:
        layout: the ``StateLayout`` the simulated module was generated from.
        params: the parameter dict of the simulation (``p``); the expressions of the
            scene are evaluated with it.

    Every expression may use configuration variables and parameters of *layout* only;
    anything else raises ``ValueError`` naming the leaf when the element is added.
    """

    def __init__(self, layout: StateLayout, params: dict | None = None):
        from geomech.codegen import StateLayout

        if not isinstance(layout, StateLayout):
            raise TypeError(f"Scene: layout must be a StateLayout, got {type(layout).__name__}")
        self.layout = layout
        self.params = dict(params or {})
        self.bodies: list[RigidBody] = []
        self.points: list[PointMass] = []
        self.rods: list[Rod] = []
        self.traces: list[Trace] = []
        self._compiled: dict[int, tuple] = {}

    # ------------------------------------------------------------------ building

    def _vector(self, expr, owner, what) -> _Evaluated | _Constant:
        from geomech.core.base.expressions import Expr
        from geomech.core.base.types import ExprType

        if expr is None or not isinstance(expr, Expr):
            return _Constant(_numeric3(expr, owner, what))
        if expr.type != ExprType.VECTOR:
            raise ValueError(
                f"{owner}: the {what} '{expr}' must be a vector, got {expr.type.name}"
            )
        _leaf_kinds(expr, self.layout, owner, what)
        return _Evaluated(expr, self.layout, "vector", owner)

    def _offset(self, offset, owner):
        from geomech.core.base.expressions import Expr

        if offset is None:
            return np.zeros(3)
        if isinstance(offset, Expr):
            from geomech.core.base.types import ExprType

            if offset.type != ExprType.VECTOR:
                raise ValueError(f"{owner}: offset '{offset}' must be a vector")
            _leaf_kinds(offset, self.layout, owner, "offset")
            if any(
                s.kind == "configuration"
                for s in (self.layout.slot(k) for k in _leaf_names(offset))
                if s is not None
            ):
                raise ValueError(
                    f"{owner}: offset '{offset}' must be constant (parameters only); it is"
                    " a body-frame vector"
                )
            from geomech.codegen import lambdify

            return np.asarray(
                lambdify(offset, self.layout)(np.zeros(self.layout.n_state), None, self.params),
                dtype=float,
            ).reshape(3)
        return _numeric3(offset, owner, "offset")

    def add(self, element):
        """Add a ``RigidBody``, ``PointMass``, ``Rod`` or ``Trace``; returns the element."""
        from geomech.core.base.types import ExprType

        if isinstance(element, RigidBody):
            owner = f"RigidBody '{element.name}'"
            R = _as_expr(element.R, owner, "R")
            if R.type != ExprType.MATRIX:
                raise ValueError(f"{owner}: R '{R}' must be a matrix (an SO3 variable)")
            _leaf_kinds(R, self.layout, owner, "attitude R")
            rot = _Evaluated(R, self.layout, "matrix", owner)
            origin = self._vector(element.origin, owner, "origin")
            offset = self._offset(element.offset, owner)
            if element.color is None:
                element.color = _DEFAULT_BODY[len(self.bodies) % len(_DEFAULT_BODY)]
            self._compiled[id(element)] = (rot, origin, offset)
            self.bodies.append(element)
        elif isinstance(element, PointMass):
            owner = f"PointMass '{element.name}'"
            self._compiled[id(element)] = (self._vector(element.position, owner, "position"),)
            if element.color is None:
                element.color = "load"
            self.points.append(element)
        elif isinstance(element, Rod):
            owner = f"Rod '{element.name}'"
            a = self._vector(element.a, owner, "end point a")
            b = self._vector(element.b, owner, "end point b")
            self._compiled[id(element)] = (a, b)
            if element.color is None:
                element.color = "cable"
            self.rods.append(element)
        elif isinstance(element, Trace):
            owner = f"Trace '{element.name}'"
            self._compiled[id(element)] = (self._vector(element.point, owner, "point"),)
            if element.color is None:
                element.color = "trace"
            self.traces.append(element)
        else:
            raise TypeError(
                "Scene.add: expected a RigidBody, PointMass, Rod or Trace, got"
                f" {type(element).__name__}"
            )
        return element

    # ------------------------------------------------------------------ evaluation

    def _states(self, traj) -> tuple[np.ndarray, np.ndarray]:
        state = np.asarray(getattr(traj, "state", traj), dtype=float)
        if state.ndim == 1:
            state = state[None, :]
        if state.ndim != 2 or state.shape[1] != self.layout.n_state:
            raise ValueError(
                f"Scene.frames: states must have {self.layout.n_state} entries per row"
                f" (the layout's n_state), got shape {state.shape}"
            )
        module = getattr(traj, "module", None)
        lay = getattr(module, "LAYOUT", None)
        if lay is not None:
            mine = self.layout.to_dict()
            if lay.get("configuration") != mine["configuration"]:
                raise ValueError(
                    "Scene.frames: the trajectory's module was generated from another"
                    " layout (configuration slots differ); build the Scene with the same"
                    " StateLayout as the module"
                )
        t = np.asarray(getattr(traj, "t", np.arange(len(state), dtype=float)), dtype=float)
        return t, state

    def frames(self, traj, indices: Sequence[int] | np.ndarray | None = None) -> SceneFrames:
        """The scene evaluated at every stored state of *traj* (a ``Trajectory`` or an
        ``(N, n_state)`` array), or at the rows *indices* of it.

        Raises ``ValueError`` when an attitude is not a rotation (``RᵀR ≠ I`` or
        ``det R ≠ 1`` beyond 1e-6) at some state.
        """
        t, states = self._states(traj)
        if indices is not None:
            idx = np.asarray(indices, dtype=int)
            t, states = t[idx], states[idx]
        p = self.params
        out = SceneFrames(t=t)
        for body in self.bodies:
            rot, origin, offset = self._compiled[id(body)]
            R = rot(states, p)
            err = np.abs(np.einsum("nki,nkj->nij", R, R) - np.eye(3)).max(initial=0.0)
            det = np.linalg.det(R) if len(R) else np.ones(0)
            if err > _ORTHO_TOL or np.any(np.abs(det - 1.0) > _ORTHO_TOL):
                raise ValueError(
                    f"Scene.frames: the attitude '{body.R}' of '{body.name}' is not a rotation"
                    f" (max |RᵀR − I| = {err:.2e}, min det = {det.min():.6f})"
                )
            pos = origin(states, p)
            center = pos + np.einsum("nij,j->ni", R, offset)
            out.bodies.append({"position": pos, "rotation": R, "center": center})
        for pm in self.points:
            (pos,) = self._compiled[id(pm)]
            out.points.append(pos(states, p))
        for rod in self.rods:
            a, b = self._compiled[id(rod)]
            out.rods.append(np.stack([a(states, p), b(states, p)], axis=1))
        for tr in self.traces:
            (pt,) = self._compiled[id(tr)]
            out.traces.append(pt(states, p))
        return out

    def body_offset(self, body: RigidBody) -> np.ndarray:
        """The numeric body-frame offset of *body*'s shape."""
        return self._compiled[id(body)][2].copy()

    def __repr__(self):
        return (
            f"Scene(bodies={len(self.bodies)}, points={len(self.points)},"
            f" rods={len(self.rods)}, traces={len(self.traces)})"
        )


def _leaf_names(expr) -> list[str]:
    out: list[str] = []

    def walk(node):
        if getattr(node, "nodes", None):
            for c in node.nodes:
                walk(c)
        else:
            out.append(str(node))

    walk(expr)
    return out
