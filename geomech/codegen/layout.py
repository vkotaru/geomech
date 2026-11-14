"""``StateLayout``: how the symbols of a system map to flat numeric arrays.

One layout is the single source of truth for the generated code and its callers:

==========================  ===================  =================  ====================
Variable                    Configuration slot   Velocity slot      Acceleration symbol
==========================  ===================  =================  ====================
``s: Scalar`` (R)           ``s`` (1)            ``ṡ`` (1)          ``s̈``
``x: Vector`` (R³)          ``x`` (3)            ``ẋ`` (3)          ``ẍ``
``q: S2``                   ``q`` (3, unit)      ``ω`` (3, ⊥ q)     ``ω̇``
``R: SO3``                  ``R`` (9, row-major) ``Ω`` (3)          ``Ω̇``
==========================  ===================  =================  ====================

``state`` holds every configuration slot in variable order (``SystemVariables`` order:
scalars, vectors, matrices), then every velocity slot in the same order. ``a`` stacks the
accelerations in the same order, ``u`` the inputs in the given order (scalars 1, vectors 3)
and ``p`` is a ``dict`` keyed by ``str(symbol)`` of every parameter.
"""

from __future__ import annotations

import keyword
import re
from dataclasses import dataclass
from typing import Any, Iterable, Sequence

import numpy as np

from geomech.core.base.expressions import (
    S2,
    SO3,
    TS2,
    TSO3,
    Expr,
    Matrix,
    Scalar,
    Vector,
    is_number_literal,
)
from geomech.core.base.types import ExprType
from geomech.core.operations.calculus import TimeDerivative, TimeIntegral, Variation
from geomech.dynamics.eom import check_variables
from geomech.dynamics.variables import SystemVariables

__all__ = ["Slot", "StateLayout", "base_symbol", "collect_params", "is_literal_leaf"]

#: Names the generated functions use for their own arguments and helpers.
_RESERVED = {"state", "a", "u", "p", "np", "ca", "out", "math"}

_SHAPES = {ExprType.SCALAR: (), ExprType.VECTOR: (3,), ExprType.MATRIX: (3, 3)}


@dataclass(frozen=True)
class Slot:
    """One symbol of a layout.

    ``array`` is the argument that holds it (``"state"``, ``"a"``, ``"u"`` or ``"p"``),
    ``start:stop`` its entries there (``None`` for a parameter), ``shape`` its numeric
    shape (``()``, ``(3,)``, ``(3, 3)``; an SO(3) point occupies 9 row-major entries) and
    ``name`` the Python identifier the generated code binds it to. ``manifold`` is
    ``"R"``, ``"R3"``, ``"S2"`` or ``"SO3"`` for the slots of a configuration variable.
    """

    symbol: Expr
    key: str
    name: str
    kind: str
    array: str
    shape: tuple
    start: int | None = None
    stop: int | None = None
    manifold: str = ""

    @property
    def size(self) -> int:
        return int(np.prod(self.shape)) if self.shape else 1

    @property
    def slice(self) -> slice:
        if self.start is None:
            raise ValueError(f"StateLayout: the parameter '{self.key}' has no slice")
        return slice(self.start, self.stop)

    def to_dict(self) -> dict:
        out: dict[str, Any] = {"symbol": self.key, "name": self.name, "kind": self.kind}
        if self.manifold:
            out["manifold"] = self.manifold
        out["shape"] = list(self.shape)
        if self.start is not None:
            out["slice"] = [self.start, self.stop]
        return out


def is_literal_leaf(node) -> bool:
    """A number, or a zero or identity constant: inlined in the generated code."""
    if is_number_literal(node):
        return True
    if isinstance(node, (Scalar, Vector, Matrix)) and node.is_zero:
        return True
    return isinstance(node, Matrix) and node.flags.is_identity


def _identifier(key: str) -> str:
    """A Python identifier for the printed symbol *key* (``\\omega_{q}`` → ``omega_q``)."""
    name = re.sub(r"[^0-9A-Za-z_]+", "_", key.replace("\\", ""))
    name = re.sub(r"_+", "_", name).strip("_")
    if not name or name[0].isdigit():
        name = "s_" + name
    if keyword.iskeyword(name) or name in _RESERVED:
        name += "_"
    return name


def base_symbol(sym):
    """The leaf under a chain of ``TimeDerivative`` nodes (``None`` for anything else)."""
    while isinstance(sym, TimeDerivative):
        sym = sym.expr
    return sym if isinstance(sym, (Scalar, Vector, Matrix)) else None


def _kinematic_symbols(var):
    """(manifold, velocity symbol, acceleration symbol) of a configuration variable."""
    if isinstance(var, S2):
        w = var.get_tangent_vector()
        return "S2", w, w.t_diff()
    if isinstance(var, SO3):
        W = var.get_tangent_vector()
        return "SO3", W, W.t_diff()
    v = var.t_diff()
    return ("R" if isinstance(var, Scalar) else "R3"), v, v.t_diff()


class StateLayout:
    """How the configuration, velocity, acceleration, input and parameter symbols of a
    system map to the flat arrays ``state``, ``a``, ``u`` and the dict ``p``.

    Args:
        variables: the ``SystemVariables`` of the system (scalars, then vectors, then
            matrices; each matrix an ``SO3``).
        inputs: input leaves (``Scalar``/``Vector``), stacked into ``u`` in this order.
        params: named constants (any ``Scalar``/``Vector``/``Matrix`` leaf that is neither
            a configuration variable nor an input, or the ``TimeDerivative`` of one for a
            prescribed rate), passed as ``p[str(symbol)]``. Numeric literals are inlined
            by the code generator and are not parameters.

    Raises ``TypeError``/``ValueError`` naming the symbol for an invalid entry and for two
    symbols that print alike (``str`` is the identity of a symbol in the engine and the key
    of ``p``).
    """

    def __init__(
        self,
        variables: SystemVariables,
        inputs: Sequence[Expr] = (),
        params: Sequence[Expr] = (),
    ):
        if not isinstance(variables, SystemVariables):
            raise TypeError(
                f"StateLayout: variables must be a SystemVariables, got {type(variables).__name__}"
            )
        check_variables(variables, caller="StateLayout")
        self.variables = variables
        self.inputs = tuple(_as_sequence(inputs, "inputs"))
        self.params = tuple(_as_sequence(params, "params"))
        self._slots: dict[str, Slot] = {}
        self._names: set[str] = set()
        self.configuration: list[Slot] = []
        self.velocity: list[Slot] = []
        self.acceleration: list[Slot] = []
        self._build_state()
        self.input_slots = self._build_inputs()
        self.param_slots = self._build_params()

    # ------------------------------------------------------------------ building

    def _unique_name(self, key: str) -> str:
        base = _identifier(key)
        name, k = base, 2
        while name in self._names:
            name, k = f"{base}_{k}", k + 1
        self._names.add(name)
        return name

    def _add(self, sym, kind, array, start, manifold="", name=None) -> Slot:
        key = str(sym)
        if key in self._slots:
            other = self._slots[key]
            raise ValueError(
                f"StateLayout: the {kind} '{key}' prints like the {other.kind} '{other.key}';"
                " two symbols with one printed name are one symbol to the engine"
            )
        shape = _SHAPES[sym.type]
        size = int(np.prod(shape)) if shape else 1
        stop = None if start is None else start + size
        slot = Slot(
            symbol=sym,
            key=key,
            name=self._unique_name(name or key),
            kind=kind,
            array=array,
            shape=shape,
            start=start,
            stop=stop,
            manifold=manifold,
        )
        self._slots[key] = slot
        return slot

    def _build_state(self):
        cfg_vars = self.config_variables
        offset = 0
        for var in cfg_vars:
            manifold = _kinematic_symbols(var)[0]
            slot = self._add(var, "configuration", "state", offset, manifold)
            self.configuration.append(slot)
            offset += slot.size
        self.n_config = offset
        acc_offset = 0
        for var, cfg in zip(cfg_vars, self.configuration):
            manifold, vel, acc = _kinematic_symbols(var)
            vname = None if manifold in ("S2", "SO3") else cfg.name + "_dot"
            vslot = self._add(vel, "velocity", "state", offset, manifold, name=vname)
            self.velocity.append(vslot)
            offset += vslot.size
            aname = vslot.name + ("_dot" if manifold in ("S2", "SO3") else "")
            if manifold in ("R", "R3"):
                aname = cfg.name + "_ddot"
            aslot = self._add(acc, "acceleration", "a", acc_offset, manifold, name=aname)
            self.acceleration.append(aslot)
            acc_offset += aslot.size
        self.n_state = offset
        self.n_accel = acc_offset

    def _build_inputs(self) -> list[Slot]:
        slots, offset = [], 0
        for u in self.inputs:
            if not isinstance(u, (Scalar, Vector)) or isinstance(u, (S2, TS2, TSO3)):
                raise TypeError(
                    f"StateLayout: input '{u}' must be a plain Scalar or Vector leaf, got"
                    f" {type(u).__name__}"
                )
            _check_size(u, "input")
            slot = self._add(u, "input", "u", offset)
            slots.append(slot)
            offset += slot.size
        self.n_inputs = offset
        return slots

    def _build_params(self) -> list[Slot]:
        slots = []
        for sym in self.params:
            base = base_symbol(sym)
            if base is None or isinstance(sym, (Variation, TimeIntegral)):
                raise TypeError(
                    f"StateLayout: parameter '{sym}' must be a Scalar, Vector or Matrix leaf"
                    f" (or the time derivative of one), got {type(sym).__name__}"
                )
            if is_literal_leaf(sym):
                raise ValueError(
                    f"StateLayout: parameter '{sym}' is a numeric literal or a zero/identity"
                    " constant; literals are inlined in the generated code, not parameters"
                )
            _check_size(base, "parameter")
            slots.append(self._add(sym, "param", "p", None))
        return slots

    # ------------------------------------------------------------------ queries

    @property
    def config_variables(self) -> list[Expr]:
        v = self.variables
        return [*v.scalars, *v.vectors, *v.matrices]

    @property
    def slices(self) -> dict[str, slice]:
        """``{str(symbol): slice}`` of every configuration and velocity slot in ``state``."""
        return {s.key: s.slice for s in self.configuration + self.velocity}

    @property
    def accel_slices(self) -> dict[str, slice]:
        """``{str(symbol): slice}`` of every acceleration in ``a``."""
        return {s.key: s.slice for s in self.acceleration}

    @property
    def input_slices(self) -> dict[str, slice]:
        """``{str(symbol): slice}`` of every input in ``u``."""
        return {s.key: s.slice for s in self.input_slots}

    @property
    def param_names(self) -> list[str]:
        return [s.key for s in self.param_slots]

    def slot(self, key) -> Slot | None:
        """The slot of a symbol or printed name, ``None`` if it is not in the layout."""
        return self._slots.get(key if isinstance(key, str) else str(key))

    def __iter__(self):
        return iter(self._slots.values())

    def __repr__(self):
        return (
            f"StateLayout(state={self.n_state}, a={self.n_accel}, u={self.n_inputs},"
            f" params={self.param_names})"
        )

    # ------------------------------------------------------------------ conversion

    def to_dict(self) -> dict:
        """The layout as plain data (lists, dicts, strings, ints), embedded in generated
        modules as ``LAYOUT``."""
        return {
            "n_state": self.n_state,
            "n_config": self.n_config,
            "n_accel": self.n_accel,
            "n_inputs": self.n_inputs,
            "configuration": [s.to_dict() for s in self.configuration],
            "velocity": [s.to_dict() for s in self.velocity],
            "acceleration": [s.to_dict() for s in self.acceleration],
            "inputs": [s.to_dict() for s in self.input_slots],
            "params": [s.to_dict() for s in self.param_slots],
        }

    def state_from_env(self, env: dict) -> np.ndarray:
        """Flat ``state`` from an ``Evaluator`` env (printed-name keys, ``(3, 1)`` columns)."""
        return self._pack(self.configuration + self.velocity, self.n_state, env)

    def accel_from_env(self, env: dict) -> np.ndarray:
        """Flat ``a`` from an ``Evaluator`` env."""
        return self._pack(self.acceleration, self.n_accel, env)

    def inputs_from_env(self, env: dict) -> np.ndarray:
        """Flat ``u`` from an ``Evaluator`` env."""
        return self._pack(self.input_slots, self.n_inputs, env)

    def params_from_env(self, env: dict) -> dict:
        """``p`` from an ``Evaluator`` env (vectors as ``(3,)`` arrays)."""
        out: dict[str, Any] = {}
        for s in self.param_slots:
            val = env[s.key]
            out[s.key] = float(val) if not s.shape else np.asarray(val, float).reshape(s.shape)
        return out

    def env_from_arrays(self, state, a=None, u=None, p=None) -> dict:
        """The ``Evaluator`` env of flat arrays: the inverse of the ``*_from_env`` methods."""
        env: dict[str, Any] = {}
        groups = [(self.configuration + self.velocity, state), (self.acceleration, a)]
        groups.append((self.input_slots, u))
        for slots, arr in groups:
            if arr is None:
                continue
            arr = np.asarray(arr, dtype=float).reshape(-1)
            for s in slots:
                env[s.key] = _column(arr[s.slice], s.shape)
        for key, val in (p or {}).items():
            env[key] = val if np.isscalar(val) else _column(np.asarray(val, float), None)
        return env

    @staticmethod
    def _pack(slots, n, env) -> np.ndarray:
        out = np.zeros(n)
        for s in slots:
            if s.key not in env:
                raise KeyError(f"StateLayout: no value for '{s.key}' in env")
            out[s.slice] = np.asarray(env[s.key], dtype=float).reshape(-1)
        return out


def _column(vals, shape):
    """Evaluator form of a slot value: float, ``(3, 1)`` column or ``(3, 3)`` matrix."""
    vals = np.asarray(vals, dtype=float)
    if shape == () or vals.size == 1:
        return float(vals.reshape(-1)[0])
    if vals.size == 3:
        return vals.reshape(3, 1)
    return vals.reshape(3, 3)


def _as_sequence(items, what):
    if isinstance(items, (Expr, str)) or not hasattr(items, "__iter__"):
        hint = f"; write {what}=[{items}]" if isinstance(items, Expr) else ""
        raise TypeError(
            f"StateLayout: {what} must be a sequence of symbols, got"
            f" {type(items).__name__} {str(items)!r}{hint}"
        )
    return list(items)


def _check_size(sym, what):
    size = tuple(getattr(sym, "size", ()) or ())
    if isinstance(sym, Vector) and size != (3,):
        raise ValueError(f"StateLayout: {what} '{sym}' has size {size}; only 3-vectors")
    if isinstance(sym, Matrix) and size != (3, 3):
        raise ValueError(f"StateLayout: {what} '{sym}' has size {size}; only 3×3 matrices")


# ---------------------------------------------------------------------------
# Parameter discovery
# ---------------------------------------------------------------------------


def _symbols(expr, out):
    """Every leaf and ``TimeDerivative``-of-leaf of *expr*, keyed by ``str``."""
    if isinstance(expr, TimeDerivative) and base_symbol(expr) is not None:
        out.setdefault(str(expr), expr)
        return
    nodes = getattr(expr, "nodes", None)
    if nodes:
        for n in nodes:
            _symbols(n, out)
    elif isinstance(expr, (Scalar, Vector, Matrix)):
        out.setdefault(str(expr), expr)


def collect_params(
    exprs: Iterable[Expr], variables: SystemVariables, inputs: Sequence[Expr] = ()
) -> list[Expr]:
    """The symbols of *exprs* that a ``StateLayout(variables, inputs)`` does not cover.

    Numeric literals and the zero/identity constants are left out (they are inlined), as
    are the configuration, velocity and acceleration symbols of *variables* and the
    *inputs*. The rest, sorted by printed name, are the parameters to pass as ``params``.
    """
    known = StateLayout(variables, inputs)
    found: dict[str, Expr] = {}
    for e in exprs:
        _symbols(e, found)
    return [
        sym
        for key, sym in sorted(found.items())
        if known.slot(key) is None and not is_literal_leaf(sym)
    ]
