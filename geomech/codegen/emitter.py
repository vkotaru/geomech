"""Translate expression trees into the body of a generated Python function.

``Emitter`` walks a tree once per node class (the same ``match``/``case`` dispatch as the
rest of the engine) and produces Python source through a *target* (``NumpyTarget`` or
``CasadiTarget``) that knows how to write each operation. Shared by both targets:

- **Leaves** are resolved through the ``StateLayout``: configuration, velocity,
  acceleration, input and parameter symbols are bound once at the top of the function,
  numeric literals and the zero/identity constants are inlined. Anything else raises
  ``ValueError`` naming the leaf (the codegen analogue of the stray-symbol check).
- **Rows**: there is one ``VECTOR`` type; every vector, row or column, is generated as the
  same 3-vector, and the row/column structure (``is_row_vector``) only decides between the
  inner and the outer product of a ``VVMul`` and the side of an ``MVMul`` (``vᵀM`` is
  ``Mᵀv``).
- **Common subexpressions**: a subtree that occurs more than once is computed once into a
  local ``_tN``. Subtrees are identified structurally: a leaf by its class and ``str``
  (the engine's identity), an operation by its class and the identities of its operands.
  The identity of each node object is computed once per emission pass (a dict keyed by
  ``id(node)``), so generation is linear in the number of distinct node objects, also on
  trees that share subtrees, whose ``str`` grows exponentially with depth. A long
  expression is also split into locals so that every generated line stays short (and
  ruff-clean).
"""

from __future__ import annotations

import json
from collections import Counter
from dataclasses import dataclass
from typing import Any, Iterable

from geomech.codegen.layout import Slot, StateLayout, base_symbol, is_literal_leaf
from geomech.core.base.expressions import TANGENT_LEAVES, Matrix, Scalar, Vector
from geomech.core.base.types import ExprType
from geomech.core.operations.addition import Add, MAdd, VAdd
from geomech.core.operations.calculus import TimeDerivative, TimeIntegral, Variation
from geomech.core.operations.geometry import Cross, Dot, Hat, Transpose, Vee, is_row_vector
from geomech.core.operations.multiplication import MMMul, Mul, MVMul, SMMul, SVMul, VVMul

__all__ = ["ADD", "ALL_KINDS", "ATOM", "MUL", "Code", "Emitter", "pystr"]

#: Length above which a generated expression is moved into a local of its own.
MAX_INLINE = 72
#: Length above which an operand of another operation is moved into a local first.
MAX_OPERAND = 30

ALL_KINDS = frozenset({"configuration", "velocity", "acceleration", "input", "param"})


def pystr(s: str) -> str:
    """A double-quoted Python string literal (the form ruff format keeps)."""
    return json.dumps(s)


#: Operator precedence of generated code: a sum, a product (``*``/``@``), an atom.
ADD, MUL, ATOM = 1, 2, 3


@dataclass(frozen=True)
class Code:
    """A generated expression: its source and its precedence (``ADD``, ``MUL``, ``ATOM``),
    which decides where it needs parentheses as an operand."""

    text: str
    atomic: bool = False
    prec: int = ADD

    def __post_init__(self):
        if self.atomic:
            object.__setattr__(self, "prec", ATOM)

    def wrap(self, minimum: int) -> str:
        """The text, parenthesized when its precedence is below *minimum*."""
        return self.text if self.prec >= minimum else f"({self.text})"

    def operand(self) -> str:
        return self.wrap(ATOM)


def _children(node):
    return getattr(node, "nodes", None) or ()


class Emitter:
    """Generate the statements that compute a set of expression trees.

    Args:
        layout: the ``StateLayout`` that resolves the leaves.
        target: a ``NumpyTarget`` or ``CasadiTarget``.
        kinds: the slot kinds the expressions may use (``ALL_KINDS`` by default); a leaf
            of another kind raises ``ValueError`` (a mass block with an input in it).
        what: a name for the expressions, used in error messages.
    """

    def __init__(
        self,
        layout: StateLayout,
        target: Any,
        kinds: Iterable[str] = ALL_KINDS,
        what: str = "expression",
    ):
        self.layout = layout
        self.target = target
        self.kinds = frozenset(kinds)
        self.what = what
        self.lines: list[str] = []
        self.used: dict[str, Slot] = {}
        self._cache: dict[int, Code] = {}
        self._counts: Counter = Counter()
        # id(node) -> (structural key, node); the node is kept so that its id is not reused
        self._keys: dict[int, tuple[int, Any]] = {}
        self._intern: dict[tuple, int] = {}
        self._n_temps = 0

    # ------------------------------------------------------------------ public

    def prepare(self, exprs) -> None:
        """Count the subtrees of every expression of the function (for CSE)."""
        for e in exprs:
            self._count(e)

    def emit(self, node) -> Code:
        """The code of *node*; statements it needs are appended to ``self.lines``."""
        key = self._key(node)
        if key in self._cache:
            return self._cache[key]
        code = self._dispatch(node)
        if _children(node) and (self._counts[key] > 1 or len(code.text) > MAX_INLINE):
            code = self.assign(code)
        self._cache[key] = code
        return code

    def assign(self, code: Code) -> Code:
        """Bind *code* to a new local and return the local's name."""
        if code.atomic and code.text.isidentifier():
            return code
        name = f"_t{self._n_temps}"
        self._n_temps += 1
        self.lines.append(f"{name} = {code.text}")
        return Code(name, atomic=True)

    def preamble(self) -> list[str]:
        """The statements that bind every slot the expressions used, in layout order."""
        order = {key: i for i, key in enumerate(s.key for s in self.layout)}
        slots = sorted(self.used.values(), key=lambda s: order[s.key])
        arrays = [a for a in ("state", "a", "u") if any(s.array == a for s in slots)]
        lines = [self.target.convert_array(a, self._array_size(a)) for a in arrays]
        lines += [f"{s.name} = {self.target.unpack(s)}" for s in slots]
        return lines

    # ------------------------------------------------------------------ counting

    def _key(self, node) -> int:
        """The structural identity of *node*: equal for equal subtrees (see the module
        docstring), computed once per node object."""
        hit = self._keys.get(id(node))
        if hit is not None:
            return hit[0]
        children = _children(node)
        if children:
            sig: tuple = (type(node), tuple(self._key(c) for c in children))
        else:
            sig = (type(node), str(node))
        key = self._intern.setdefault(sig, len(self._intern))
        self._keys[id(node)] = (key, node)
        return key

    def _count(self, node):
        key = self._key(node)
        self._counts[key] += 1
        if self._counts[key] == 1:
            for child in _children(node):
                self._count(child)

    def _array_size(self, array):
        return {"state": self.layout.n_state, "a": self.layout.n_accel}.get(
            array, self.layout.n_inputs
        )

    # ------------------------------------------------------------------ dispatch

    def _dispatch(self, node) -> Code:
        t = self.target
        match node:
            case Add() | VAdd() | MAdd():
                return self._sum(node)
            case Mul() | SVMul() | SMMul():
                return t.scale(self._op(node.left), self._op(node.right))
            case MVMul() if is_row_vector(node.left):
                return t.matvec(t.transpose_matrix(self._op(node.right)), self._op(node.left))
            case MVMul() | MMMul():
                return t.matmul(self._op(node.left), self._op(node.right))
            case VVMul() if node.type == ExprType.SCALAR:
                return t.dot(self._op(node.left), self._op(node.right))
            case VVMul():
                return t.outer(self._op(node.left), self._op(node.right))
            case Dot():
                return t.dot(self._op(node.left), self._op(node.right))
            case Cross():
                return t.cross(self._op(node.left), self._op(node.right))
            case Hat():
                return t.call("_hat", self._op(node.expr))
            case Vee():
                return t.call("_vee", self._op(node.expr))
            case Transpose():
                return self._transpose(node)
            case TimeDerivative():
                return self._time_derivative(node)
            case Variation() | TimeIntegral():
                raise ValueError(
                    f"codegen: {self.what} contains the {type(node).__name__} '{node}'; only"
                    " configuration, velocity, acceleration, input and parameter symbols can"
                    " be generated (a variation or an integral has no numeric value here)"
                )
            case Scalar() | Vector() | Matrix():
                return self._leaf(node)
            case _:
                raise TypeError(
                    f"codegen: no translation for the node class {type(node).__name__} '{node}'"
                )

    def _op(self, node) -> Code:
        """The code of an operand, moved into a local when it is long (short lines)."""
        code = self.emit(node)
        if not code.atomic and len(code.text) > MAX_OPERAND:
            code = self.assign(code)
        return code

    def _sum(self, node) -> Code:
        """``a + b + …``; a long sum is accumulated in pieces that each fit on a line."""
        codes = [self._op(n) for n in node.nodes]
        parts = [codes[0].wrap(ADD)] + [c.wrap(MUL) for c in codes[1:]]
        text = parts[0]
        for part in parts[1:]:
            if len(text) + 3 + len(part) > MAX_INLINE:
                text = self.assign(Code(text)).text
            text = f"{text} + {part}"
        return Code(text)

    def _transpose(self, node) -> Code:
        inner = self._op(node.expr)
        if node.expr.type == ExprType.MATRIX:
            return self.target.transpose_matrix(inner)
        return inner  # a scalar, or the row vᵀ of a column: one 3-vector either way

    # ------------------------------------------------------------------ leaves

    def _time_derivative(self, node) -> Code:
        base = base_symbol(node)
        if base is not None and self.layout.slot(node) is not None:
            return self._slot(node)
        if base is not None:
            raise ValueError(
                f"codegen: {self.what} contains '{node}', the time derivative of the"
                f" leaf '{base}', which is neither a velocity nor an acceleration of"
                " the layout (list a prescribed rate in params)"
            )
        raise ValueError(
            f"codegen: {self.what} contains '{node}', a TimeDerivative of the non-leaf"
            f" expression '{node.expr}'; only the time derivative of a configuration"
            " variable or of a velocity can be generated (expand the derivative first)"
        )

    def _leaf(self, node) -> Code:
        if is_literal_leaf(node):
            return self.target.literal(node)
        if self.layout.slot(node) is None and isinstance(node, TANGENT_LEAVES):
            raise ValueError(
                f"codegen: {self.what} contains '{node}', a tangent or variation vector that"
                " is not a velocity of the layout; only configuration, velocity,"
                " acceleration, input and parameter symbols can be generated (a variation"
                " vector such as ξ or η, the Variation of a manifold point, has no numeric"
                " value here, and is not a parameter)"
            )
        if self.layout.slot(node) is None:
            raise ValueError(
                f"codegen: {self.what} contains the symbol '{node}' ({type(node).__name__}),"
                " which is not a configuration, velocity, acceleration, input or parameter"
                " of the layout; add it to params (a constant) or inputs"
            )
        return self._slot(node)

    def _slot(self, node) -> Code:
        slot = self.layout.slot(node)
        assert slot is not None
        if slot.symbol.type != node.type:
            raise ValueError(
                f"codegen: '{node}' is a {node.type.name} in the {self.what} but the layout's"
                f" {slot.kind} '{slot.key}' is a {slot.symbol.type.name}"
            )
        if slot.kind not in self.kinds:
            allowed = ", ".join(sorted(self.kinds))
            raise ValueError(
                f"codegen: {self.what} contains the {slot.kind} '{slot.key}'; it may only"
                f" depend on: {allowed}"
            )
        self.used[slot.key] = slot
        return Code(slot.name, atomic=True)
