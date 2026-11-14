"""The NumPy target: vectors are ``(3,)`` arrays, matrices ``(3, 3)``, scalars floats.

=====================  ==========================================  =========================
Node                   NumPy                                       Note
=====================  ==========================================  =========================
``Add/VAdd/MAdd``      ``a + b + …``
``Mul/SVMul/SMMul``    ``a * s``                                   scalar times anything
``MVMul``              ``M @ v``                                   a row ``vᵀ M`` is ``Mᵀ v``
``MMMul``              ``A @ B``
``VVMul``              inner ``a @ b``, outer ``np.outer(a, b)``   by ``is_row_vector``
``Dot``                ``a @ b``
``Cross``              ``np.cross(a, b)``
``Hat`` / ``Vee``      ``_hat(v)`` / ``_vee(M)``                   helpers emitted once
``Transpose``          ``M.T``; of a vector, the vector itself
zero / identity        ``0.0``, ``np.zeros(3)``, ``np.zeros((3, 3))``, ``np.eye(3)``
=====================  ==========================================  =========================
"""

from __future__ import annotations

from geomech.codegen.emitter import MUL, Code, pystr
from geomech.codegen.layout import Slot
from geomech.core.base.expressions import Matrix, Vector

__all__ = ["NumpyTarget"]

_HELPERS = """
def _arr(x, n, name):
    x = np.asarray(x, dtype=float).reshape(-1)
    if x.size != n:
        raise ValueError(f"{name} must have {n} entries, got {x.size}")
    return x


def _vec(x):
    return np.asarray(x, dtype=float).reshape(3)


def _mat(x):
    return np.asarray(x, dtype=float).reshape(3, 3)


def _hat(v):
    return np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])


def _vee(M):
    return np.array([M[2, 1], M[0, 2], M[1, 0]])


def _check_params(p, names):
    p = {} if p is None else p
    missing = [k for k in names if k not in p]
    if missing:
        raise KeyError(f"missing parameter {missing[0]!r} (expected {list(names)})")
    extra = [k for k in p if k not in names]
    if extra:
        raise ValueError(f"unknown parameter {extra[0]!r} (expected {list(names)})")
"""


def _number(value) -> str:
    if isinstance(value, bool):
        value = int(value)
    if isinstance(value, int):
        return repr(value)
    return repr(float(value))


class NumpyTarget:
    """Writes each operation of the translation table as NumPy source."""

    name = "numpy"
    imports: tuple[str, ...] = ("import numpy as np",)
    helpers = _HELPERS

    # ------------------------------------------------------------------ leaves

    def literal(self, node) -> Code:
        if isinstance(node, Matrix):
            return Code("np.zeros((3, 3))" if node.is_zero else "np.eye(3)", atomic=True)
        if isinstance(node, Vector):
            return Code("np.zeros(3)", atomic=True)
        if node.is_zero:
            return Code("0.0", atomic=True)
        return Code(_number(node.value), atomic=True)

    def convert_array(self, array: str, n: int) -> str:
        return f'{array} = _arr({array}, {n}, "{array}")'

    def unpack(self, slot: Slot) -> str:
        if slot.array == "p":
            conv = {(): "", (3,): "_vec", (3, 3): "_mat"}[slot.shape]
            return f"{conv}(p[{pystr(slot.key)}])" if conv else f"p[{pystr(slot.key)}]"
        if slot.shape == ():
            return f"{slot.array}[{slot.start}]"
        if slot.shape == (3, 3):
            return f"{slot.array}[{slot.start}:{slot.stop}].reshape(3, 3)"
        return f"{slot.array}[{slot.start}:{slot.stop}]"

    # ------------------------------------------------------------------ operations

    def scale(self, a: Code, s: Code) -> Code:
        return Code(f"{a.wrap(MUL)} * {s.operand()}", prec=MUL)

    def matmul(self, a: Code, b: Code) -> Code:
        return Code(f"{a.wrap(MUL)} @ {b.operand()}", prec=MUL)

    matvec = matmul
    dot = matmul

    def outer(self, a: Code, b: Code) -> Code:
        return Code(f"np.outer({a.text}, {b.text})", atomic=True)

    def cross(self, a: Code, b: Code) -> Code:
        return Code(f"np.cross({a.text}, {b.text})", atomic=True)

    def call(self, fn: str, a: Code) -> Code:
        return Code(f"{fn}({a.text})", atomic=True)

    def transpose_matrix(self, a: Code) -> Code:
        return Code(f"{a.operand()}.T", atomic=True)

    # ------------------------------------------------------------------ assembly

    def flatten_matrix(self, a: Code) -> Code:
        """Row-major 9-vector of a 3×3 matrix."""
        return Code(f"{a.operand()}.reshape(9)", atomic=True)

    def check_params(self, names_const: str) -> str:
        return f"_check_params(p, {names_const})"

    def stack(self, n: int, parts) -> list[str]:
        """Statements building the flat vector ``out`` of length *n* from
        ``(start, size, code)`` parts; the last statement returns it."""
        lines = [f"out = np.zeros({n})"]
        for start, size, code in parts:
            lines.append(f"out[{_index(start, size)}] = {code.text}")
        return lines + ["return out"]

    def block_matrix(self, row_spans, col_spans, blocks) -> list[str]:
        """Statements building and returning the matrix ``out`` whose block rows and
        columns are the ``(start, size)`` spans, from ``{(i, j): code}``."""
        n_rows = sum(n for _, n in row_spans)
        n_cols = sum(n for _, n in col_spans)
        lines = [f"out = np.zeros(({n_rows}, {n_cols}))"]
        for (i, j), code in sorted(blocks.items()):
            ri, ci = _index(*row_spans[i]), _index(*col_spans[j])
            lines.append(f"out[{ri}, {ci}] = {code.text}")
        return lines + ["return out"]


def _index(start: int, size: int) -> str:
    return str(start) if size == 1 else f"{start}:{start + size}"
