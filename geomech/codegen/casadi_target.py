"""The CasADi target: the NumPy translation table written with ``casadi`` operations.

Vectors are ``3×1`` columns, matrices ``3×3``; ``ca.mtimes`` for matrix products,
``ca.dot``, ``ca.cross``, ``ca.vertcat``/``ca.horzcat`` for ``_hat``/``_vee`` and the block
assembly. The generated functions accept NumPy arrays or CasADi expressions (``SX``,
``MX``, ``DM``) for ``state``, ``a``, ``u`` and the values of ``p``, and return CasADi
objects (``DM`` for numeric arguments), so ``ca.Function`` and ``ca.jacobian`` apply to
them directly. The generated module imports ``casadi`` (extra ``geomech[casadi]``);
``geomech`` itself does not.
"""

from __future__ import annotations

from geomech.codegen.emitter import Code
from geomech.codegen.layout import Slot
from geomech.codegen.numpy_target import NumpyTarget

__all__ = ["CasadiTarget", "require_casadi"]

_HELPERS = """
def _is_ca(x):
    return isinstance(x, (ca.SX, ca.MX, ca.DM))


def _arr(x, n, name):
    if not _is_ca(x):
        x = ca.DM(np.asarray(x, dtype=float).reshape(-1, 1))
    if x.numel() != n:
        raise ValueError(f"{name} must have {n} entries, got {x.numel()}")
    return ca.reshape(x, n, 1)


def _vec(x):
    return _arr(x, 3, "parameter")


def _mat(x):
    if _is_ca(x) and x.shape == (3, 3):
        return x
    if _is_ca(x):
        return ca.reshape(x, 3, 3).T
    return ca.DM(np.asarray(x, dtype=float).reshape(3, 3))


def _hat(v):
    return ca.vertcat(
        ca.horzcat(0, -v[2], v[1]),
        ca.horzcat(v[2], 0, -v[0]),
        ca.horzcat(-v[1], v[0], 0),
    )


def _vee(M):
    return ca.vertcat(M[2, 1], M[0, 2], M[1, 0])


def _check_params(p, names):
    p = {} if p is None else p
    missing = [k for k in names if k not in p]
    if missing:
        raise KeyError(f"missing parameter {missing[0]!r} (expected {list(names)})")
    extra = [k for k in p if k not in names]
    if extra:
        raise ValueError(f"unknown parameter {extra[0]!r} (expected {list(names)})")
"""


def require_casadi():
    """Import ``casadi`` or raise ``ImportError`` saying how to install it."""
    try:
        import casadi
    except ImportError as err:  # pragma: no cover - depends on the environment
        raise ImportError(
            "geomech.codegen: the CasADi target needs the optional dependency 'casadi';"
            " install it with `pip install geomech[casadi]` (or `pip install casadi`)"
        ) from err
    return casadi


class CasadiTarget(NumpyTarget):
    """Writes each operation of the translation table as CasADi source."""

    name = "casadi"
    imports = ("import casadi as ca", "import numpy as np")
    helpers = _HELPERS

    def literal(self, node) -> Code:
        code = super().literal(node).text
        code = {
            "np.zeros((3, 3))": "ca.DM.zeros(3, 3)",
            "np.eye(3)": "ca.DM.eye(3)",
            "np.zeros(3)": "ca.DM.zeros(3, 1)",
        }.get(code, code)
        return Code(code, atomic=True)

    def unpack(self, slot: Slot) -> str:
        if slot.array == "p":
            return super().unpack(slot)
        if slot.shape == (3, 3):
            return f"ca.reshape({slot.array}[{slot.start}:{slot.stop}], 3, 3).T"
        return super().unpack(slot)

    def matmul(self, a: Code, b: Code) -> Code:
        return Code(f"ca.mtimes({a.text}, {b.text})", atomic=True)

    matvec = matmul

    def dot(self, a: Code, b: Code) -> Code:
        return Code(f"ca.dot({a.text}, {b.text})", atomic=True)

    def outer(self, a: Code, b: Code) -> Code:
        return Code(f"ca.mtimes({a.text}, {b.operand()}.T)", atomic=True)

    def cross(self, a: Code, b: Code) -> Code:
        return Code(f"ca.cross({a.text}, {b.text})", atomic=True)

    def flatten_matrix(self, a: Code) -> Code:
        return Code(f"ca.reshape({a.operand()}.T, 9, 1)", atomic=True)

    def stack(self, n: int, parts) -> list[str]:
        if not parts:
            return [f"return ca.DM.zeros({n}, 1)"]
        lines, items, pos = [], [], 0
        for start, size, code in parts:
            if start > pos:
                items.append(f"ca.DM.zeros({start - pos}, 1)")
            items.append(code.text)
            pos = start + size
        if pos < n:
            items.append(f"ca.DM.zeros({n - pos}, 1)")
        lines.append("return ca.vertcat(")
        lines += [f"    {item}," for item in items]
        return lines + [")"]

    def block_matrix(self, row_spans, col_spans, blocks) -> list[str]:
        """``ca.vertcat`` of one ``ca.horzcat`` per block row; a missing block is a zero
        ``DM``, and a scalar row's 3-vector block is transposed to ``1×3``."""
        if not col_spans or not row_spans:
            n_rows, n_cols = sum(n for _, n in row_spans), sum(n for _, n in col_spans)
            return [f"return ca.DM.zeros({n_rows}, {n_cols})"]
        lines = []
        for i, (_, nr) in enumerate(row_spans):
            items = []
            for j, (_, nc) in enumerate(col_spans):
                code = blocks.get((i, j))
                if code is None:
                    items.append(f"ca.DM.zeros({nr}, {nc})")
                elif nr == 1 and nc == 3:
                    items.append(f"{code.operand()}.T")
                else:
                    items.append(code.text)
            lines.append(f"_r{i} = ca.horzcat(")
            lines += [f"    {item}," for item in items]
            lines.append(")")
        rows = ", ".join(f"_r{i}" for i in range(len(row_spans)))
        return lines + [f"return ca.vertcat({rows})"]
