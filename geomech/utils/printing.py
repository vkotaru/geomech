"""Printing utilities for expressions and equations of motion.

print_eom(eqns)     — display the EOM dict as LaTeX equations
print_tree(expr)    — display the expression tree structure for debugging
render_eom(eqns)    — compile EOM to PDF via pdflatex
"""

from __future__ import annotations

import logging
import os
import subprocess
import tempfile

logger = logging.getLogger("geomech")

# ---------------------------------------------------------------------------
# EOM printer
# ---------------------------------------------------------------------------


def _eom_line(var_vec, eqn):
    """``\\int ξ \\cdot \\left( E \\right) \\, dt = 0`` for one equation, through ``to_latex``."""
    rhs = _strip_outer_parens(to_latex(eqn))
    return "\\int " + to_latex(var_vec) + " \\cdot \\left(" + rhs + "\\right) \\, dt = 0"


def print_eom(eqns):
    """Print equations of motion as LaTeX integral equations (via ``to_latex``).

    *eqns* is the dict returned by ``compute_eom``:
        {str(variation_vector): (variation_vector, equation)}
    """
    for key, (var_vec, eqn) in eqns.items():
        print(_eom_line(var_vec, eqn))


# ---------------------------------------------------------------------------
# Tree printer
# ---------------------------------------------------------------------------


def print_tree(expr, indent=0, label="", style="indent"):
    """Print the expression tree structure for debugging.

    Each node shows its class name; each leaf shows its name, value and the flags the
    engine keys on (constant, zero, unit-norm, symmetric, skew-symmetric, orthogonal,
    identity, numeric, manifold) and, for a tangent or variation vector, its manifold
    point (``of=S2 'q'``). When two different point objects with the same name occur in
    the tree, each is numbered in order of first appearance (``of=S2 'q' #1``,
    ``S2('q' #2, ...)``), so they are told apart without printing memory addresses.
    The last line is the ``repr`` of the expression.

    Parameters
    ----------
    style : str
        'indent' (default) — one node per line, children indented below their parent.
        'topdown' — compact graphical top-down tree with / and \\ branches; readable
        for small expressions only.
    """
    if style == "topdown":
        lines, _, _ = _render_topdown(expr)
        print("\n".join(lines))
    else:
        lines = []
        _build_tree(expr, lines, indent, label, _point_labels(expr))
        print("\n".join(lines))
    print()
    print(f"expr: {repr(expr)}")


def tree_str(expr, indent=0, label="", style="topdown"):
    """Return the expression tree as a string (without printing).

    The default here stays ``'topdown'`` (the compact form); ``print_tree`` defaults to
    ``'indent'``, which shows every leaf's flags.
    """
    if style == "topdown":
        lines, _, _ = _render_topdown(expr)
        return "\n".join(lines)
    else:
        lines = []
        _build_tree(expr, lines, indent, label, _point_labels(expr))
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Human-readable repr
# ---------------------------------------------------------------------------


def _leaf_prefix(expr):
    """Type prefix of a leaf in ``repr`` and the compact tree: ``v:``, ``M:``, ``S2:``, ``SO3:``.

    Scalars have none, so ``repr(Scalar("m")) != repr(Vector("m"))``.
    """
    from geomech.core.base.expressions import S2, SO3, Matrix, Vector

    if isinstance(expr, S2):
        return "S2:"
    if isinstance(expr, SO3):
        return "SO3:"
    if isinstance(expr, Vector):
        return "v:"
    if isinstance(expr, Matrix):
        return "M:"
    return ""


def _repr_leaf(expr, decorate=""):
    """``repr`` of a leaf: type prefix, name (with *decorate* appended), ``*`` if constant.

    A numeric leaf (``Number(2)``, ``Scalar("0.5", value=0.5)``) prints its name alone:
    it is constant by construction.
    """
    name = getattr(expr, "name", None)
    if name is None:
        value = getattr(expr, "value", None)
        return str(value) if value is not None else "?"
    if expr.is_numeric:
        return name
    star = "*" if expr.is_constant else ""
    return _leaf_prefix(expr) + name + decorate + star


def repr_str(expr):
    """Return a human-readable string for the expression (used by ``__repr__``).

    Leaves carry their type and constness: ``m`` (scalar), ``v:x`` (vector), ``M:J``
    (matrix), ``S2:q``, ``SO3:R``, and a trailing ``*`` marks a constant (``M:J*``,
    ``g*``); numbers print bare (``(2)``). Products are joined with ``·``, so the ``*``
    is never a multiplication: ``m*·g*·<v:e3*, v:x>``.
    """
    from geomech.core.operations.addition import Add, MAdd, VAdd
    from geomech.core.operations.calculus import TimeDerivative, TimeIntegral, Variation
    from geomech.core.operations.geometry import Cross, Dot, Hat, Transpose, Vee
    from geomech.core.operations.multiplication import MMMul, Mul, MVMul, SMMul, SVMul, VVMul

    nodes = getattr(expr, "nodes", None)

    # Leaf
    if nodes is None or len(nodes) == 0:
        return _repr_leaf(expr)

    # N-ary addition
    if isinstance(expr, (Add, VAdd, MAdd)):
        parts = []
        for i, child in enumerate(nodes):
            s = repr_str(child)
            if i > 0 and not s.startswith("-"):
                parts.append(" + ")
            elif i > 0:
                parts.append(" ")
            parts.append(s)
        return "(" + "".join(parts) + ")"

    # Binary multiplication — check for negative scalar factor
    if isinstance(expr, (Mul, SVMul, SMMul, MVMul, MMMul, VVMul)):
        left = repr_str(expr.left)
        right = repr_str(expr.right)
        # Detect multiplying by -1
        lval = getattr(expr.left, "value", None)
        rval = getattr(expr.right, "value", None)
        if isinstance(lval, (int, float)) and lval == -1:
            return f"-{right}"
        if isinstance(rval, (int, float)) and rval == -1:
            return f"-{left}"
        return f"{left}·{right}"

    # Dot, Cross
    if isinstance(expr, Dot):
        return f"<{repr_str(expr.left)}, {repr_str(expr.right)}>"
    if isinstance(expr, Cross):
        return f"cross({repr_str(expr.left)}, {repr_str(expr.right)})"

    # Hat, Vee, Transpose
    if isinstance(expr, Hat):
        return f"hat({repr_str(expr.expr)})"
    if isinstance(expr, Vee):
        return f"vee({repr_str(expr.expr)})"
    if isinstance(expr, Transpose):
        inner = repr_str(expr.expr)
        if getattr(expr.expr, "nodes", None) and not inner.startswith(("(", "<")):
            inner = f"({inner})"
        return f"{inner}'"

    # Calculus
    if isinstance(expr, Variation):
        return f"δ({repr_str(expr.expr)})"
    if isinstance(expr, TimeDerivative):
        inner = expr.expr
        name = getattr(inner, "name", None)
        # For simple leaf names, put a combining dot above: v:x → v:ẋ
        if not getattr(inner, "nodes", None) and name and name.isalnum() and len(name) <= 10:
            return _repr_leaf(inner, decorate="\u0307")
        return f"d/dt({repr_str(inner)})"
    if isinstance(expr, TimeIntegral):
        return f"∫({repr_str(expr.expr)})dt"

    # Fallback
    return str(expr)


_SHORT_NAMES = {
    "TimeDerivative": "d/dt",
    "TimeIntegral": "∫dt",
    "Variation": "δ",
    "SVMul": "S·V",
    "SMMul": "S·M",
    "MVMul": "M·V",
    "MMMul": "M·M",
    "VVMul": "V·V",
}


def _node_label(expr, compact=False):
    """Get the display label for a node."""
    cls_name = type(expr).__name__
    nodes = getattr(expr, "nodes", None)
    if nodes is None or len(nodes) == 0:
        if compact:
            return _compact_leaf(expr)
        detail = _leaf_detail(expr)
        return f"{cls_name}({detail})"
    if compact:
        return _SHORT_NAMES.get(cls_name, cls_name)
    return cls_name


def _compact_leaf(expr):
    """Short leaf label: just name/value + type prefix + flag markers."""
    name = getattr(expr, "name", None)
    value = getattr(expr, "value", None)
    flags = getattr(expr, "flags", None)

    # Type prefix (scalars get none)
    prefix = _leaf_prefix(expr)

    # Display name
    display = name if name is not None else str(value) if value is not None else "?"

    # Flag suffix: * = constant, numeric values shown inline
    suffix = ""
    if flags:
        is_const = getattr(flags, "is_constant", False)
        is_num = getattr(flags, "is_numeric", False)
        if is_num and value is not None:
            # Pure numeric: just show the value
            display = str(value)
            if is_const:
                suffix = "*"
        elif is_const:
            suffix = "*"

    return prefix + display + suffix


def _render_topdown(expr, gap=3):
    """Render expression as top-down tree with / and \\ branches.

    Returns (lines, root_col, width) where lines is a list of strings,
    root_col is the column of the root node center, and width is the
    total width of the rendered block.
    """
    label = _node_label(expr, compact=True)
    nodes = getattr(expr, "nodes", None)

    # Leaf
    if nodes is None or len(nodes) == 0:
        return [label], len(label) // 2, len(label)

    # Render children
    children = [_render_topdown(child, gap) for child in nodes]

    # Unary node — pipe connector
    if len(children) == 1:
        c_lines, c_root, c_width = children[0]
        label_start = max(0, c_root - len(label) // 2)
        width = max(c_width, label_start + len(label))
        lines = []
        lines.append(" " * label_start + label)
        lines.append(" " * c_root + "|")
        for line in c_lines:
            lines.append(line.ljust(width))
        return lines, label_start + len(label) // 2, width

    # Multiple children — place side by side with gap
    child_blocks = []
    total_width = 0
    for i, (c_lines, c_root, c_width) in enumerate(children):
        if i > 0:
            total_width += gap
        child_blocks.append((c_lines, c_root, c_width, total_width))
        total_width += c_width

    first_root = child_blocks[0][3] + child_blocks[0][1]
    last_root = child_blocks[-1][3] + child_blocks[-1][1]
    center = (first_root + last_root) // 2
    label_start = center - len(label) // 2

    # Shift everything right if label goes negative
    if label_start < 0:
        shift = -label_start
        label_start = 0
        child_blocks = [(cl, cr, cw, co + shift) for cl, cr, cw, co in child_blocks]
        total_width += shift
        first_root += shift
        last_root += shift

    width = max(total_width, label_start + len(label))
    parent_root = label_start + len(label) // 2

    lines = []
    # Parent label
    lines.append(" " * label_start + label)

    # Draw branch lines from parent down to each child
    # For 2 children, draw diagonal lines; for N, single connector row
    child_roots = [co + cr for (_, cr, _, co) in child_blocks]

    # Single connector row with / | \ pointing to each child
    connector = [" "] * width
    for cr in child_roots:
        if cr < parent_root:
            connector[cr] = "/"
        elif cr > parent_root:
            connector[cr] = "\\"
        else:
            connector[cr] = "|"
    lines.append("".join(connector).rstrip())

    # Merge child lines row by row
    max_child_lines = max(len(cb[0]) for cb in child_blocks)
    for row in range(max_child_lines):
        line = [" "] * width
        for c_lines_list, c_root, c_width, c_offset in child_blocks:
            if row < len(c_lines_list):
                text = c_lines_list[row]
                for j, ch in enumerate(text):
                    pos = c_offset + j
                    if pos < width and ch != " ":
                        line[pos] = ch
        lines.append("".join(line).rstrip())

    return lines, parent_root, width


def _point_labels(expr):
    """``{id(point): "#k"}`` for manifold points that share a printed name with another.

    Collects the ``S2``/``SO3`` leaves of *expr* and the points of its tangent and
    variation vectors in order of first appearance. A name used by one object gets no
    label; a name used by several distinct objects numbers them ``#1``, ``#2``, … .
    """
    from geomech.core.base.expressions import S2, SO3, TS2, TSO3

    by_name = {}

    def add(point):
        objs = by_name.setdefault(str(point), [])
        if not any(o is point for o in objs):
            objs.append(point)

    def walk(node):
        nodes = getattr(node, "nodes", None)
        if nodes:
            for n in nodes:
                walk(n)
        elif isinstance(node, (S2, SO3)):
            add(node)
        elif isinstance(node, TS2) and node.S2 is not None:
            add(node.S2)
        elif isinstance(node, TSO3) and node.SO3 is not None:
            add(node.SO3)

    walk(expr)
    return {
        id(o): f"#{k}" for objs in by_name.values() if len(objs) > 1 for k, o in enumerate(objs, 1)
    }


def _build_tree(expr, lines, indent, label, points=None):
    """Recursively build tree lines; *points* is the map of ``_point_labels``."""
    prefix = "  " * indent
    connector = "├─ " if indent > 0 else ""
    tag = f"{label}: " if label else ""

    cls_name = type(expr).__name__
    nodes = getattr(expr, "nodes", None)

    # Leaf node — show name/value
    if nodes is None or len(nodes) == 0:
        detail = _leaf_detail(expr, points)
        lines.append(f"{prefix}{connector}{tag}{cls_name}({detail})")
        return

    # Internal node
    lines.append(f"{prefix}{connector}{tag}{cls_name}")

    for i, child in enumerate(nodes):
        child_label = _child_label(expr, i)
        _build_tree(child, lines, indent + 1, child_label, points)


_LEAF_FLAGS = (
    "is_constant",
    "is_zero",
    "is_ones",
    "is_identity",
    "is_unit_norm",
    "is_symmetric",
    "is_skew_symmetric",
    "is_orthogonal",
    "is_numeric",
    "is_manifold",
)


def _leaf_detail(expr, points=None):
    """Format the detail string for a leaf node.

    Name, value, every flag that is set, and for a ``TS2``/``TSO3`` the manifold point it
    belongs to (``of=SO3 'R'``). *points* (from ``_point_labels``) numbers distinct point
    objects that share a name (``of=SO3 'R' #2``).
    """
    from geomech.core.base.expressions import TS2, TSO3

    points = points or {}
    parts = []
    name = getattr(expr, "name", None)
    if name is not None:
        mark = points.get(id(expr))
        parts.append(f"'{name}'" + (f" {mark}" if mark else ""))
    value = getattr(expr, "value", None)
    if value is not None and not _is_default_value(value):
        parts.append(f"value={value}")
    flags = getattr(expr, "flags", None)
    if flags is not None:
        active = [f for f in _LEAF_FLAGS if getattr(flags, f, False)]
        if active:
            parts.append(", ".join(active))
    parent = None
    if isinstance(expr, TS2):
        parent = expr.S2
    elif isinstance(expr, TSO3):
        parent = expr.SO3
    if parent is not None:
        mark = points.get(id(parent))
        parts.append(f"of={type(parent).__name__} '{parent}'" + (f" {mark}" if mark else ""))
    return ", ".join(parts)


def _is_default_value(value):
    """Check if a value is the default (array of None)."""
    try:
        import numpy as np

        if isinstance(value, np.ndarray) and all(v is None for v in value.flat):
            return True
    except (ImportError, TypeError):
        pass
    return False


def _child_label(expr, index):
    """Return a descriptive label for the i-th child of expr."""
    from geomech.core.operations.geometry import Cross, Dot
    from geomech.core.operations.multiplication import MMMul, Mul, MVMul, SMMul, SVMul, VVMul

    # Binary ops: label left/right
    if isinstance(expr, (Mul, SVMul, SMMul, MVMul, MMMul, VVMul, Dot, Cross)):
        return "L" if index == 0 else "R"

    # Unary ops: label as expr
    if hasattr(expr, "expr") and len(getattr(expr, "nodes", [])) == 1:
        return ""

    # N-ary: label with index
    return str(index)


# ---------------------------------------------------------------------------
# LaTeX conversion
# ---------------------------------------------------------------------------


def _wrap(s):
    return "\\left(" + s + "\\right)"


def _numeric_value(expr):
    """The number of a numeric leaf (``Number(2)``, ``Scalar("0.5", value=0.5)``), else None."""
    if getattr(expr, "nodes", None):
        return None
    val = getattr(expr, "value", None)
    if isinstance(val, (int, float)) and not isinstance(val, bool):
        return val
    return None


def _format_number(val):
    r"""``0.5 → \frac{1}{2}``, ``2.0 → 2``; ``val`` is non-negative."""
    if val == 0.5:
        return r"\frac{1}{2}"
    if isinstance(val, float) and val.is_integer():
        return str(int(val))
    return str(val)


def _coefficient(factors):
    """Split scalar *factors* into ``(negative, coefficient_str, symbolic_factors)``.

    The numeric factors are multiplied into one coefficient, printed first; a coefficient
    of magnitude one prints as nothing (its sign is returned separately).
    """
    coeff = 1
    symbolic = []
    for f in factors:
        val = _numeric_value(f)
        if val is None:
            symbolic.append(f)
        else:
            coeff *= val
    neg = coeff < 0
    mag = abs(coeff)
    return neg, ("" if mag == 1 else _format_number(mag)), symbolic


def _scalar_factors(expr):
    """The factors of a ``Mul`` chain, left to right (a non-``Mul`` is one factor)."""
    from geomech.core.operations.multiplication import Mul

    if isinstance(expr, Mul):
        return _scalar_factors(expr.left) + _scalar_factors(expr.right)
    return [expr]


def _scalar_factor_latex(f, alone):
    """LaTeX of one symbolic scalar factor; a ``Dot`` is wrapped unless it stands alone."""
    from geomech.core.operations.geometry import Dot

    s = to_latex(f)
    if isinstance(f, Dot) and not alone:
        return _wrap(s)
    return s


def _product_latex(scalars, base=None):
    """``[-]coefficient symbol symbol … base`` for scalar *scalars* times *base*.

    *base* is the LaTeX of the vector or matrix factor (``None`` for a pure scalar
    product). A ``Dot`` factor is wrapped whenever anything else multiplies it, so that
    ``(b·c) a`` does not read as ``b·(c a)``. A factor whose LaTeX starts with ``-`` (a
    negated transpose, ``-a^T b``) is wrapped unless it leads an unsigned product, so that
    ``x (-a^T b)`` does not read as ``x - a^T b`` and no ``--`` appears.
    """
    factors = []
    for s in scalars:
        factors.extend(_scalar_factors(s))
    neg, coeff, symbolic = _coefficient(factors)
    alone = base is None and not coeff and len(symbolic) == 1
    parts = [coeff] if coeff else []
    for f in symbolic:
        parts.append(_signed_factor(_scalar_factor_latex(f, alone), parts or neg))
    if base is not None:
        parts.append(_signed_factor(base, parts or neg))
    body = " ".join(parts) if parts else "1"
    return "-" + body if neg else body


# Operands that are wrapped in \left( \right) when they appear inside a Cross, a matrix
# product or a Transpose: products and cross products whose juxtaposition would
# otherwise change the meaning (``a × b × c``, ``R a × b``, ``R a^T``).
def _needs_parens_in_cross(e):
    from geomech.core.operations.geometry import Cross
    from geomech.core.operations.multiplication import MVMul, SVMul, VVMul

    return isinstance(e, (Cross, MVMul, SVMul, VVMul))


def _needs_parens_in_product(e):
    from geomech.core.operations.geometry import Cross
    from geomech.core.operations.multiplication import SMMul, SVMul

    return isinstance(e, (Cross, SVMul, SMMul))


def _is_transpose_atom(e):
    r"""Operands that take ``^T`` directly: a leaf (without its own superscript), ``\dot x``,
    ``\ddot x`` of a leaf, and ``\widehat{·}``."""
    from geomech.core.operations.calculus import TimeDerivative
    from geomech.core.operations.geometry import Hat

    if not getattr(e, "nodes", None):
        return "^" not in to_latex(e)
    if isinstance(e, Hat):
        return True
    if isinstance(e, TimeDerivative):
        inner = e.expr.expr if isinstance(e.expr, TimeDerivative) else e.expr
        return not getattr(inner, "nodes", None)
    return False


def _signed_factor(s, inside):
    """Wrap the LaTeX *s* of a factor if it starts with ``-`` and is *inside* a product
    (not its first factor, or under a sign): ``x \\left(-a^T b\\right)``, not ``x -a^T b``."""
    return _wrap(s) if inside and s.startswith("-") else s


def _operand(e, needs_parens, inside=False):
    s = to_latex(e)
    if needs_parens(e):
        return _wrap(s)
    return _signed_factor(s, inside)


def to_latex(expr):
    r"""Convert an expression tree to a LaTeX string.

    Returns a string suitable for use with IPython.display.Math() or inside a LaTeX
    equation environment. The output is precedence-aware: composite operands of a cross
    product, of a matrix product and of a transpose are wrapped in ``\left( \right)``
    (``a \times \left(b \times c\right)``, ``R^T \left(a \times b\right)``,
    ``\left(R a\right)^T``), a ``Dot`` multiplied by anything is wrapped
    (``\left(b \cdot c\right) a``), and the numeric coefficient of a product is printed
    once, first (``-\frac{1}{2} m v``). Sums are always wrapped.
    """
    from geomech.core.operations.addition import Add, MAdd, VAdd
    from geomech.core.operations.calculus import TimeDerivative, TimeIntegral, Variation
    from geomech.core.operations.geometry import Cross, Dot, Hat, Transpose, Vee
    from geomech.core.operations.multiplication import MMMul, Mul, MVMul, SMMul, SVMul, VVMul

    nodes = getattr(expr, "nodes", None)

    # Leaf
    if nodes is None or len(nodes) == 0:
        name = getattr(expr, "name", None)
        if name is not None:
            val = _numeric_value(expr)
            if val is not None:
                body = _format_number(abs(val))
                return "-" + body if val < 0 else body
            return name
        return "?"

    # N-ary addition
    if isinstance(expr, (Add, VAdd, MAdd)):
        parts = []
        for i, child in enumerate(nodes):
            s = to_latex(child)
            if i > 0 and not s.startswith("-"):
                parts.append(" + ")
            elif i > 0:
                parts.append(" ")
            parts.append(s)
        return _wrap("".join(parts))

    # Scalar * Scalar: one coefficient first, then the symbols
    if isinstance(expr, Mul):
        return _product_latex([expr])

    # Scalar multiples of a vector/matrix: flatten the chain, coefficient first
    if isinstance(expr, (SVMul, SMMul)):
        scalars = []
        base = expr
        while isinstance(base, type(expr)):
            scalars.append(base.right)
            base = base.left
        base_str = to_latex(base)
        neg, coeff, symbolic = _coefficient([f for s in scalars for f in _scalar_factors(s)])
        if (coeff or symbolic) and isinstance(base, Cross):
            base_str = _wrap(base_str)  # m (a × b), not m a × b
        return _product_latex(scalars, base_str)

    # Matrix/vector products (associative: nested products need no parentheses)
    if isinstance(expr, (MVMul, MMMul, VVMul)):
        left = _operand(expr.left, _needs_parens_in_product)
        right = _operand(expr.right, _needs_parens_in_product, inside=True)
        return f"{left} {right}"

    # Dot product
    if isinstance(expr, Dot):
        right = _signed_factor(to_latex(expr.right), True)
        return f"{to_latex(expr.left)} \\cdot {right}"

    # Cross product
    if isinstance(expr, Cross):
        left = _operand(expr.left, _needs_parens_in_cross)
        right = _operand(expr.right, _needs_parens_in_cross, inside=True)
        return f"{left} \\times {right}"

    # Hat (skew-symmetric)
    if isinstance(expr, Hat):
        return f"\\widehat{{{to_latex(expr.expr)}}}"

    # Vee
    if isinstance(expr, Vee):
        return f"\\left({to_latex(expr.expr)}\\right)^\\vee"

    # Transpose: x^T for an atom, \left(...\right)^T for anything composite
    if isinstance(expr, Transpose):
        if isinstance(expr.expr, (SVMul, SMMul)):
            # (s v)^T = s v^T: keep the coefficient outside the transpose
            scalars, base = [], expr.expr
            while isinstance(base, (SVMul, SMMul)):
                scalars.append(base.right)
                base = base.left
            return _product_latex(scalars, to_latex(Transpose(base)))
        inner = to_latex(expr.expr)
        if _is_transpose_atom(expr.expr) and not inner.startswith("-"):
            return f"{inner}^T"
        return f"{_wrap(inner)}^T"

    # Time derivative: \dot / \ddot of a leaf, d/dt(...) of a composite
    if isinstance(expr, TimeDerivative):
        inner = expr.expr
        if isinstance(inner, TimeDerivative) and not getattr(inner.expr, "nodes", None):
            return f"\\ddot{{{to_latex(inner.expr)}}}"
        if not getattr(inner, "nodes", None) or (
            isinstance(inner, Variation) and not getattr(inner.expr, "nodes", None)
        ):
            return f"\\dot{{{to_latex(inner)}}}"
        return f"\\frac{{d}}{{dt}}{_wrap(to_latex(inner))}"

    # Variation
    if isinstance(expr, Variation):
        if getattr(expr.expr, "nodes", None):
            return f"\\delta{_wrap(to_latex(expr.expr))}"
        return f"\\delta {to_latex(expr.expr)}"

    # Time integral
    if isinstance(expr, TimeIntegral):
        return f"\\int {to_latex(expr.expr)} \\, dt"

    # Fallback
    return str(expr)


def display_latex(expr):
    """Display an expression as rendered LaTeX in a Jupyter notebook."""
    from IPython.display import Math, display

    display(Math(to_latex(expr)))


def _strip_outer_parens(s):
    """Remove one ``\\left( ... \\right)`` pair if it wraps the whole string.

    ``\\left(a\\right) \\left(b\\right)`` starts and ends with a parenthesis that are not
    a pair; stripping those would unbalance the LaTeX, so the pair must match.
    """
    open_tok, close_tok = "\\left(", "\\right)"
    if not (s.startswith(open_tok) and s.endswith(close_tok)):
        return s
    depth = 0
    i = 0
    while i < len(s):
        if s.startswith("\\left", i):
            depth += 1
            i += len("\\left")
        elif s.startswith("\\right", i):
            depth -= 1
            i += len("\\right")
            if depth == 0 and i + 1 < len(s):
                return s  # the first group closes before the end
        else:
            i += 1
    return s[len(open_tok) : -len(close_tok)]


def display_eom(eqns):
    """Display equations of motion as rendered LaTeX in a Jupyter notebook."""
    from IPython.display import Math, display

    for key, (var_vec, eqn) in eqns.items():
        latex_str = _strip_outer_parens(to_latex(eqn)) + " = 0"
        display(Math(latex_str))


def _latex_str_key(s):
    """Convert a raw str key (from EOM dicts) to cleaner LaTeX.

    Handles patterns like \\frac{d}{dt}(\\frac{d}{dt}(x)) → \\ddot{x}
    and \\frac{d}{dt}(\\omega_{q}) → \\dot{\\omega_{q}}.
    """
    import re

    # ddot: \frac{d}{dt}(\frac{d}{dt}(X)) → \ddot{X}
    m = re.match(r"^\\frac\{d\}\{dt\}\(\\frac\{d\}\{dt\}\((.+)\)\)$", s)
    if m:
        return f"\\ddot{{{m.group(1)}}}"
    # dot: \frac{d}{dt}(X) → \dot{X}
    m = re.match(r"^\\frac\{d\}\{dt\}\((.+)\)$", s)
    if m:
        return f"\\dot{{{m.group(1)}}}"
    return s


def display_standard_form(sf):
    """Display standard form as M*ddq + G*u + f = 0 in a Jupyter notebook."""
    from IPython.display import Math, display

    for key, eq in sf.items():
        parts = []
        # M * accel terms
        for accel, coeff in eq.M.items():
            coeff_str = to_latex(coeff)
            accel_str = _latex_str_key(accel)
            parts.append(f"{coeff_str} {accel_str}")
        # G * input terms
        for inp, coeff in eq.G.items():
            g_str = to_latex(coeff)
            if g_str in ("I", "1"):
                parts.append(inp)
            else:
                parts.append(f"{g_str} {inp}")
        # f term
        f_str = to_latex(eq.f)
        if f_str not in ("0v", "0"):
            # f is one more summand of the equation: drop its own outer parentheses
            parts.append(_strip_outer_parens(f_str))
        # Join with sign-aware concatenation
        if not parts:
            equation = "0 = 0"
        else:
            result = parts[0]
            for p in parts[1:]:
                if p.startswith("-"):
                    result += f" {p}"
                else:
                    result += f" + {p}"
            equation = result + " = 0"
        display(Math(f"{_latex_str_key(key)}: \\quad {equation}"))


# ---------------------------------------------------------------------------
# PDF renderer
# ---------------------------------------------------------------------------


def eom_to_latex(eqns):
    """Return a complete LaTeX document string for the EOM dict (equations via ``to_latex``)."""
    lines = [
        r"\documentclass[12pt]{article}",
        r"\usepackage{amsmath,amssymb}",
        r"\usepackage[margin=1in]{geometry}",
        r"\begin{document}",
        r"\section*{Equations of Motion}",
    ]
    for key, (var_vec, eqn) in eqns.items():
        lines.append(r"\begin{equation}")
        lines.append(_eom_line(var_vec, eqn))
        lines.append(r"\end{equation}")
    lines.append(r"\end{document}")
    return "\n".join(lines)


def render_eom(eqns, output="eom.pdf", open_pdf=True):
    """Compile EOM to PDF via pdflatex.

    Parameters
    ----------
    eqns : dict
        The dict returned by ``compute_eom``.
    output : str
        Output PDF path (default: ``eom.pdf`` in current directory).
    open_pdf : bool
        If True, open the PDF after compilation (macOS ``open``).

    Returns the PDF path, or ``None`` when ``pdflatex`` fails; the failure (with the end
    of the ``pdflatex`` log) is a ``WARNING`` on the ``geomech`` logger, the success an
    ``INFO`` line.
    """
    tex_src = eom_to_latex(eqns)
    output = os.path.abspath(output)
    pdf_name = os.path.splitext(os.path.basename(output))[0]

    with tempfile.TemporaryDirectory() as tmpdir:
        tex_path = os.path.join(tmpdir, pdf_name + ".tex")
        with open(tex_path, "w") as f:
            f.write(tex_src)

        result = subprocess.run(
            ["pdflatex", "-interaction=nonstopmode", "-output-directory", tmpdir, tex_path],
            capture_output=True,
            text=True,
        )
        pdf_tmp = os.path.join(tmpdir, pdf_name + ".pdf")
        if result.returncode != 0 or not os.path.exists(pdf_tmp):
            logger.warning(
                "render_eom: pdflatex failed:\n%s",
                result.stdout[-2000:] if len(result.stdout) > 2000 else result.stdout,
            )
            return None

        import shutil

        shutil.copy2(pdf_tmp, output)

    logger.info("render_eom: PDF written to %s", output)
    if open_pdf:
        subprocess.run(["open", output])
    return output
