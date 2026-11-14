"""Structural invariants of expression trees: ``check_tree(expr)``.

The constructors validate most of these when a node is built; ``check_tree`` re-checks
a whole tree after the fact, because a tree can also be assembled by hand
(``Add([Add(a, b), c])`` does not flatten), mutated (``node.nodes[0] = ...``), or produced
by a rewrite that bypasses a constructor guard. Each violation raises ``InvariantError``
naming the invariant, the offending subtree and its path from the root (``root.R.0``:
the right child's first addend; ``L``/``R`` for binary nodes, the index for n-ary nodes,
``x`` for the operand of a unary node).

Checks, per node class:

operand types
    ``Add``/``VAdd``/``MAdd`` children are SCALAR/VECTOR/MATRIX; ``Mul`` (S, S);
    ``SVMul`` (V, S); ``SMMul`` (M, S); ``MVMul`` (M, V) or a row times a matrix
    (row V, M); ``MMMul`` (M, M); ``VVMul`` (V, V); ``Dot``/``Cross`` (V, V); ``Hat`` V;
    ``Vee`` M; every node has a type and the number of children its class takes.
nested sums
    no ``Add`` directly under an ``Add`` (same for ``VAdd``, ``MAdd``): the constructors
    flatten, so a nested sum means the tree was built around them.
row vectors
    ``VVMul`` has exactly one row operand (row × column is the inner product, column ×
    row the outer product); ``Dot``, ``Cross`` and ``Hat`` take columns only; a ``VAdd``
    does not mix rows and columns (zero addends are exempt).
calculus
    no ``TimeDerivative``/``Variation`` directly over an ``S2``/``SO3`` leaf (their
    kinematics are eager: ``q.t_diff()`` is ``ω × q``), no ``Variation`` of a
    ``Variation``, and no ``Variation`` of a variation vector (``ξ``, ``η``).
tangent vectors
    every ``TS2``/``TSO3`` leaf has its manifold point set (``.S2``/``.SO3``) and is
    named like that point's tangent or variation vector.
numeric scalars
    a ``Scalar`` has a ``value`` exactly when it is flagged numeric, the value is a
    real number, and a scalar printed like a number in parentheses (``(-1)``,
    ``(0.5)``, the form ``Number`` prints) carries that value: a named
    ``Scalar("(-1)")`` without one would print like the number and be treated as a
    symbol.
"""

from __future__ import annotations

import math
import numbers

from geomech.utils.errors import InvariantError

__all__ = ["check_tree", "InvariantError"]


def _short(expr, limit=120):
    s = str(expr)
    return s if len(s) <= limit else s[: limit - 3] + "..."


def _fail(path, expr, what):
    raise InvariantError(f"check_tree: {what} at {path}: {type(expr).__name__} '{_short(expr)}'")


def _type_name(t):
    return getattr(t, "name", str(t))


def _numeric_name_value(name):
    """The number a parenthesised scalar name like ``(-1)`` or ``(0.5)`` stands for, else None.

    ``Number(-1)`` prints ``(-1)``; a bare ``2`` (``Number("2")``) is a named constant by
    design and is not checked.
    """
    if not (isinstance(name, str) and name.startswith("(") and name.endswith(")")):
        return None
    try:
        return float(name[1:-1])
    except ValueError:
        return None


def check_tree(expr):
    """Raise ``InvariantError`` if *expr* violates a structural invariant; else return None.

    See the module docstring for the list of checks. The walk is iterative, so deep
    trees do not hit the recursion limit; shared subtrees are checked once.
    """
    from geomech.core.base.expressions import S2, SO3, TS2, TSO3, Expr, Scalar
    from geomech.core.base.types import ExprType
    from geomech.core.operations.addition import Add, MAdd, VAdd
    from geomech.core.operations.calculus import TimeDerivative, TimeIntegral, Variation
    from geomech.core.operations.geometry import Cross, Dot, Hat, Transpose, Vee, is_row_vector
    from geomech.core.operations.multiplication import MMMul, Mul, MVMul, SMMul, SVMul, VVMul

    S, V, M = ExprType.SCALAR, ExprType.VECTOR, ExprType.MATRIX
    nary = {Add: S, VAdd: V, MAdd: M}
    binary = {
        Mul: ((S, S),),
        SVMul: ((V, S),),
        SMMul: ((M, S),),
        MVMul: ((M, V), (V, M)),
        MMMul: ((M, M),),
        VVMul: ((V, V),),
        Dot: ((V, V),),
        Cross: ((V, V),),
    }
    unary = {Hat: V, Vee: M, Transpose: None, TimeDerivative: None, Variation: None}
    unary[TimeIntegral] = None

    seen = set()
    stack = [(expr, "root")]
    while stack:
        node, path = stack.pop()
        if id(node) in seen:
            continue
        seen.add(id(node))
        if not isinstance(node, Expr):
            raise InvariantError(
                f"check_tree: {path} is a {type(node).__name__} ({node!r}), not an expression"
            )
        cls = type(node)
        children = getattr(node, "nodes", None)
        for i, child in enumerate(children or ()):
            if not isinstance(child, Expr):
                raise InvariantError(
                    f"check_tree: operand {i} of {cls.__name__} at {path} is a"
                    f" {type(child).__name__} ({child!r}), not an expression"
                )
        if isinstance(node, tuple(unary)) and (not children or len(children) != 1):
            _fail(path, node, f"{cls.__name__} needs exactly one operand")
        node_type = getattr(node, "type", None)
        if node_type not in (S, V, M):
            _fail(path, node, f"node has no expression type (type is {node_type!r})")

        if cls in nary or isinstance(node, (Add, VAdd, MAdd)):
            want = nary[next(c for c in nary if isinstance(node, c))]
            if not children:
                _fail(path, node, "empty sum")
            for i, child in enumerate(children):
                if getattr(child, "type", None) != want:
                    _fail(
                        f"{path}.{i}",
                        child,
                        f"{cls.__name__} addend has type"
                        f" {_type_name(getattr(child, 'type', None))}, expected {want.name}",
                    )
                if type(child) is cls:
                    _fail(
                        f"{path}.{i}", child, f"{cls.__name__} nested directly in a {cls.__name__}"
                    )
            if cls is VAdd:
                orient = {i: is_row_vector(c) for i, c in enumerate(children) if not c.is_zero}
                if len(set(orient.values())) > 1:
                    _fail(path, node, "VAdd mixes row and column vectors")
            for i, child in enumerate(children):
                stack.append((child, f"{path}.{i}"))
            continue

        if isinstance(node, tuple(binary)):
            kinds = next(v for c, v in binary.items() if isinstance(node, c))
            if not children or len(children) != 2:
                _fail(path, node, f"{cls.__name__} needs 2 operands, has {len(children or ())}")
            l, r = children
            pair = (getattr(l, "type", None), getattr(r, "type", None))
            if pair not in kinds:
                _fail(
                    path,
                    node,
                    f"{cls.__name__} operand types {_type_name(pair[0])} x {_type_name(pair[1])}"
                    f" (expected {' or '.join(a.name + ' x ' + b.name for a, b in kinds)})",
                )
            if isinstance(node, VVMul) and is_row_vector(l) == is_row_vector(r):
                both = "rows" if is_row_vector(l) else "columns"
                _fail(path, node, f"VVMul needs exactly one row operand; both operands are {both}")
            if isinstance(node, MVMul) and pair == (V, M) and not is_row_vector(l):
                _fail(f"{path}.L", l, "a vector times a matrix needs a row vector on the left")
            if isinstance(node, MVMul) and pair == (M, V) and is_row_vector(r):
                _fail(f"{path}.R", r, "a matrix times a row vector is not defined")
            if isinstance(node, (Dot, Cross)):
                for side, child in (("L", l), ("R", r)):
                    if is_row_vector(child):
                        _fail(
                            f"{path}.{side}",
                            child,
                            f"{cls.__name__} operand is a row vector; {cls.__name__} takes"
                            " columns (use the column, or VVMul for rowᵀ·column)",
                        )
            stack.append((r, f"{path}.R"))
            stack.append((l, f"{path}.L"))
            continue

        if isinstance(node, tuple(unary)):
            if not children or len(children) != 1 or children[0] is None:
                _fail(path, node, f"{cls.__name__} needs exactly one operand")
            inner = children[0]
            want = unary[next(c for c in unary if isinstance(node, c))]
            if want is not None and getattr(inner, "type", None) != want:
                _fail(
                    f"{path}.x",
                    inner,
                    f"{cls.__name__} operand has type"
                    f" {_type_name(getattr(inner, 'type', None))}, expected {want.name}",
                )
            if isinstance(node, Hat) and is_row_vector(inner):
                _fail(f"{path}.x", inner, "Hat operand is a row vector; Hat takes a column")
            if isinstance(node, (TimeDerivative, Variation)) and isinstance(inner, (S2, SO3)):
                _fail(
                    path,
                    node,
                    f"{cls.__name__} directly over the manifold point '{inner}' (its kinematics"
                    " are eager: use .t_diff()/.delta())",
                )
            if isinstance(node, Variation):
                if isinstance(inner, Variation):
                    _fail(path, node, "variation of a variation")
                if isinstance(inner, (TS2, TSO3)):
                    parent = inner.S2 if isinstance(inner, TS2) else inner.SO3
                    if (
                        parent is not None
                        and inner.name == parent.manifold_info.variation_vector_name
                    ):
                        _fail(path, node, "variation of a variation vector")
            stack.append((inner, f"{path}.x"))
            continue

        # ---- leaves ----
        if children:
            _fail(path, node, f"unknown node class {cls.__name__} with children")
        if isinstance(node, (TS2, TSO3)):
            parent = node.S2 if isinstance(node, TS2) else node.SO3
            kind = "S2" if isinstance(node, TS2) else "SO3"
            if parent is None:
                _fail(path, node, f"{cls.__name__} leaf has no {kind} point (.{kind} is None)")
            info = parent.manifold_info
            if node.name not in (info.tangent_vector_name, info.variation_vector_name):
                _fail(
                    path,
                    node,
                    f"{cls.__name__} leaf is named neither like the tangent vector"
                    f" '{info.tangent_vector_name}' nor like the variation vector"
                    f" '{info.variation_vector_name}' of its point '{parent}'",
                )
        if isinstance(node, Scalar):
            value = node.value
            if (value is not None) != bool(node.is_numeric):
                _fail(
                    path,
                    node,
                    f"Scalar value {value!r} inconsistent with is_numeric={node.is_numeric}",
                )
            if value is not None:
                if isinstance(value, bool) or not isinstance(value, numbers.Real):
                    _fail(path, node, f"numeric Scalar has a non-real value {value!r}")
            printed = _numeric_name_value(node.name)
            if printed is not None:
                if value is None:
                    _fail(
                        path,
                        node,
                        "Scalar is printed like a number but has no value; it would be"
                        " treated as a symbol (use Number(...) for a numeric constant)",
                    )
                elif not math.isclose(float(value), printed, rel_tol=1e-9, abs_tol=0.0):
                    _fail(path, node, f"Scalar printed '{node.name}' has value {value!r}")
    return None
