"""Indexed families and finite sums with a symbolic number of terms.

An indexed sum ``Σ_{j=1}^{i} l_j q_j`` is kept *outside* the expression tree, as a list
of terms. Each term is a summand (an ordinary ``geomech`` expression whose leaves carry
the index in their names, ``q_{j}``) together with its *bound* indices and the
*order constraints* between them (``1 ≤ j ≤ i ≤ n``). The engine never sees a sum
node: the derivation runs on each summand with the existing rewrites, and only the
bookkeeping of indices happens here.

Index sets
----------
Every bound index runs over ``1, …, n``; ``n`` is the symbolic number of members of
every family. A term's constraints are pairs ``a ≤ b`` or ``a < b`` between *atoms*:
index names, the size ``n`` and integers. ``1 ≤ j`` and ``j ≤ n`` hold for every index
and are not stored. Two facts are used and nothing else:

- substitution: renaming an index in the summand and in the constraints together
  leaves the sum unchanged (the leaves are renamed through their names, ``q_{j}`` →
  ``q_{k}``);
- the Kronecker rule: the coefficient of ``ξ_p`` in ``Σ_{j ∈ S} c_j · ξ_j`` is
  ``c_p`` when ``p ∈ S`` and zero otherwise, so the bound index ``j`` is replaced by
  the free index ``p`` and the condition ``p ∈ S`` becomes constraints on the remaining
  indices (``Σ_i [p ≤ i] f = Σ_{i ≥ p} f``).

A sum can be expanded for a given ``n`` (``instantiate``), which gives the ordinary
expression the same model has with ``n`` members; the tests check the derived
equations through it.
"""

from __future__ import annotations

import itertools
import re
from dataclasses import dataclass, field

from geomech.core.base import (
    S1,
    S2,
    SO3,
    TS1,
    TS2,
    TSO3,
    ExprType,
    Matrix,
    Scalar,
    Vector,
)
from geomech.core.base.constants import zero_of
from geomech.core.base.expressions import DeferredOperand, Expr, Number, manifold_point_of
from geomech.core.transformations.utils import addends, map_children, sum_of

__all__ = [
    "Index",
    "Size",
    "n",
    "Family",
    "Sum",
    "IndexedExpr",
    "Term",
    "indices",
    "rename_indices",
    "as_indexed",
]

#: the symbolic number of members of every family
SIZE = "n"

_FRESH = itertools.count(1)


def _fresh() -> str:
    """A bound-index name no user index uses (``ι1``, ``ι2``, …)."""
    return f"ι{next(_FRESH)}"


# ---------------------------------------------------------------------------
# Indices and bounds
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Index:
    """A summation or family index, ``i = Index("i")``.

    ``i - 1`` and ``i + 1`` are strict bounds: ``Sum(j, 1, i - 1, …)`` runs over
    ``j < i`` and ``Sum(j, i + 1, n, …)`` over ``j > i``.
    """

    name: str

    def __post_init__(self):
        if not isinstance(self.name, str) or not re.fullmatch(r"[A-Za-z][A-Za-z0-9]*", self.name):
            raise ValueError(f"Index: name must be letters and digits, got {self.name!r}")
        if self.name == SIZE:
            raise ValueError(f"Index: '{SIZE}' is reserved for the number of members")

    def __sub__(self, k):
        if k != 1:
            raise ValueError(f"Index: only {self.name} - 1 is supported, got {self.name} - {k}")
        return _Strict(self.name, below=True)

    def __add__(self, k):
        if k != 1:
            raise ValueError(f"Index: only {self.name} + 1 is supported, got {self.name} + {k}")
        return _Strict(self.name, below=False)

    def __str__(self):
        return self.name


@dataclass(frozen=True)
class _Strict:
    """``i - 1`` (``below=True``) or ``i + 1`` as a summation bound."""

    name: str
    below: bool


@dataclass(frozen=True)
class Size:
    """The symbolic number of members ``n`` (there is one, shared by every family)."""

    def __str__(self):
        return SIZE


n = Size()


def indices(names: str) -> list[Index]:
    """``i, j, k = indices("i j k")``."""
    return [Index(s) for s in names.split()]


def _atom(x):
    """A bound as a constraint atom: an index name, ``"n"`` or an int."""
    if isinstance(x, Index):
        return x.name
    if isinstance(x, Size):
        return SIZE
    if isinstance(x, bool) or not isinstance(x, int):
        raise TypeError(f"Sum: a bound must be an Index, n, i ± 1 or an int, got {x!r}")
    if x < 1:
        raise ValueError(f"Sum: an integer bound must be at least 1, got {x}")
    return x


# ---------------------------------------------------------------------------
# Constraints
# ---------------------------------------------------------------------------

#: one constraint: (a, b, strict) means a < b when strict, else a ≤ b
Constraint = tuple


def _sub_atom(a, mapping):
    return mapping.get(a, a) if isinstance(a, str) else a


def _value(a, env, size):
    if isinstance(a, int):
        return a
    if a == SIZE:
        return size
    return env[a]


def normalize(cons, bound=()) -> frozenset | None:
    """*cons* without the constraints every index satisfies; ``None`` if one is false.

    Every index lies in ``1, …, n``, so ``1 ≤ x``, ``x ≤ n`` and ``a ≤ a`` are dropped,
    and ``a < a``, ``x < 1`` and ``n < x`` make the set empty. Two integers are compared;
    an integer against ``n`` (``3 ≤ n``) is kept, as it holds only for some n.
    """
    out = set()
    for a, b, strict in cons:
        if isinstance(a, int) and isinstance(b, int):
            if (a < b) if strict else (a <= b):
                continue
            return None
        if a == b:
            if strict:
                return None
            continue
        # 1 ≤ x and x ≤ n hold for every index x (and 1 ≤ n); c ≤ n for an int c ≥ 2
        # depends on n and is kept
        if not strict and (a == 1 or (b == SIZE and isinstance(a, str))):
            continue
        if strict and (b == 1 or (a == SIZE and isinstance(b, str))):
            return None
        out.add((a, b, strict))
    return frozenset(out)


def _cyclic(cons) -> bool:
    """True when the constraints force ``a < a`` through a chain (an empty sum)."""
    graph: dict = {}
    for a, b, strict in cons:
        graph.setdefault(a, []).append((b, strict))
    for start in graph:
        stack = [(start, False)]
        seen = set()
        while stack:
            node, strict = stack.pop()
            for nxt, s in graph.get(node, ()):
                st = strict or s
                if nxt == start and st:
                    return True
                if (nxt, st) not in seen:
                    seen.add((nxt, st))
                    stack.append((nxt, st))
    return False


# ---------------------------------------------------------------------------
# Families and renaming
# ---------------------------------------------------------------------------

_KINDS = (Scalar, Vector, Matrix, S1, S2, SO3)
_AUTO_ATTR = {"Manifold", "UnitNorm", "OrthogonalMatrix"}
_INDEX_IN_NAME = re.compile(r"_\{([^{}]+)\}")


@dataclass(frozen=True)
class Family:
    """Indexed symbols ``q_1, …, q_n`` of one kind: ``q = Family(S2, "q")``, ``q[i]``.

    *kind* is ``Scalar``, ``Vector``, ``Matrix``, ``S1``, ``S2`` or ``SO3``; *attr* is
    passed to every member (``attr=["Constant"]`` for masses and lengths), *axis* to the
    members of an ``S1`` family. ``q[i]`` is the leaf ``q_{i}``, ``q[3]`` the leaf
    ``q_{3}``.
    """

    kind: type
    name: str
    attr: tuple = ()
    axis: Expr | None = field(default=None, compare=False)

    def __init__(self, kind, name, attr=None, axis=None):
        if kind not in _KINDS:
            raise TypeError(
                f"Family: kind must be one of {', '.join(k.__name__ for k in _KINDS)}, got"
                f" {kind!r}"
            )
        if not isinstance(name, str) or "{" in name or "}" in name or "*" in name:
            raise ValueError(f"Family: invalid name {name!r}")
        object.__setattr__(self, "kind", kind)
        object.__setattr__(self, "name", name)
        object.__setattr__(self, "attr", tuple(attr or ()))
        object.__setattr__(self, "axis", axis)
        _FAMILY_NAMES.add(name)

    def leaf_name(self, key) -> str:
        return f"{self.name}_{{{key}}}"

    def __getitem__(self, idx):
        if isinstance(idx, Index):
            key = idx.name
        elif isinstance(idx, int) and not isinstance(idx, bool) and idx >= 1:
            key = str(idx)
        elif isinstance(idx, str):  # internal: an index name
            key = idx
        else:
            raise TypeError(f"Family {self.name}: index must be an Index or an int ≥ 1")
        name = self.leaf_name(key)
        attr = list(self.attr)
        if self.kind is S1:
            return S1(name, axis=self.axis, attr=attr)
        return self.kind(name, attr=attr)

    def member_key(self, leaf) -> str | None:
        """The index of *leaf* in this family (``"j"``, ``"3"``), else ``None``."""
        name = getattr(leaf, "name", None)
        if not isinstance(name, str) or not isinstance(leaf, self.kind):
            return None
        if self.kind in (Scalar, Vector, Matrix) and isinstance(leaf, (S1, S2, SO3)):
            return None
        m = re.fullmatch(re.escape(self.name) + r"_\{([^{}]+)\}", name)
        return m.group(1) if m else None

    def __str__(self):
        return self.name


#: names of every Family created; ``<name>_{k}`` is reserved for its members
_FAMILY_NAMES: set = set()
_MEMBER = re.compile(r"(.+)_\{([^{}]+)\}")


def member_of(name):
    """``(family name, key)`` when *name* is ``<family>_{key}`` of a Family, else None."""
    if not isinstance(name, str):
        return None
    m = _MEMBER.fullmatch(name)
    if m and m.group(1) in _FAMILY_NAMES:
        return m.group(1), m.group(2)
    return None


def _sub_name(name: str, mapping) -> str:
    """*name* with its index renamed, for a family member ``q_{j}``; other names as is."""
    hit = member_of(name)
    if hit is None or hit[1] not in mapping:
        return name
    return f"{hit[0]}_{{{mapping[hit[1]]}}}"


def _rename_leaf(leaf, mapping, memo):
    key = id(leaf)
    if key in memo:
        return memo[key]
    out = _rename_leaf_uncached(leaf, mapping, memo)
    memo[key] = out
    return out


def _rename_leaf_uncached(leaf, mapping, memo):
    if isinstance(leaf, (TS1, TS2, TSO3)):
        point = manifold_point_of(leaf)
        if point is None:
            new = _sub_name(leaf.name, mapping)
            return leaf if new == leaf.name else type(leaf)(new)
        new_point = _rename_leaf(point, mapping, memo)
        if new_point is point:
            return leaf
        info = point.manifold_info
        if leaf.name == info.tangent_vector_name:
            return new_point.get_tangent_vector()
        if leaf.name == info.variation_vector_name:
            return new_point.get_variation_vector()
        raise ValueError(f"rename_indices: '{leaf}' is neither tangent nor variation of '{point}'")
    name = getattr(leaf, "name", None)
    if not isinstance(name, str) or getattr(leaf, "is_numeric", False):
        return leaf
    new = _sub_name(name, mapping)
    if new == name:
        return leaf
    if isinstance(leaf, (S1, S2, SO3)):
        attr = [a for a in (leaf.attr or []) if a not in _AUTO_ATTR]
        if isinstance(leaf, S1):
            return S1(new, axis=_rename_leaf(leaf.axis, mapping, memo), attr=attr)
        return type(leaf)(new, attr=attr)
    attr = list(leaf.attr or [])
    if type(leaf) is Scalar:
        return Scalar(new, attr=attr)
    if type(leaf) is Vector:
        return Vector(new, size=leaf.size, attr=attr)
    if type(leaf) is Matrix:
        return Matrix(new, size=leaf.size, attr=attr)
    raise TypeError(f"rename_indices: cannot rename the {type(leaf).__name__} leaf '{leaf}'")


def rename_indices(expr: Expr, mapping: dict) -> Expr:
    """*expr* with every index ``k`` of *mapping* replaced in the leaf names.

    ``q_{j}`` becomes ``q_{mapping[j]}``, and the tangent and variation leaves follow
    their point (``ω_{q_{j}}`` → ``ω_{q_{k}}``). Values may be names or ints.
    """
    mapping = {k: str(v) for k, v in mapping.items() if str(v) != k}
    if not mapping:
        return expr
    memo: dict = {}
    cache: dict = {}

    def walk(e):
        key = id(e)
        if key in cache:
            return cache[key]
        if getattr(e, "nodes", None):
            out = map_children(e, walk, "rename_indices")
        else:
            out = _rename_leaf(e, mapping, memo)
        cache[key] = out
        return out

    return walk(expr)


def leaf_index_names(expr: Expr) -> set[str]:
    """The index names that occur in the leaf names of *expr*."""
    out: set[str] = set()

    def walk(e):
        nodes = getattr(e, "nodes", None)
        if nodes:
            for c in nodes:
                walk(c)
        else:
            hit = member_of(getattr(e, "name", None))
            if hit is not None:
                out.add(hit[1])

    walk(expr)
    return out


# ---------------------------------------------------------------------------
# Terms and indexed expressions
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Term:
    """``Σ_{bound ∈ cons} body``: a summand with its bound indices and constraints."""

    bound: tuple
    cons: frozenset
    body: Expr

    def rename(self, mapping: dict) -> Term:
        """The same term with indices renamed (bound or free) in the body and constraints."""
        mapping = {k: v for k, v in mapping.items() if v != k}
        if not mapping:
            return self
        bound = tuple(mapping.get(b, b) for b in self.bound)
        cons = frozenset(
            (_sub_atom(a, mapping), _sub_atom(b, mapping), s) for a, b, s in self.cons
        )
        return Term(bound, cons, rename_indices(self.body, mapping))

    def freshen(self) -> Term:
        """The same term with new names for its bound indices (no capture when combined)."""
        return self.rename({b: _fresh() for b in self.bound})

    def substitute(self, index: str, value) -> Term | None:
        """``index := value`` (an index name or an int); the index leaves the bound set.

        Returns ``None`` when the constraints become false (an empty sum).
        """
        body = rename_indices(self.body, {index: value})
        cons = frozenset(
            (value if a == index else a, value if b == index else b, s) for a, b, s in self.cons
        )
        if isinstance(value, int) and index in self.bound:
            # the bound index ranged over 1..n: pinned to c, the term needs c ≤ n
            cons = cons | {(value, SIZE, False)}
        kept = normalize(cons)
        if kept is None:
            return None
        bound = tuple(b for b in self.bound if b != index)
        return Term(bound, kept, body)._eliminate_equalities()

    def _eliminate_equalities(self) -> Term | None:
        """Replace a bound index ``b`` with ``a ≤ b ≤ a`` by ``a``; ``None`` if empty."""
        term = self
        while True:
            if _cyclic(term.cons):
                return None
            hit = None
            for b in term.bound:
                for a, c, s in term.cons:
                    if s or b not in (a, c):
                        continue
                    other = c if a == b else a
                    if (b, other, False) in term.cons and (other, b, False) in term.cons:
                        hit = (b, other)
                        break
                if hit:
                    break
            if hit is None:
                return term
            b, other = hit
            nxt = term.substitute(b, other)
            if nxt is None:
                return None
            term = nxt

    def satisfied(self, env: dict, size: int) -> bool:
        for a, b, strict in self.cons:
            va, vb = _value(a, env, size), _value(b, env, size)
            if (va >= vb) if strict else (va > vb):
                return False
        return True


def _type_of(expr) -> ExprType:
    return expr.type


class IndexedExpr(DeferredOperand):
    """A sum of terms ``Σ_{bound ∈ cons} body`` with a common expression type.

    Built with ``Sum``; supports ``+``, ``-``, ``*`` (with expressions, numbers and other
    indexed expressions), ``dot``, ``cross``, ``t_diff`` and ``T``. A product of two sums
    is the sum over the pairs of their bound indices.
    """

    __slots__ = ("terms", "type")

    def __init__(self, terms, expr_type):
        self.terms = tuple(terms)
        self.type = expr_type

    # -- construction ------------------------------------------------------

    @staticmethod
    def of(expr) -> IndexedExpr:
        return as_indexed(expr)

    def _map(self, fn) -> IndexedExpr:
        terms = [Term(t.bound, t.cons, fn(t.body)) for t in self.terms]
        typ = _type_of(terms[0].body) if terms else self.type
        return IndexedExpr(terms, typ)

    def _pair(self, other, fn) -> IndexedExpr:
        other = as_indexed(other)
        terms = []
        for a in self.terms:
            for b in other.terms:
                b = b.freshen()
                cons = normalize(a.cons | b.cons)
                if cons is None:
                    continue
                terms.append(Term(a.bound + b.bound, cons, fn(a.body, b.body)))
        if not terms:
            probe = fn(_probe(self.type), _probe(other.type))
            return IndexedExpr((), probe.type)
        return IndexedExpr(terms, _type_of(terms[0].body))

    # -- arithmetic --------------------------------------------------------

    def __add__(self, other):
        if isinstance(other, (int, float)) and other == 0:
            return self
        other = as_indexed(other)
        if other.type != self.type:
            raise TypeError(f"Sum: cannot add a {other.type.name} to a {self.type.name} sum")
        return IndexedExpr(self.terms + other.terms, self.type)

    def __radd__(self, other):
        if isinstance(other, (int, float)) and other == 0:
            return self
        return as_indexed(other) + self

    def __neg__(self):
        return self._map(lambda b: b * -1)

    def __sub__(self, other):
        return self + (-as_indexed(other))

    def __rsub__(self, other):
        return as_indexed(other) + (-self)

    def __mul__(self, other):
        if isinstance(other, (int, float)) and not isinstance(other, bool):
            return self._map(lambda b: b * other)
        return self._pair(other, lambda a, b: a * b)

    def __rmul__(self, other):
        if isinstance(other, (int, float)) and not isinstance(other, bool):
            return self._map(lambda b: Number(other) * b)
        return as_indexed(other)._pair(self, lambda a, b: a * b)

    def dot(self, other):
        from geomech.core.operations import Dot

        return self._pair(other, Dot)

    def cross(self, other):
        from geomech.core.operations import Cross

        return self._pair(other, Cross)

    def t_diff(self):
        return self._map(lambda b: b.t_diff())

    def T(self):
        return self._map(lambda b: b.T())

    # -- evaluation --------------------------------------------------------

    def instantiate(self, size: int, env: dict | None = None) -> Expr:
        """The expression for ``n = size``, with the free indices of *env* set to ints."""
        return instantiate_terms(self.terms, self.type, size, env or {})

    def __repr__(self):
        return f"IndexedExpr({len(self.terms)} terms, {self.type.name})"

    def __str__(self):
        from geomech.indexed.printing import terms_str

        return terms_str(self.terms)


def _probe(expr_type):
    return {
        ExprType.SCALAR: Scalar("probe"),
        ExprType.VECTOR: Vector("probe"),
        ExprType.MATRIX: Matrix("probe"),
    }[expr_type]


def as_indexed(expr) -> IndexedExpr:
    """*expr* as an indexed expression (a plain expression is a term with no bound index)."""
    if isinstance(expr, IndexedExpr):
        return expr
    if isinstance(expr, (int, float)) and not isinstance(expr, bool):
        expr = Number(expr)
    if not isinstance(expr, Expr):
        raise TypeError(f"expected a geomech expression or a Sum, got {type(expr).__name__}")
    return IndexedExpr((Term((), frozenset(), expr),), expr.type)


def Sum(index: Index, lower, upper, body) -> IndexedExpr:
    """``Σ_{index = lower}^{upper} body``.

    *lower* and *upper* are ``1`` (or another int), ``n``, an index ``i``, or ``i - 1`` /
    ``i + 1`` for a strict bound. The body may itself contain sums, whose bounds may refer
    to *index*: ``Sum(i, 1, n, Sum(j, 1, i, l[j] * q[j]))``.
    """
    if not isinstance(index, Index):
        raise TypeError(f"Sum: the summation index must be an Index, got {index!r}")
    body = as_indexed(body)
    terms = []
    for t in body.terms:
        f = _fresh()
        t = t.rename({index.name: f})
        extra = []
        for bound, is_lower in ((lower, True), (upper, False)):
            if isinstance(bound, _Strict):
                atom = bound.name
                if atom == index.name:
                    raise ValueError(f"Sum: a bound of {index} refers to {index} itself")
                if is_lower and not bound.below:
                    extra.append((atom, f, True))
                elif not is_lower and bound.below:
                    extra.append((f, atom, True))
                else:
                    raise ValueError(
                        f"Sum: use {atom} + 1 as a lower bound and {atom} - 1 as an upper bound"
                    )
            else:
                atom = _atom(bound)
                if atom == index.name:
                    raise ValueError(f"Sum: a bound of {index} refers to {index} itself")
                extra.append((atom, f, False) if is_lower else (f, atom, False))
        cons = normalize(t.cons | frozenset(extra))
        if cons is None or _cyclic(cons):
            continue
        terms.append(Term((f,) + t.bound, cons, t.body))
    return IndexedExpr(terms, body.type)


def instantiate_terms(terms, expr_type, size: int, env: dict) -> Expr:
    """The terms expanded for ``n = size``: an ordinary expression."""
    if not isinstance(size, int) or size < 1:
        raise ValueError(f"instantiate: the size must be an int ≥ 1, got {size!r}")
    parts = []
    for t in terms:
        for values in itertools.product(range(1, size + 1), repeat=len(t.bound)):
            local = dict(env)
            local.update(zip(t.bound, values))
            try:
                ok = t.satisfied(local, size)
            except KeyError as err:
                raise ValueError(
                    f"instantiate: the free index {err.args[0]} has no value; pass it in env"
                ) from None
            if not ok:
                continue
            mapping = {k: v for k, v in local.items() if isinstance(k, str)}
            body = rename_indices(t.body, mapping)
            parts.extend(addends(body))
    if not parts:
        return zero_of(expr_type)
    return sum_of(parts, expr_type)


def _free_names(term: Term):
    names = set()
    for a, b, _ in term.cons:
        for x in (a, b):
            if isinstance(x, str) and x != SIZE and x not in term.bound:
                names.add(x)
    return names
