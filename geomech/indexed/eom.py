"""Equations of motion for a symbolic number of bodies.

``compute_indexed_eom(L, dW, variables)`` derives, for each family of configuration
variables (``q_1, …, q_n``), the equation of its *i-th* member, with ``i`` and ``n``
symbolic; and for each ordinary configuration variable its equation, with the sums over
the families left as sums.

Each term ``Σ_{bound ∈ cons} body`` of the Lagrangian and of the virtual work is taken
through the stages of ``compute_eom`` on its own: the variation of the body, the
manifold rules, ``full_simplify``, integration by parts against the variation vectors of
the members that occur in the body, ``expand`` and the extraction of each coefficient.
The variation and the time derivative commute with a sum of fixed length, and Hamilton's
principle holds for each member of a family, so the coefficient of ``ξ_i`` is the sum
of the coefficients of ``ξ_j`` over the terms with ``j := i`` (``Term.substitute``):
when ``j`` is bound by the term, the constraints on ``j`` become constraints on the
other indices with ``i`` in its place; when the member is ``q_3``, the coefficient
counts for ``i = 3`` only. The terms of one equation are then merged (bound indices
renamed to a canonical order, and each addend to its least printed form under the
renamings that leave the constraints unchanged) and simplified.

``IndexedEquation.instantiate(size, at)`` expands an equation for given ``n`` and ``i``;
``IndexedEOM.instantiate(size)`` gives the ``compute_eom``-shaped dictionary of the
model with ``n = size``, which ``check_eom`` verifies against the unrolled Lagrangian.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field

from geomech.core.base import S1, S2, SO3, TS1, TS2, TSO3, Matrix, Scalar, Vector
from geomech.core.base.constants import zero_of
from geomech.core.base.expressions import manifold_point_of
from geomech.core.math.extract import same_node
from geomech.core.math.ibp import _integrate_by_parts
from geomech.core.transformations.expansion import expand
from geomech.core.transformations.manifold_rules import apply_manifold_rules
from geomech.core.transformations.simplification import full_simplify
from geomech.core.transformations.utils import addends, sum_of
from geomech.dynamics.eom import _separate_variations
from geomech.dynamics.variables import SystemVariables
from geomech.indexed.core import (
    Family,
    IndexedExpr,
    Term,
    as_indexed,
    instantiate_terms,
    member_of,
    normalize,
    rename_indices,
)

__all__ = ["IndexedVariables", "IndexedEquation", "IndexedEOM", "compute_indexed_eom", "FREE"]

#: the name of the free index of a family equation inside its terms
FREE = "ι0"

_CANON = [f"ν{k}" for k in range(1, 10)]


@dataclass
class IndexedVariables:
    """Configuration variables: families (``q_1, …, q_n``) and ordinary variables."""

    families: list = field(default_factory=list)
    scalars: list = field(default_factory=list)
    vectors: list = field(default_factory=list)
    matrices: list = field(default_factory=list)

    def __post_init__(self):
        for fam in self.families:
            if not isinstance(fam, Family):
                raise TypeError(f"IndexedVariables: families entries must be Family, got {fam!r}")
            if fam.attr and "Constant" in fam.attr:
                raise ValueError(f"IndexedVariables: the family '{fam}' is constant")
        names = [f.name for f in self.families]
        if len(set(names)) != len(names):
            raise ValueError("IndexedVariables: two families have the same name")
        for var in self._globals():
            hit = member_of(getattr(var, "name", None))
            if hit is not None:
                raise ValueError(
                    f"IndexedVariables: '{var}' is named like member {hit[1]} of the family"
                    f" '{hit[0]}'; give the variable another name or index the family"
                )

    def members(self, size: int) -> SystemVariables:
        """The ``SystemVariables`` of the model with ``n = size``."""
        scalars, vectors, matrices = list(self.scalars), list(self.vectors), list(self.matrices)
        for fam in self.families:
            leaves = [fam[c] for c in range(1, size + 1)]
            if fam.kind is Scalar:
                scalars += leaves
            elif fam.kind in (Matrix, SO3):
                matrices += leaves
            else:
                vectors += leaves
        return SystemVariables(scalars=scalars, vectors=vectors, matrices=matrices)

    def _globals(self):
        return list(self.scalars) + list(self.vectors) + list(self.matrices)


def _variation(var):
    return (
        var.delta()
        if isinstance(var, Scalar) and not isinstance(var, TS1)
        else (var.get_variation_vector())
    )


@dataclass
class IndexedEquation:
    """The equation ``Σ terms = 0`` of one configuration variable or family.

    For a family, ``variation`` is the variation vector of its member ``FREE``
    (``ξ_{q_{ι0}}``) and the terms have ``FREE`` as their free index.
    """

    name: str
    family: Family | None
    variation: object
    terms: tuple
    type: object

    def instantiate(self, size: int, at: int | None = None):
        """The equation for ``n = size`` (and member ``at`` of a family)."""
        if self.family is not None:
            if at is None:
                raise ValueError(f"instantiate: give the member of '{self.name}' with at=")
            if not 1 <= at <= size:
                raise ValueError(f"instantiate: at must be in 1..{size}, got {at}")
            env = {FREE: at}
        else:
            env = {}
        return instantiate_terms(self.terms, self.type, size, env)

    def to_latex(self, free: str = "i") -> str:
        from geomech.indexed.printing import terms_latex

        return terms_latex(self.terms, {FREE: free}) + " = 0"

    def __str__(self):
        from geomech.indexed.printing import terms_str

        return terms_str(self.terms, {FREE: "i"}) + " = 0"


class IndexedEOM(dict):
    """``{name: IndexedEquation}``, families by family name, variables by variation name."""

    def __init__(self, equations, variables):
        super().__init__(equations)
        self.variables = variables

    def instantiate(self, size: int) -> dict:
        """The ``compute_eom`` dictionary of the model with ``n = size``."""
        out = {}
        for eq in self.values():
            if eq.family is None:
                out[str(eq.variation)] = (eq.variation, eq.instantiate(size))
            else:
                for c in range(1, size + 1):
                    vec = _variation(eq.family[c])
                    out[str(vec)] = (vec, eq.instantiate(size, c))
        return out


def _members_in(expr, variables: IndexedVariables):
    """``{(family, key)}`` of the family members whose leaves occur in *expr*."""
    found = set()

    def walk(e):
        nodes = getattr(e, "nodes", None)
        if nodes:
            for c in nodes:
                walk(c)
            return
        leaf = manifold_point_of(e) if isinstance(e, (TS1, TS2, TSO3)) else e
        if leaf is None:
            return
        for fam in variables.families:
            key = fam.member_key(leaf)
            if key is not None:
                found.add((fam, key))

    walk(expr)
    return found


def _check_bound(expr, term: Term) -> None:
    """Raise when a family member in *expr* has an index that no Sum binds."""

    def walk(e):
        nodes = getattr(e, "nodes", None)
        if nodes:
            for c in nodes:
                walk(c)
            return
        leaf = manifold_point_of(e) if isinstance(e, (TS1, TS2, TSO3)) else e
        hit = member_of(getattr(leaf, "name", None))
        if hit is None:
            return
        fam, key = hit
        if key not in term.bound and not key.isdigit():
            raise ValueError(
                f"compute_indexed_eom: the index '{key}' of {fam}_{{{key}}} is not bound by a"
                " Sum; write every family member inside the Sum over its index"
            )

    walk(expr)


def _is_zero(expr) -> bool:
    return bool(getattr(expr, "is_zero", False))


def _derive_term(term: Term, expr, variables: IndexedVariables, is_variation: bool):
    """The coefficients of the variation vectors in one term: ``[(owner, Term)]``."""
    _check_bound(term.body, term)
    members = sorted(_members_in(expr, variables), key=lambda fk: (fk[0].name, fk[1]))
    owners = []  # (variable leaf, family or None, key)
    for var in variables._globals():
        owners.append((var, None, None))
    for fam, key in members:
        if key not in term.bound and not key.isdigit():
            raise ValueError(
                f"compute_indexed_eom: the index '{key}' of {fam.leaf_name(key)} is not bound"
                " by a Sum; write every family member inside the Sum over its index"
            )
        owners.append((fam[key], fam, key))
    scalars = [v for v, _, _ in owners if isinstance(v, Scalar) and not isinstance(v, TS1)]
    matrices = [v for v, _, _ in owners if isinstance(v, (Matrix, SO3))]
    vectors = [v for v, _, _ in owners if v not in scalars and v not in matrices]
    X = SystemVariables(scalars=scalars, vectors=vectors, matrices=matrices)
    vecs = [_variation(v) for v, _, _ in owners]
    dots = [v.t_diff() for v in vecs]

    dS = apply_manifold_rules(expr, X)
    dS = full_simplify(dS)
    dS = _integrate_by_parts(dS, dots)
    dS = expand(dS)
    for part in addends(dS):
        if _is_zero(part):
            continue
        if not any(_contains(part, v) for v in vecs):
            raise ValueError(
                "compute_indexed_eom: a term of δS holds no variation of a configuration"
                f" variable: '{part}'; list the variable (or its family) in IndexedVariables"
                " or mark the symbol 'Constant'"
            )
    eom = _separate_variations(dS, vecs)
    out = []
    for (var, fam, key), vec in zip(owners, vecs):
        coeff = eom[str(vec)][1]
        if _is_zero(coeff):
            continue
        t = Term(term.bound, term.cons, coeff)
        if fam is None:
            out.append(((str(vec), None, vec), t))
            continue
        if key.isdigit():
            c = int(key)
            body = rename_indices(coeff, {})
            cons = normalize(term.cons | {(FREE, c, False), (c, FREE, False)})
            if cons is None:
                continue
            t = Term(term.bound, cons, body)
        else:
            t = t.substitute(key, FREE)
            if t is None:
                continue
        out.append(((fam.name, fam, _variation(fam[FREE])), t))
    return out


def _contains(expr, query) -> bool:
    if same_node(expr, query):
        return True
    return any(_contains(c, query) for c in getattr(expr, "nodes", ()) or ())


def combine_terms(terms, expr_type):
    """Merge terms equal up to renaming their bound indices, and simplify each summand."""
    groups: dict = {}
    for t in terms:
        k = len(t.bound)
        if k > len(_CANON):
            raise ValueError(f"combine_terms: more than {len(_CANON)} bound indices")
        best_key, auts = None, []
        for perm in itertools.permutations(t.bound):
            mp = dict(zip(perm, _CANON))
            cons = frozenset(
                (
                    mp.get(a, a) if isinstance(a, str) else a,
                    mp.get(b, b) if isinstance(b, str) else b,
                    s,
                )
                for a, b, s in t.cons
            )
            key = (k, tuple(sorted(repr(c) for c in cons)))
            if best_key is None or key < best_key:
                best_key, auts = key, [(mp, cons)]
            elif key == best_key:
                auts.append((mp, cons))
        cons = auts[0][1]
        group = groups.setdefault(best_key, (tuple(_CANON[:k]), cons, []))
        for part in addends(t.body):
            if _is_zero(part):
                continue
            renamed = [rename_indices(part, mp) for mp, _ in auts]
            group[2].append(min(renamed, key=str))
    out = []
    for bound, cons, parts in groups.values():
        if not parts:
            continue
        body = full_simplify(sum_of(parts, expr_type))
        if _is_zero(body):
            continue
        out.append(Term(bound, cons, body))
    out.sort(key=lambda t: (len(t.bound), str(t.body)))
    return tuple(out)


def compute_indexed_eom(lagrangian, inf_work, variables: IndexedVariables) -> IndexedEOM:
    """Equations of motion of a model with families of ``n`` members, ``n`` symbolic.

    *lagrangian* and *inf_work* are expressions or ``Sum``s (``inf_work`` may be ``None``
    or 0). Returns an ``IndexedEOM`` with one ``IndexedEquation`` per family (the equation
    of its i-th member) and per ordinary configuration variable.
    """
    if not isinstance(variables, IndexedVariables):
        raise TypeError("compute_indexed_eom: variables must be IndexedVariables")
    if not variables.families and not variables._globals():
        raise ValueError("compute_indexed_eom: no configuration variables")
    L = as_indexed(lagrangian)
    if L.type.name != "SCALAR":
        raise TypeError("compute_indexed_eom: the Lagrangian must be a scalar")
    jobs = [(t, t.body.delta(), True) for t in L.terms]
    if inf_work is not None and not (isinstance(inf_work, (int, float)) and inf_work == 0):
        W = as_indexed(inf_work)
        if W.type.name != "SCALAR":
            raise TypeError("compute_indexed_eom: the virtual work must be a scalar")
        jobs += [(t, t.body, False) for t in W.terms]

    collected: dict = {}
    for term, expr, is_var in jobs:
        for (name, fam, vec), t in _derive_term(term, expr, variables, is_var):
            entry = collected.setdefault(name, (fam, vec, []))
            entry[2].append(t)

    equations = {}
    for var in variables._globals():
        vec = _variation(var)
        fam, _, terms = collected.get(str(vec), (None, vec, []))
        equations[str(vec)] = IndexedEquation(
            str(vec), None, vec, combine_terms(terms, vec.type), vec.type
        )
    for fam in variables.families:
        vec = _variation(fam[FREE])
        _, _, terms = collected.get(fam.name, (fam, vec, []))
        equations[fam.name] = IndexedEquation(
            fam.name, fam, vec, combine_terms(terms, vec.type), vec.type
        )
    return IndexedEOM(equations, variables)


def unroll(expr, size: int):
    """An expression or ``Sum`` expanded for ``n = size``."""
    if isinstance(expr, IndexedExpr):
        return expr.instantiate(size)
    if expr is None:
        return None
    return expr


__all__ += ["unroll", "combine_terms"]
_ = (S1, S2, Vector, zero_of)
