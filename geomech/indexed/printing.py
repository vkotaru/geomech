"""Printing of indexed terms: nested sums with their bounds, in text and LaTeX.

The bound indices of a term are printed outermost first, those that index a
non-constant leaf (``ω̇_k``) before those that index only constants (``m_j``), so that
the constants collect inside, as in ``Σ_k (Σ_{j ≥ max(i,k)} m_j) l_i l_k (…)``. A
constraint between two bound indices bounds the inner one; a constraint with no bound
index is an Iverson bracket ``[i ≤ 3]``.
"""

from __future__ import annotations

import re

from geomech.indexed.core import SIZE, Term, rename_indices

_LETTERS = ["j", "k", "r", "s", "t", "u", "v", "w"]
_INDEX_IN_NAME = re.compile(r"_\{([^{}]+)\}")


def _dynamic(term: Term, b: str) -> bool:
    """True when the bound index *b* indexes a non-constant leaf of the body."""
    found = False

    def walk(e):
        nonlocal found
        nodes = getattr(e, "nodes", None)
        if nodes:
            for c in nodes:
                walk(c)
            return
        name = getattr(e, "name", None)
        if isinstance(name, str) and b in _INDEX_IN_NAME.findall(name):
            if not getattr(e, "is_constant", False):
                found = True

    walk(term.body)
    return found


def _order(term: Term) -> list:
    return sorted(term.bound, key=lambda b: (not _dynamic(term, b), term.bound.index(b)))


def _layout(term: Term, names: dict):
    """``(order, bounds, guards)`` with every index renamed through *names*."""
    order = _order(term)
    names = dict(names)
    letters = [x for x in _LETTERS if x not in names.values()]
    for b, letter in zip(order, letters):
        names[b] = letter

    def atom(a):
        return str(a) if isinstance(a, int) else (SIZE if a == SIZE else names.get(a, a))

    bounds, guards = [], []
    placed: set = set()
    for b in order:
        lows, highs = [], []
        for a, c, strict in term.cons:
            if c == b and (a not in term.bound or a in placed):
                lows.append(atom(a) + ("+1" if strict else ""))
            if a == b and (c not in term.bound or c in placed):
                highs.append(atom(c) + ("-1" if strict else ""))
        bounds.append((names[b], sorted(lows) or ["1"], sorted(highs) or [SIZE]))
        placed.add(b)
    for a, c, strict in term.cons:
        if a not in term.bound and c not in term.bound:
            guards.append((atom(a), "<" if strict else "≤", atom(c)))
    return order, names, bounds, sorted(guards)


def term_str(term: Term, names: dict) -> str:
    order, names, bounds, guards = _layout(term, names)
    out = ""
    for g in guards:
        out += f"[{g[0]} {g[1]} {g[2]}] "
    for letter, lows, highs in bounds:
        lo = lows[0] if len(lows) == 1 else f"max({', '.join(lows)})"
        hi = highs[0] if len(highs) == 1 else f"min({', '.join(highs)})"
        out += f"Σ_{{{letter}={lo}}}^{{{hi}}} "
    body = rename_indices(term.body, names)
    return out + (f"({body})" if order or guards else str(body))


def terms_str(terms, names: dict | None = None) -> str:
    if not terms:
        return "0"
    return " + ".join(term_str(t, names or {}) for t in terms)


def term_latex(term: Term, names: dict) -> str:
    from geomech.utils.printing import to_latex

    order, names, bounds, guards = _layout(term, names)
    out = ""
    rel = {"≤": r"\le", "<": "<"}
    for g in guards:
        out += rf"[{g[0]} {rel[g[1]]} {g[2]}]\, "
    for letter, lows, highs in bounds:
        lo = lows[0] if len(lows) == 1 else rf"\max({', '.join(lows)})"
        hi = highs[0] if len(highs) == 1 else rf"\min({', '.join(highs)})"
        out += rf"\sum_{{{letter}={lo}}}^{{{hi}}} "
    body = to_latex(rename_indices(term.body, names))
    return out + (rf"\left({body}\right)" if order or guards else body)


def terms_latex(terms, names: dict | None = None) -> str:
    if not terms:
        return "0"
    return " + ".join(term_latex(t, names or {}) for t in terms)
