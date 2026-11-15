"""Sums over a symbolic number of bodies (``geomech.indexed``).

The generic equation of the i-th body is checked three ways: against Lee et al.
eq. 5.17 for the n-link spherical chain (for every n tested, not only n = 2), with
``check_eom`` on the unrolled Lagrangian for n = 1, …, 6 for every model, and by
asserting the index structure of the derived sums (which terms carry which bounds).
"""

from __future__ import annotations

import functools

import numpy as np
import pytest

from geomech import S2, SO3, Dot, Scalar, Vector, check_eom
from geomech.indexed import (
    FREE,
    Family,
    Index,
    IndexedEquation,
    IndexedVariables,
    Sum,
    Term,
    compute_indexed_eom,
    indices,
    n,
    rename_indices,
)
from geomech.indexed.core import normalize
from geomech.indexed.systems import MODELS, spherical_chain
from geomech.numeric import Evaluator, hat, sample_state

SIZES = (1, 2, 3, 4, 5, 6)
TOL = 1e-5


@functools.cache
def derived(name):
    model = MODELS[name]()
    return model, model.eom()


# ---------------------------------------------------------------------------
# Sums and index sets
# ---------------------------------------------------------------------------


def _family(kind=Scalar, name="a", attr=("Constant",)):
    return Family(kind, name, attr=list(attr))


def test_sum_instantiates_to_members():
    (i,) = indices("i")
    a = _family()
    assert str(Sum(i, 1, n, a[i]).instantiate(3)) == "(a_{1}+a_{2}+a_{3})"


def test_nested_sum_with_dependent_bound():
    """``Σ_{i=1}^{n} Σ_{j=1}^{i} a_j b_i`` has n(n+1)/2 terms."""
    i, j = indices("i j")
    a, b = _family(name="a"), _family(name="b")
    s = Sum(i, 1, n, Sum(j, 1, i, a[j] * b[i]))
    for size in (1, 2, 3, 4):
        parts = str(s.instantiate(size)).strip("()").split("+")
        assert len(parts) == size * (size + 1) // 2
    assert "a_{3}*b_{2}" not in str(s.instantiate(3))
    assert "a_{2}*b_{3}" in str(s.instantiate(3))


def test_strict_bounds():
    """``Σ_{j=1}^{i-1}`` runs over j < i and ``Σ_{j=i+1}^{n}`` over j > i."""
    i, j = indices("i j")
    a, b = _family(name="a"), _family(name="b")
    below = str(Sum(i, 1, n, Sum(j, 1, i - 1, a[j] * b[i])).instantiate(3))
    assert below == "(a_{1}*b_{2}+a_{1}*b_{3}+a_{2}*b_{3})"
    above = str(Sum(i, 1, n, Sum(j, i + 1, n, a[j] * b[i])).instantiate(3))
    assert above == "(a_{2}*b_{1}+a_{3}*b_{1}+a_{3}*b_{2})"


def test_product_of_sums_uses_distinct_indices():
    """``(Σ_j a_j)(Σ_j a_j)`` is a double sum: its two factors get different indices."""
    (j,) = indices("j")
    a = _family()
    s = Sum(j, 1, n, a[j])
    sq = s * s
    (term,) = sq.terms
    assert len(term.bound) == 2 and len(set(term.bound)) == 2
    assert len(str(sq.instantiate(2)).split("+")) == 4


def test_expression_times_sum():
    """A plain expression on the left of ``*`` and ``+`` defers to the sum."""
    (j,) = indices("j")
    a = _family()
    g = Scalar("g", attr=["Constant"])
    assert str((g * Sum(j, 1, n, a[j])).instantiate(2)) == "(g*a_{1}+g*a_{2})"
    assert str((g + Sum(j, 1, n, a[j])).instantiate(1)) == "(g+a_{1})"


def test_normalize_drops_trivial_and_detects_empty():
    assert normalize({(1, "j", False), ("j", "n", False), ("j", "j", False)}) == frozenset()
    assert normalize({("j", "j", True)}) is None
    assert normalize({("j", 1, True)}) is None
    assert normalize({(2, 3, False)}) == frozenset()
    assert normalize({(3, 2, False)}) is None


def test_kronecker_substitution():
    """The coefficient of ξ_p in ``Σ_{j ≤ i} c_j ξ_j`` is ``Σ_{i ≥ p} c_p``."""
    a = _family()
    t = Term(("i", "j"), frozenset({("j", "i", False)}), a["j"])
    s = t.substitute("j", "p")
    assert s.bound == ("i",)
    assert s.cons == frozenset({("p", "i", False)})
    assert str(s.body) == "a_{p}"


def test_equality_constraint_eliminates_index():
    """``p ≤ i ≤ p`` fixes i = p, and ``p < p`` empties the sum."""
    a = _family()
    t = Term(("i",), frozenset({("p", "i", False), ("i", "p", False)}), a["i"])
    s = t._eliminate_equalities()
    assert s.bound == () and str(s.body) == "a_{p}"
    assert Term(("i",), frozenset({("i", "p", True)}), a["i"]).substitute("i", "p") is None


def test_rename_follows_tangent_leaves():
    q = Family(S2, "q")
    w = q["j"].get_tangent_vector()
    out = rename_indices(Dot(w, q["j"]) + Dot(q["j"].get_variation_vector(), q["k"]), {"j": "k"})
    assert str(out) == "(Dot(\\omega_{q_{k}},q_{k})+Dot(\\xi_{q_{k}},q_{k}))"
    R = Family(SO3, "R")
    assert str(rename_indices(R["j"].get_tangent_vector(), {"j": 3})) == "\\Omega_{R_{3}}"


def test_index_errors():
    with pytest.raises(ValueError, match="reserved"):
        Index("n")
    (i,) = indices("i")
    with pytest.raises(ValueError, match="itself"):
        Sum(i, 1, i, _family()[i])
    with pytest.raises(TypeError, match="kind"):
        Family(int, "x")


def test_unbound_member_raises():
    """A family member outside the Sum over its index has no meaning; it is reported."""
    (i,) = indices("i")
    q = Family(S2, "q")
    L = Dot(q[i].get_tangent_vector(), q[i].get_tangent_vector())
    with pytest.raises(ValueError, match="not bound by a Sum"):
        compute_indexed_eom(L, None, IndexedVariables(families=[q]))


def test_fixed_member_counts_only_for_its_index():
    """A force on q_1 alone appears in the equation of member 1 only."""
    (i,) = indices("i")
    q = Family(S2, "q")
    m = Family(Scalar, "m", attr=["Constant"])
    F = Vector("F")
    w = q[i].get_tangent_vector()
    L = Sum(i, 1, n, m[i] * Dot(w, w) * 0.5)
    dW = Dot(q[1].get_variation_vector(), F)
    eom = compute_indexed_eom(L, dW, IndexedVariables(families=[q]))
    eq = eom["q"]
    assert "F" in str(eq.instantiate(3, 1))
    assert "F" not in str(eq.instantiate(3, 2))


# ---------------------------------------------------------------------------
# Every model, every n: check_eom on the unrolled Lagrangian
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", sorted(MODELS))
@pytest.mark.parametrize("size", SIZES)
def test_generic_equations_check_eom(name, size):
    model, eom = derived(name)
    L, dW, X, inputs = model.unroll(size)
    err = check_eom(L, dW, X, eom.instantiate(size), inputs=inputs)
    assert err < TOL, f"{name}, n = {size}: check_eom error {err:.2e}"


@pytest.mark.parametrize("name", sorted(MODELS))
def test_dropping_any_term_is_detected(name):
    """``check_eom`` is sensitive to every derived term: removing one fails it."""
    model, eom = derived(name)
    size = 3
    L, dW, X, inputs = model.unroll(size)
    for key, eq in eom.items():
        for k in range(len(eq.terms)):
            cut = dict(eom)
            cut[key] = IndexedEquation(
                eq.name, eq.family, eq.variation, eq.terms[:k] + eq.terms[k + 1 :], eq.type
            )
            mutated = type(eom)(cut, eom.variables)
            err = check_eom(L, dW, X, mutated.instantiate(size), inputs=inputs)
            assert err > 1e-3, f"{name}: dropping term {k} of {key} went unnoticed ({err:.1e})"


# ---------------------------------------------------------------------------
# The n-link spherical chain against Lee et al. eq. 5.17, for every n
# ---------------------------------------------------------------------------


def _chain_env(size, rng):
    model, _ = derived("spherical_chain")
    _, _, X, inputs = model.unroll(size)
    env = sample_state(X, inputs, rng)
    env["g"] = 9.81
    env["e3"] = np.array([[0.0], [0.0], [1.0]])
    for k in range(1, size + 1):
        env[f"m_{{{k}}}"] = float(rng.uniform(0.5, 2.0))
        env[f"l_{{{k}}}"] = float(rng.uniform(0.5, 1.5))
    return env


def _eq_5_17(env, size):
    """Rows of Lee et al. eq. 5.17 with constant ``m_ij = M_ij I``, ``M_ij = l_i l_j
    Σ_{k ≥ max(i,j)} m_k``, and ``U = Σ_k m_k g e₃·x_k``:

    ``M_ii ω̇_i + Σ_{j≠i} M_ij S(q_i)ᵀ S(q_j) ω̇_j − Σ_{j≠i} M_ij |ω_j|² S(q_i) q_j
    + g l_i (Σ_{k≥i} m_k) S(q_i) e₃``.
    """

    def vec(key):
        return np.asarray(env[key], float).reshape(3)

    q = [vec(f"q_{{{k}}}") for k in range(1, size + 1)]
    w = [vec(f"\\omega_{{q_{{{k}}}}}") for k in range(1, size + 1)]
    wd = [vec(f"\\frac{{d}}{{dt}}(\\omega_{{q_{{{k}}}}})") for k in range(1, size + 1)]
    m = [env[f"m_{{{k}}}"] for k in range(1, size + 1)]
    l_ = [env[f"l_{{{k}}}"] for k in range(1, size + 1)]
    e3 = vec("e3")

    def M(a, b):
        return l_[a] * l_[b] * sum(m[max(a, b) :])

    rows = []
    for a in range(size):
        row = M(a, a) * wd[a] + env["g"] * l_[a] * sum(m[a:]) * hat(q[a]) @ e3
        for b in range(size):
            if b != a:
                row = row + M(a, b) * hat(q[a]).T @ hat(q[b]) @ wd[b]
                row = row - M(a, b) * (w[b] @ w[b]) * hat(q[a]) @ q[b]
        rows.append(row)
    return q, rows


@pytest.mark.parametrize("size", SIZES)
def test_spherical_chain_matches_lee_5_17(size):
    """The generic chain equation, at member k, is ``E_k = −Row_k + f_k`` on the tangent
    plane of q_k, where ``Row_k`` is Lee et al. eq. 5.17."""
    _, eom = derived("spherical_chain")
    eq = eom["q"]
    rng = np.random.default_rng(size)
    for _ in range(3):
        env = _chain_env(size, rng)
        ev = Evaluator(env=env, rng=rng)
        q, rows = _eq_5_17(env, size)
        for k in range(1, size + 1):
            got = np.asarray(ev(eq.instantiate(size, k)), float).reshape(3)
            ref = -rows[k - 1] + np.asarray(env[f"f_{{{k}}}"], float).reshape(3)
            P = np.eye(3) - np.outer(q[k - 1], q[k - 1])
            err = np.max(np.abs(P @ (got - ref))) / max(1.0, np.max(np.abs(P @ ref)))
            assert err < 1e-9, f"n = {size}, member {k}: {got} vs {ref}"


# ---------------------------------------------------------------------------
# Structure of the derived sums
# ---------------------------------------------------------------------------


def _shapes(eq):
    """``{(number of bound indices, constraints)}`` of an equation's terms, with the
    bound indices numbered in the canonical order and the free index written ``i``."""
    out = set()
    for t in eq.terms:
        names = {b: f"b{k}" for k, b in enumerate(t.bound, 1)}
        names[FREE] = "i"
        cons = frozenset(
            (
                names.get(a, a) if isinstance(a, str) else a,
                names.get(b, b) if isinstance(b, str) else b,
                s,
            )
            for a, b, s in t.cons
        )
        out.add((len(t.bound), cons))
    return out


def test_spherical_chain_structure():
    """Input, gravity ``Σ_{k ≥ i} m_k``, and inertia ``Σ_j Σ_{k ≥ max(i, j)} m_k``: the
    mass sums of M_ij."""
    _, eom = derived("spherical_chain")
    eq = eom["q"]
    assert _shapes(eq) == {
        (0, frozenset()),
        (1, frozenset({("i", "b1", False)})),
        (2, frozenset({("i", "b2", False), ("b1", "b2", False)})),
    } or _shapes(eq) == {
        (0, frozenset()),
        (1, frozenset({("i", "b1", False)})),
        (2, frozenset({("i", "b1", False), ("b2", "b1", False)})),
    }
    text = str(eq)
    assert "Σ_{k=max(i, j)}^{n}" in text and "Σ_{j=i}^{n}" in text
    assert text.startswith("f_{i} + ")


def test_rigid_chain_structure():
    """Joints are ``Σ_{j<i}``: the masses below link i enter through strict bounds."""
    _, eom = derived("rigid_chain")
    text = str(eom["R"])
    assert "Σ_{j=i+1}^{n}" in text and "Σ_{j=1}^{i-1}" in text
    assert "Σ_{k=max(i+1, j+1)}^{n}" in text


def test_team_attitude_equation_is_eulers():
    """Each quadrotor's attitude decouples: ``−J_i Ω̇_i − Ω_i × J_iΩ_i + M_i = 0``."""
    for name in ("multi_quadrotor_point_mass", "multi_quadrotor_rigid_load"):
        _, eom = derived(name)
        assert str(eom["R"]) == (
            "(J_{i}*\\frac{d}{dt}(\\Omega_{R_{i}})*(-1)"
            "+Cross(\\Omega_{R_{i}},J_{i}*\\Omega_{R_{i}})*(-1)+M_{i}) = 0"
        )


def test_flexible_cable_quadrotor_translation():
    """The quadrotor's translation sees the whole cable: ``−(m_Q + Σ_k m_k) ẍ_Q`` appears
    as ``−m_Q ẍ_Q`` plus a sum over every link mass."""
    _, eom = derived("quadrotor_flexible_cable")
    eq = eom["\\delta{x_Q}"]
    text = str(eq)
    assert "\\frac{d}{dt}(\\frac{d}{dt}(x_Q))*m_Q*(-1)" in text
    assert "Σ_{j=1}^{n} (" in text and "m_{j}" in text


def test_latex_has_nested_sums():
    _, eom = derived("spherical_chain")
    tex = eom["q"].to_latex()
    assert r"\sum_{k=\max(i, j)}^{n}" in tex and tex.endswith("= 0")


def test_spherical_chain_model_is_reusable():
    """A model can be derived twice; the equations print the same."""
    a = spherical_chain().eom()["q"]
    b = spherical_chain().eom()["q"]
    assert str(a) == str(b)


@pytest.mark.parametrize("name", sorted(MODELS))
def test_instantiated_equations_feed_the_downstream_tools(name):
    """``eom.instantiate(n)`` is a ``compute_eom`` dictionary: standard form, VBL and code
    generation take it unchanged."""
    from geomech import linearize, to_standard_form
    from geomech.codegen.layout import StateLayout, collect_params
    from geomech.codegen.module import generate_module

    model, eom = derived(name)
    _, _, X, inputs = model.unroll(3)
    eqs = eom.instantiate(3)
    sf = to_standard_form(eqs, X, inputs)
    assert set(sf) == set(eqs)
    lin = linearize(eqs, X, inputs)
    assert lin.A and lin.B
    src = generate_module(sf, StateLayout(X, inputs, collect_params(sf, X, inputs)))
    assert "def " in src


# ---------------------------------------------------------------------------
# Findings of the independent review (each reproduced a wrong result before the fix)
# ---------------------------------------------------------------------------


def test_index_pinned_to_an_integer_keeps_the_size_condition():
    """``Σ_{k=2}^{2}`` is empty for n = 1: the derived equation must vanish there."""
    j, k = indices("j k")
    x = Family(Scalar, "x")
    m = Family(Scalar, "m", attr=["Constant"])
    L = Sum(j, 1, n, Sum(k, 2, 2, 0.5 * m[k] * x[j].t_diff() * x[j].t_diff()))
    assert str(L.instantiate(1)) == "0"
    eq = compute_indexed_eom(L, None, IndexedVariables(families=[x]))["x"]
    assert getattr(eq.instantiate(1, 1), "is_zero", False)
    assert "m_{2}" in str(eq.instantiate(2, 1))


def test_normalize_keeps_integer_against_size():
    assert normalize({(3, "n", False)}) == frozenset({(3, "n", False)})
    assert normalize({("n", 3, True)}) == frozenset({("n", 3, True)})
    assert normalize({(1, "n", False)}) == frozenset()


def test_only_family_members_are_renamed():
    """A constant ``e_{z}`` that is not a family member keeps its name inside Sum(z, …)."""
    (z,) = indices("z")
    p = Family(Vector, "p")
    e = Vector("e_{z}", attr=["Constant"])
    assert str(Sum(z, 1, n, p[z].dot(e)).instantiate(2)) == "(Dot(p_{1},e_{z})+Dot(p_{2},e_{z}))"


def test_variable_named_like_a_member_is_rejected():
    x = Family(Scalar, "x")
    with pytest.raises(ValueError, match="named like member"):
        IndexedVariables(families=[x], scalars=[Scalar("x_{1}")])


def test_unbound_constant_member_is_rejected():
    (i,) = indices("i")
    k = Index("k")
    q = Family(S2, "q")
    m = Family(Scalar, "m", attr=["Constant"])
    w = q[i].get_tangent_vector()
    L = Sum(i, 1, n, m[k] * Dot(w, w))
    with pytest.raises(ValueError, match="'k' of m_\\{k\\} is not bound"):
        compute_indexed_eom(L, None, IndexedVariables(families=[q]))
