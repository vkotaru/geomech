"""Phase 5b: constants module, like-term keys, literal spelling, frozen nodes, D2 identity.

Fix plan items 5.5-5.7 and decisions D2 (a) and D3 (a). Every equation-level change is
checked numerically (``geomech.numeric.evaluate``), not only by its printed form.
"""

from __future__ import annotations

import copy
import dataclasses
import pickle

import numpy as np
import pytest

import geomech
from geomech import (
    S2,
    SO3,
    Add,
    Cross,
    Dot,
    Matrix,
    MMMul,
    Mul,
    MVMul,
    Number,
    One,
    Scalar,
    SMMul,
    SVMul,
    SystemVariables,
    TimeDerivative,
    Transpose,
    Variation,
    Vector,
    VVMul,
    Zero,
    ZeroMatrix,
    ZeroVector,
    compute_eom,
    evaluate,
    full_simplify,
)
from geomech.core.base.constants import HALF, NEG_ONE, NEG_TWO, negate, zero_like, zero_of
from geomech.core.base.types import ExprFlags, ExprType
from geomech.core.math.extract import extract_linear_coeff, has_algebraic, same_node
from geomech.core.transformations.simplification import term_key

# ---------------------------------------------------------------------------
# 5.5 constants
# ---------------------------------------------------------------------------


def test_constant_literals_print_like_number():
    assert str(NEG_ONE) == str(Number(-1)) == "(-1)"
    assert str(NEG_TWO) == str(Number(-2)) == "(-2)"
    assert str(HALF) == str(Number(0.5)) == "(0.5)"
    assert NEG_ONE.value == -1 and HALF.value == 0.5


@pytest.mark.parametrize(
    "expr",
    [Scalar("m"), Vector("v"), Vector("v").T(), Matrix("A")],
    ids=["scalar", "vector", "row", "matrix"],
)
def test_negate_is_the_same_tree_as_unary_minus(expr):
    assert str(negate(expr)) == str(-expr)
    assert type(negate(expr)) is type(-expr)


def test_negate_rejects_a_non_expression():
    with pytest.raises(TypeError, match="negate"):
        negate(3)


def test_zero_like():
    assert zero_like(Scalar("m")) is Zero
    assert zero_like(Vector("v")) is ZeroVector
    assert zero_like(Matrix("A")) is ZeroMatrix
    assert zero_of(ExprType.VECTOR) is ZeroVector
    assert zero_of(None) is Zero


def test_constant_symbol_variations_are_the_zero_leaves():
    """Phase 4 follow-up 5: δ and d/dt of a constant vector/matrix are ZeroVector/ZeroMatrix
    (printed 0v / 0, as conventions.md lists the zeros); a non-default size is kept."""
    e3, J = Vector("e3", attr=["Constant"]), Matrix("J", attr=["Constant"])
    assert e3.delta() is ZeroVector and e3.t_diff() is ZeroVector
    assert J.delta() is ZeroMatrix and J.t_diff() is ZeroMatrix
    v6 = Vector("c", size=(6,), attr=["Constant"])
    assert v6.delta().size == (6,) and v6.delta().is_zero


def test_public_helper_names():
    """has_algebraic, same_node, check_variables, wrap_numeric are public now."""
    from geomech.core.operations.multiplication import _wrap_numeric, wrap_numeric
    from geomech.dynamics.eom import check_variables

    x = Vector("x")
    assert has_algebraic(Dot(x, Vector("f")), x)
    assert same_node(x, Vector("x")) and not same_node(x, Scalar("x"))
    assert str(wrap_numeric(2)) == "(2)" and _wrap_numeric is wrap_numeric
    check_variables(SystemVariables(vectors=[x]))


# ---------------------------------------------------------------------------
# like terms: linear grouping by term_key, inner products commute
# ---------------------------------------------------------------------------


def test_inner_products_merge_either_way_round():
    """Phase 5a follow-up 1: aᵀb s − bᵀa s is zero, aᵀb s + bᵀa s is 2 aᵀb s."""
    a, b, s = Vector("a"), Vector("b"), Scalar("s")
    ab, ba = a.T() * b, b.T() * a
    assert full_simplify(ab * s - ba * s) is Zero
    total = ab * s + ba * s
    merged = full_simplify(total)
    assert not isinstance(merged, Add)
    env = {}
    assert np.isclose(evaluate(merged, env), evaluate(total, env))


def test_term_key_symmetries():
    a, b, c = Vector("a"), Vector("b"), Vector("c")
    m, k = Scalar("m"), Scalar("k")
    assert term_key(Dot(a, b)) == term_key(Dot(b, a))
    assert term_key(Mul(Mul(m, k), m)) == term_key(Mul(m, Mul(m, k)))
    assert term_key(VVMul(a.T(), b)) == term_key(VVMul(b.T(), a))
    # order still matters for every other product
    assert term_key(Cross(a, b)) != term_key(Cross(b, a))
    assert term_key(MMMul(Matrix("A"), Matrix("B"))) != term_key(MMMul(Matrix("B"), Matrix("A")))
    assert term_key(VVMul(a, b.T())) != term_key(VVMul(b, a.T()))  # outer products
    # class and name both count; the key is unambiguous for names that concatenate
    assert term_key(Scalar("x")) != term_key(Vector("x"))
    assert term_key(Mul(Scalar("ab"), c.T() * c)) != term_key(Mul(Scalar("a"), c.T() * c))


def test_like_terms_group_in_one_pass_and_keep_order():
    """A long sum of repeated terms is merged by key; groups keep first-term order and the
    value is unchanged."""
    syms = [Scalar(f"s{i}") for i in range(40)]
    terms = [syms[i % 7] * syms[(i * 3) % 11] for i in range(200)]
    total = Add(*terms)
    out = full_simplify(total)
    env = {}
    assert np.isclose(evaluate(out, env), evaluate(total, env))
    assert len(out.nodes) <= 77


# ---------------------------------------------------------------------------
# phase 5a advisory: literal spelling
# ---------------------------------------------------------------------------


def test_literal_spellings_are_unified():
    x, v = Scalar("x"), Vector("v")
    assert str(full_simplify(Number(2.0))) == "(2)"
    assert str(full_simplify(x * Scalar("0.5", value=0.5))) == str(full_simplify(x * 0.5))
    assert str(full_simplify(v * Number(2.0))) == str(full_simplify(v * 2))
    # Zero and One keep their names; a named scalar with a value is a symbol (D1)
    assert "m" in str(full_simplify(Scalar("m", value=2) * x))
    assert str(full_simplify(Scalar("2", value=2) * x)) == str(full_simplify(2 * x))


# ---------------------------------------------------------------------------
# D3: frozen nodes
# ---------------------------------------------------------------------------


def _sample_tree():
    q, R, x = S2("q"), SO3("R"), Vector("x")
    m = Scalar("m", attr=["Constant"])
    return Dot(q.t_diff(), MVMul(R, x)) * m + Number(2)


@pytest.mark.parametrize(
    "mutate",
    [
        lambda e, x, w: setattr(e, "nodes", ()),
        lambda e, x, w: setattr(x, "name", "y"),
        lambda e, x, w: setattr(x.flags, "is_constant", True),
        lambda e, x, w: setattr(w, "S2", None),
        lambda e, x, w: setattr(x, "anything", 1),
        lambda e, x, w: delattr(x, "name"),
    ],
    ids=["nodes", "name", "flags", "tangent-point", "new-attribute", "delete"],
)
def test_nodes_are_frozen(mutate):
    e, x, w = _sample_tree(), Vector("x"), S2("q").get_tangent_vector()
    before = (str(e), hash(e), str(x), x.flags)
    with pytest.raises(dataclasses.FrozenInstanceError):
        mutate(e, x, w)
    assert (str(e), hash(e), str(x), x.flags) == before


def test_children_are_tuples_and_flags_frozen():
    e = _sample_tree()
    stack = [e]
    while stack:
        node = stack.pop()
        nodes = getattr(node, "nodes", None)
        if nodes is not None and not isinstance(node, (Scalar, Vector, Matrix)):
            assert isinstance(nodes, tuple), type(node)
            stack.extend(nodes)
    flags = ExprFlags.from_attr_list(["Constant"])
    assert dataclasses.replace(flags, is_zero=True).is_zero and not flags.is_zero


@pytest.mark.parametrize(
    "roundtrip",
    [lambda o: pickle.loads(pickle.dumps(o)), copy.copy, copy.deepcopy],
    ids=["pickle", "copy", "deepcopy"],
)
def test_pickle_and_copy_of_a_frozen_tree(roundtrip):
    e = _sample_tree()
    w = S2("q").get_tangent_vector()
    e2, w2 = roundtrip(e), roundtrip(w)
    assert str(e2) == str(e) and e2 == e and hash(e2) == hash(e)
    assert isinstance(e2.nodes, tuple)
    assert str(w2.S2) == "q" and w2.flags == w.flags
    with pytest.raises(dataclasses.FrozenInstanceError):
        e2.nodes = ()


def test_unpickling_list_operands_gives_tuples():
    """A tree pickled before the freeze held its operands in a list."""
    e = Add(Scalar("a"), Scalar("b"))
    state = dict(e.__dict__, nodes=list(e.nodes))
    clone = Add.__new__(Add)
    clone.__setstate__(state)
    assert isinstance(clone.nodes, tuple) and str(clone) == str(e)


def test_shared_constants_are_safe_to_share():
    """NEG_ONE sits in many trees; nothing can change it under them."""
    t1, t2 = negate(Scalar("a")), negate(Vector("b"))
    assert t1.right is NEG_ONE and t2.right is NEG_ONE
    with pytest.raises(dataclasses.FrozenInstanceError):
        NEG_ONE.value = 1


def test_compute_eom_unchanged_by_the_freeze():
    """A pinned system still derives (the string is checked in test_manifold_pipelines)."""
    q = S2("q")
    m, g, l = (Scalar(n, attr=["Constant"]) for n in "mgl")
    e3 = Vector("e3", attr=["Constant"])
    omega = q.get_tangent_vector()
    L = HALF * m * l * l * Dot(omega, omega) - m * g * l * Dot(q, e3)
    eom = compute_eom(L, Zero, SystemVariables(vectors=[q]))
    assert geomech.check_eom(L, Zero, SystemVariables(vectors=[q]), eom) < 1e-6


# ---------------------------------------------------------------------------
# small structure checks
# ---------------------------------------------------------------------------


def test_one_is_the_removed_factor_marker():
    """``_remove_factor`` tests ``is One`` instead of ``str(x) == "1"``."""
    from geomech.core.math.ibp import _remove_factor

    a, x = Scalar("a"), Scalar("x")
    assert _remove_factor(x, x) is One
    assert str(_remove_factor(Mul(a, x), x)) == "a"


def test_extract_through_a_shared_literal():
    x, k = Vector("x"), Scalar("k")
    expr = Dot(x.delta(), SVMul(Vector("f"), Mul(k, NEG_ONE)))
    coeff = extract_linear_coeff(expr, x.delta())
    assert np.allclose(evaluate(coeff, {}), evaluate(SVMul(Vector("f"), Mul(k, NEG_ONE)), {}))


def test_transpose_and_scaled_matrix_types():
    A = Matrix("A")
    assert Transpose(A).type is ExprType.MATRIX
    assert SMMul(A, 2).type is ExprType.MATRIX
    assert TimeDerivative(Vector("x")).type is ExprType.VECTOR
    assert Variation(Scalar("x")).type is ExprType.SCALAR


# ---------------------------------------------------------------------------
# D2 (a): separator in the product family, type guard in has(), calculus identity
# ---------------------------------------------------------------------------


def test_products_print_a_separator_and_no_longer_collide():
    """Engine audit F13: Mul(ab, c) == Mul(a, bc) was True (both printed 'abc')."""
    ab, c, a, bc = Scalar("ab"), Scalar("c"), Scalar("a"), Scalar("bc")
    assert str(Mul(ab, c)) == "ab*c" and str(Mul(a, bc)) == "a*bc"
    assert Mul(ab, c) != Mul(a, bc)
    AB, C, A, BC = Matrix("AB"), Matrix("C"), Matrix("A"), Matrix("BC")
    assert MMMul(AB, C) != MMMul(A, BC)
    v, w = Vector("v"), Vector("w")
    assert str(SVMul(v, Scalar("m"))) == "v*m"
    assert str(MVMul(Matrix("J"), v)) == "J*v"
    assert str(SMMul(Matrix("J"), 2)) == "J*(2)"
    assert str(VVMul(v.T(), w)) == "(v)'*w" and str(VVMul(v, w.T())) == "v*(w)'"


def test_repr_keeps_its_own_product_dot():
    """``str`` uses '*', ``repr`` keeps '·' (and '*' as the constant marker)."""
    m, x = Scalar("m", attr=["Constant"]), Vector("x")
    assert str(x * m) == "x*m"
    assert repr(x * m) == "v:x·m*"


def test_calculus_nodes_compare_the_operand_class():
    x, X = Scalar("x"), Vector("x")
    assert Variation(x) != Variation(X)
    assert TimeDerivative(Scalar("m")) != TimeDerivative(Vector("m"))
    assert TimeDerivative(Variation(x)) != TimeDerivative(Variation(X))
    assert Variation(x) == Variation(Scalar("x"))
    assert hash(Variation(x)) == hash(Variation(Scalar("x")))


def test_has_is_type_guarded():
    m_s, m_v, w = Scalar("m"), Vector("m"), Vector("w")
    assert not m_s.has(m_v)
    assert not Dot(m_v, w).has(m_s)
    assert Dot(m_v, w).has(Vector("m"))
    assert Mul(m_s, Scalar("k")).has(Scalar("m"))


def test_scalar_and_vector_variables_with_one_name_extract_apart():
    """Engine audit F13 reproducer: a scalar x next to a vector x made extraction fail."""
    x, X, k, F = Scalar("x"), Vector("x"), Scalar("k"), Vector("F")
    expr = Add(Mul(Variation(x), k), Dot(Variation(X), F))
    assert str(full_simplify(extract_linear_coeff(expr, Variation(x)))) == "k"
    assert str(full_simplify(extract_linear_coeff(expr, Variation(X)))) == "F"


# ---------------------------------------------------------------------------
# The separator is reserved (review F2), Transpose compares its operand class (F3),
# dispatch tables match the exact class (advisory)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "make",
    [
        lambda n: Scalar(n),
        lambda n: Scalar(name=n),
        lambda n: Vector(n),
        lambda n: Matrix(n),
        lambda n: S2(n),
        lambda n: SO3(n),
        lambda n: geomech.SkewSymmMatrix(n),
        lambda n: geomech.TS2(n),
        lambda n: geomech.TSO3(n),
        lambda n: Number(n),
    ],
)
def test_a_leaf_name_may_not_contain_the_product_separator(make):
    with pytest.raises(ValueError, match=r"'m\*l'.*reserved as the product separator"):
        make("m*l")


def test_a_lumped_symbol_without_the_separator_stays_apart():
    m, l, g = Scalar("m"), Scalar("l"), Scalar("g")
    assert Mul(Scalar("ml"), g) != Mul(m, Mul(l, g))
    assert str(Vector("x^\\ast")) == "x^\\ast"


def test_transpose_compares_the_operand_class():
    xv, xm = Vector("x"), Matrix("x")
    assert Transpose(xv) != Transpose(xm)
    assert Variation(Transpose(xv)) != Variation(Transpose(xm))
    assert TimeDerivative(Transpose(xm)) != TimeDerivative(Transpose(xv))
    assert Transpose(xm) == Transpose(Matrix("x"))
    assert hash(Transpose(xm)) == hash(Transpose(Matrix("x")))
    assert Variation(Transpose(xv)) == Variation(Transpose(Vector("x")))


def test_dispatch_tables_match_the_exact_class_only():
    """A subclass of a node class is not dispatched: simplify, pull and to_latex treat
    it like a leaf. No such subclass exists; a new one needs its own table entries."""
    from geomech import pull, simplify, to_latex

    class MyDot(Dot):
        pass

    v, w = Vector("v"), Vector("w")
    scaled = MyDot(SVMul(v, Number(2)), w)
    assert str(simplify(Dot(SVMul(v, Number(2)), w))) == "Dot(v,w)*(2)"
    assert simplify(scaled) is scaled
    assert pull(scaled) is scaled
    assert to_latex(scaled) == str(scaled)
