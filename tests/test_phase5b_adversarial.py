"""Adversarial review of phase 5b (structure, frozen nodes, typing, product separator).

Strict xfail tests are confirmed findings (``reason="F<n>: ..."``); they pass once the
finding is fixed. The plain tests pin what held: the separator equalities, the class-aware
calculus equality, the inner-product merge (numerically, against outer products and
non-symmetric bilinear forms), frozen nodes, pickling, and the corrected forms of two
negative assertions that the separator made vacuous (F1).
"""

from __future__ import annotations

import copy
import dataclasses
import os
import pickle
import sys

import numpy as np
import pytest

from geomech import (
    SO3,
    Add,
    Cross,
    Dot,
    Hat,
    MAdd,
    Matrix,
    MMMul,
    Mul,
    MVMul,
    Number,
    Scalar,
    SMMul,
    SVMul,
    TimeDerivative,
    Transpose,
    VAdd,
    Variation,
    Vector,
    VVMul,
    full_simplify,
    simplify,
)
from geomech.core.base.types import ExprType
from geomech.core.math.ibp import _t_diff_product
from geomech.core.transformations.simplification import _collect_like_terms, term_key
from geomech.utils.printing import _latex_str_key

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))


def _source(rel):
    with open(os.path.join(ROOT, rel), encoding="utf-8") as fh:
        return fh.read()


# ---------------------------------------------------------------------------
# A small numeric evaluator (scalars, 3-vectors, 3x3 matrices; no calculus)
# ---------------------------------------------------------------------------


def _value(expr, env, rng):
    nodes = getattr(expr, "nodes", None)
    if not nodes:
        key = (type(expr).__name__, str(expr))
        if expr.is_zero:
            return {ExprType.SCALAR: 0.0, ExprType.VECTOR: np.zeros(3)}.get(
                expr.type, np.zeros((3, 3))
            )
        if getattr(expr, "value", None) is not None and expr.type == ExprType.SCALAR:
            return float(expr.value)
        if key not in env:
            env[key] = {
                ExprType.SCALAR: lambda: float(rng.normal()),
                ExprType.VECTOR: lambda: rng.normal(size=3),
            }.get(expr.type, lambda: rng.normal(size=(3, 3)))()
        return env[key]
    v = [_value(n, env, rng) for n in nodes]
    match expr:
        case Add() | VAdd() | MAdd():
            return sum(v[1:], v[0])
        case Mul() | SVMul() | SMMul():
            return v[0] * v[1]
        case Dot():
            return float(v[0] @ v[1])
        case Cross():
            return np.cross(v[0], v[1])
        case Transpose():
            return v[0].T if np.ndim(v[0]) == 2 else v[0]
        case VVMul():
            return float(v[0] @ v[1]) if expr.type == ExprType.SCALAR else np.outer(v[0], v[1])
        case MVMul() | MMMul():
            return v[0] @ v[1]
        case Hat():
            a = v[0]
            return np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    raise NotImplementedError(type(expr).__name__)


def _same_value(a, b, seeds=3):
    for s in range(seeds):
        env, rng = {}, np.random.default_rng(s)
        va, vb = _value(a, env, rng), _value(b, env, rng)
        if not np.allclose(va, vb, rtol=1e-10, atol=1e-10):
            return False
    return True


# ---------------------------------------------------------------------------
# F1: negative assertions the separator made vacuous (lost coverage)
# ---------------------------------------------------------------------------


def test_f1_manifold_pipelines_negative_check_uses_the_separator():
    src = _source("tests/test_manifold_pipelines.py")
    assert '"(R)\'R" not in' not in src


def test_f1_vbl_negative_check_uses_the_separator():
    src = _source("tests/test_vbl.py")
    assert '"(-1)(-1)" not in' not in src
    assert '"(-1.0)(-1" not in' not in src
    assert '"(-1)(-1.0" not in' not in src


def test_f1_the_vacuous_patterns_cannot_match_any_more():
    """The old patterns never occur in a product's printed form now; the new ones do."""
    R = SO3("R")
    assert "(R)'R" not in str(MMMul(Transpose(R), R))
    assert "(R)'*R" in str(MMMul(Transpose(R), R))
    x = Scalar("x")
    chained = Mul(Mul(x, Number(-1)), Number(-1))
    assert "(-1)(-1)" not in str(chained)
    assert "(-1)*(-1)" in str(chained)


def test_f1_rigid_pendulum_f_has_no_RtR_with_the_separator():
    """What test_f_no_R_transpose_R meant, written with the separator: it holds."""
    from geomech import to_standard_form
    from tests.test_manifold_pipelines import _rigid_pendulum_com

    eom, variables, inputs, _ = _rigid_pendulum_com()
    key = list(eom.keys())[0]
    f = str(to_standard_form(eom, variables, inputs)[key].f)
    assert "(R)'*R" not in f


@pytest.mark.parametrize(
    "name",
    [
        "point_mass_r3",
        "spring_mass_r",
        "free_rigid_body",
        "rigid_pendulum_pivot",
        "spherical_pendulum_s2",
        "quadrotor",
        "double_rigid_pendulum",
        "elastic_spherical_pendulum",
    ],
)
def test_f1_vbl_blocks_have_no_unfolded_numbers_with_the_separator(name):
    """What test_no_unfolded_numeric_factors meant, written with the separator: it holds."""
    from tests.test_vbl import _all_block_strings, _derive

    _, _, lin = _derive(name)
    for text in _all_block_strings(lin):
        assert "(-1)*(-1)" not in text, text
        assert "(-1.0)*(-1" not in text, text
        assert "(-1)*(-1.0" not in text, text


# ---------------------------------------------------------------------------
# F2: a leaf name containing the separator reopens the F13 collision
# ---------------------------------------------------------------------------


def test_f2_separator_in_a_name_does_not_collide_with_a_product():
    m, l, g = Scalar("m"), Scalar("l"), Scalar("g")
    try:
        lumped = Scalar("m*l")
    except ValueError:
        return  # fixed by rejecting the separator in a symbol name
    assert Mul(lumped, g) != Mul(m, Mul(l, g))


def test_f2_t_diff_product_keeps_a_lumped_symbol_apart():
    x, y = Vector("x"), Vector("y")
    m, l = Scalar("m", attr=["Constant"]), Scalar("l", attr=["Constant"])
    try:
        lumped = Scalar("m*l", attr=["Constant"])
    except ValueError:
        return
    f = Dot(SVMul(x, Mul(m, l)), y)
    g = Dot(SVMul(x, lumped), y)
    out = _t_diff_product(Mul(f, g))
    # the product rule of two different factors has two terms
    assert isinstance(out, Add) and len(out.nodes) == 2, str(out)


def test_f2_like_terms_keep_a_lumped_symbol_apart():
    """term_key writes a leaf as Class#len#name, so like-term merging is not fooled."""
    m, l, g = Scalar("m"), Scalar("l"), Scalar("g")
    # the constructor now rejects '*' in a name (F2); build the leaf around it to keep
    # checking that term_key does not rely on that
    lumped = object.__new__(Scalar)
    lumped.__setstate__(vars(Scalar("ml")) | {"name": "m*l"})
    assert term_key(Mul(lumped, g)) != term_key(Mul(m, Mul(l, g)))
    out = full_simplify(Add(Mul(lumped, g), Mul(Mul(m, Mul(l, g)), Number(-1))))
    assert isinstance(out, Add)  # not cancelled to 0


# ---------------------------------------------------------------------------
# F3: Transpose is type-preserving but compares no operand class
# ---------------------------------------------------------------------------


def test_f3_transpose_compares_the_operand_class():
    xv, xm = Vector("x"), Matrix("x")
    assert Transpose(xv).type != Transpose(xm).type
    assert Transpose(xv) != Transpose(xm)
    assert Variation(Transpose(xv)) != Variation(Transpose(xm))


# ---------------------------------------------------------------------------
# F4: the engine fuzzer's shrinker lost its fresh-leaf candidates (tuple operands)
# ---------------------------------------------------------------------------


def test_f4_engine_fuzz_shrinker_offers_fresh_leaf_candidates():
    sys.path.insert(0, os.path.join(ROOT, "specs", "reviews", "scripts"))
    try:
        import engine_fuzz
    finally:
        sys.path.pop(0)
    a, b, v = Scalar("a"), Scalar("b"), Vector("v")
    cands = {str(c) for c in engine_fuzz.shrink_candidates(Dot(SVMul(v, Mul(a, b)), v))}
    assert "Dot(z_v,v)" in cands


# ---------------------------------------------------------------------------
# F5: display_standard_form turns every '*' of a key into a space
# ---------------------------------------------------------------------------


def test_f5_latex_key_keeps_a_star_in_a_name():
    assert _latex_str_key("u^*") == "u^*"


def test_f5_latex_key_product_and_derivatives_still_work():
    assert _latex_str_key("R*e3*f") == "R e3 f"
    assert _latex_str_key("\\frac{d}{dt}(\\omega_{q})") == "\\dot{\\omega_{q}}"
    assert _latex_str_key("\\frac{d}{dt}(\\frac{d}{dt}(x))") == "\\ddot{x}"


# ---------------------------------------------------------------------------
# What held: equality semantics
# ---------------------------------------------------------------------------


def test_trees_built_two_ways_are_still_equal():
    a, b = Scalar("a"), Scalar("b")
    v, w, M = Vector("v"), Vector("w"), Matrix("M")
    assert a * b == Mul(a, b) and hash(a * b) == hash(Mul(a, b))
    assert v * a == SVMul(v, a) == a * v
    assert M * v == MVMul(M, v)
    assert -a == a * -1 and -v == v * -1 and -M == M * -1
    assert v.T() * w == VVMul(Transpose(v), w)
    assert Variation(a) == Variation(Scalar("a"))
    assert TimeDerivative(Variation(v)) == TimeDerivative(Variation(Vector("v")))


def test_f13_pairs_are_unequal_and_calculus_compares_operand_class():
    ab, c, a, bc = Scalar("ab"), Scalar("c"), Scalar("a"), Scalar("bc")
    assert Mul(ab, c) != Mul(a, bc)
    assert MMMul(Matrix("AB"), Matrix("C")) != MMMul(Matrix("A"), Matrix("BC"))
    x, X = Scalar("x"), Vector("x")
    assert Variation(x) != Variation(X)
    assert TimeDerivative(TimeDerivative(x)) != TimeDerivative(TimeDerivative(X))
    assert not Dot(X, Vector("w")).has(x)
    assert {Variation(x): 1}.get(Variation(X)) is None


def test_reassociated_products_print_alike_and_have_one_value():
    """The separator does not mark grouping: (ab)c and a(bc) are equal, and equal in value."""
    a, b, c = Scalar("a"), Scalar("b"), Scalar("c")
    A, B, v = Matrix("A"), Matrix("B"), Vector("v")
    assert Mul(Mul(a, b), c) == Mul(a, Mul(b, c))
    assert MVMul(MMMul(A, B), v) == MVMul(A, MVMul(B, v))
    assert _same_value(MVMul(MMMul(A, B), v), MVMul(A, MVMul(B, v)))


# ---------------------------------------------------------------------------
# What held: the inner-product merge
# ---------------------------------------------------------------------------


def test_inner_product_merge_is_value_preserving_and_restricted():
    a, b, s = Vector("a"), Vector("b"), Scalar("s")
    J = Matrix("J", attr=["SymmetricMatrix"])
    M = Matrix("M")
    # aᵀb s − bᵀa s → 0
    e = Add(Mul(VVMul(a.T(), b), s), Mul(Mul(VVMul(b.T(), a), s), Number(-1)))
    assert str(simplify(e)) == "0"
    # (Ma)ᵀb = bᵀ(Ma): a genuine inner product of two columns, merged
    e = Add(VVMul(Transpose(MVMul(M, a)), b), VVMul(Transpose(b), MVMul(M, a)))
    out = simplify(e)
    assert not isinstance(out, Add) and _same_value(out, e)
    # aᵀMb vs bᵀMa: not equal for a general M, not merged
    e = Add(VVMul(a.T(), MVMul(M, b)), Mul(VVMul(b.T(), MVMul(M, a)), Number(-1)))
    out = simplify(e)
    assert isinstance(out, Add) and _same_value(out, e)
    # outer products a bᵀ vs b aᵀ: not merged
    e = MAdd(VVMul(a, b.T()), SMMul(VVMul(b, a.T()), -1))
    out = simplify(e)
    assert isinstance(out, MAdd) and _same_value(out, e)
    # symmetric J: aᵀJb = bᵀJa holds but is not recognised (conservative, no value change)
    e = Add(VVMul(a.T(), MVMul(J, b)), Mul(VVMul(b.T(), MVMul(J, a)), Number(-1)))
    out = simplify(e)
    assert isinstance(out, Add) and _same_value(out, e)


def test_like_term_grouping_random_sums_keep_their_value():
    """Seeded sums of inner/outer products, Dot, Mul and Cross with numeric coefficients."""
    rng = np.random.default_rng(5)
    vecs = [Vector(n) for n in "abc"]
    mats = [Matrix("M"), Matrix("N")]
    s = Scalar("s")
    nums = [Number(2), Number(-1), Number(0.5), Number(-2)]
    for _ in range(300):
        kind = rng.integers(3)
        terms = []
        for _ in range(rng.integers(2, 6)):
            a, b = (vecs[i] for i in rng.integers(3, size=2))
            k = nums[rng.integers(len(nums))]
            if kind == 0:
                base = [
                    VVMul(a.T(), b),
                    VVMul(b.T(), a),
                    Dot(a, b),
                    Mul(s, Dot(b, a)),
                    VVMul(a.T(), MVMul(mats[0], b)),
                ][rng.integers(5)]
                terms.append(Mul(base, k))
            elif kind == 1:
                base = [Cross(a, b), Cross(b, a), MVMul(mats[rng.integers(2)], a)][rng.integers(3)]
                terms.append(SVMul(base, k))
            else:
                base = [VVMul(a, b.T()), VVMul(b, a.T()), MMMul(mats[0], mats[1])][rng.integers(3)]
                terms.append(SMMul(base, k))
        expr_type = terms[0].type
        out = _collect_like_terms(list(terms), expr_type)
        build = {ExprType.SCALAR: Add, ExprType.VECTOR: VAdd, ExprType.MATRIX: MAdd}[expr_type]
        before = build(*terms)
        if not out:
            for seed in range(2):
                env, r = {}, np.random.default_rng(seed)
                assert np.allclose(_value(before, env, r), 0, atol=1e-10)
            continue
        after = out[0] if len(out) == 1 else build(*out)
        assert _same_value(before, after), (str(before), str(after))


# ---------------------------------------------------------------------------
# What held: frozen nodes, pickling, shared constants
# ---------------------------------------------------------------------------


def test_no_mutation_path_on_nodes_and_flags():
    R = SO3("R")
    e = Dot(MVMul(R, Vector("v")), Cross(R.get_tangent_vector(), Vector("w")))
    for target, attr, val in [
        (e, "nodes", ()),
        (e.left, "nodes", ()),
        (R, "name", "S"),
        (R, "flags", None),
        (R, "manifold_info", None),
        (R.get_tangent_vector(), "SO3", None),
        (VVMul(Vector("a").T(), Vector("b")), "_result_type", ExprType.MATRIX),
    ]:
        with pytest.raises(dataclasses.FrozenInstanceError):
            setattr(target, attr, val)
    with pytest.raises(dataclasses.FrozenInstanceError):
        R.flags.is_orthogonal = False
    with pytest.raises(TypeError):
        e.nodes[0] = Vector("z")  # tuple


def test_pickle_and_deepcopy_keep_identity_hash_and_type():
    R = SO3("R")
    e = Add(
        Mul(VVMul(Vector("a").T(), Vector("b")), Scalar("s")),
        Dot(R.get_tangent_vector(), Vector("w")),
    )
    for clone in (pickle.loads(pickle.dumps(e)), copy.deepcopy(e)):
        assert clone == e and hash(clone) == hash(e) and str(clone) == str(e)
        assert isinstance(clone.nodes, tuple)
        assert clone.nodes[0].left.type == ExprType.SCALAR  # the VVMul result type survives
        with pytest.raises(dataclasses.FrozenInstanceError):
            clone.nodes = ()
        assert {e: 1}[clone] == 1


def test_simplify_does_not_alias_shared_constants_into_wrong_values():
    from geomech.core.base.constants import HALF, NEG_ONE, negate

    x = Scalar("x")
    a = negate(x)
    b = negate(Scalar("y"))
    assert a.right is NEG_ONE and b.right is NEG_ONE
    assert str(full_simplify(Add(a, x))) == "0"
    assert str(NEG_ONE) == "(-1)" and str(HALF) == "(0.5)"
