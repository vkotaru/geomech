"""Adversarial review of phase 5a (simplifier completeness, decision D1).

Strict xfails (``reason="F<n>: ..."``) are confirmed findings of the review;
every other test pins a property that held:

- every rule firing is sound along curves: the two sides of each ``rule()`` call made by
  ``expand``/``pull``/``vector_rules``/``simplify`` agree as full (t, ε) Taylor jets on the
  constraint manifold, including under ``TimeDerivative``/``Variation`` and nested below
  ``Dot``/``MVMul``/``Hat``/``Transpose``;
- like-term merging keeps the order of matrix products, cross products and row forms;
- the hypothesis of ``simplify.transpose_skew`` is checked (products of skew matrices,
  derivatives of skew matrices, symmetric + skew sums are left alone);
- decision D1: a named scalar with a value is a symbol in every pass, in LaTeX and in
  standard form; ``evaluate`` uses the value; equal values do not merge;
- reviewer Lagrangians (non-symmetric and skew matrices, row forms, named values, S²,
  SO(3), products, scalar configurations) agree with the numerical Euler-Lagrange check and their
  linearizations with finite differences.
"""

from __future__ import annotations

import collections
import importlib
import importlib.util
import os
import zlib

import numpy as np
import pytest

import geomech.utils.trace as _trace
from geomech import (
    S2,
    SO3,
    Add,
    Cross,
    Dot,
    Hat,
    IdentityMatrix,
    MAdd,
    Matrix,
    MMMul,
    Mul,
    MVMul,
    Number,
    Scalar,
    SMMul,
    SVMul,
    SystemVariables,
    TimeDerivative,
    Transpose,
    VAdd,
    Variation,
    Vector,
    Vee,
    VVMul,
    compute_eom,
    evaluate,
    expand,
    full_simplify,
    linearize,
    pull,
    simplify,
    to_latex,
    to_standard_form,
    vector_rules,
)
from geomech.numeric import Evaluator, check_eom, finite_difference_check
from tests.vbl_systems import HALF

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FUZZ_PATH = os.path.join(REPO, "scripts", "engine_fuzz.py")

C = ["Constant"]


def _sc(n, **kw):
    return Scalar(n, attr=C, **kw)


def _vc(n):
    return Vector(n, attr=C)


def _mc(n, *extra):
    return Matrix(n, attr=C + list(extra))


class _NamedEvaluator(Evaluator):
    """Leaves drawn from a generator seeded by (name, seed): independent of traversal order."""

    def __init__(self, seed):
        super().__init__(env={})
        self.seed = seed

    def _reseed(self, key):
        self.rng = np.random.default_rng([zlib.crc32(key.encode()), self.seed])

    def leaf(self, node):
        if str(node) not in self.env:
            self._reseed(str(node))
        return super().leaf(node)

    def _derived_leaf(self, node):
        if str(node) not in self.env:
            self._reseed("d:" + str(node))
        return super()._derived_leaf(node)


def _close(e1, e2, seeds=range(4)):
    for seed in seeds:
        ev = _NamedEvaluator(seed)
        x = np.asarray(ev(e1), float)
        y = np.asarray(ev(e2), float)
        if x.shape != y.shape:
            x, y = x.reshape(-1), y.reshape(-1)
        if not np.allclose(x, y, atol=1e-12, rtol=1e-12):
            return False
    return True


# ---------------------------------------------------------------------------
# 1. Like-term merging keeps non-commutative order
# ---------------------------------------------------------------------------

a, b, c = Vector("a"), Vector("b"), Vector("c")
x, y, m = Scalar("x"), Scalar("y"), Scalar("m")
A, B = Matrix("A"), Matrix("B")
K = Matrix("K", attr=["SkewSymmetricMatrix"])
Js = Matrix("J", attr=["SymmetricMatrix"])

NON_COMMUTING = {
    # AB s - BA s is not 0
    "matrix_product_order": MAdd(SMMul(MMMul(A, B), m), SMMul(SMMul(MMMul(B, A), m), Number(-1))),
    # (a×b) m + (b×a) m is 0, not 2 (a×b) m
    "cross_order": VAdd(SVMul(Cross(a, b), m), SVMul(Cross(b, a), m)),
    # (a×b)·c m + m (b×a)·c is 0
    "triple_product_order": Add(Mul(Dot(Cross(a, b), c), m), Mul(m, Dot(Cross(b, a), c))),
    # aᵀ A b m - m bᵀ A a is not 0 for a non-symmetric A
    "row_form_order": Add(
        Mul(VVMul(Transpose(a), MVMul(A, b)), m),
        Mul(Mul(m, VVMul(Transpose(b), MVMul(A, a))), Number(-1)),
    ),
    # A a m - Aᵀ a m is not 0
    "matrix_vs_transpose": VAdd(SVMul(MVMul(A, a), m), SVMul(MVMul(Transpose(A), a), Mul(m, -1))),
    # â b̂ m - b̂ â m is not 0
    "hat_product_order": MAdd(
        SMMul(MMMul(Hat(a), Hat(b)), m), SMMul(MMMul(Hat(b), Hat(a)), Mul(m, -1))
    ),
    # outer products: a bᵀ m - b aᵀ m is not 0
    "outer_order": MAdd(
        SMMul(VVMul(a, Transpose(b)), m), SMMul(VVMul(b, Transpose(a)), Mul(m, -1))
    ),
}


@pytest.mark.parametrize("name", list(NON_COMMUTING))
def test_like_terms_keep_noncommutative_order(name):
    e = NON_COMMUTING[name]
    out = full_simplify(e)
    assert _close(e, out), f"{name}: {e} -> {out}"


def test_commutative_merges_are_made():
    assert (
        str(full_simplify(Add(Mul(Mul(x, Number(2)), y), Mul(Mul(x, y), Number(3))))) == "x*y*(5)"
    )
    assert str(full_simplify(Add(Mul(Dot(a, b), m), Mul(m, Dot(b, a))))) == "Dot(a,b)*m*(2)"
    e = VAdd(SVMul(Cross(a, b), Mul(m, x)), SVMul(SVMul(Cross(a, b), x), Mul(m, Number(-1))))
    assert str(full_simplify(e)) == "0v"


# ---------------------------------------------------------------------------
# 2. Hypotheses of the new simplify rules
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "inner",
    [
        MMMul(K, K),
        MMMul(Hat(a), Hat(b)),
        MAdd(Js, K),
        MMMul(Js, K),
        TimeDerivative(K),
    ],
    ids=["KK", "hat_hat", "J_plus_K", "JK", "dK"],
)
def test_transpose_skew_needs_a_skew_operand(inner):
    out = simplify(Transpose(inner))
    # no rule may turn these into -inner: none of them is skew-symmetric
    assert str(out) != str(SMMul(inner, Number(-1)))
    if not isinstance(inner, TimeDerivative):
        assert _close(Transpose(inner), out)


@pytest.mark.parametrize("inner", [K, MAdd(K, Hat(a)), MAdd(K, K), SMMul(K, m)])
def test_transpose_skew_preserves_value(inner):
    out = full_simplify(Transpose(inner))
    assert _close(Transpose(inner), out)
    assert _close(out, SMMul(inner, Number(-1)))


def test_mmmul_reassociate_and_combine_preserve_value():
    R = SO3("R")
    for e in [
        MMMul(MMMul(A, B), K),
        MMMul(Transpose(R), MMMul(R, A)),
        MMMul(MMMul(A, Transpose(R)), MMMul(R, B)),
        MMMul(Transpose(IdentityMatrix), MMMul(IdentityMatrix, A)),
        MMMul(Transpose(A), MMMul(A, B)),  # generic A: AᵀA is not I
    ]:
        out = full_simplify(e)
        assert _close(e, out), (e, out)
    assert "A" in str(full_simplify(MMMul(Transpose(A), MMMul(A, B))))


def test_vee_over_madd_is_linear_for_any_matrix():
    for e in [Vee(MAdd(A, B)), Vee(MAdd(A, Transpose(A), Hat(a))), Vee(MAdd(SMMul(A, m), K))]:
        out = full_simplify(expand(e))
        assert _close(e, out), (e, out)


@pytest.mark.parametrize(
    "e",
    [
        MVMul(Transpose(SVMul(a, m)), SMMul(A, x)),
        VVMul(Transpose(a), SVMul(b, m)),
        VVMul(SVMul(a, m), Transpose(SVMul(b, x))),
        Transpose(SVMul(a, m)),
        Transpose(SMMul(A, m)),
        Hat(SVMul(a, m)),
        Vee(SMMul(A, m)),
        Dot(Hat(SVMul(a, m)) * b, MVMul(Transpose(SMMul(A, x)), c)),
    ],
)
def test_new_pull_rules_preserve_value(e):
    assert _close(e, pull(e))
    assert _close(e, full_simplify(e))


# ---------------------------------------------------------------------------
# 3. Every rule firing is sound along curves (jets), including inside calculus nodes
# ---------------------------------------------------------------------------


def _load_fuzz():
    spec = importlib.util.spec_from_file_location("_engine_fuzz_for_review", FUZZ_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_RULE_MODULES = [
    "geomech.core.transformations.vector_identities",
    "geomech.core.transformations.pulling",
    "geomech.core.transformations.simplification",
    "geomech.core.transformations.expansion",
]


@pytest.fixture
def captured_rules(monkeypatch):
    fired = []

    def rec(name, before, after):
        fired.append((name, before, after))
        return after

    for name in _RULE_MODULES:
        monkeypatch.setattr(importlib.import_module(name), "rule", rec)
    monkeypatch.setattr(_trace, "_TRACE_RULES", True)  # the *_fold rules log only when on
    return fired


def _calculus_cases(L):
    q, w, xi = L.q, L.w, L.xi
    return [
        TimeDerivative(Dot(L.a, MVMul(L.M, Cross(q, Cross(w, q))))),
        TimeDerivative(TimeDerivative(Dot(L.c, Cross(q, Cross(w, q))))),
        TimeDerivative(MVMul(Hat(Cross(q, Cross(xi, q))), L.a)),
        Variation(Dot(L.a, Cross(q, Cross(w, q)))),
        TimeDerivative(Variation(Dot(Cross(w, q), L.b))),
        TimeDerivative(Dot(Cross(q, Cross(TimeDerivative(w), q)), L.a)),
        TimeDerivative(Dot(q, q) * L.s1),
        Variation(Dot(L.a, MVMul(MMMul(Transpose(L.R), L.R), L.b))),
        TimeDerivative(Vee(MAdd(MMMul(Transpose(L.R), MMMul(L.R, Hat(L.Om))), L.M))),
        TimeDerivative(Transpose(MAdd(L.K, Hat(Cross(q, Cross(w, q)))))),
        TimeDerivative(Dot(MVMul(Transpose(SVMul(L.a, L.s1)), SMMul(L.M, L.s2)), L.b)),
        Variation(Dot(Cross(q, Cross(w, q)), Cross(q, Cross(w, q)))),
        TimeDerivative(Dot(MVMul(L.R, w), MVMul(L.R, Cross(q, Cross(w, q))))),
    ]


def test_every_rule_firing_is_sound_along_curves(captured_rules):
    ef = _load_fuzz()
    L = ef.Leaves()
    rng = np.random.default_rng(5)
    gen = ef.Gen(rng, L)
    cases = _calculus_cases(L) + [gen.any(int(rng.integers(1, 5))) for _ in range(150)]
    counts, bad = collections.Counter(), []
    for i, expr in enumerate(cases):
        for fn in (expand, pull, vector_rules, simplify, full_simplify):
            captured_rules.clear()
            try:
                fn(expr)
            except ValueError:
                continue  # F1: a calculus node rebuilt over a manifold leaf
            ev = ef.JetEvaluator(np.random.default_rng(1000 + i), mode="consistent", env={})
            for name, before, after in list(captured_rules):
                counts[name] += 1
                try:
                    jb, ja = ev(before), ev(after)
                except (ef.Unsupported, ef.ShapeMismatch, ValueError, KeyError):
                    continue
                cb, ca = jb.c, ja.c
                if cb.shape != ca.shape:
                    if cb.size != ca.size:
                        assert not cb.any() and not ca.any(), (name, before, after)
                        continue
                    cb, ca = cb.reshape(-1), ca.reshape(-1)
                scale = max(1.0, float(np.abs(cb).max()), float(np.abs(ca).max()))
                err = float(np.abs(cb - ca).max()) / scale
                if err > 1e-6:  # jnormalize truncation is ~1e-7 at worst
                    bad.append((name, err, str(before)[:120], str(after)[:120]))
    assert not bad, bad[:5]
    # the phase-5a rules were exercised
    for name in (
        "pull.transpose_scaled",
        "pull.hat_scaled",
        "pull.vee_scaled",
        "pull.row_mvmul_scaled",
        "pull.vvmul_scaled",
        "expand.vee_over_madd",
        "simplify.mmmul_combine",
        "vector.s2_tangent_recovery",
    ):
        assert counts[name] > 0, name


# ---------------------------------------------------------------------------
# 4. Decision D1: named scalars with a value are symbols
# ---------------------------------------------------------------------------


def test_d1_named_scalars_stay_symbolic_in_every_pass():
    mv = _sc("m", value=2)
    lv = _sc("l", value=1)
    zv = _sc("z", value=0)
    v = Vector("v")
    e = VAdd(SVMul(v, Mul(mv, lv)), SVMul(v, zv), SVMul(v, Mul(Number(2), Number(0.5))))
    for fn in (expand, pull, vector_rules, simplify, full_simplify):
        out = fn(e)
        assert "m" in str(out) and "l" in str(out) and "z" in str(out), (fn.__name__, out)
        assert _close(e, out)
    assert "m" in to_latex(full_simplify(e)) and "l" in to_latex(full_simplify(e))


def test_d1_equal_values_do_not_merge_but_equal_names_do():
    p, r = _sc("p", value=2), _sc("r", value=2)
    assert str(full_simplify(p * x + r * x)) == "(p*x+r*x)"
    assert str(full_simplify(p - r)) == "(p+r*(-1))"
    assert str(full_simplify(p * x - x * p)) == "0"
    assert str(full_simplify(p * x + x * Number(2))) == "(p*x+x*(2))"
    # evaluate uses the values
    assert evaluate(full_simplify(p * x + r * x), {"x": 1.5}) == pytest.approx(6.0)
    assert evaluate(full_simplify(p - r)) == pytest.approx(0.0)


def test_d1_named_value_in_pipeline_standard_form_and_linearization():
    mv, lv, gv = _sc("m", value=2), _sc("l", value=1), _sc("g", value=9.81)
    e3 = _vc("e3")
    q, f = S2("q"), Vector("f")
    w, xi = q.get_tangent_vector(), q.get_variation_vector()
    L = HALF * mv * lv * lv * Dot(w, w) - mv * gv * lv * Dot(q, e3)
    var = SystemVariables(vectors=[q])
    eom = compute_eom(L, Dot(xi, f), var)
    eq = eom[str(xi)][1]
    assert str(eq) == "(\\frac{d}{dt}(\\omega_{q})*m*l*l*(-1)+Cross(q,e3)*m*g*l*(-1)+f)"
    assert check_eom(L, Dot(xi, f), var, eom, inputs=[f]) < 1e-6
    sf = to_standard_form(eom, var, [f])[str(xi)]
    (Mblock,) = sf.M.values()
    assert str(Mblock) == "I*m*l*l*(-1)" and str(sf.f) == "Cross(q,e3)*m*g*l*(-1)"
    assert "m" in to_latex(eq) and "g" in to_latex(eq)
    lin = linearize(eom, var, [f])
    assert finite_difference_check(lin, eom, var, [f]) < 1e-6
    blocks = [str(bl) for row in lin.A.values() for bl in row.values()]
    assert any("m*g*l" in s for s in blocks)


# ---------------------------------------------------------------------------
# 5. Printing / equality
# ---------------------------------------------------------------------------


def test_minus_one_forms():
    n_sym = Scalar("(-1)")  # a symbol named "(-1)", no value: not a literal
    assert str(full_simplify(x * Number(-1.0))) == str(full_simplify(x * Number(-1))) == "x*(-1)"
    # the symbol is not folded with the literal and not dropped
    out = full_simplify(x * n_sym * Number(-1))
    assert str(out) == "x*(-1)*(-1)"
    assert repr(x * n_sym) == "x·(-1)"  # printed as a name, not as a minus sign


# ---------------------------------------------------------------------------
# 6. Reviewer Lagrangians: numerical check and linearization
# ---------------------------------------------------------------------------


def _s2(name="q", fname="f"):
    q, f = S2(name), Vector(fname)
    return q, q.get_tangent_vector(), q.get_variation_vector(), f


def _so3(name="R", mname="M"):
    R, M = SO3(name), Vector(mname)
    return R, R.get_tangent_vector(), R.get_variation_vector(), M


def sys_s_skew_transpose():
    q, w, xi, f = _s2()
    K, c, g, e3 = _mc("K", "SkewSymmetricMatrix"), _vc("c"), _sc("g"), _vc("e3")
    L = HALF * Dot(w, w) + Dot(MVMul(Transpose(K), w), c) + Dot(MVMul(Transpose(K), q), w)
    return L - g * Dot(q, e3), Dot(xi, f), SystemVariables(vectors=[q]), [f]


def sys_s_identity_transpose():
    q, w, xi, f = _s2()
    c = _vc("c")
    L = HALF * Dot(w, MVMul(Transpose(IdentityMatrix), w))
    L = L + Dot(c, MVMul(MMMul(Transpose(IdentityMatrix), IdentityMatrix), w))
    return L, Dot(xi, f), SystemVariables(vectors=[q]), [f]


def sys_s2s2_double_pendulum():
    q1, w1, x1, f1 = _s2("q1", "f1")
    q2, w2, x2, f2 = _s2("q2", "f2")
    m1, m2, l1, l2, g, e3 = _sc("m1"), _sc("m2"), _sc("l1"), _sc("l2"), _sc("g"), _vc("e3")
    v1 = q1.t_diff() * l1
    v2 = v1 + q2.t_diff() * l2
    L = HALF * m1 * Dot(v1, v1) + HALF * m2 * Dot(v2, v2) - m1 * g * Dot(q1 * l1, e3)
    L = L - m2 * g * Dot(q1 * l1 + q2 * l2, e3)
    return L, Dot(x1, f1) + Dot(x2, f2), SystemVariables(vectors=[q1, q2]), [f1, f2]


def sys_s2s2_skew_coupling():
    q1, w1, x1, f1 = _s2("q1", "f1")
    q2, w2, x2, f2 = _s2("q2", "f2")
    K, k = _mc("K", "SkewSymmetricMatrix"), _sc("k")
    L = HALF * Dot(w1, w1) + HALF * Dot(w2, w2) + k * Dot(w1, K * w2)
    L = L + Dot(q1, Transpose(K) * q2)
    return L, Dot(x1, f1) + Dot(x2, f2), SystemVariables(vectors=[q1, q2]), [f1, f2]


def sys_s2_qdot_skew_nested():
    q, w, xi, f = _s2()
    K, A = _mc("K", "SkewSymmetricMatrix"), _mc("A")
    qd = q.t_diff()
    L = HALF * Dot(qd, qd) + Dot(qd, MVMul(MMMul(A, Transpose(K + K)), q))
    return L, Dot(xi, f), SystemVariables(vectors=[q]), [f]


def sys_s_like_terms_scaled():
    q, w, xi, f = _s2()
    a, b, k, m = _vc("a"), _vc("b"), _sc("k"), _sc("m")
    L = HALF * Dot(w, w) + k * m * Dot(Cross(a, w), b) + m * k * Dot(Cross(w, a), b)
    L = L + k * Dot(Cross(a, w) * m, b) * 2
    return L, Dot(xi, f), SystemVariables(vectors=[q]), [f]


def sys_so3_RtR_J():
    R, Om, eta, M = _so3()
    J = _mc("J", "SymmetricMatrix")
    L = HALF * Dot(Om, MVMul(MMMul(Transpose(R), MMMul(R, J)), Om))
    return L, Dot(eta, M), SystemVariables(matrices=[R]), [M]


def sys_so3_skew_sandwich():
    R, Om, eta, M = _so3()
    J, K, a, b = _mc("J", "SymmetricMatrix"), _mc("K", "SkewSymmetricMatrix"), _vc("a"), _vc("b")
    L = HALF * Dot(Om, J * Om) + Dot(a, MVMul(MMMul(MMMul(Transpose(R), K), R), b))
    return L, Dot(eta, M), SystemVariables(matrices=[R]), [M]


def sys_so3_hat_transpose_J():
    R, Om, eta, M = _so3()
    J, a, b = _mc("J", "SymmetricMatrix"), _vc("a"), _vc("b")
    chain = MMMul(MMMul(Transpose(Hat(Om)), J), Hat(Om))
    L = HALF * Dot(Om, J * Om) + Dot(a, MVMul(chain, b))
    return L, Dot(eta, M), SystemVariables(matrices=[R]), [M]


def sys_so3_scaled_hat():
    R, Om, eta, M = _so3()
    J, k, a, b = _mc("J", "SymmetricMatrix"), _sc("k"), _vc("a"), _vc("b")
    L = HALF * Dot(Om, J * Om) + Dot(a, MVMul(Hat(Om * k), MVMul(R, b))) + Dot(Vee(Hat(Om) * k), a)
    return L, Dot(eta, M), SystemVariables(matrices=[R]), [M]


def sys_so3_row_scaled():
    R, Om, eta, M = _so3()
    J, k, a, b = _mc("J", "SymmetricMatrix"), _sc("k"), _vc("a"), _vc("b")
    L = HALF * Dot(Om, J * Om) + VVMul(Transpose(Om * k), MVMul(Transpose(R), a))
    L = L + VVMul(Transpose(a), (R * b) * k)
    return L, Dot(eta, M), SystemVariables(matrices=[R]), [M]


def sys_so3so3_coupling():
    R1, O1, e1, M1 = _so3("R1", "M1")
    R2, O2, e2, M2 = _so3("R2", "M2")
    J1, J2 = _mc("J1", "SymmetricMatrix"), _mc("J2", "SymmetricMatrix")
    k, a, b = _sc("k"), _vc("a"), _vc("b")
    L = HALF * Dot(O1, J1 * O1) + HALF * Dot(O2, J2 * O2) + k * Dot(R1 * O1, R2 * O2)
    L = L + Dot(a, MVMul(MMMul(Transpose(R1), R2), b))
    return L, Dot(e1, M1) + Dot(e2, M2), SystemVariables(matrices=[R1, R2]), [M1, M2]


def sys_se3_offset():
    R, Om, eta, M = _so3()
    x, F = Vector("x"), Vector("F")
    J, m, g, e3, rho = _mc("J", "SymmetricMatrix"), _sc("m"), _sc("g"), _vc("e3"), _vc("rho")
    xd = x.t_diff()
    L = HALF * m * Dot(xd, xd) + HALF * Dot(Om, J * Om) - m * g * Dot(e3, x)
    L = L + m * Dot(xd, MVMul(MMMul(R, Hat(Om)), rho))
    return L, Dot(x.delta(), F) + Dot(eta, M), SystemVariables(vectors=[x], matrices=[R]), [F, M]


def sys_r2_coupled_named():
    s1, s2, u1, u2 = Scalar("s1"), Scalar("s2"), Scalar("u1"), Scalar("u2")
    k, kv = _sc("k"), _sc("kv", value=2)
    d1, d2 = s1.t_diff(), s2.t_diff()
    L = HALF * d1 * d1 + HALF * d2 * d2 + k * d1 * d2 * s1 - HALF * kv * s1 * s2
    L = L + d2 * s1 * k - k * s1 * d2
    return L, s1.delta() * u1 + s2.delta() * u2, SystemVariables(scalars=[s1, s2]), [u1, u2]


def sys_r_so3_mixed():
    s, fs = Scalar("s"), Scalar("f_s")
    R, Om, eta, M = _so3()
    J, a = _mc("J", "SymmetricMatrix"), _vc("a")
    sd = s.t_diff()
    L = HALF * s * Dot(Om, J * Om) + HALF * sd * sd + sd * Dot(a, Om)
    return L, s.delta() * fs + Dot(eta, M), SystemVariables(scalars=[s], matrices=[R]), [fs, M]


def sys_r3_row_forms_skew():
    x, F = Vector("x"), Vector("F")
    Mm, K, A = _mc("Mm", "SymmetricMatrix"), _mc("K", "SkewSymmetricMatrix"), _mc("A")
    xd = x.t_diff()
    L = HALF * VVMul(Transpose(xd), Mm * xd) + VVMul(Transpose(x), K * xd)
    L = L + HALF * Dot(xd, A * xd) + Dot(x, MVMul(Transpose(K), xd))
    return L, Dot(x.delta(), F), SystemVariables(vectors=[x]), [F]


def sys_r3_s2_so3_mixed():
    x, F = Vector("x"), Vector("F")
    q, w, xi, f = _s2()
    R, Om, eta, M = _so3()
    m, l, J = _sc("m"), _sc("l"), _mc("J", "SymmetricMatrix")
    v = x.t_diff() + MVMul(R, Cross(q, w)) * l
    L = HALF * m * Dot(v, v) + HALF * Dot(Om, J * Om) + HALF * Dot(w, w)
    dW = Dot(x.delta(), F) + Dot(xi, f) + Dot(eta, M)
    return L, dW, SystemVariables(vectors=[x, q], matrices=[R]), [F, f, M]


REVIEW_SYSTEMS = [
    sys_s_skew_transpose,
    sys_s_identity_transpose,
    sys_s2s2_double_pendulum,
    sys_s2s2_skew_coupling,
    sys_s2_qdot_skew_nested,
    sys_s_like_terms_scaled,
    sys_so3_RtR_J,
    sys_so3_skew_sandwich,
    sys_so3_hat_transpose_J,
    sys_so3_scaled_hat,
    sys_so3_row_scaled,
    sys_so3so3_coupling,
    sys_se3_offset,
    sys_r2_coupled_named,
    sys_r_so3_mixed,
    sys_r3_row_forms_skew,
    sys_r3_s2_so3_mixed,
]


@pytest.mark.parametrize("build", REVIEW_SYSTEMS, ids=lambda f: f.__name__[4:])
def test_review_system_matches_oracle_and_linearizes(build):
    L, dW, var, inputs = build()
    eom = compute_eom(L, dW, var)
    assert check_eom(L, dW, var, eom, inputs=inputs) < 1e-6
    lin = linearize(eom, var, inputs)
    assert finite_difference_check(lin, eom, var, inputs) < 1e-6
    to_standard_form(eom, var, inputs)


def test_like_terms_merge_in_pipeline():
    """The three copies of the k m (a×ω)·b term merge to one (4 k m)."""
    L, dW, var, _ = sys_s_like_terms_scaled()
    eom = compute_eom(L, dW, var)
    (eq,) = [e for _, e in eom.values()]
    assert str(eq) == "(\\frac{d}{dt}(\\omega_{q})*(-1)+Cross(\\omega_{q},Cross(a,b))*k*m*(-4)+f)"


# ---------------------------------------------------------------------------
# Findings (were strict xfails; fixed, see the fix log of the review)
# ---------------------------------------------------------------------------

_q = S2("q")


@pytest.mark.parametrize(
    "e",
    [
        TimeDerivative(SVMul(_q, Dot(_q, _q))),
        Variation(SVMul(_q, Dot(_q, _q))),
        TimeDerivative(VAdd(_q, Cross(_q, Cross(_q, _q)))),
    ],
    ids=["dt_q_qq", "delta_q_qq", "dt_q_plus_qxqxq"],
)
def test_calculus_node_over_reduced_manifold_leaf_does_not_raise(e):
    full_simplify(e)


_R = SO3("R")


@pytest.mark.parametrize(
    "e, eager",
    [
        # d/dt(q·1) and δ(q·1) raised before 5a already (via simplify's zero/one elimination)
        (TimeDerivative(SVMul(_q, Number(1))), _q.t_diff()),
        (Variation(SVMul(_q, Number(1))), _q.delta()),
        (TimeDerivative(SVMul(_q, Dot(_q, _q))), _q.t_diff()),
        (Variation(SVMul(_q, Dot(_q, _q))), _q.delta()),
        (TimeDerivative(MVMul(MMMul(Transpose(_R), _R), _q)), _q.t_diff()),
        (Variation(MVMul(IdentityMatrix, _q)), _q.delta()),
        (TimeDerivative(SMMul(_R, Number(1))), _R.t_diff()),
        (Variation(MMMul(IdentityMatrix, _R)), _R.delta()),
    ],
    ids=[
        "dt_q_one",
        "delta_q_one",
        "dt_q_qq",
        "delta_q_qq",
        "dt_RtR_q",
        "delta_I_q",
        "dt_R_one",
        "delta_I_R",
    ],
)
def test_calculus_node_over_reduced_manifold_leaf_is_eager(e, eager):
    """The rebuilt calculus node over a bare S2/SO3 leaf is the leaf's own kinematics
    (ω×q, ξ×q, R·Hat(Ω), R·Hat(η)), never the node the constructor rejects."""
    for fn in (expand, pull, vector_rules, simplify, full_simplify):
        out = fn(e)
        assert "frac{d}{dt}(q)" not in str(out) and "delta{q}" not in str(out), (fn, out)
        assert "frac{d}{dt}(R)" not in str(out) and "delta{R}" not in str(out), (fn, out)
    assert str(full_simplify(e)) == str(full_simplify(eager))


def test_integer_valued_float_normalised_for_vectors_and_matrices():
    v = Vector("v")
    assert str(full_simplify(v * Number(2.0))) == str(full_simplify(v * Number(2)))
    assert str(full_simplify(A * Number(-1.0))) == str(full_simplify(A * Number(-1)))


def test_huge_integer_literal_folds():
    assert "x" in str(full_simplify(Number(10**400) * x))


def test_huge_integer_literal_is_a_literal_and_folds_exactly():
    from geomech.core.base.expressions import is_number_literal

    assert is_number_literal(Number(10**400))
    assert not is_number_literal(Scalar("m", value=10**400))
    assert str(full_simplify(Number(10**400) * x * Number(10**400))) == f"x*({10**800})"


def test_repr_of_named_constant_with_value_keeps_marker():
    assert repr(_sc("m", value=2)) == repr(_sc("m"))


def test_repr_marks_named_values_but_not_literals():
    assert repr(_sc("m", value=2)) == "m*"
    assert repr(Scalar("m", value=2)) == "m*"  # a value makes it constant in time
    assert repr(Number(2)) == "(2)" and repr(Scalar("0.5", value=0.5)) == "0.5"
    assert repr(Vector("v") * _sc("m", value=2)) == "v:v·m*"


def _fuzz_inserted_path():
    """The directory engine_fuzz.py puts first on sys.path (its line 55 expression)."""
    with open(FUZZ_PATH, encoding="utf-8") as fh:
        src = fh.read()
    line = next(ln for ln in src.splitlines() if ln.startswith("sys.path.insert(0, "))
    expr = line[len("sys.path.insert(0, ") : -1]
    return eval(expr, {"os": os, "__file__": FUZZ_PATH})  # noqa: S307


def test_engine_fuzz_imports_its_own_tree():
    assert os.path.abspath(_fuzz_inserted_path()) == REPO
