#!/usr/bin/env python
"""Property-based numerical fuzzing of the geomech symbolic engine.

Every value-preserving transformation (expand, pull, vector_rules, simplify,
full_simplify, collect, apply_manifold_rules, the delta()/t_diff() product
rules, ibp._t_diff_product) is applied to thousands of random, well-typed
expression trees and the result is compared numerically with the input.
extract_linear_coeff is checked against the finite difference
expr(query=v) - expr(query=0) == coeff . v.

Numeric model
-------------
Every leaf is a curve x(t, eps) represented as a truncated bivariate Taylor
series ("jet") of order KT in t and KE in eps.  Products of jets are truncated
Cauchy products, ``TimeDerivative`` is the exact coefficient shift d/dt and
``Variation`` the exact shift d/d(eps), so calculus nodes are evaluated exactly
(no finite differences).  Manifold leaves are sampled as curves that satisfy
their constraints to all orders:

  S2  q      : q = p / |p| for a random jet p; omega := q x dq/dt, xi := q x dq/deps
               so that dq/dt = omega x q and dq/deps = xi x q hold exactly.
  SO3 R      : R = R0 (I - A)^-1 (I + A) (Cayley) with A = hat(d), d a nilpotent jet;
               Omega := vee(R^T dR/dt), eta := vee(R^T dR/deps).
  UnitNorm v : normalised jet.   Symmetric J : A + A^T.   Orthogonal Q : Cayley.
  Constant   : only the (0, 0) coefficient is non-zero.

In "adversarial" mode the flagged leaves violate their declared properties
(omega not orthogonal to q, |u| != 1, R not orthogonal), which shows which
rules rely on which declared property.

Usage
-----
    source venv/bin/activate
    python specs/reviews/scripts/engine_fuzz.py [--cases N] [--seed S] [--quick]

Exit status is 0; findings are printed and written as JSON next to the script
(engine_fuzz_findings.json) unless --no-json is given.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import signal
import sys
import time
import traceback

import numpy as np

import geomech as _geomech
from geomech import (
    S2,
    SO3,
    TS2,
    TSO3,
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
    One,
    Scalar,
    SMMul,
    SVMul,
    SystemVariables,
    TimeDerivative,
    TimeIntegral,
    Transpose,
    VAdd,
    Variation,
    Vector,
    Vee,
    VVMul,
    Zero,
    ZeroMatrix,
    ZeroVector,
    collect,
    expand,
    extract_linear_coeff,
    full_simplify,
    pull,
    simplify,
    vector_rules,
)
from geomech.core.base.types import ExprType
from geomech.core.math.ibp import _t_diff_product
from geomech.core.transformations.manifold_rules import apply_manifold_rules
from geomech.utils.errors import (
    AlgebraicError,
    ExpressionMismatchError,
    NonLinearError,
    SizeMismatchError,
    UndefinedCaseError,
)

DOCUMENTED_ERRORS = (
    ExpressionMismatchError,
    SizeMismatchError,
    AlgebraicError,
    NonLinearError,
    UndefinedCaseError,
)
# Raised on purpose by the engine but not part of utils/errors.py
SEMI_DOCUMENTED_ERRORS = (NotImplementedError, ValueError)

RTOL = 1e-8
CASE_TIMEOUT_S = 10

# ---------------------------------------------------------------------------
# Jets: truncated bivariate Taylor series in (t, eps)
# ---------------------------------------------------------------------------

KT, KE = 5, 2  # derived leaves (omega = q x dq/dt) lose the top t-order, so trees with
# nested d/dt depth <= KT-1 (omega/Omega/xi/eta count as one level) are evaluated exactly
DEGS = [(p, q) for p in range(KT + 1) for q in range(KE + 1)]
NDEG = len(DEGS)
DIDX = {d: i for i, d in enumerate(DEGS)}

_I, _J, _K = [], [], []
for i, (p1, q1) in enumerate(DEGS):
    for j, (p2, q2) in enumerate(DEGS):
        if p1 + p2 <= KT and q1 + q2 <= KE:
            _I.append(i)
            _J.append(j)
            _K.append(DIDX[(p1 + p2, q1 + q2)])
_I, _J = np.array(_I), np.array(_J)
_S = np.zeros((len(_K), NDEG))
_S[np.arange(len(_K)), _K] = 1.0
_DT = np.zeros((NDEG, NDEG))
_DE = np.zeros((NDEG, NDEG))
for i, (p, q) in enumerate(DEGS):
    if p + 1 <= KT:
        _DT[DIDX[(p + 1, q)], i] = p + 1
    if q + 1 <= KE:
        _DE[DIDX[(p, q + 1)], i] = q + 1


class ShapeMismatch(Exception):
    """Row/column (or scalar/matrix) shape confusion between operands."""


class Jet:
    __slots__ = ("c",)

    def __init__(self, c):
        self.c = np.asarray(c, dtype=float)

    def _same_shape(self, o):
        """Operands of + must agree in shape; an identically-zero operand (e.g. the engine's
        Transpose(ZeroVector)) is accepted in either orientation."""
        if self.c.shape != o.c.shape:
            if not o.c.any():
                return Jet(self.c.copy())
            if not self.c.any():
                return Jet(o.c.copy())
            raise ShapeMismatch(f"{self.c.shape[:2]} vs {o.c.shape[:2]}")
        return None

    @property
    def shape(self):
        return self.c.shape[:2]

    def val(self):
        return self.c[:, :, 0]

    def __add__(self, o):
        z = self._same_shape(o)
        if z is not None:
            return z
        return Jet(self.c + o.c)

    def __sub__(self, o):
        z = self._same_shape(o)
        if z is not None:
            return Jet(-z.c) if not self.c.any() and o.c.any() else z
        return Jet(self.c - o.c)

    def __neg__(self):
        return Jet(-self.c)

    def T(self):
        return Jet(self.c.transpose(1, 0, 2))

    def dt(self):
        return Jet(self.c @ _DT)

    def de(self):
        return Jet(self.c @ _DE)

    def col(self):
        """Make a vector jet a column (3,1)."""
        if self.shape == (1, 3):
            return self.T()
        return self


def jconst(arr):
    arr = np.atleast_2d(np.asarray(arr, dtype=float))
    c = np.zeros(arr.shape + (NDEG,))
    c[:, :, 0] = arr
    return Jet(c)


def jmatmul(A, B):
    prod = np.einsum("ikn,kjn->ijn", A.c[:, :, _I], B.c[:, :, _J])
    return Jet(prod @ _S)


def jsmul(s, B):
    """scalar jet * any jet."""
    prod = s.c[0, 0, _I][None, None, :] * B.c[:, :, _J]
    return Jet(prod @ _S)


def jhat(v):
    v = v.col()
    c = v.c
    H = np.zeros((3, 3, NDEG))
    H[0, 1] = -c[2, 0]
    H[0, 2] = c[1, 0]
    H[1, 0] = c[2, 0]
    H[1, 2] = -c[0, 0]
    H[2, 0] = -c[1, 0]
    H[2, 1] = c[0, 0]
    return Jet(H)


def jvee(M):
    c = M.c
    return Jet(np.stack([c[2, 1], c[0, 2], c[1, 0]])[:, None, :])


def jcross(a, b):
    return jmatmul(jhat(a), b.col())


def jdot(a, b):
    return jmatmul(a.col().T(), b.col())


def jrsqrt(x):
    """x^(-1/2) for a scalar jet via Taylor expansion around the (0,0) coefficient."""
    x0 = float(x.c[0, 0, 0])
    d = Jet(x.c.copy())
    d.c[0, 0, 0] = 0.0
    out = jconst([[x0**-0.5]])
    dk = jconst([[1.0]])
    coeff = x0**-0.5
    for k in range(1, KT + KE + 1):
        coeff *= (-0.5 - (k - 1)) / x0
        dk = jsmul(dk, d)
        out = out + Jet(dk.c * coeff / math.factorial(k))
    return out


def jnormalize(p):
    return jsmul(jrsqrt(jdot(p, p)), p)


def rand_jet(rng, shape, decay=0.6, constant=False):
    c = rng.normal(size=shape + (NDEG,))
    if constant:
        c[:, :, 1:] = 0.0
    else:
        for i, (p, q) in enumerate(DEGS):
            c[:, :, i] *= decay ** (p + q)
    return Jet(c)


def rand_rot(rng):
    Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    if np.linalg.det(Q) < 0:
        Q[:, 0] = -Q[:, 0]
    return Q


def rand_orthogonal_jet(rng):
    """R(t, eps) in SO(3) to all orders: R0 (I - A)^-1 (I + A), A = hat(d), d nilpotent."""
    d = rand_jet(rng, (3, 1))
    d.c[:, :, 0] = 0.0
    A = jhat(d)
    Id = jconst(np.eye(3))
    S = Id
    Ak = Id
    for _ in range(KT + KE):
        Ak = jmatmul(Ak, A)
        S = S + Ak
    return jmatmul(jconst(rand_rot(rng)), jmatmul(S, Id + A))


# ---------------------------------------------------------------------------
# Jet evaluator for expression trees
# ---------------------------------------------------------------------------


class Unsupported(Exception):
    pass


class JetEvaluator:
    """Evaluate an expression tree to a Jet.

    mode="consistent": every flagged leaf satisfies its declared property to all
    orders.  mode="adversarial": flagged leaves are random (properties violated).
    ``overrides`` maps str(node) -> Jet and takes precedence over structural
    evaluation (used to set the query of extract_linear_coeff to v or 0).
    """

    def __init__(self, rng, mode="consistent", env=None, overrides=None, gauge="full"):
        self.rng = rng
        self.mode = mode
        self.env = {} if env is None else env
        self.overrides = {} if overrides is None else overrides
        # gauge="tangential": Variation(omega) for an S2 tangent vector omega is evaluated as
        # the projection of d(omega)/d(eps) onto the tangent plane of q (drops the component
        # along q).  This was the representative that manifold_rules substituted before
        # fix plan 1.1 (xi_dot - omega x xi); the exact rule is now xi_dot - 2 omega x xi,
        # so the reference for the manifold_rules campaign is gauge="full".
        self.gauge = gauge

    # -- leaves -----------------------------------------------------------
    def _sample_s2(self, node):
        key = str(node)
        if key in self.env:
            return self.env[key]
        q = jnormalize(rand_jet(self.rng, (3, 1)))
        self.env[key] = q
        wname = node.manifold_info.tangent_vector_name
        xname = node.manifold_info.variation_vector_name
        if self.mode == "consistent":
            self.env[wname] = jcross(q, q.dt())
            self.env[xname] = jcross(q, q.de())
        else:
            self.env[wname] = rand_jet(self.rng, (3, 1))
            self.env[xname] = rand_jet(self.rng, (3, 1))
        return q

    def _sample_so3(self, node):
        key = str(node)
        if key in self.env:
            return self.env[key]
        if self.mode == "consistent":
            R = rand_orthogonal_jet(self.rng)
        else:
            R = rand_jet(self.rng, (3, 3))
        self.env[key] = R
        wname = node.manifold_info.tangent_vector_name
        xname = node.manifold_info.variation_vector_name
        self.env[wname] = jvee(jmatmul(R.T(), R.dt()))
        self.env[xname] = jvee(jmatmul(R.T(), R.de()))
        return R

    def leaf(self, node):
        if isinstance(node, S2):
            return self._sample_s2(node)
        if isinstance(node, SO3):
            return self._sample_so3(node)
        if isinstance(node, TS2) and node.S2 is not None:
            self._sample_s2(node.S2)
            return self.env[str(node)]
        if isinstance(node, TSO3) and node.SO3 is not None:
            self._sample_so3(node.SO3)
            return self.env[str(node)]
        if node.is_zero:
            if isinstance(node, Matrix):
                return jconst(np.zeros((3, 3)))
            if isinstance(node, Vector):
                return jconst(np.zeros((3, 1)))
            return jconst([[0.0]])
        if isinstance(node, Matrix) and node.flags.is_identity:
            return jconst(np.eye(3))
        key = str(node)
        if key in self.env:
            return self.env[key]
        rng = self.rng
        consistent = self.mode == "consistent"
        if isinstance(node, Scalar):
            if node.value is not None:
                val = jconst([[float(node.value)]])
            else:
                val = rand_jet(rng, (1, 1), constant=node.is_constant)
        elif isinstance(node, Matrix):
            if node.is_orthogonal and consistent:
                val = rand_orthogonal_jet(rng)
            else:
                val = rand_jet(rng, (3, 3), constant=node.is_constant)
                if node.is_symmetric and consistent:
                    val = val + val.T()
                if node.is_skew_symmetric and consistent:
                    val = val - val.T()
        elif isinstance(node, Vector):
            val = rand_jet(rng, (3, 1), constant=node.is_constant)
            if node.is_unit_norm and consistent:
                val = jnormalize(val)
        else:
            raise Unsupported(f"unknown leaf {type(node)} {key}")
        self.env[key] = val
        return val

    # -- tree ----------------------------------------------------------------
    def __call__(self, node):
        key = str(node)
        if self.overrides and key in self.overrides:
            return self.overrides[key]
        if isinstance(node, (Add, VAdd, MAdd)):
            out = self(node.nodes[0])
            for n in node.nodes[1:]:
                out = out + self(n)
            return out
        if isinstance(node, Mul):
            return jsmul(self(node.left), self(node.right))
        if isinstance(node, (SVMul, SMMul)):
            return jsmul(self(node.right), self(node.left))
        if isinstance(node, (MVMul, MMMul, VVMul)):
            return jmatmul(self(node.left), self(node.right))
        if isinstance(node, Dot):
            return jdot(self(node.left), self(node.right))
        if isinstance(node, Cross):
            return jcross(self(node.left), self(node.right))
        if isinstance(node, Hat):
            return jhat(self(node.expr))
        if isinstance(node, Vee):
            return jvee(self(node.expr))
        if isinstance(node, Transpose):
            return self(node.expr).T()
        if isinstance(node, TimeDerivative):
            return self(node.expr).dt()
        if isinstance(node, Variation):
            inner = node.expr
            full = self(inner).de()
            if (
                self.gauge == "tangential"
                and isinstance(inner, TS2)
                and inner.S2 is not None
                and inner.name == inner.S2.manifold_info.tangent_vector_name
            ):
                q = self(inner.S2)
                return full - jsmul(jdot(q, full), q)
            return full
        if isinstance(node, TimeIntegral):
            raise Unsupported("TimeIntegral cannot be evaluated pointwise")
        return self.leaf(node)


def values_close(a, b, rtol=RTOL):
    """Compare the (t=0, eps=0) values of two jets (shape-insensitive for vectors)."""
    va, vb = np.asarray(a.val()), np.asarray(b.val())
    if va.shape != vb.shape:
        if va.size == vb.size:
            va, vb = va.reshape(-1), vb.reshape(-1)
        else:
            return False, float("inf")
    scale = max(1.0, float(np.max(np.abs(va))), float(np.max(np.abs(vb))))
    err = float(np.max(np.abs(va - vb))) / scale
    return err <= rtol, err


# ---------------------------------------------------------------------------
# Random expression generator
# ---------------------------------------------------------------------------


def contains_variation(expr):
    """True if expr holds a Variation node or a variation vector (xi/eta) leaf."""
    if isinstance(expr, Variation):
        return True
    if isinstance(expr, TS2) and expr.S2 is not None:
        if expr.name == expr.S2.manifold_info.variation_vector_name:
            return True
    if isinstance(expr, TSO3) and expr.SO3 is not None:
        if expr.name == expr.SO3.manifold_info.variation_vector_name:
            return True
    for n in getattr(expr, "nodes", None) or []:
        if contains_variation(n):
            return True
    return False


def contains_calculus(expr):
    if isinstance(expr, (Variation, TimeDerivative, TimeIntegral)):
        return True
    return any(contains_calculus(n) for n in getattr(expr, "nodes", None) or [])


def tree_size(expr):
    return 1 + sum(tree_size(n) for n in getattr(expr, "nodes", None) or [])


def tdepth(expr):
    """Nested d/dt depth; tangent/variation vectors of manifold leaves count as one level
    because their jets are derived from the manifold point's jet."""
    if isinstance(expr, (TS2, TSO3)):
        return 1
    base = 1 if isinstance(expr, TimeDerivative) else 0
    kids = getattr(expr, "nodes", None) or []
    return base + max((tdepth(n) for n in kids), default=0)


def low_orders_close(a, b, atol=1e-10):
    """Compare all coefficients that are exact for derived leaves (p <= KT-2, q <= KE-1)."""
    idx = [i for i, (p, q) in enumerate(DEGS) if p <= KT - 2 and q <= KE - 1]
    return np.allclose(a.c[:, :, idx], b.c[:, :, idx], atol=atol)


class Leaves:
    def __init__(self):
        self.s1, self.s2 = Scalar("s_1"), Scalar("s_2")
        self.m = Scalar("m", attr=["Constant"])
        self.kval = Scalar("k", value=2.5)
        self.q = S2("q")
        self.p = S2("p")
        self.w = self.q.get_tangent_vector()
        self.xi = self.q.get_variation_vector()
        self.wp = self.p.get_tangent_vector()
        self.xip = self.p.get_variation_vector()
        self.a, self.b, self.c = Vector("a"), Vector("b"), Vector("c")
        self.e3 = Vector("e_3", attr=["Constant"])
        self.u = Vector("u", attr=["UnitNorm"])
        self.M, self.N = Matrix("M"), Matrix("N")
        self.J = Matrix("J", attr=["SymmetricMatrix", "Constant"])
        self.K = Matrix("K", attr=["SkewSymmetricMatrix"])
        self.Q = Matrix("Q", attr=["OrthogonalMatrix"])
        self.R = SO3("R")
        self.Om = self.R.get_tangent_vector()
        self.eta = self.R.get_variation_vector()
        from geomech import I as Ident

        self.Ident = Ident

        self.scalars = [
            self.s1,
            self.s2,
            self.m,
            self.kval,
            Number(2),
            Number(-1),
            Number(0.5),
            Zero,
            One,
        ]
        self.vectors = [
            self.a,
            self.b,
            self.c,
            self.e3,
            self.u,
            self.q,
            self.w,
            self.xi,
            ZeroVector,
        ]
        self.matrices = [self.M, self.N, self.J, self.K, self.Q, self.R, self.Ident, ZeroMatrix]
        self.variation_leaves = {self.xi, self.xip, self.eta}


class Gen:
    """Typed random tree generator.

    opts:
      calc      - allow TimeDerivative/Variation nodes (on leaves and subtrees)
      variation - allow Variation nodes (requires calc)
      motifs    - add S2/SO3 structural motifs that exercise vector_rules/simplify
      rowvec    - allow Transpose(MVMul(Transpose(v), M)) row-vector motif
      second_s2 - include the second S2 point p
    """

    def __init__(
        self, rng, L: Leaves, calc=True, variation=True, motifs=True, rowvec=True, second_s2=True
    ):
        self.rng = rng
        self.L = L
        self.calc = calc
        self.variation = variation and calc
        self.motifs = motifs
        self.rowvec = rowvec
        self.second_s2 = second_s2

    def pick(self, seq):
        return seq[self.rng.integers(len(seq))]

    def chance(self, p):
        return self.rng.random() < p

    # -- leaves ------------------------------------------------------------
    def scalar_leaf(self):
        s = self.pick(self.L.scalars)
        return self.maybe_calc(s)

    def vector_leaf(self):
        pool = list(self.L.vectors)
        if self.second_s2:
            pool += [self.L.p, self.L.wp, self.L.xip]
        pool += [self.L.Om, self.L.eta]
        v = self.pick(pool)
        return self.maybe_calc(v)

    def matrix_leaf(self):
        return self.maybe_calc(self.pick(self.L.matrices))

    def maybe_calc(self, leaf):
        if not self.calc or leaf.is_constant or self.chance(0.65):
            return leaf
        if isinstance(leaf, (S2, SO3)):
            # manifold leaves are eager: use their kinematic trees
            return leaf.t_diff() if self.chance(0.5) else leaf.delta()
        if leaf in self.L.variation_leaves:
            return TimeDerivative(leaf)
        r = self.rng.random()
        if r < 0.45:
            return TimeDerivative(leaf)
        if r < 0.6:
            return TimeDerivative(TimeDerivative(leaf))
        if self.variation:
            if r < 0.85:
                return Variation(leaf)
            return TimeDerivative(Variation(leaf))
        return TimeDerivative(leaf)

    def wrap_calc(self, e):
        """Wrap a subtree in TimeDerivative or Variation when legal."""
        if self.variation and not contains_variation(e) and self.chance(0.5):
            return Variation(e)
        return TimeDerivative(e)

    def plain_vector(self, depth):
        """A vector subtree that is not headed by Transpose (VVMul operand)."""
        for _ in range(6):
            v = self.vector(depth)
            if not isinstance(v, Transpose):
                return v
        return self.pick([self.L.a, self.L.b, self.L.c])

    # -- typed trees ---------------------------------------------------------
    def scalar(self, depth):
        if depth <= 0 or self.chance(0.2):
            return self.scalar_leaf()
        L = self.L
        r = self.rng.random()
        d = depth - 1
        if r < 0.18:
            return Add(*[self.scalar(d) for _ in range(self.rng.integers(2, 4))])
        if r < 0.36:
            return Mul(self.scalar(d), self.scalar(d))
        if r < 0.62:
            return Dot(self.vector(d), self.vector(d))
        if r < 0.70:
            return VVMul(Transpose(self.plain_vector(d)), self.plain_vector(d))
        if r < 0.78 and self.calc:
            return self.wrap_calc(self.scalar(d))
        if self.motifs:
            motifs = [
                lambda: Dot(L.q, L.xi),
                lambda: Dot(L.w, L.q),
                lambda: Dot(L.q, L.q),
                lambda: Dot(L.u, L.u),
                lambda: Dot(L.q, Cross(L.w, L.q)),
                lambda: Dot(Cross(L.w, L.q), L.q),
                lambda: Dot(L.q, TimeDerivative(L.w)),
                lambda: Dot(self.vector(d), Cross(self.vector(d), self.vector(d))),
                lambda: Dot(Cross(L.a, L.b), L.b),
                lambda: Dot(L.a, Cross(L.b, L.b)),
                lambda: Dot(MVMul(L.R, self.vector(d)), MVMul(L.R, self.vector(d))),
                lambda: Dot(L.p, L.xi),  # xi of q against p: rule must not fire
                lambda: Dot(SVMul(L.q, self.scalar(d)), SVMul(L.q, self.scalar(d))),
                lambda: Mul(Dot(L.q, L.q), Dot(L.w, L.q)),
            ]
            return self.pick(motifs)()
        return Dot(self.vector(d), self.vector(d))

    def vector(self, depth):
        if depth <= 0 or self.chance(0.2):
            return self.vector_leaf()
        L = self.L
        r = self.rng.random()
        d = depth - 1
        if r < 0.15:
            return VAdd(*[self.vector(d) for _ in range(self.rng.integers(2, 4))])
        if r < 0.30:
            return SVMul(self.vector(d), self.scalar(d))
        if r < 0.45:
            return MVMul(self.matrix(d), self.vector(d))
        if r < 0.62:
            return Cross(self.vector(d), self.vector(d))
        if r < 0.68:
            return Vee(self.matrix(d))
        if r < 0.74 and self.calc:
            return self.wrap_calc(self.vector(d))
        if r < 0.78 and self.rowvec:
            return Transpose(MVMul(Transpose(self.vector(d)), self.matrix(d)))
        if self.motifs:
            motifs = [
                lambda: Cross(L.w, L.q),
                lambda: Cross(L.xi, L.q),
                lambda: Cross(L.q, Cross(L.w, L.q)),
                lambda: Cross(L.q, Cross(L.q, L.w)),
                lambda: Cross(Cross(L.w, L.q), L.q),
                lambda: Cross(Cross(L.q, L.xi), L.q),
                lambda: Cross(L.w, Cross(L.w, L.q)),
                lambda: Cross(L.xi, Cross(L.xi, L.q)),
                lambda: Cross(L.q, Cross(self.vector(d), L.q)),  # generic a: not tangent
                lambda: Cross(L.q, Cross(L.q, self.vector(d))),
                lambda: Cross(L.q, Cross(Cross(self.vector(d), L.q), L.q)),
                lambda: Cross(L.q, Cross(TimeDerivative(L.w), L.q)),
                lambda: Cross(L.q, Cross(Variation(L.w), L.q)),
                lambda: Cross(L.q, Cross(L.wp, L.q)),  # tangent of p, not q
                lambda: Cross(L.p, Cross(L.w, L.p)),
                lambda: Cross(L.u, Cross(L.a, L.u)),  # unit-norm but not S2
                lambda: Cross(self.vector(d), SVMul(self.vector(d), self.scalar(d))),
                lambda: Cross(L.a, SVMul(L.a, self.scalar(d))),
                lambda: MVMul(Hat(self.vector(d)), self.vector(d)),
                lambda: MVMul(MMMul(Transpose(L.R), L.R), self.vector(d)),
                lambda: MVMul(Transpose(L.R), MVMul(L.R, self.vector(d))),
                lambda: MVMul(L.R, MVMul(Transpose(L.R), self.vector(d))),
                lambda: MVMul(Transpose(L.Q), MVMul(L.Q, self.vector(d))),
                lambda: MVMul(MMMul(self.matrix(d), self.matrix(d)), self.vector(d)),
                lambda: MVMul(L.Ident, self.vector(d)),
                lambda: Vee(Hat(self.vector(d))),
                lambda: Vee(MAdd(Hat(self.vector(d)), Hat(self.vector(d)))),
                lambda: MVMul(Transpose(Hat(self.vector(d))), self.vector(d)),
                lambda: Cross(Cross(L.q, Cross(L.w, L.q)), L.w),
                lambda: SVMul(L.q, Dot(L.q, self.vector(d))),
                lambda: Cross(L.q, Cross(L.a, L.a)),
            ]
            return self.pick(motifs)()
        return Cross(self.vector(d), self.vector(d))

    def matrix(self, depth):
        if depth <= 0 or self.chance(0.25):
            return self.matrix_leaf()
        L = self.L
        r = self.rng.random()
        d = depth - 1
        if r < 0.15:
            return MAdd(*[self.matrix(d) for _ in range(self.rng.integers(2, 4))])
        if r < 0.28:
            return SMMul(self.matrix(d), self.scalar(d))
        if r < 0.45:
            return MMMul(self.matrix(d), self.matrix(d))
        if r < 0.60:
            return Hat(self.vector(d))
        if r < 0.70:
            return Transpose(self.matrix(d))
        if r < 0.78:
            return VVMul(self.plain_vector(d), Transpose(self.plain_vector(d)))
        if r < 0.84 and self.calc:
            return self.wrap_calc(self.matrix(d))
        if self.motifs:
            motifs = [
                lambda: MMMul(Transpose(L.R), L.R),
                lambda: MMMul(L.R, Transpose(L.R)),
                lambda: MMMul(Transpose(L.Q), L.Q),
                lambda: MMMul(Transpose(L.M), L.M),  # generic: must not become I
                lambda: Transpose(Hat(self.vector(d))),
                lambda: Transpose(MMMul(self.matrix(d), self.matrix(d))),
                lambda: Transpose(Transpose(self.matrix(d))),
                lambda: Transpose(L.J),
                lambda: Transpose(L.K),
                lambda: Transpose(SMMul(self.matrix(d), self.scalar(d))),
                lambda: Transpose(VVMul(self.plain_vector(d), Transpose(self.plain_vector(d)))),
                lambda: Hat(Vee(self.matrix(d))),  # generic M: not skew
                lambda: Hat(Vee(L.K)),
                lambda: Hat(Vee(Hat(self.vector(d)))),
                lambda: Hat(Vee(MMMul(Transpose(L.R), L.R.t_diff()))),
                lambda: Hat(VAdd(self.vector(d), self.vector(d))),
                lambda: MMMul(L.Ident, self.matrix(d)),
                lambda: MMMul(SMMul(L.R, self.scalar(d)), SMMul(Transpose(L.R), self.scalar(d))),
                lambda: Transpose(MAdd(self.matrix(d), self.matrix(d))),
                lambda: MMMul(Transpose(L.R), Transpose(Transpose(L.R))),
            ]
            return self.pick(motifs)()
        return MMMul(self.matrix(d), self.matrix(d))

    def any(self, depth):
        r = self.rng.random()
        if r < 0.4:
            return self.scalar(depth)
        if r < 0.8:
            return self.vector(depth)
        return self.matrix(depth)


# ---------------------------------------------------------------------------
# Tree surgery for shrinking
# ---------------------------------------------------------------------------


def rebuild(expr, new_nodes):
    cls = type(expr)
    if isinstance(expr, (Add, VAdd, MAdd)):
        return cls(*new_nodes)
    if isinstance(expr, (Mul, SVMul, SMMul, MVMul, MMMul, VVMul, Dot, Cross)):
        return cls(new_nodes[0], new_nodes[1])
    if isinstance(expr, (Hat, Vee, Transpose, Variation, TimeDerivative, TimeIntegral)):
        return cls(new_nodes[0])
    return expr


_FRESH = {
    ExprType.SCALAR: [Scalar("z_s")],
    ExprType.VECTOR: [Vector("z_v")],
    ExprType.MATRIX: [Matrix("Z_m")],
}


def shrink_candidates(expr):
    """Simpler trees of the same type: same-typed descendants, then nodes with one
    child replaced by a fresh leaf (keeping a Transpose wrapper when present)."""
    out = []
    nodes = getattr(expr, "nodes", None) or []
    etype = expr.type
    stack = list(nodes)
    seen = set()
    while stack:
        n = stack.pop(0)
        if getattr(n, "type", None) == etype and str(n) not in seen:
            seen.add(str(n))
            out.append(n)
        stack.extend(getattr(n, "nodes", None) or [])
    if isinstance(expr, (Add, VAdd, MAdd)) and len(nodes) > 2:
        for i in range(len(nodes)):
            out.append(rebuild(expr, nodes[:i] + nodes[i + 1 :]))
    for i, child in enumerate(nodes):
        ctype = getattr(child, "type", None)
        if ctype is None or ctype not in _FRESH:
            continue
        if isinstance(child, Transpose):
            inner_t = getattr(child.expr, "type", None)
            if inner_t not in _FRESH:
                continue
            repl = Transpose(_FRESH[inner_t][0])
        else:
            repl = _FRESH[ctype][0]
        if str(repl) == str(child):
            continue
        try:
            out.append(rebuild(expr, nodes[:i] + [repl] + nodes[i + 1 :]))
        except Exception:
            pass
    return out


def shrink(expr, fails, budget=400):
    """Greedy shrink: repeatedly move to the smallest failing candidate."""
    tries = 0
    while tries < budget:
        improved = False
        cands = sorted(shrink_candidates(expr), key=tree_size)
        for cand in cands:
            tries += 1
            try:
                ok = fails(cand)
            except Exception:
                ok = False
            if ok:
                expr = cand
                improved = True
                break
            if tries >= budget:
                break
        if not improved:
            break
    return expr


# ---------------------------------------------------------------------------
# Harness bookkeeping
# ---------------------------------------------------------------------------


class CaseTimeout(Exception):
    pass


def _alarm(signum, frame):
    raise CaseTimeout()


class Stats:
    def __init__(self, name):
        self.name = name
        self.cases = 0
        self.ok = 0
        self.unsupported = 0
        self.documented_errors = {}
        self.findings = []
        self.seen = set()

    def add_finding(self, kind, key, **info):
        if key in self.seen:
            return
        self.seen.add(key)
        self.findings.append(dict(kind=kind, transform=self.name, **info))


ALL = {}


def stats(name):
    if name not in ALL:
        ALL[name] = Stats(name)
    return ALL[name]


def classify_exception(exc):
    if isinstance(exc, DOCUMENTED_ERRORS):
        return "documented"
    if isinstance(exc, SEMI_DOCUMENTED_ERRORS):
        return "semi"
    return "crash"


def run_transform_case(name, transform, expr, rng, mode="consistent", shrinker=True, gauge="full"):
    """Apply *transform* to expr and compare values. Records findings."""
    st = stats(name)
    st.cases += 1
    if tdepth(expr) > KT - 1:
        st.unsupported += 1
        return
    env = {}
    mk = lambda: JetEvaluator(rng, mode=mode, env=env, gauge=gauge)
    ev = mk()
    try:
        before = ev(expr)
    except (Unsupported, ShapeMismatch):
        st.unsupported += 1
        return
    signal.setitimer(signal.ITIMER_REAL, CASE_TIMEOUT_S)
    try:
        out = transform(expr)
    except CaseTimeout:
        signal.setitimer(signal.ITIMER_REAL, 0)
        st.add_finding("timeout", f"timeout:{name}", expr=str(expr), repr=repr(expr))
        return
    except Exception as exc:
        signal.setitimer(signal.ITIMER_REAL, 0)
        c = classify_exception(exc)
        if c == "documented":
            st.documented_errors[type(exc).__name__] = (
                st.documented_errors.get(type(exc).__name__, 0) + 1
            )
        exc_type = type(exc)
        key = f"{c}:{exc_type.__name__}:{traceback.extract_tb(exc.__traceback__)[-1].lineno}"

        def fails_exc(e):
            try:
                transform(e)
            except exc_type:
                return True
            except Exception:
                return False
            return False

        small = shrink(expr, fails_exc) if shrinker else expr
        tb = traceback.extract_tb(exc.__traceback__)[-1]
        st.add_finding(
            {
                "crash": "crash",
                "semi": "semi-documented-exception",
                "documented": "documented-error-on-legal-input",
            }[c],
            key,
            expr=str(small),
            repr=repr(small),
            error=f"{type(exc).__name__}: {exc}",
            where=f"{os.path.relpath(tb.filename)}:{tb.lineno}",
        )
        return
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
    try:
        after = mk()(out)
    except Unsupported:
        st.unsupported += 1
        return
    except ShapeMismatch as exc:

        def fails_shape(e):
            try:
                mk()(e)
                mk()(transform(e))
            except ShapeMismatch:
                return True
            except Exception:
                return False
            return False

        small = shrink(expr, fails_shape) if shrinker else expr
        st.add_finding(
            "row-column-confusion",
            f"shape:{name}:{str(transform(small))[:120]}",
            expr=str(small),
            repr=repr(small),
            output=str(transform(small)),
            error=f"operand shapes {exc} in the output tree",
            mode=mode,
        )
        return
    ok, err = values_close(before, after)
    if ok:
        st.ok += 1
        return

    def fails_val(e):
        b = mk()(e)
        o = transform(e)
        a = mk()(o)
        return not values_close(b, a)[0]

    small = shrink(expr, fails_val) if shrinker else expr
    small_out = transform(small)
    ev3 = mk()
    b, a = ev3(small), ev3(small_out)
    st.add_finding(
        "wrong-value" if mode == "consistent" else "flag-dependent",
        f"val:{name}:{str(small_out)[:120]}",
        expr=str(small),
        repr=repr(small),
        output=str(small_out),
        before=np.round(b.val(), 6).tolist(),
        after=np.round(a.val(), 6).tolist(),
        rel_err=err,
        mode=mode,
    )


# ---------------------------------------------------------------------------
# Fuzz campaigns
# ---------------------------------------------------------------------------


def fuzz_value_preserving(n_cases, seed, mode="consistent"):
    rng = np.random.default_rng(seed)
    L = Leaves()
    transforms = [
        ("expand", expand),
        ("pull", pull),
        ("vector_rules", vector_rules),
        ("simplify", simplify),
        ("full_simplify", full_simplify),
    ]
    suffix = "" if mode == "consistent" else f"[{mode}]"
    gen = Gen(rng, L)
    for i in range(n_cases):
        depth = int(rng.integers(1, 5))
        expr = gen.any(depth)
        for name, fn in transforms:
            run_transform_case(name + suffix, fn, expr, rng, mode=mode)
        # full_simplify convergence / idempotence check
        if mode == "consistent":
            st = stats("full_simplify-convergence")
            st.cases += 1
            try:
                out = full_simplify(expr)
                again = simplify(out)
                if str(again) != str(out):
                    st.add_finding(
                        "cosmetic",
                        f"conv:{str(out)[:80]}",
                        expr=str(expr),
                        repr=repr(expr),
                        output=str(out),
                        note=(
                            "simplify(full_simplify(e)) != full_simplify(e):"
                            " not a fixed point after max_iter"
                        ),
                    )
                else:
                    st.ok += 1
            except Exception:
                pass


def _delta_of_node(e):
    # e is Variation(x): compare with x.delta(); other roots (shrink candidates) are left as is
    return e.expr.delta() if isinstance(e, Variation) else e


def _tdiff_of_node(e):
    return e.expr.t_diff() if isinstance(e, TimeDerivative) else e


def fuzz_delta_tdiff2(n_cases, seed):
    rng = np.random.default_rng(seed)
    L = Leaves()
    gen = Gen(rng, L, calc=True, variation=False)
    for i in range(n_cases):
        depth = int(rng.integers(1, 5))
        expr = gen.any(depth)
        if contains_variation(expr):
            continue
        run_transform_case("delta-rule", _delta_of_node, Variation(expr), rng)
        run_transform_case("t_diff-rule", _tdiff_of_node, TimeDerivative(expr), rng)


def fuzz_manifold_rules(n_cases, seed):
    rng = np.random.default_rng(seed)
    L = Leaves()
    gen = Gen(rng, L, calc=True, variation=True)
    variables = SystemVariables(vectors=[L.q, L.p], matrices=[L.R])

    def tr(e):
        return apply_manifold_rules(e, variables)

    for i in range(n_cases):
        depth = int(rng.integers(1, 4))
        expr = gen.any(depth)
        # make sure Variation(omega)/Variation(Omega) occur often
        r = rng.random()
        if r < 0.3:
            expr = (
                Dot(Variation(L.w), gen.vector(depth))
                if rng.random() < 0.5
                else Cross(Variation(L.Om), gen.vector(depth))
            )
        elif r < 0.45:
            expr = (
                TimeDerivative(Variation(L.w))
                if rng.random() < 0.5
                else VAdd(Variation(L.wp), gen.vector(depth))
            )
        # reference: the exact delta(omega) = xi_dot - 2 omega x xi (fix plan 1.1); the
        # "tangential" gauge is the pre-fix formula and is no longer a valid reference.
        run_transform_case("manifold_rules[full delta]", tr, expr, rng, gauge="full")


def fuzz_collect(n_cases, seed):
    rng = np.random.default_rng(seed)
    L = Leaves()
    gen = Gen(rng, L, calc=True, variation=True)
    targets = [
        L.xi,
        TimeDerivative(L.xi),
        L.eta,
        TimeDerivative(L.eta),
        L.a,
        Variation(L.a),
        TimeDerivative(Variation(L.a)),
    ]
    for i in range(n_cases):
        vec = targets[rng.integers(len(targets))]
        depth = int(rng.integers(1, 4))
        d = depth
        # scalar trees that contain the target somewhere inside Dot/Cross/MVMul
        builders = [
            lambda: Dot(gen.vector(d), vec),
            lambda: Dot(vec, gen.vector(d)),
            lambda: Dot(MVMul(gen.matrix(d), vec), Cross(gen.vector(d), gen.vector(d))),
            lambda: Dot(MVMul(gen.matrix(d), gen.vector(d)), Cross(vec, gen.vector(d))),
            lambda: Dot(MVMul(gen.matrix(d), gen.vector(d)), Cross(gen.vector(d), vec)),
            lambda: Dot(Cross(gen.vector(d), gen.vector(d)), MVMul(gen.matrix(d), vec)),
            lambda: Dot(MVMul(L.R, vec), MVMul(L.R, gen.vector(d))),
            lambda: Dot(MVMul(L.R, gen.vector(d)), MVMul(L.R, Cross(vec, gen.vector(d)))),
            lambda: Dot(MVMul(gen.matrix(d), vec), MVMul(gen.matrix(d), gen.vector(d))),
            lambda: Dot(MVMul(gen.matrix(d), Cross(vec, gen.vector(d))), gen.vector(d)),
            lambda: Dot(MVMul(gen.matrix(d), Cross(gen.vector(d), vec)), gen.vector(d)),
            lambda: Dot(gen.vector(d), MVMul(gen.matrix(d), Cross(gen.vector(d), vec))),
            lambda: Dot(gen.vector(d), Cross(vec, gen.vector(d))),
            lambda: Dot(gen.vector(d), Cross(gen.vector(d), vec)),
            lambda: Dot(gen.vector(d), Cross(Cross(vec, gen.vector(d)), gen.vector(d))),
            lambda: Dot(gen.vector(d), Cross(gen.vector(d), Cross(gen.vector(d), vec))),
            lambda: Dot(Cross(Cross(gen.vector(d), vec), gen.vector(d)), gen.vector(d)),
            lambda: Dot(MVMul(gen.matrix(d), MVMul(gen.matrix(d), vec)), gen.vector(d)),
            lambda: Dot(MVMul(gen.matrix(d), SVMul(vec, gen.scalar(d))), gen.vector(d)),
            lambda: Dot(
                gen.vector(d),
                MVMul(gen.matrix(d), Cross(gen.vector(d), Cross(vec, gen.vector(d)))),
            ),
            lambda: Mul(
                gen.scalar(d), Dot(MVMul(gen.matrix(d), vec), Cross(gen.vector(d), gen.vector(d)))
            ),
            lambda: Add(
                Dot(gen.vector(d), vec),
                Mul(
                    Dot(MVMul(gen.matrix(d), gen.vector(d)), Cross(gen.vector(d), vec)),
                    gen.scalar(d),
                ),
            ),
        ]

        # Target inside the *matrix* factor of an MVMul (phase-1 review F1): Hat(v),
        # v cᵀ, M Hat(v), on either side, with and without a matrix on the other side.
        # Transposition cannot expose it; collect must raise AlgebraicError (a documented
        # error), never recurse.
        def _col(v):  # VVMul(v, wᵀ) needs a column vector on the left
            while isinstance(v, Transpose):
                v = v.expr
            return v

        matrix_factor_builders = [
            lambda: Dot(MVMul(Hat(vec), gen.vector(d)), gen.vector(d)),
            lambda: Dot(gen.vector(d), MVMul(Hat(vec), gen.vector(d))),
            lambda: Dot(MVMul(VVMul(vec, Transpose(gen.vector(d))), gen.vector(d)), gen.vector(d)),
            lambda: Dot(
                gen.vector(d), MVMul(VVMul(_col(gen.vector(d)), Transpose(vec)), gen.vector(d))
            ),
            lambda: Dot(MVMul(MMMul(gen.matrix(d), Hat(vec)), gen.vector(d)), gen.vector(d)),
            lambda: Dot(MVMul(Hat(vec), gen.vector(d)), MVMul(gen.matrix(d), gen.vector(d))),
            lambda: Dot(MVMul(gen.matrix(d), gen.vector(d)), MVMul(Hat(vec), gen.vector(d))),
            lambda: Dot(MVMul(Hat(vec), gen.vector(d)), Cross(gen.vector(d), gen.vector(d))),
            lambda: Dot(Cross(gen.vector(d), gen.vector(d)), MVMul(Hat(vec), gen.vector(d))),
            lambda: Dot(MVMul(gen.matrix(d), MVMul(Hat(vec), gen.vector(d))), gen.vector(d)),
        ]
        if rng.random() < 0.25:
            expr = matrix_factor_builders[rng.integers(len(matrix_factor_builders))]()
            label = "collect[target inside matrix factor]"
        else:
            expr = builders[rng.integers(len(builders))]()
            linear = _count_occurrences(expr, vec) == 1
            label = (
                "collect[target occurs once]" if linear else "collect[target occurs several times]"
            )
        run_transform_case(label, lambda e, v=vec: collect(e, v), expr, rng)


def fuzz_t_diff_product(n_cases, seed):
    rng = np.random.default_rng(seed)
    L = Leaves()
    gen = Gen(rng, L, calc=True, variation=False)
    factors_pool = [L.s1, L.s2, L.m, L.kval, Number(2), One]
    for i in range(n_cases):
        n = int(rng.integers(1, 5))
        fs = []
        for _ in range(n):
            r = rng.random()
            if r < 0.6:
                fs.append(factors_pool[rng.integers(len(factors_pool))])
            elif r < 0.8:
                fs.append(Dot(gen.vector(1), gen.vector(1)))
            else:
                fs.append(gen.scalar(2))
        if contains_variation(Add(*fs)) if len(fs) > 1 else contains_variation(fs[0]):
            continue
        chain = fs[0]
        for f in fs[1:]:
            chain = Mul(chain, f) if rng.random() < 0.7 else Mul(f, chain)
        run_transform_case(
            "ibp._t_diff_product", lambda e: _t_diff_product(e.expr), TimeDerivative(chain), rng
        )


# -- extract_linear_coeff -------------------------------------------------------


class LinGen:
    """Trees that are linear (degree exactly 1) in a query leaf."""

    def __init__(self, rng, gen: Gen, query):
        self.rng = rng
        self.gen = gen
        self.query = query
        self.qtype = query.type

    def chance(self, p):
        return self.rng.random() < p

    def _free(self, make, fallback):
        """A generic subtree that does not contain the query (retry, then fallback)."""
        for _ in range(8):
            e = make()
            if not e.has(self.query):
                return e
        return fallback

    def gs(self, d):
        return self._free(lambda: self.gen.scalar(d), Scalar("z_s"))

    def gv(self, d):
        return self._free(lambda: self.gen.vector(d), Vector("z_v"))

    def gpv(self, d):
        return self._free(lambda: self.gen.plain_vector(d), Vector("z_v"))

    def gm(self, d):
        return self._free(lambda: self.gen.matrix(d), Matrix("Z_m"))

    def plain_lin_vector(self, d):
        for _ in range(6):
            v = self.lin_vector(d)
            if not isinstance(v, Transpose):
                return v
        return self.lin_vector(0)

    def lin_vector(self, d):
        """Vector expression linear in query."""
        q = self.query
        if d <= 0:
            if self.qtype == ExprType.VECTOR:
                return q
            if self.qtype == ExprType.SCALAR:
                return SVMul(self.gv(1), q)
            return MVMul(q, self.gv(1))
        r = self.rng.random()
        if r < 0.15:
            return VAdd(self.lin_vector(d - 1), self.gv(d - 1))
        if r < 0.25:
            return VAdd(self.gv(d - 1), self.lin_vector(d - 1), self.lin_vector(d - 1))
        if r < 0.40:
            return SVMul(self.lin_vector(d - 1), self.gs(d - 1))
        if r < 0.50:
            return SVMul(self.gv(d - 1), self.lin_scalar(d - 1))
        if r < 0.62:
            return MVMul(self.gm(d - 1), self.lin_vector(d - 1))
        if r < 0.70:
            return MVMul(self.lin_matrix(d - 1), self.gv(d - 1))
        if r < 0.80:
            return Cross(self.gv(d - 1), self.lin_vector(d - 1))
        if r < 0.90:
            return Cross(self.lin_vector(d - 1), self.gv(d - 1))
        if r < 0.95:
            return Vee(self.lin_matrix(d - 1))
        return Transpose(Transpose(self.lin_vector(d - 1)))

    def lin_scalar(self, d):
        q = self.query
        if d <= 0:
            if self.qtype == ExprType.SCALAR:
                return q
            if self.qtype == ExprType.VECTOR:
                return Dot(self.gv(1), q) if self.chance(0.5) else Dot(q, self.gv(1))
            return Dot(self.gv(1), MVMul(q, self.gv(1)))
        r = self.rng.random()
        if r < 0.2:
            return Add(self.lin_scalar(d - 1), self.gs(d - 1))
        if r < 0.3:
            return Add(self.gs(d - 1), self.lin_scalar(d - 1), self.lin_scalar(d - 1))
        if r < 0.45:
            return Mul(self.lin_scalar(d - 1), self.gs(d - 1))
        if r < 0.55:
            return Mul(self.gs(d - 1), self.lin_scalar(d - 1))
        if r < 0.75:
            return Dot(self.gv(d - 1), self.lin_vector(d - 1))
        if r < 0.90:
            return Dot(self.lin_vector(d - 1), self.gv(d - 1))
        if r < 0.95:
            return VVMul(Transpose(self.plain_lin_vector(d - 1)), self.gpv(d - 1))
        return VVMul(Transpose(self.gpv(d - 1)), self.plain_lin_vector(d - 1))

    def lin_matrix(self, d):
        q = self.query
        if d <= 0:
            if self.qtype == ExprType.MATRIX:
                return q
            if self.qtype == ExprType.SCALAR:
                return SMMul(self.gm(1), q)
            return Hat(q)
        r = self.rng.random()
        if r < 0.2:
            return MAdd(self.lin_matrix(d - 1), self.gm(d - 1))
        if r < 0.35:
            return SMMul(self.lin_matrix(d - 1), self.gs(d - 1))
        if r < 0.45:
            return SMMul(self.gm(d - 1), self.lin_scalar(d - 1))
        if r < 0.60:
            return MMMul(self.gm(d - 1), self.lin_matrix(d - 1))
        if r < 0.75:
            return MMMul(self.lin_matrix(d - 1), self.gm(d - 1))
        if r < 0.85:
            return Hat(self.lin_vector(d - 1))
        if r < 0.92:
            return Transpose(self.lin_matrix(d - 1))
        if r < 0.96:
            return VVMul(self.plain_lin_vector(d - 1), Transpose(self.gpv(d - 1)))
        return VVMul(self.gpv(d - 1), Transpose(self.plain_lin_vector(d - 1)))

    def nonlinear(self, d):
        """Trees of degree 2 in the query (should raise NonLinearError)."""
        r = self.rng.random()
        if self.qtype == ExprType.VECTOR:
            if r < 0.25:
                return Dot(self.lin_vector(d), self.lin_vector(d))
            if r < 0.5:
                return SVMul(self.lin_vector(d), self.lin_scalar(d))
            if r < 0.7:
                return MVMul(self.lin_matrix(d), self.lin_vector(d))
            if r < 0.85:
                return Cross(self.lin_vector(d), self.lin_vector(d))
            return Mul(self.lin_scalar(d), self.lin_scalar(d))
        if self.qtype == ExprType.SCALAR:
            if r < 0.4:
                return Mul(self.lin_scalar(d), self.lin_scalar(d))
            if r < 0.7:
                return SVMul(self.lin_vector(d), self.lin_scalar(d))
            return Dot(self.lin_vector(d), self.lin_vector(d))
        return MMMul(self.lin_matrix(d), self.lin_matrix(d))


def _coeff_times(coeff, qval, qtype, etype):
    """Numeric coeff . query according to the documented convention."""
    c = coeff
    if qtype == ExprType.SCALAR:
        return jsmul(qval, c)  # c * s for any etype
    if qtype == ExprType.VECTOR:
        if etype == ExprType.SCALAR:
            return jdot(c, qval)  # expr = query^T coeff
        if etype == ExprType.VECTOR:
            return jmatmul(c, qval.col())  # expr = C q
        return None
    if qtype == ExprType.MATRIX:
        return None
    return None


def fuzz_extract(n_cases, seed):
    rng = np.random.default_rng(seed)
    L = Leaves()
    gen = Gen(rng, L, calc=True, variation=True, motifs=False, rowvec=False)
    queries = [L.xi, L.eta, L.a, L.s1, L.M, Variation(L.b), Variation(L.s2)]
    for i in range(n_cases):
        query = queries[rng.integers(len(queries))]
        qtype = query.type
        lg = LinGen(rng, gen, query)
        depth = int(rng.integers(0, 4))
        r = rng.random()
        if r < 0.45:
            expr = lg.lin_scalar(depth)
        elif r < 0.85:
            expr = lg.lin_vector(depth)
        else:
            expr = lg.lin_matrix(depth)
        etype = expr.type
        name = f"extract[{qtype.name.lower()} query, {etype.name.lower()} expr]"
        st = stats(name)
        st.cases += 1
        qkey = str(query)
        # the query must not also occur through calculus wrapping elsewhere -- it does not by
        # construction (gen trees never contain the query leaf except via LinGen)
        if _contains_str(expr, qkey) != "linear-ok" or tdepth(expr) > KT - 1:
            st.unsupported += 1
            continue
        shape = {ExprType.SCALAR: (1, 1), ExprType.VECTOR: (3, 1), ExprType.MATRIX: (3, 3)}[qtype]
        v = rand_jet(rng, shape)
        zero = jconst(np.zeros(shape))
        env = {}
        try:
            ev_v = JetEvaluator(rng, env=env, overrides={qkey: v})
            e_v = ev_v(expr)
            e_0 = JetEvaluator(rng, env=env, overrides={qkey: zero})(expr)
        except Unsupported:
            st.unsupported += 1
            continue
        signal.setitimer(signal.ITIMER_REAL, CASE_TIMEOUT_S)
        try:
            coeff = extract_linear_coeff(expr, query)
        except CaseTimeout:
            st.add_finding("timeout", f"timeout:{name}", expr=str(expr), repr=repr(expr))
            continue
        except Exception as exc:
            signal.setitimer(signal.ITIMER_REAL, 0)
            c = classify_exception(exc)
            if c == "documented":
                st.documented_errors[type(exc).__name__] = (
                    st.documented_errors.get(type(exc).__name__, 0) + 1
                )
                # a linear-by-construction tree must not raise NonLinearError
                if isinstance(exc, NonLinearError):
                    small = shrink(
                        expr,
                        lambda e: _raises(lambda: extract_linear_coeff(e, query), NonLinearError),
                    )
                    st.add_finding(
                        "false-nonlinear",
                        f"falsenl:{name}:{type(exc).__name__}",
                        expr=str(small),
                        repr=repr(small),
                        query=qkey,
                        error=f"{type(exc).__name__}: {exc}",
                    )
                elif isinstance(exc, AlgebraicError):
                    small = shrink(
                        expr,
                        lambda e: _raises(lambda: extract_linear_coeff(e, query), AlgebraicError),
                    )
                    st.add_finding(
                        "unsupported-linear-term",
                        f"alg:{name}:{str(exc)[:60]}",
                        expr=str(small),
                        repr=repr(small),
                        query=qkey,
                        error=f"{type(exc).__name__}: {exc}",
                    )
                continue
            tb = traceback.extract_tb(exc.__traceback__)[-1]
            exc_type = type(exc)
            small = shrink(
                expr, lambda e: _raises(lambda: extract_linear_coeff(e, query), exc_type)
            )
            st.add_finding(
                "crash" if c == "crash" else "semi-documented-exception",
                f"{c}:{name}:{exc_type.__name__}:{tb.lineno}",
                expr=str(small),
                repr=repr(small),
                query=qkey,
                error=f"{type(exc).__name__}: {exc}",
                where=f"{os.path.relpath(tb.filename)}:{tb.lineno}",
            )
            continue
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
        if coeff.has(query):
            st.add_finding(
                "coefficient-contains-query",
                f"hasq:{name}:{str(coeff)[:60]}",
                expr=str(expr),
                repr=repr(expr),
                query=qkey,
                output=str(coeff),
            )
            continue
        if qtype == ExprType.MATRIX or (qtype == ExprType.VECTOR and etype == ExprType.MATRIX):
            # no numeric convention for these cells: check only that a term depending on the
            # query is never reported as a zero coefficient
            depends = not values_close(e_v, e_0)[0]
            if coeff.is_zero and depends:

                def fails_zero(e):
                    try:
                        cc = extract_linear_coeff(e, query)
                    except Exception:
                        return False
                    d_ = JetEvaluator(rng, env=env, overrides={qkey: v})(e) - JetEvaluator(
                        rng, env=env, overrides={qkey: zero}
                    )(e)
                    return cc.is_zero and not values_close(d_, Jet(np.zeros_like(d_.c)))[0]

                small = shrink(expr, fails_zero)
                st.add_finding(
                    "silent-zero",
                    f"zero:{name}:{str(small)[:80]}",
                    expr=str(small),
                    repr=repr(small),
                    query=qkey,
                    output=str(extract_linear_coeff(small, query)),
                    note=(
                        "expression depends on the query but the returned coefficient"
                        " is a zero constant"
                    ),
                )
            else:
                st.ok += 1
            continue
        try:
            cval = JetEvaluator(rng, env=env, overrides={qkey: zero})(coeff)
            pred = _coeff_times(cval, v, qtype, etype)
        except Unsupported:
            st.unsupported += 1
            continue
        diff = e_v - e_0
        ok, err = values_close(diff, pred)
        if ok:
            st.ok += 1
            continue

        def fails_val(e):
            try:
                cc = extract_linear_coeff(e, query)
                if cc.has(query):
                    return False
                d_ = JetEvaluator(rng, env=env, overrides={qkey: v})(e) - JetEvaluator(
                    rng, env=env, overrides={qkey: zero}
                )(e)
                p_ = _coeff_times(
                    JetEvaluator(rng, env=env, overrides={qkey: zero})(cc), v, qtype, e.type
                )
                return p_ is not None and not values_close(d_, p_)[0]
            except Exception:
                return False

        small = shrink(expr, fails_val)
        small_coeff = extract_linear_coeff(small, query)
        st.add_finding(
            "wrong-value",
            f"val:{name}:{str(small_coeff)[:80]}",
            expr=str(small),
            repr=repr(small),
            query=qkey,
            output=str(small_coeff),
            rel_err=err,
        )


def fuzz_extract_nonlinear(n_cases, seed):
    rng = np.random.default_rng(seed)
    L = Leaves()
    gen = Gen(rng, L, calc=True, variation=True, motifs=False, rowvec=False)
    queries = [L.xi, L.eta, L.a, L.s1, L.M]
    for i in range(n_cases):
        query = queries[rng.integers(len(queries))]
        lg = LinGen(rng, gen, query)
        depth = int(rng.integers(0, 3))
        expr = lg.nonlinear(depth)
        qtype, etype = query.type, expr.type
        name = f"extract-nonlinear[{qtype.name.lower()} query, {etype.name.lower()} expr]"
        st = stats(name)
        st.cases += 1
        qkey = str(query)
        signal.setitimer(signal.ITIMER_REAL, CASE_TIMEOUT_S)
        try:
            coeff = extract_linear_coeff(expr, query)
        except CaseTimeout:
            st.add_finding("timeout", f"timeout:{name}", expr=str(expr), repr=repr(expr))
            continue
        except NonLinearError:
            st.ok += 1
            continue
        except Exception as exc:
            c = classify_exception(exc)
            if c == "documented":
                st.documented_errors[type(exc).__name__] = (
                    st.documented_errors.get(type(exc).__name__, 0) + 1
                )
                continue
            tb = traceback.extract_tb(exc.__traceback__)[-1]
            st.add_finding(
                "crash" if c == "crash" else "semi-documented-exception",
                f"{c}:{name}:{type(exc).__name__}:{tb.lineno}",
                expr=str(expr),
                repr=repr(expr),
                query=qkey,
                error=f"{type(exc).__name__}: {exc}",
                where=f"{os.path.relpath(tb.filename)}:{tb.lineno}",
            )
            continue
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
        # returned something for a quadratic tree: is it at least numerically right?
        shape = {ExprType.SCALAR: (1, 1), ExprType.VECTOR: (3, 1), ExprType.MATRIX: (3, 3)}[qtype]
        env = {}
        try:
            v = rand_jet(rng, shape)
            zero = jconst(np.zeros(shape))
            e_v = JetEvaluator(rng, env=env, overrides={qkey: v})(expr)
            e_0 = JetEvaluator(rng, env=env, overrides={qkey: zero})(expr)
            e_2v = JetEvaluator(rng, env=env, overrides={qkey: Jet(2 * v.c)})(expr)
            quad = not values_close(e_2v - e_0, Jet(2 * (e_v - e_0).c))[0]
        except Unsupported:
            st.unsupported += 1
            continue
        if not quad:
            # degree-2 tree that happens to be identically zero or linear (e.g. Cross(q,q))
            st.ok += 1
            continue

        def fails_nl(e):
            try:
                extract_linear_coeff(e, query)
                a0 = JetEvaluator(rng, env=env, overrides={qkey: zero})(e)
                a1 = JetEvaluator(rng, env=env, overrides={qkey: v})(e)
                a2 = JetEvaluator(rng, env=env, overrides={qkey: Jet(2 * v.c)})(e)
                return not values_close(a2 - a0, Jet(2 * (a1 - a0).c))[0]
            except Exception:
                return False

        small = shrink(expr, fails_nl)
        try:
            small_coeff = str(extract_linear_coeff(small, query))
        except Exception:
            small_coeff = str(coeff)
        st.add_finding(
            "missed-nonlinear",
            f"nl:{name}:{small_coeff[:80]}",
            expr=str(small),
            repr=repr(small),
            query=qkey,
            output=small_coeff,
            note="query occurs with degree 2 but no NonLinearError was raised",
        )


def _count_occurrences(expr, target):
    if str(expr) == str(target) and type(expr) is type(target):
        return 1
    return sum(_count_occurrences(n, target) for n in getattr(expr, "nodes", None) or [])


def _raises(fn, exc_type):
    try:
        fn()
    except exc_type:
        return True
    except Exception:
        return False
    return False


def _contains_str(expr, key):
    """'linear-ok' if the query leaf appears only as itself (not wrapped in calculus)."""

    def walk(e, inside_calc):
        if str(e) == key:
            return "bad" if inside_calc else "ok"
        res = "none"
        for n in getattr(e, "nodes", None) or []:
            r = walk(n, inside_calc or isinstance(e, (TimeDerivative, Variation, TimeIntegral)))
            if r == "bad":
                return "bad"
            if r == "ok":
                res = "ok"
        return res

    return "linear-ok" if walk(expr, False) == "ok" else "bad"


# ---------------------------------------------------------------------------
# Self-test of the jet evaluator against tests/numeric_helpers.Evaluator
# ---------------------------------------------------------------------------


def selftest(seed=1):
    rng = np.random.default_rng(seed)
    L = Leaves()
    ev = JetEvaluator(rng)
    # manifold constraints hold to all orders
    q = ev(L.q)
    w, xi = ev(L.w), ev(L.xi)
    assert low_orders_close(jdot(q, q), jconst([[1.0]]))
    assert low_orders_close(jdot(q, w), jconst([[0.0]]))
    assert low_orders_close(jdot(q, xi), jconst([[0.0]]))
    assert low_orders_close(q.dt(), jcross(w, q))
    assert low_orders_close(q.de(), jcross(xi, q))
    R = ev(L.R)
    assert low_orders_close(jmatmul(R.T(), R), jconst(np.eye(3)))
    Om, eta = ev(L.Om), ev(L.eta)
    assert low_orders_close(R.dt(), jmatmul(R, jhat(Om)))
    assert low_orders_close(R.de(), jmatmul(R, jhat(eta)))
    # delta(omega) = xi_dot - omega x xi (manifold_rules formula) and
    # delta(Omega) = eta_dot + Omega x eta
    full_ok = values_close(w.de(), xi.dt() - jcross(w, xi))[0]
    two_ok = values_close(w.de(), xi.dt() - Jet(2 * jcross(w, xi).c))[0]
    tang = w.de() - jsmul(jdot(q, w.de()), q)
    tang_ok = values_close(tang, xi.dt() - jcross(w, xi))[0]
    print(
        f"self-test S2: delta(omega) == xi_dot - omega x xi: {full_ok};"
        f" == xi_dot - 2 omega x xi: {two_ok}; tangential part == xi_dot - omega x xi: {tang_ok}"
    )
    assert values_close(Om.de(), eta.dt() + jcross(Om, eta))[0]
    # a unit-norm flagged vector, a symmetric and an orthogonal matrix
    u, J, Q = ev(L.u), ev(L.J), ev(L.Q)
    assert low_orders_close(jdot(u, u), jconst([[1.0]]))
    assert low_orders_close(J, J.T())
    assert low_orders_close(jmatmul(Q.T(), Q), jconst(np.eye(3)))
    # agreement with the repository Evaluator on calculus-free trees
    sys.path.insert(0, os.path.join(REPO, "tests"))
    from numeric_helpers import Evaluator as RefEvaluator

    gen = Gen(rng, L, calc=False, motifs=True)
    n = 0
    for _ in range(300):
        e = gen.any(int(rng.integers(1, 4)))
        env = {}
        jev = JetEvaluator(rng, env=env)
        try:
            jv = jev(e)
        except Unsupported:
            continue
        ref_env = {
            k: v.val() if v.shape != (1, 1) else float(v.val()[0, 0]) for k, v in env.items()
        }

        # the reference evaluator keys d/dt(leaf) and delta(leaf) chains by their printed name
        def export_derived(node):
            if isinstance(node, (TimeDerivative, Variation)):
                val = jev(node)
                ref_env[str(node)] = val.val() if val.shape != (1, 1) else float(val.val()[0, 0])
            for n in getattr(node, "nodes", None) or []:
                export_derived(n)

        export_derived(e)
        # the reference evaluator keys S2/SO3 leaves by name and draws its own values for
        # anything missing; everything in env is already there.
        ref = RefEvaluator(env=ref_env, rng=np.random.default_rng(0))
        try:
            rv = ref(e)
        except Exception:
            continue
        rv = np.atleast_2d(np.asarray(rv, dtype=float))
        ok, err = values_close(jv, jconst(rv))
        assert ok, (str(e), jv.val(), rv, err)
        n += 1
    return n


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

REPO = os.path.dirname(os.path.dirname(os.path.abspath(_geomech.__file__)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", type=int, default=3000, help="cases per campaign")
    ap.add_argument("--seed", type=int, default=20261113)
    ap.add_argument("--quick", action="store_true", help="200 cases per campaign")
    ap.add_argument("--no-json", action="store_true")
    ap.add_argument("--json", default=None, help="findings JSON path")
    args = ap.parse_args()
    n = 200 if args.quick else args.cases
    signal.signal(signal.SIGALRM, _alarm)
    sys.setrecursionlimit(4000)

    t0 = time.time()
    n_self = selftest(args.seed)
    print(
        f"self-test: jet evaluator agrees with tests/numeric_helpers.Evaluator on {n_self} trees"
    )

    campaigns = [
        ("value-preserving (consistent)", lambda: fuzz_value_preserving(n, args.seed)),
        (
            "value-preserving (adversarial)",
            lambda: fuzz_value_preserving(max(n // 3, 50), args.seed + 1, mode="adversarial"),
        ),
        ("delta/t_diff rules", lambda: fuzz_delta_tdiff2(n, args.seed + 2)),
        ("manifold_rules", lambda: fuzz_manifold_rules(n, args.seed + 3)),
        ("collect", lambda: fuzz_collect(n, args.seed + 4)),
        ("ibp._t_diff_product", lambda: fuzz_t_diff_product(n, args.seed + 5)),
        ("extract_linear_coeff", lambda: fuzz_extract(n, args.seed + 6)),
        ("extract nonlinear", lambda: fuzz_extract_nonlinear(n, args.seed + 7)),
    ]
    for label, fn in campaigns:
        t1 = time.time()
        fn()
        print(f"campaign {label}: {time.time() - t1:.1f}s")

    print(f"\ntotal time {time.time() - t0:.1f}s\n")
    print(
        f"{'transform':48s} {'cases':>7s} {'ok':>7s} {'unsup':>6s} {'findings':>8s}"
        "  documented errors"
    )
    for name, st in ALL.items():
        print(
            f"{name:48s} {st.cases:7d} {st.ok:7d} {st.unsupported:6d} {len(st.findings):8d}"
            f"  {st.documented_errors}"
        )

    findings = [f for st in ALL.values() for f in st.findings]
    print(f"\n{len(findings)} distinct findings")
    for i, f in enumerate(findings, 1):
        print(f"\n--- finding {i}: [{f['kind']}] in {f['transform']}")
        for k, v in f.items():
            if k in ("kind", "transform"):
                continue
            print(f"  {k}: {v}")
    if not args.no_json:
        path = args.json or os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "engine_fuzz_findings.json"
        )
        with open(path, "w") as fh:
            json.dump(
                {
                    "stats": {
                        k: dict(
                            cases=v.cases,
                            ok=v.ok,
                            unsupported=v.unsupported,
                            documented_errors=v.documented_errors,
                        )
                        for k, v in ALL.items()
                    },
                    "findings": findings,
                },
                fh,
                indent=1,
                default=str,
            )
        print(f"\nfindings written to {path}")


if __name__ == "__main__":
    main()
