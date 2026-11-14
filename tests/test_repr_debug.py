"""Debug printers show what the engine keys on (fix plan 3.6; debuggability audit F6).

``repr`` prefixes non-scalar leaves with their type and marks constants with ``*``, and
joins products with ``·`` so that ``*`` is never a multiplication; ``print_tree`` defaults
to the indented style, whose leaves list every flag and the manifold point of a tangent or
variation vector. None of this touches ``str``, the equality key.
"""

from __future__ import annotations

import re

import pytest

from geomech import (
    S2,
    SO3,
    Cross,
    Dot,
    Hat,
    IdentityMatrix,
    Matrix,
    Number,
    Scalar,
    SkewSymmMatrix,
    Transpose,
    Vector,
    ZeroVector,
)
from geomech.utils.printing import print_tree, tree_str


class TestReprLeaves:
    @pytest.mark.parametrize(
        "leaf, text",
        [
            (Scalar("m"), "m"),
            (Scalar("m", attr=["Constant"]), "m*"),
            (Vector("x"), "v:x"),
            (Vector("e3", attr=["Constant"]), "v:e3*"),
            (Matrix("J"), "M:J"),
            (Matrix("J", attr=["Constant", "SymmetricMatrix"]), "M:J*"),
            (SkewSymmMatrix("S"), "M:S"),
            (S2("q"), "S2:q"),
            (SO3("R"), "SO3:R"),
            (S2("q").get_tangent_vector(), "v:\\omega_{q}"),
            (SO3("R").get_variation_vector(), "v:\\eta_{R}"),
            (Number(2), "(2)"),
            (Scalar("0.5", value=0.5, attr=["Constant"]), "0.5"),
            (ZeroVector, "v:0v*"),
            (IdentityMatrix, "M:I*"),
        ],
    )
    def test_leaf(self, leaf, text):
        assert repr(leaf) == text

    def test_same_name_different_type_differ(self):
        """``repr(Scalar('m')) == repr(Vector('m'))`` before phase 3."""
        assert repr(Scalar("m")) != repr(Vector("m"))
        assert repr(Vector("R")) != repr(SO3("R"))

    def test_constness_is_visible(self):
        assert repr(Scalar("m")) != repr(Scalar("m", attr=["Constant"]))


class TestReprComposites:
    def test_lagrangian(self):
        m, g = Scalar("m", attr=["Constant"]), Scalar("g", attr=["Constant"])
        e3, x = Vector("e3", attr=["Constant"]), Vector("x")
        half = Scalar("0.5", value=0.5, attr=["Constant"])
        v = x.t_diff()
        L = half * m * Dot(v, v) - m * g * Dot(e3, x)
        xd = "x\u0307"  # x + combining dot above
        assert repr(L) == f"(0.5·m*·<v:{xd}, v:{xd}> -m*·g*·<v:e3*, v:x>)"

    def test_no_multiplication_star(self):
        m = Scalar("m", attr=["Constant"])
        x = Vector("x")
        text = repr(m * Dot(x, x) * 2)
        # every * is a constant marker directly after a name
        assert all(re.match(r"\w", text[i - 1]) for i, ch in enumerate(text) if ch == "*")

    def test_rigid_body(self):
        R = SO3("R")
        Om = R.get_tangent_vector()
        J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
        assert repr(Dot(Om, J * Om)) == "<v:\\Omega_{R}, M:J*·v:\\Omega_{R}>"
        assert repr(Transpose(R) * Vector("x")) == "SO3:R'·v:x"
        assert repr(Transpose(R * Vector("x"))) == "(SO3:R·v:x)'"

    def test_sphere(self):
        q = S2("q")
        assert repr(q.t_diff()) == "cross(v:\\omega_{q}, S2:q)"

    def test_negation_and_hat(self):
        x = Vector("x")
        assert repr(-x) == "-v:x"
        assert repr(Hat(x)) == "hat(v:x)"
        assert repr(Scalar("y").t_diff()) == "y\u0307"

    def test_str_is_unchanged(self):
        """``str`` is the equality key and is not touched by the repr changes."""
        R = SO3("R")
        J = Matrix("J", attr=["Constant"])
        Om = R.get_tangent_vector()
        assert str(Dot(Om, J * Om)) == "Dot(\\Omega_{R},J*\\Omega_{R})"
        assert str(Transpose(R) * Vector("x")) == "(R)'*x"


class TestPrintTree:
    def test_default_style_is_indent(self, capsys):
        x = Vector("x")
        print_tree(Dot(x, x))
        out = capsys.readouterr().out
        assert out.splitlines()[0] == "Dot"
        assert "├─ L: Vector('x')" in out
        assert out.rstrip().endswith("expr: <v:x, v:x>")

    def test_topdown_still_available(self, capsys):
        x = Vector("x")
        print_tree(Dot(x, x), style="topdown")
        out = capsys.readouterr().out
        assert "v:x" in out and "/" in out

    def test_tree_str_default_is_unchanged(self):
        assert "v:x" in tree_str(Vector("x"))

    def test_leaf_flags(self):
        R = SO3("R")
        text = tree_str(R, style="indent")
        assert text == "SO3('R', is_orthogonal, is_manifold)"
        text = tree_str(SkewSymmMatrix("S"), style="indent")
        assert text == "SkewSymmMatrix('S', is_skew_symmetric)"
        assert "is_identity" in tree_str(IdentityMatrix, style="indent")
        assert "is_unit_norm" in tree_str(S2("q"), style="indent")

    def test_tangent_vector_shows_its_parent(self):
        q = S2("q")
        om = q.get_tangent_vector()
        text = tree_str(om, style="indent")
        assert text == "TS2('\\omega_{q}', of=S2 'q')"  # deterministic: no address (review F6)
        R = SO3("R")
        eta = R.get_variation_vector()
        assert "of=SO3 'R'" in tree_str(eta, style="indent")
        assert "@0x" not in tree_str(eta, style="indent")

    def test_two_parents_with_the_same_name_are_told_apart(self):
        """Debuggability audit C6: two ``S2("q")`` compare equal but are different objects."""
        q1, q2 = S2("q"), S2("q")
        t1 = tree_str(Cross(q1.get_tangent_vector(), q1), style="indent")
        t2 = tree_str(Cross(q2.get_tangent_vector(), q1), style="indent")
        assert q1.get_tangent_vector() == q2.get_tangent_vector()
        assert t1 != t2
