"""Precedence-aware ``to_latex`` (fix plan 3.5; dynamics audit F8, debuggability audit F6).

Before phase 3, ``to_latex`` never parenthesised a composite operand, so different trees
printed the same text (``a \\times b \\times c`` for both ``a×(b×c)`` and ``(a×b)×c``),
and ``print_eom``/``eom_to_latex``/``render_eom`` emitted ``str(eqn)`` inside LaTeX math.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from geomech import (
    SO3,
    Cross,
    Dot,
    Matrix,
    Mul,
    MVMul,
    Scalar,
    SVMul,
    Transpose,
    Vector,
    compute_eom,
)
from geomech.utils.printing import _strip_outer_parens, eom_to_latex, print_eom, to_latex

sys.path.insert(0, str(Path(__file__).parent))
import vbl_systems  # noqa: E402

a, b, c, v = Vector("a"), Vector("b"), Vector("c"), Vector("v")
R = SO3("R")
J = Matrix("J")
s, m = Scalar("s"), Scalar("m")


# Every row of the dynamics audit's F8 table, as (tree, tree it used to be confused with).
AMBIGUOUS_PAIRS = [
    pytest.param(Cross(a, Cross(b, c)), Cross(Cross(a, b), c), id="cross-associativity"),
    pytest.param(
        MVMul(Transpose(R), Cross(a, b)),
        Cross(MVMul(Transpose(R), a), b),
        id="Rt(a x b) vs (Rt a) x b",
    ),
    pytest.param(MVMul(R, Cross(a, b)), Cross(MVMul(R, a), b), id="R(a x b) vs (R a) x b"),
    pytest.param(Transpose(Cross(a, b)), Cross(a, Transpose(b)), id="(a x b)^T vs a x b^T"),
    pytest.param(Transpose(MVMul(R, a)), MVMul(R, Transpose(Transpose(a))), id="(R a)^T vs R a"),
    pytest.param(SVMul(a, Dot(b, c)), Dot(b, SVMul(c, Scalar("a_s"))), id="(b.c) a vs b.(c a)"),
    pytest.param(Cross(a, SVMul(b, -1)), SVMul(Cross(a, b), -1), id="a x (-b) vs -(a x b)"),
    pytest.param(
        Mul(Dot(a, b), Dot(b, c)), Dot(a, SVMul(b, Dot(b, c))), id="(a.b)(b.c) vs a.((b.c) b)"
    ),
]


@pytest.mark.parametrize("left, right", AMBIGUOUS_PAIRS)
def test_ambiguous_pairs_print_differently(left, right):
    assert str(left) != str(right)
    assert to_latex(left) != to_latex(right), (to_latex(left), to_latex(right))


@pytest.mark.parametrize(
    "expr, text",
    [
        (Cross(a, Cross(b, c)), r"a \times \left(b \times c\right)"),
        (Cross(Cross(a, b), c), r"\left(a \times b\right) \times c"),
        (MVMul(Transpose(R), Cross(a, b)), r"R^T \left(a \times b\right)"),
        (Cross(MVMul(Transpose(R), a), b), r"\left(R^T a\right) \times b"),
        (Cross(MVMul(R, a), b), r"\left(R a\right) \times b"),
        (Transpose(Cross(a, b)), r"\left(a \times b\right)^T"),
        (Transpose(MVMul(R, a)), r"\left(R a\right)^T"),
        (SVMul(a, Dot(b, c)), r"\left(b \cdot c\right) a"),
        (Cross(a, SVMul(b, -1)), r"a \times \left(-b\right)"),
        (SVMul(Cross(a, b), -1), r"-a \times b"),
        (SVMul(Cross(a, b), m), r"m \left(a \times b\right)"),
        (Mul(Mul(s, s), -0.5), r"-\frac{1}{2} s s"),
        (Mul(s, 2.5), r"2.5 s"),
        (Mul(Mul(s, 2), Mul(m, -3)), r"-6 s m"),
        (Mul(Dot(a, b), -1), r"-a \cdot b"),
        (Mul(Dot(a, b), Dot(b, c)), r"\left(a \cdot b\right) \left(b \cdot c\right)"),
        (Transpose(v) * s, r"s v^T"),
        (Transpose(R) * v, r"R^T v"),
        (MVMul(J, MVMul(R, a)), r"J R a"),  # associative: no parentheses
        (MVMul(J, SVMul(a, m)), r"J \left(m a\right)"),
        (Transpose(Transpose(R) * J), r"\left(R^T J\right)^T"),
        (a.t_diff(), r"\dot{a}"),
        (a.t_diff().t_diff(), r"\ddot{a}"),
    ],
)
def test_exact_text(expr, text):
    assert to_latex(expr) == text


def test_composite_time_derivative_and_variation_are_wrapped():
    from geomech import TimeDerivative, Variation

    assert to_latex(TimeDerivative(MVMul(J, a))) == r"\frac{d}{dt}\left(J a\right)"
    assert to_latex(Variation(MVMul(J, a))) == r"\delta\left(J a\right)"
    assert to_latex(Variation(a)) == r"\delta a"
    assert to_latex(Transpose(Variation(a))) == r"\left(\delta a\right)^T"


def test_transpose_never_doubles_a_superscript():
    from geomech import Vee

    assert "^T^T" not in to_latex(Transpose(Transpose(v) * s))
    assert to_latex(Transpose(Vee(J))) == r"\left(\left(J\right)^\vee\right)^T"
    assert to_latex(Transpose(Vector("x^*"))) == r"\left(x^*\right)^T"


def test_strip_outer_parens_only_strips_a_matching_pair():
    assert _strip_outer_parens(r"\left(a + b\right)") == "a + b"
    two_groups = r"\left(a\right) \left(b\right)"
    assert _strip_outer_parens(two_groups) == two_groups


# ---------------------------------------------------------------------------
# print_eom / eom_to_latex / render_eom use to_latex
# ---------------------------------------------------------------------------


def _double_rigid_pendulum_eom():
    sysm = vbl_systems.double_rigid_pendulum()
    return compute_eom(sysm.L, sysm.dW, sysm.variables)


def test_print_eom_uses_to_latex(capsys):
    eom = _double_rigid_pendulum_eom()
    print_eom(eom)
    out = capsys.readouterr().out.splitlines()
    assert len(out) == 2
    for line, (var_vec, eqn) in zip(out, eom.values()):
        assert line.startswith(r"\int " + to_latex(var_vec) + r" \cdot \left(")
        assert line.endswith(r"\right) \, dt = 0")
        assert _strip_outer_parens(to_latex(eqn)) in line
        assert "Cross(" not in line and "Dot(" not in line  # no str() form


def test_eom_to_latex_uses_to_latex():
    eom = _double_rigid_pendulum_eom()
    doc = eom_to_latex(eom)
    assert "Cross(" not in doc and "Hat(" not in doc and "^T^T" not in doc
    for var_vec, eqn in eom.values():
        assert _strip_outer_parens(to_latex(eqn)) in doc


def test_double_rigid_pendulum_row_is_unambiguous():
    """The audit's real example: ``-m1 Ω×ρ×Ω×ρ`` now shows its nesting."""
    eom = _double_rigid_pendulum_eom()
    row = to_latex(eom[r"\eta_{R1}"][1])
    assert (
        r"-m1 \left(\Omega_{R1} \times \left(\rho_1 \times \left(\Omega_{R1} \times"
        r" \rho_1\right)\right)\right)" in row
    )
    assert r"R1^T R2 \left(\Omega_{R2} \times \rho_2\right)" in row


@pytest.mark.skipif(shutil.which("pdflatex") is None, reason="pdflatex not installed")
@pytest.mark.parametrize("build", ["double_rigid_pendulum", "quadrotor_with_load"])
def test_eom_document_compiles(tmp_path, build):
    sysm = getattr(vbl_systems, build)()
    eom = compute_eom(sysm.L, sysm.dW, sysm.variables)
    tex = tmp_path / "eom.tex"
    tex.write_text(eom_to_latex(eom))
    result = subprocess.run(
        ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", str(tex)],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=120,
    )
    errors = [line for line in result.stdout.splitlines() if line.startswith("!")]
    assert result.returncode == 0, errors
    assert (tmp_path / "eom.pdf").exists()
