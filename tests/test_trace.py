"""Phase 4.1: ``EOMTrace`` / ``LinearizationTrace`` keep every pipeline stage.

The traces must (a) not change any result, (b) record the stages in the order the code
runs them, (c) let a debugger reproduce each stage from the previous one with the public
functions, and (d) keep the stages before a failure.
"""

from __future__ import annotations

import pytest

import geomech
from geomech import (
    SO3,
    AlgebraicError,
    Dot,
    EOMTrace,
    Hat,
    LinearizationTrace,
    Matrix,
    MMMul,
    SystemVariables,
    Vector,
    Vee,
    apply_manifold_rules,
    collect,
    compute_eom,
    expand,
    extract_linear_coeff,
    full_simplify,
    integrate_by_parts,
    linearize,
)
from geomech.core.math.ibp import _apply_ibp
from tests.vbl_systems import (
    point_mass_r3,
    quadrotor_with_load,
    rigid_pendulum_pivot,
    spherical_pendulum_s2,
)


def _strs(eom):
    return {k: str(eq) for k, (_, eq) in eom.items()}


def test_exports():
    assert geomech.EOMTrace is EOMTrace
    assert geomech.LinearizationTrace is LinearizationTrace
    from geomech.core.transformations.manifold_rules import apply_manifold_rules as amr

    assert apply_manifold_rules is amr
    assert geomech.core.transformations.apply_manifold_rules is amr


@pytest.mark.parametrize("build", [point_mass_r3, rigid_pendulum_pivot, quadrotor_with_load])
def test_trace_does_not_change_the_result(build):
    s = build()
    plain = compute_eom(s.L, s.dW, s.variables)
    tr = EOMTrace()
    traced = compute_eom(s.L, s.dW, s.variables, trace=tr)
    assert _strs(plain) == _strs(traced)
    assert {k: str(v) for k, v in tr.equations.items()} == _strs(plain)


def test_trace_is_keyword_only():
    s = point_mass_r3()
    with pytest.raises(TypeError):
        compute_eom(s.L, s.dW, s.variables, (), EOMTrace())


def test_stage_order_rigid_pendulum():
    s = rigid_pendulum_pivot()
    tr = EOMTrace()
    compute_eom(s.L, s.dW, s.variables, trace=tr)
    d_eta = "\\frac{d}{dt}(\\eta_{R})"
    assert tr.stage_names() == [
        "dL",
        "dW",
        "dS",
        "manifold_rules",
        "simplify",
        "ibp.prepare",
        f"ibp.collect[{d_eta}]",
        f"ibp.integrate[{d_eta}]",
        f"ibp.simplify[{d_eta}]",
        "ibp",
        "expand",
        "raw[\\eta_{R}]",
        "equation[\\eta_{R}]",
    ]
    assert [str(v) for v in tr.variation_vectors] == ["\\eta_{R}"]
    assert [str(t) for t in tr.ibp_targets] == [d_eta]
    assert set(tr.ibp[d_eta]) == {"collect", "integrate", "simplify"}


def test_each_stage_reproduces_from_the_previous_one():
    """The debugger's recipe of the audit's scenario (b), now with only public names."""
    s = rigid_pendulum_pivot()
    tr = EOMTrace()
    compute_eom(s.L, s.dW, s.variables, trace=tr)
    assert str(tr.dL) == str(s.L.delta())
    assert str(tr.dS) == str(s.L.delta() + s.dW)
    assert str(tr.after_manifold_rules) == str(apply_manifold_rules(tr.dS, s.variables))
    assert str(tr.after_simplify) == str(full_simplify(tr.after_manifold_rules))
    assert str(tr.after_ibp) == str(integrate_by_parts(tr.after_simplify, tr.ibp_targets))
    (target,) = tr.ibp_targets
    steps = tr.ibp[str(target)]
    assert str(steps["collect"]) == str(collect(tr.after_ibp_prepare, target))
    assert str(steps["integrate"]) == str(_apply_ibp(steps["collect"], target))
    assert str(steps["simplify"]) == str(full_simplify(steps["integrate"]))
    assert str(tr.after_expand) == str(expand(tr.after_ibp))
    (vec,) = tr.variation_vectors
    raw = tr.raw_equations[str(vec)]
    assert str(raw) == str(extract_linear_coeff(tr.after_expand, vec))
    assert str(tr.equations[str(vec)]) == str(full_simplify(raw))


def test_getitem_accepts_stage_and_attribute_names():
    s = point_mass_r3()
    tr = EOMTrace()
    compute_eom(s.L, s.dW, s.variables, trace=tr)
    assert tr["expand"] is tr.after_expand
    assert tr["after_expand"] is tr.after_expand
    assert tr["simplify"] is tr["after_simplify"]
    assert tr["equations"] is tr.equations
    with pytest.raises(KeyError, match="no stage 'nope'"):
        tr["nope"]


def test_summary_lists_every_stage_with_node_counts(capsys):
    s = spherical_pendulum_s2()
    tr = EOMTrace()
    compute_eom(s.L, s.dW, s.variables, trace=tr)
    text = tr.format_summary()
    lines = text.splitlines()
    assert len(lines) == len(tr.stages)
    for line, name in zip(lines, tr.stage_names()):
        assert line.startswith(name)
        assert " nodes  " in line
    # the repr is cut at 200 characters
    assert all(len(line.split(" nodes  ", 1)[1]) <= 200 for line in lines)
    tr.summary()
    assert capsys.readouterr().out.strip() == text
    assert str(tr) == text


def test_show_prints_the_indent_tree(capsys):
    s = point_mass_r3()
    tr = EOMTrace()
    compute_eom(s.L, s.dW, s.variables, trace=tr)
    tr.show("equation[\\delta{x}]")
    out = capsys.readouterr().out
    assert (
        "VAdd" in out and "├─" in out and out.rstrip().endswith(repr(tr.equations["\\delta{x}"]))
    )
    tr.show("raw_equations")
    assert "--- raw_equations[\\delta{x}]" in capsys.readouterr().out


def test_trace_keeps_the_stages_before_a_failure():
    """Audit C3: collect has no rule for Dot(w, Vee(Hat(Ω) A)); the trace shows how far
    the pipeline got and the term that kept its d/dt(η)."""
    R = SO3("R")
    Om = R.get_tangent_vector()
    A = Matrix("A", attr=["Constant"])
    w = Vector("w", attr=["Constant"])
    L = Dot(w, Vee(MMMul(Hat(Om), A)))
    tr = EOMTrace()
    with pytest.raises(AlgebraicError):
        compute_eom(L, None, SystemVariables(matrices=[R]), trace=tr)
    assert tr.stage_names()[-1] == "expand"
    assert tr.equations == {} and tr.raw_equations == {}
    assert "\\frac{d}{dt}(\\eta_{R})" in str(tr.after_expand)


def test_linearization_trace_spherical_pendulum():
    s = spherical_pendulum_s2()
    eom = compute_eom(s.L, s.dW, s.variables)
    plain = linearize(eom, s.variables, s.inputs)
    tr = LinearizationTrace()
    lin = linearize(eom, s.variables, s.inputs, trace=tr)
    for blocks in ("M", "A", "B"):
        assert {
            r: {k: str(v) for k, v in row.items()} for r, row in getattr(lin, blocks).items()
        } == {r: {k: str(v) for k, v in row.items()} for r, row in getattr(plain, blocks).items()}
    label = "\\xi_{q}"
    assert set(tr.dE) == {label}
    assert str(tr.dE[label]) == str(expand(full_simplify(eom[label][1].delta())))
    # every bucket has a raw coefficient and a final block, which is the one in M/A/B
    assert set(tr.buckets[label]) == set(tr.raw_blocks[label]) == set(tr.blocks[label])
    for sym, block in tr.blocks[label].items():
        target = lin.M if sym in lin.M[label] else lin.A if sym in lin.A[label] else lin.B
        assert str(target[label][sym]) == str(block)
    # the mass block is the negated, simplified raw coefficient
    dvel = "\\frac{d}{dt}(\\delta{\\omega_{q}})"
    raw = tr.raw_blocks[label][dvel]
    assert str(lin.M[label][dvel]) == str(full_simplify(full_simplify(raw) * -1))
    assert tr.stage_names()[0] == f"dE[{label}]"
    assert f"block[{label}][{dvel}]" in tr.stage_names()


def test_linearize_trace_is_keyword_only():
    s = point_mass_r3()
    eom = compute_eom(s.L, s.dW, s.variables)
    with pytest.raises(TypeError):
        linearize(eom, s.variables, s.inputs, (), LinearizationTrace())
