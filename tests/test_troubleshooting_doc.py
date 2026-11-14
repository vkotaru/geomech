"""Every error quoted in ``docs/troubleshooting.md`` is produced, word for word, by the code.

``tests/test_phase4_adversarial.py`` holds one trigger per row of the "Typed errors"
tables (``TRIGGERS``). This file adds the rows it left out (``linearize(..., params=[l])``)
and the third ``AlgebraicError`` block, and asserts that no row of the tables is left
without a check, so a new row cannot be added without its trigger.
"""

from __future__ import annotations

import re
from pathlib import Path

from geomech import (
    S2,
    Dot,
    EOMTrace,
    Scalar,
    SystemVariables,
    Vector,
    compute_eom,
    linearize,
)
from tests.test_phase4_adversarial import DOC, TRIGGERS, _msg

ROOT = Path(__file__).resolve().parent.parent
TEXT = (ROOT / "docs/troubleshooting.md").read_text(encoding="utf-8")


def _pendulum_with_length():
    m, g = Scalar("m", attr=["Constant"]), Scalar("g", attr=["Constant"])
    l = Scalar("l")
    e3 = Vector("e3", attr=["Constant"])
    q, f = S2("q"), Vector("f")
    om, xi = q.get_tangent_vector(), q.get_variation_vector()
    L = Scalar("0.5", value=0.5) * m * l * l * Dot(om, om) - m * g * l * Dot(q, e3)
    return L, Dot(xi, f), q, l, f


def _linearize_with_params():
    L, dW, q, l, f = _pendulum_with_length()
    variables = SystemVariables(vectors=[q])
    eom = compute_eom(L, dW, variables, params=[l])
    linearize(eom, variables, [f], params=[l])


def _second_order_lagrangian():
    x = Vector("x")
    compute_eom(
        Dot(x.t_diff().t_diff(), x.t_diff().t_diff()) * 0.5, None, SystemVariables([], [x])
    )


EXTRA = {
    "equations with a prescribed param, `linearize(..., params=[l])`": _linearize_with_params,
}

# rows checked elsewhere: the warning's message ends in a repr placeholder
# (test_troubleshooting_full_simplify_warning_prefix)
CHECKED_ELSEWHERE = {"`full_simplify(expr, max_iter=1)` at the cap"}


def _block(label):
    return TEXT.split(label, 1)[1].split("```", 2)[1]


def test_linearize_params_row_matches_the_code():
    label = "equations with a prescribed param, `linearize(..., params=[l])`"
    assert label in DOC
    assert _msg(EXTRA[label]) == DOC[label]


def test_every_typed_error_row_has_a_trigger():
    missing = sorted(set(DOC) - set(TRIGGERS) - set(EXTRA) - CHECKED_ELSEWHERE)
    assert not missing, f"rows of troubleshooting.md without a trigger: {missing}"


def test_second_order_lagrangian_block_matches_the_code():
    expected = "".join(_block("**3. A Lagrangian in accelerations**").split())
    got = "".join(_msg(_second_order_lagrangian).split())
    assert got == expected


def test_second_order_lagrangian_points_to_a_recorded_stage():
    x = Vector("x")
    tr = EOMTrace()
    try:
        compute_eom(
            Dot(x.t_diff().t_diff(), x.t_diff().t_diff()) * 0.5,
            None,
            SystemVariables(vectors=[x]),
            trace=tr,
        )
    except Exception as err:  # noqa: BLE001
        stages = re.findall(r"stage '([^']*)'", str(err))
    assert stages and all(s in tr.stage_names() for s in stages)
    assert tr.failed_at == "raw[\\delta{x}]" and tr.error == "AlgebraicError"


def test_section_links_resolve():
    """Every ``(#anchor)`` link of the page points to one of its headings."""
    anchors = set()
    for heading in re.findall(r"^#+ (.+)$", TEXT, flags=re.MULTILINE):
        slug = re.sub(r"[^\w\- ]", "", heading.lower()).strip().replace(" ", "-")
        anchors.add(slug)
    links = set(re.findall(r"\]\(#([^)]+)\)", TEXT))
    assert links and links <= anchors, sorted(links - anchors)
