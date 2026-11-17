"""The model pages in ``docs/models/`` show what the library derives.

* The math between ``<!-- model:energy -->`` and ``<!-- model:eom -->`` markers is what
  ``scripts/model_pages.py`` generates from the library today.
* The bounds stated in the checks table hold.
* The code block of each page builds the library system: a hand-written block builds
  the same ``L`` and ``δW``; a block that calls ``get(...)`` gets the right system. On a
  page with a symbolic number of bodies n, the block builds the ``geomech.indexed`` model
  and the ``<!-- model:indexed -->`` block is its generic derivation.
* ``docs/models/index.md`` lists every page.
"""

import re
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import model_pages  # noqa: E402

from geomech.numeric import check_eom  # noqa: E402
from geomech.systems import get, list_systems  # noqa: E402

PAGES = model_pages.existing()
IDS = [p.name for p in PAGES]


def _page(page):
    return (model_pages.PAGES / f"{page.name}.md").read_text()


def _block(text, block):
    m = re.search(rf"<!-- model:{block} -->\n(.*?)<!-- /model:{block} -->", text, re.S)
    assert m, f"no <!-- model:{block} --> block"
    return m.group(1)


@pytest.fixture(scope="module")
def derived():
    cache = {}

    def get_derived(page):
        if page.name not in cache:
            system = page.system()
            cache[page.name] = (system, system.derive())
        return cache[page.name]

    return get_derived


@pytest.mark.parametrize("page", PAGES, ids=IDS)
def test_generated_math_matches_the_library(page, derived):
    system, eom = derived(page)
    text = _page(page)
    stale = (
        f"docs/models/{page.name}.md is out of date: run "
        f"python scripts/model_pages.py --only {page.name}"
    )
    for block, make in (("energy", model_pages.energy_block), ("eom", model_pages.eom_block)):
        assert _block(text, block) == make(system, eom), stale
    if page.indexed:
        assert _block(text, "indexed") == model_pages.indexed_block(page), stale
    else:
        assert "<!-- model:indexed -->" not in text


@pytest.mark.parametrize("page", PAGES, ids=IDS)
def test_checks_table_bounds_hold(page, derived):
    system, eom = derived(page)
    bounds = [
        10.0 ** int(b)
        for b in re.findall(r"below \$10\^\{(-?\d+)\}\$", _block(_page(page), "checks"))
    ]
    err = check_eom(system.L, system.dW, system.variables, eom, inputs=system.inputs)
    _, drift = model_pages.energy_drift(system, eom, page)
    measured = [err, drift]
    if page.indexed:
        # the generic equations at n = size, and the generic L against the library's
        measured += list(model_pages.indexed_checks(page))
    assert len(bounds) == len(measured), "one stated bound per measured check"
    for value, bound in zip(measured, bounds):
        assert value < bound


@pytest.mark.parametrize("page", PAGES, ids=IDS)
def test_code_block_builds_the_library_system(page):
    code = re.search(r"## Code\n(?:[^#`].*?\n)*?\n?```python\n(.*?)```", _page(page), re.S)
    assert code, f"docs/models/{page.name}.md has no python block under '## Code'"
    ns: dict = {}
    exec(code.group(1), ns)
    system = page.system()
    if page.indexed:
        model = model_pages.indexed_model(page)
        assert str(ns["L"]) == str(model.lagrangian)
        assert str(ns["dW"]) == str(model.inf_work)
        assert set(ns["eom"]) == set(model.eom())
        assert ns["system"].name == page.name
        assert len(ns["system"].x0) == len(system.x0), "the block builds the page's n"
    elif "L" in ns:
        assert str(ns["L"]) == str(system.L)
        assert str(ns["dW"]) == str(system.dW)
        assert set(ns["eom"]) == set(system.derive())
    else:
        assert ns["system"].name == page.name, "the block builds L, dW, eom or a system"


def test_index_lists_every_page():
    assert (model_pages.PAGES / "index.md").read_text() == model_pages.index_page(), (
        "docs/models/index.md is out of date: run python scripts/model_pages.py"
    )


def test_every_model_is_a_library_system():
    names = {info.name for info in list_systems()}
    assert {p.name for p in model_pages.MODELS} <= names


def test_render_rejects_a_page_without_a_block():
    with pytest.raises(ValueError, match="eom"):
        model_pages.render(
            "<!-- model:energy -->\nx<!-- /model:energy -->\n", {"energy": "", "eom": ""}
        )


def test_hover_carries_the_weight():
    system = get("multi_quadrotor_point_mass")
    layout = system.layout()
    u = model_pages.hover(system, layout)
    p = system.params
    total = p["m_L"] + p["m1"] + p["m2"] + p["m3"]
    thrusts = [u[s] for k, s in layout.input_slices.items() if k.startswith("f")]
    assert np.isclose(np.sum(thrusts), total * p["g"])
