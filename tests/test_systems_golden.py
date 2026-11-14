"""The stored reference output of every system (spec ``specs/systems-library.md``,
verification item 3).

``tests/golden/<name>.json`` holds each equation and standard-form block of a system as a
printed string and as values at fixed stored states (``scripts/regenerate_golden.py``
documents the format and the tolerance). Each test re-derives one system: a changed value
fails as a change of the equations; a changed string with every value unchanged fails
with the instruction to regenerate the files.
"""

from __future__ import annotations

import copy
import importlib.util
import json

import pytest

from tests.test_systems import ROOT


def _script():
    spec = importlib.util.spec_from_file_location(
        "regenerate_golden", ROOT / "scripts" / "regenerate_golden.py"
    )
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


SCRIPT = _script()
ENTRIES = SCRIPT.entries()


def test_golden_file_set_is_complete():
    stored = sorted(p.stem for p in SCRIPT.OUT.glob("*.json"))
    assert stored == sorted(stem for stem, _, _ in ENTRIES), (
        f"tests/golden/ does not match geomech.systems; run `{SCRIPT.REGENERATE}`"
    )
    # every parametric system is stored at two sizes
    assert len(ENTRIES) == len({name for _, name, _ in ENTRIES}) + 4


@pytest.mark.parametrize("stem, name, kwargs", ENTRIES, ids=[e[0] for e in ENTRIES])
def test_golden(stem, name, kwargs):
    values, strings = SCRIPT.check_one(stem, name, kwargs)
    assert not values, "the equations changed:\n" + "\n".join(values)
    assert not strings, (
        "the printed form changed while every value agrees:\n"
        + "\n".join(strings)
        + f"\nCheck the numbers, then regenerate the files with `{SCRIPT.REGENERATE}`."
    )


def _stored(stem):
    return json.loads((SCRIPT.OUT / f"{stem}.json").read_text())


def test_golden_files_store_states_values_and_strings():
    data = _stored("quadrotor_suspended_load")
    assert data["rtol"] == SCRIPT.RTOL and len(data["states"]) == SCRIPT.N_STATES
    labels = list(data["blocks"])
    assert any(lab.startswith("eom[") for lab in labels)
    assert any(lab.startswith("M[") for lab in labels)
    assert any(lab.startswith("G[") for lab in labels)
    for block in data["blocks"].values():
        assert block["string"] and len(block["values"]) == SCRIPT.N_STATES


def test_a_changed_value_is_reported_as_a_change_of_the_equations():
    stem, name, kwargs = ENTRIES[0]
    stored = _stored(stem)
    label = next(iter(stored["blocks"]))
    bad = copy.deepcopy(stored)
    bad["blocks"][label]["values"][0][0] += 1e-3
    _, exprs = SCRIPT.derive(name, kwargs)
    values, strings = SCRIPT.compare(bad, exprs)
    assert values and label in values[0] and "value changed" in values[0]


def test_a_changed_string_alone_asks_for_regeneration():
    stem, name, kwargs = ENTRIES[0]
    stored = _stored(stem)
    label = next(iter(stored["blocks"]))
    bad = copy.deepcopy(stored)
    bad["blocks"][label]["string"] = "an older printing"
    _, exprs = SCRIPT.derive(name, kwargs)
    values, strings = SCRIPT.compare(bad, exprs)
    assert not values and strings and "an older printing" in strings[0]


def test_check_mode_passes_on_the_checked_in_files(capsys):
    assert SCRIPT.main(["--check"]) == 0
    assert "checked" in capsys.readouterr().out


def test_a_changed_default_parameter_or_x0_is_reported():
    stored = _stored("quadrotor")
    s, _ = SCRIPT.derive("quadrotor", {})
    assert SCRIPT.compare_defaults(stored, s) == []
    bad = copy.deepcopy(stored)
    bad["params"]["e3"] = [0.0, 0.0, -1.0]
    bad["x0"][0] += 1e-3
    problems = SCRIPT.compare_defaults(bad, s)
    assert any("params['e3']" in m for m in problems) and any("x0" in m for m in problems)


def test_a_hand_edited_tolerance_is_reported():
    stem, name, kwargs = ENTRIES[0]
    bad = copy.deepcopy(_stored(stem))
    bad["rtol"] = 1e9
    _, exprs = SCRIPT.derive(name, kwargs)
    values, _ = SCRIPT.compare(bad, exprs)
    assert values and "rtol" in values[0]
