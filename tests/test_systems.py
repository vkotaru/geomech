"""The API of ``geomech.systems`` (spec ``specs/systems-library.md``): the registry, the
``System`` record, the conventions every system follows, and the docs table (item 5)."""

from __future__ import annotations

import dataclasses
import re
from pathlib import Path

import numpy as np
import pytest

import geomech.systems as systems
from geomech import S2, SO3, Matrix, Scalar, Vector
from geomech.core.base.expressions import TS2, TSO3, is_number_literal
from geomech.core.operations import TimeDerivative, Variation
from geomech.systems import MAX_SIZE, Reference, System, get, list_systems, systems_table

ROOT = Path(__file__).resolve().parents[1]

#: The table of the spec, in order: name → dof at the builder's defaults.
EXPECTED = {
    "point_mass": 3,
    "spring_mass": 1,
    "planar_pendulum": 1,
    "double_planar_pendulum": 2,
    "spherical_pendulum": 2,
    "particle_on_rotating_sphere": 2,
    "double_spherical_pendulum": 4,
    "n_link_spherical_chain": 6,
    "elastic_spherical_pendulum": 3,
    "cart_spherical_pendulum": 3,
    "free_rigid_body": 3,
    "heavy_top": 3,
    "rigid_pendulum": 3,
    "double_rigid_pendulum": 6,
    "n_link_rigid_chain": 9,
    "two_rigid_bodies_spherical_joint": 9,
    "rigid_body": 6,
    "spacecraft_reaction_wheels": 6,
    "spacecraft_reaction_wheels_s1": 6,
    "biped_3d": 6,
    "quadrotor": 6,
    "quadrotor_offset_com": 6,
    "quadrotor_suspended_load": 8,
    "quadrotor_elastic_load": 9,
    "quadrotor_flexible_cable": 12,
    "multi_quadrotor_point_mass": 18,
    "multi_quadrotor_rigid_load": 21,
}
PARAMETRIC = {
    "n_link_spherical_chain": "n_links",
    "n_link_rigid_chain": "n_links",
    "quadrotor_flexible_cable": "n_links",
    "multi_quadrotor_point_mass": "n_quads",
    "multi_quadrotor_rigid_load": "n_quads",
}
NAMES = list(EXPECTED)


def _cases():
    out = [(n, {}) for n in NAMES]
    out += [(n, {p: size}) for size in (1, 2) for n, p in PARAMETRIC.items()]
    return out


CASES = _cases()
IDS = [n + ("" if not kw else f"-{next(iter(kw.values()))}") for n, kw in CASES]


def test_list_systems_is_the_spec_table_in_order():
    infos = list_systems()
    assert [i.name for i in infos] == NAMES
    assert {i.name: i.dof for i in infos} == EXPECTED
    assert {i.name: i.parameters for i in infos if i.parameters} == {
        k: (v,) for k, v in PARAMETRIC.items()
    }
    assert list_systems() == infos  # stable
    for info in infos:
        assert info.summary.endswith(".") and len(info.summary) > 20, info


def test_every_name_is_a_builder_and_get_rebuilds_it():
    for name in NAMES:
        assert getattr(systems, name).__name__ == name
        assert name in systems.__all__
        s = get(name)
        assert isinstance(s, System) and s.name == name


def test_get_unknown_name_lists_the_systems():
    with pytest.raises(KeyError, match="no system named 'pendulum'.*spherical_pendulum"):
        get("pendulum")


@pytest.mark.parametrize("name, param", list(PARAMETRIC.items()))
def test_parametric_builders_take_keyword_arguments_only(name, param):
    build = getattr(systems, name)
    with pytest.raises(TypeError):
        build(2)
    for bad in (0, -1, 2.0, True, "2", MAX_SIZE + 1, 10**6):
        with pytest.raises(ValueError, match=f"{param} must be an int from 1 to {MAX_SIZE}"):
            build(**{param: bad})
    assert build(**{param: 1}).options == {param: 1}
    assert build(**{param: MAX_SIZE}).options == {param: MAX_SIZE}
    small, large = build(**{param: 2}), build(**{param: 4})
    assert small.options == {param: 2} and large.options == {param: 4}
    assert large.dof > small.dof


@pytest.mark.parametrize("name", [n for n in NAMES if n not in PARAMETRIC])
def test_plain_builders_take_no_argument(name):
    with pytest.raises(TypeError):
        get(name, n=2)


def test_system_is_frozen():
    s = get("spherical_pendulum")
    with pytest.raises(dataclasses.FrozenInstanceError):
        s.L = None  # type: ignore[misc]
    assert not s.x0.flags.writeable


@pytest.mark.parametrize("name", NAMES)
def test_system_shares_no_mutable_state(name):
    """params, symbols and options are read-only mappings, every array is read-only and
    owned by the system: mutating one system cannot change another."""
    s = get(name)
    key = next(iter(s.params))
    with pytest.raises(TypeError):
        s.params[key] = 0.0  # type: ignore[index]
    with pytest.raises(TypeError):
        s.symbols["new"] = 0  # type: ignore[index]
    with pytest.raises(TypeError):
        s.options["n"] = 0  # type: ignore[index]
    for k, val in s.params.items():
        if isinstance(val, np.ndarray):
            assert not val.flags.writeable, k
            with pytest.raises(ValueError):
                val[...] = 0.0
            assert not np.shares_memory(val, get(name).params[k]), k
        else:
            assert isinstance(val, float)
    for val in s.symbols.values():
        assert not isinstance(val, (list, dict, np.ndarray))
    p = dict(s.params)  # the documented way to change a value
    p[key] = 123.0
    assert get(name).params[key] is not p[key]


@pytest.mark.parametrize("name", NAMES)
def test_system_deep_copies_and_its_mappings_pickle(name):
    """Read-only mappings must not break copy.deepcopy or pickling (multiprocessing,
    joblib): the copy is equal in value, independent of the original, and still read-only."""
    import copy
    import pickle

    s = get(name)
    c = copy.deepcopy(s)
    assert c.name == s.name and str(c.L) == str(s.L) and np.array_equal(c.x0, s.x0)
    for k, val in s.params.items():
        assert np.array_equal(c.params[k], val), k
    with pytest.raises(TypeError, match="read-only; copy first"):
        c.params["m"] = 1.0  # type: ignore[index]
    for mapping in (s.params, s.symbols, s.options):
        back = pickle.loads(pickle.dumps(mapping))
        assert back.keys() == mapping.keys() and type(back) is type(mapping)
        with pytest.raises(TypeError):
            back["x"] = 0  # type: ignore[index]


def test_read_only_error_says_how_to_change_a_value():
    s = get("spherical_pendulum")
    for write in (
        lambda: s.params.update(m=2.0),
        lambda: s.params.pop("m"),
        lambda: s.params.setdefault("z", 1.0),
        lambda: s.params.clear(),
    ):
        with pytest.raises(TypeError, match=r"p = dict\(sys.params\)"):
            write()


def test_module_constants_are_read_only():
    from geomech.systems import _base, aerial, rigid

    for arr in (_base.E3, _base.J_A, _base.J_B, aerial.J_QUAD, rigid.J_LINK1, rigid.J_LINK2):
        with pytest.raises(ValueError):
            arr[0] = -1.0


#: Inertias given about a pivot, not about the centre of mass:
#: system → {inertia: (mass key, function of params → centre of mass from the pivot)}.
PIVOT_INERTIAS = {
    "heavy_top": {"J": ("m", lambda p: p["l"] * np.array([0.0, 0.0, 1.0]))},
    "rigid_pendulum": {"J": ("m", lambda p: np.asarray(p["\\rho"]))},
}


@pytest.mark.parametrize("name, kwargs", CASES, ids=IDS)
def test_default_parameters_are_physical(name, kwargs):
    """e₃ = (0, 0, 1); e₁ a horizontal unit vector; masses, lengths, g, stiffness positive;
    every inertia (shifted to the centre of mass where it is given about a pivot) symmetric
    positive definite with principal moments satisfying the triangle inequality."""
    s = get(name, **kwargs)
    p = s.params
    if "e3" in p:
        assert np.array_equal(p["e3"], [0.0, 0.0, 1.0])
    if "e1" in p:
        assert np.array_equal(p["e1"], [1.0, 0.0, 0.0])
    for k, val in p.items():
        if np.isscalar(val) and k[0] in "mlgkML":
            assert val > 0, (k, val)
    if "g" in p:
        assert p["g"] == 9.81
    shifts = PIVOT_INERTIAS.get(name, {})
    inertias = [k for k in p if k.startswith("J")]
    for k in inertias:
        J = np.asarray(p[k], dtype=float)
        if k in shifts:
            mkey, com = shifts[k]
            r = com(p)
            J = J - p[mkey] * ((r @ r) * np.eye(3) - np.outer(r, r))
        assert np.allclose(J, J.T), k
        ev = np.linalg.eigvalsh(J)
        assert ev[0] > 0, f"{name}: {k} about the centre of mass has eigenvalues {ev}"
        assert ev[0] + ev[1] >= ev[2] - 1e-12, f"{name}: {k} violates the triangle inequality"
    if name in ("free_rigid_body", "rigid_pendulum", "heavy_top", "quadrotor"):
        assert inertias, name


def _leaves(expr, out):
    nodes = getattr(expr, "nodes", None)
    if isinstance(expr, (TimeDerivative, Variation)) or not nodes:
        out.append(expr)
        return
    for n in nodes:
        _leaves(n, out)


@pytest.mark.parametrize("name, kwargs", CASES, ids=IDS)
def test_symbols_follow_the_conventions(name, kwargs):
    """Unique printed names; every parameter is a constant leaf (no prescribed params);
    inputs are non-constant leaves; the default params are exactly the layout's."""
    s = get(name, **kwargs)
    leaves: list = []
    _leaves(s.L, leaves)
    _leaves(s.dW, leaves)
    seen: dict = {}
    for leaf in leaves:
        if is_number_literal(leaf):
            continue
        key = str(leaf)
        seen.setdefault(key, repr(leaf))
        assert seen[key] == repr(leaf), f"two symbols print as {key!r}"
    layout = s.layout()
    for p in layout.params:
        assert isinstance(p, (Scalar, Vector, Matrix)) and p.is_constant, p
        assert not isinstance(p, (S2, SO3, TS2, TSO3))
    for u in s.inputs:
        assert not u.is_constant
    assert set(s.params) == set(layout.param_names)
    for key, val in s.params.items():
        shape = layout.slot(key).shape
        assert np.shape(val) == shape, (key, np.shape(val), shape)


@pytest.mark.parametrize("name, kwargs", CASES, ids=IDS)
def test_x0_is_on_the_manifold(name, kwargs):
    s = get(name, **kwargs)
    layout = s.layout()
    vals = layout.unpack(s.x0)
    assert s.x0.shape == (layout.n_state,)
    for slot in layout.configuration:
        v = vals[slot.name]
        if slot.manifold == "S2":
            w = vals[f"omega_{slot.name}"]
            assert abs(np.linalg.norm(v) - 1) < 1e-12 and abs(v @ w) < 1e-12
            assert np.linalg.norm(w) > 0.1
        if slot.manifold == "SO3":
            assert np.allclose(v.T @ v, np.eye(3), atol=1e-12) and np.linalg.det(v) > 0


@pytest.mark.parametrize("name", NAMES)
def test_reference_is_documented(name):
    ref = get(name).reference
    assert isinstance(ref, Reference)
    bib = (ROOT / "docs" / "references.bib").read_text()
    keys = set(re.findall(r"@\w+\{(\w+),", bib))
    assert ref.source in keys | {"elementary", "check_eom"}
    if ref.source == "lee2018global":
        assert ref.section, f"{name}: a Lee et al. reference names its section"
    if ref.equation:  # a compared equation: see tests/test_systems_book.py
        assert ref.source == "lee2018global"
    assert ref.note
    doc = getattr(systems, name).__doc__
    assert doc and "L" in doc and "δW" in doc, f"{name}: state L and δW in the docstring"


def test_docs_table_matches_list_systems():
    page = (ROOT / "docs" / "systems.md").read_text()
    m = re.search(r"<!-- systems-table -->\n(.*?)<!-- /systems-table -->", page, re.S)
    assert m, "docs/systems.md has no <!-- systems-table --> block"
    assert m.group(1) == systems_table(), (
        "docs/systems.md is out of date: paste geomech.systems.systems_table() between the"
        " <!-- systems-table --> markers"
    )


def test_api_reference_documents_every_builder():
    page = (ROOT / "docs" / "api" / "systems.rst").read_text()
    for mod in ("spheres", "rigid", "aerial"):
        assert f"geomech.systems.{mod}" in page


def test_skill_mentions_the_library():
    a = (ROOT / "skills" / "geomech" / "SKILL.md").read_bytes()
    b = (ROOT / ".claude" / "skills" / "geomech" / "SKILL.md").read_bytes()
    assert a == b
    assert b"geomech.systems" in a and b"get(" in a
