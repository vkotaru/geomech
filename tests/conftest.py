"""Suite-wide fixtures.

``check_pipeline_trees`` (session, autouse) runs ``geomech.check_tree`` on every tree that
goes into or comes out of a successful ``compute_eom``, ``to_standard_form`` or
``linearize`` call anywhere in the suite: the Lagrangian and virtual work, every
equation, every standard-form block and every linearization block. A test that expects
the pipeline to raise is unaffected (the check runs only after a successful call), and a
pipeline that builds a tree violating an invariant fails the test that ran it with
``InvariantError`` naming the subtree.

Test modules bind the pipeline functions at import (``from geomech import compute_eom``),
so the fixture replaces the bound name in every loaded ``tests.*`` module, on the
``geomech`` packages and in the defining modules (so ``geomech.linearize is
geomech.dynamics.vbl.linearize`` still holds) for the session.
"""

from __future__ import annotations

import functools
import sys

import pytest

import geomech
import geomech.dynamics
import geomech.dynamics.eom
import geomech.dynamics.standard_form
import geomech.dynamics.vbl
from geomech import check_tree

CHECKED_CALLS = {"compute_eom": 0, "to_standard_form": 0, "linearize": 0}


def _check_all(*trees):
    for tree in trees:
        if tree is not None and hasattr(tree, "type"):
            check_tree(tree)


def _wrap_compute_eom(func):
    @functools.wraps(func)
    def wrapper(lagrangian, inf_work, variables, *args, **kwargs):
        out = func(lagrangian, inf_work, variables, *args, **kwargs)
        _check_all(lagrangian, inf_work)
        for _, eq in out.values():
            check_tree(eq)
        CHECKED_CALLS["compute_eom"] += 1
        return out

    return wrapper


def _wrap_to_standard_form(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        out = func(*args, **kwargs)
        for sf in out.values():
            _check_all(*sf.M.values(), sf.f, *sf.G.values())
        CHECKED_CALLS["to_standard_form"] += 1
        return out

    return wrapper


def _wrap_linearize(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        out = func(*args, **kwargs)
        for blocks in (out.M, out.A, out.B):
            for row in blocks.values():
                _check_all(*row.values())
        CHECKED_CALLS["linearize"] += 1
        return out

    return wrapper


_WRAPPERS = {
    "compute_eom": (geomech.dynamics.eom.compute_eom, _wrap_compute_eom),
    "to_standard_form": (geomech.dynamics.standard_form.to_standard_form, _wrap_to_standard_form),
    "linearize": (geomech.dynamics.vbl.linearize, _wrap_linearize),
}


@pytest.fixture(scope="session", autouse=True)
def check_pipeline_trees():
    patch = pytest.MonkeyPatch()
    wrapped = {name: wrap(orig) for name, (orig, wrap) in _WRAPPERS.items()}
    targets = [geomech, geomech.dynamics]
    targets += [geomech.dynamics.eom, geomech.dynamics.standard_form, geomech.dynamics.vbl]
    targets += [m for name, m in list(sys.modules.items()) if name.startswith("tests.") and m]
    for module in targets:
        for name, (orig, _) in _WRAPPERS.items():
            if getattr(module, name, None) is orig:
                patch.setattr(module, name, wrapped[name])
    yield CHECKED_CALLS
    patch.undo()
