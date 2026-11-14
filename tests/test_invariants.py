"""Phase 4.4: ``check_tree`` and the configuration-variable leaf check.

Each invariant gets a tree that violates it (built around the constructor guards by
passing a list to ``Add`` or by mutating ``nodes``, the two ways the review found to get
past them) and a check that the message names the invariant, the subtree and its path.
"""

from __future__ import annotations

import pytest

import geomech
from geomech import (
    S2,
    SO3,
    TS2,
    Add,
    Cross,
    Dot,
    Hat,
    InvariantError,
    Number,
    Scalar,
    SystemVariables,
    TimeDerivative,
    Transpose,
    VAdd,
    Variation,
    Vector,
    VVMul,
    check_tree,
    compute_eom,
    getScalars,
    linearize,
    to_standard_form,
)
from tests import conftest
from tests.vbl_systems import HALF, quadrotor_with_load, spherical_pendulum_s2


def _vectors():
    return Vector("a"), Vector("b")


def test_exports():
    assert geomech.check_tree is check_tree
    assert geomech.InvariantError is InvariantError
    assert issubclass(InvariantError, Exception)


def test_pipeline_outputs_pass():
    s = quadrotor_with_load()
    eom = compute_eom(s.L, s.dW, s.variables)
    for tree in (s.L, s.dW, *(eq for _, eq in eom.values())):
        assert check_tree(tree) is None
    for sf in to_standard_form(eom, s.variables, s.inputs).values():
        for block in (*sf.M.values(), sf.f, *sf.G.values()):
            check_tree(block)
    lin = linearize(eom, s.variables, s.inputs)
    for blocks in (lin.M, lin.A, lin.B):
        for row in blocks.values():
            for block in row.values():
                check_tree(block)


def test_conftest_fixture_checks_pipeline_calls(check_pipeline_trees):
    before = dict(check_pipeline_trees)
    s = spherical_pendulum_s2()
    eom = compute_eom(s.L, s.dW, s.variables)
    to_standard_form(eom, s.variables, s.inputs)
    linearize(eom, s.variables, s.inputs)
    for name in ("compute_eom", "to_standard_form", "linearize"):
        assert check_pipeline_trees[name] == before[name] + 1


# ---------------------------------------------------------------------------
# one test per invariant
# ---------------------------------------------------------------------------


def test_nested_same_class_sum():
    a, b, c = getScalars("a b c")
    nested = Add([Add(a, b), c])  # a list operand is not flattened
    with pytest.raises(InvariantError, match=r"Add nested directly in a Add at root\.0: Add"):
        check_tree(nested)
    check_tree(Add(Add(a, b), c))  # flattened by the constructor


def test_operand_type_of_a_sum():
    a, b = getScalars("a b")
    node = Add(a, b)
    node.nodes[1] = Vector("x")
    with pytest.raises(
        InvariantError, match=r"Add addend has type VECTOR, expected SCALAR at root\.1"
    ):
        check_tree(node)


def test_operand_types_of_a_product():
    a, b = _vectors()
    node = Cross(a, b)
    node.nodes[0] = Scalar("m")
    with pytest.raises(InvariantError, match=r"Cross operand types SCALAR x VECTOR"):
        check_tree(node)


def test_vvmul_needs_exactly_one_row():
    a, b = _vectors()
    node = VVMul(Transpose(a), b)
    check_tree(node)
    node.nodes[0] = a  # two columns
    with pytest.raises(InvariantError, match="exactly one row operand; both operands are columns"):
        check_tree(node)


@pytest.mark.parametrize(
    "build, path",
    [
        (lambda a, b: Dot(Transpose(a), b), "root.L"),
        (lambda a, b: Cross(a, Transpose(b)), "root.R"),
    ],
    ids=["Dot", "Cross"],
)
def test_dot_and_cross_take_columns(build, path):
    """Phase 3 review F2: rows were accepted where only columns make sense."""
    a, b = _vectors()
    with pytest.raises(InvariantError, match=rf"operand is a row vector.* at {path}: Transpose"):
        check_tree(build(a, b))


def test_hat_takes_a_column():
    a, _ = _vectors()
    with pytest.raises(InvariantError, match=r"Hat operand is a row vector.*root\.x"):
        check_tree(Hat(Transpose(a)))


def test_vadd_does_not_mix_rows_and_columns():
    a, b = _vectors()
    with pytest.raises(InvariantError, match="VAdd mixes row and column vectors at root"):
        check_tree(Transpose(a) + b)
    check_tree(Transpose(a) + Transpose(b))
    check_tree(VAdd(Transpose(a), geomech.ZeroVector))  # a zero addend fits either


def test_no_calculus_over_a_manifold_point():
    q = S2("q")
    node = TimeDerivative(Vector("x"))
    node.nodes[0] = q
    with pytest.raises(
        InvariantError, match="TimeDerivative directly over the manifold point 'q'"
    ):
        check_tree(node)
    node = Variation(Vector("x"))
    node.nodes[0] = SO3("R")
    with pytest.raises(InvariantError, match="Variation directly over the manifold point 'R'"):
        check_tree(node)


def test_no_variation_of_a_variation_vector():
    q = S2("q")
    node = Variation(q.get_tangent_vector())
    node.nodes[0] = q.get_variation_vector()
    with pytest.raises(InvariantError, match="variation of a variation vector"):
        check_tree(node)


def test_tangent_vector_needs_its_point():
    with pytest.raises(InvariantError, match=r"TS2 leaf has no S2 point \(\.S2 is None\)"):
        check_tree(Dot(TS2("w"), Vector("a")))


def test_tangent_vector_named_like_its_point():
    q = S2("q")
    with pytest.raises(InvariantError, match="named neither like the tangent vector"):
        check_tree(TS2("w", S2=q))
    check_tree(q.get_tangent_vector())
    check_tree(q.get_variation_vector())


def test_numeric_scalars_are_consistent():
    check_tree(Number(-1) * Scalar("m"))
    check_tree(HALF)
    with pytest.raises(InvariantError, match="printed like a number but has no value"):
        check_tree(Scalar("(-1)") * Scalar("m"))
    bad = Scalar("k", value=2)
    bad.flags.is_numeric = False
    with pytest.raises(InvariantError, match="inconsistent with is_numeric=False"):
        check_tree(bad)
    with pytest.raises(InvariantError, match=r"printed '\(2\)' has value 3"):
        check_tree(Scalar("(2)", value=3, attr=["Constant"]))


def test_unary_needs_its_operand():
    node = Transpose(Vector("a"))
    node.nodes = []
    with pytest.raises(InvariantError, match="Transpose needs exactly one operand"):
        check_tree(node)


def test_not_an_expression():
    node = Add(*getScalars("a b"))
    node.nodes[0] = 3
    with pytest.raises(InvariantError, match=r"operand 0 of Add at root is a int \(3\)"):
        check_tree(node)


def test_error_names_the_deep_path():
    a, b = _vectors()
    m = Scalar("m")
    tree = m * Dot(a, b) + m * Dot(Transpose(a), b)
    with pytest.raises(InvariantError, match=r"at root\.1\.R\.L: Transpose '\(a\)''"):
        check_tree(tree)


# ---------------------------------------------------------------------------
# GEOMECH_DEBUG runs check_tree at pipeline entry
# ---------------------------------------------------------------------------


def _nested_add_lagrangian():
    x = Scalar("x")
    v = x.t_diff()
    (m,) = getScalars("m", attr=["Constant"])
    return Add([Add(HALF * m * v * v, x * x)]), x


def test_debug_mode_checks_the_input(monkeypatch):
    L, x = _nested_add_lagrangian()
    variables = SystemVariables(scalars=[x])
    monkeypatch.delenv("GEOMECH_DEBUG", raising=False)
    unchecked = conftest._WRAPPERS["compute_eom"][0]  # not wrapped by the suite fixture
    unchecked(L, None, variables)  # fine without the switch
    monkeypatch.setenv("GEOMECH_DEBUG", "1")
    with pytest.raises(InvariantError, match="Add nested directly in a Add"):
        compute_eom(L, None, variables)


def test_debug_mode_checks_linearize_and_standard_form_inputs(monkeypatch):
    s = spherical_pendulum_s2()
    eom = compute_eom(s.L, s.dW, s.variables)
    key = next(iter(eom))
    vec, eq = eom[key]
    bad = {key: (vec, VAdd([VAdd(eq, Vector("a"))]))}
    monkeypatch.setenv("GEOMECH_DEBUG", "1")
    with pytest.raises(InvariantError):
        linearize(bad, s.variables, s.inputs)
    with pytest.raises(InvariantError):
        to_standard_form(bad, s.variables, s.inputs)


# ---------------------------------------------------------------------------
# Phase 2 F-2: a configuration variable is a single leaf
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "entry",
    [
        lambda: Vector("x") + Vector("y"),
        lambda: Scalar("m") * Vector("z"),
        lambda: Cross(Vector("a"), Vector("b")),
    ],
    ids=["sum", "scaled", "cross"],
)
def test_composite_configuration_variable_rejected(entry):
    var = entry()
    with pytest.raises(TypeError, match="expression, not a single symbol"):
        compute_eom(Dot(var, var), None, SystemVariables(vectors=[var]))
