"""Adversarial review of phase 4 (debuggability).

Confirmed findings are strict xfails (``reason="F<n>: ..."``); the passing tests pin
what held: trace fidelity, failure behaviour, zero cost when off, the rule-name sync
test, ``check_tree`` on bad and on fuzzed trees, the conftest fixture, error chaining,
and every error text quoted in ``docs/troubleshooting.md``.
"""

from __future__ import annotations

import ast
import logging
import pickle
import re
import subprocess
import sys
import types
import warnings
from pathlib import Path

import pytest

import geomech
from geomech import (
    S2,
    SO3,
    TS2,
    Add,
    AlgebraicError,
    Cross,
    Dot,
    EOMTrace,
    Hat,
    InvariantError,
    LinearizationTrace,
    Matrix,
    MVMul,
    NonLinearError,
    Number,
    Scalar,
    SystemVariables,
    TimeDerivative,
    Transpose,
    Variation,
    Vector,
    Vee,
    check_tree,
    compute_eom,
    extract_linear_coeff,
    full_simplify,
    linearize,
    set_trace,
    to_standard_form,
)
from geomech.utils import trace as trace_mod
from tests.vbl_systems import (
    point_mass_r3,
    quadrotor_with_load,
    rigid_pendulum_pivot,
    spherical_pendulum_s2,
)

ROOT = Path(__file__).resolve().parent.parent


def _strs(eom):
    return {k: str(eq) for k, (_, eq) in eom.items()}


@pytest.fixture
def rules_log():
    """Capture ``geomech.rules`` lines with tracing on (no stderr handler involved)."""
    lines = []

    class _H(logging.Handler):
        def emit(self, record):
            lines.append(record.getMessage())

    handler = _H()
    log = logging.getLogger("geomech.rules")
    log.addHandler(handler)
    set_trace(True)
    yield lines
    set_trace(False)
    log.removeHandler(handler)


# ---------------------------------------------------------------------------
# Held: trace fidelity
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "build", [point_mass_r3, rigid_pendulum_pivot, spherical_pendulum_s2, quadrotor_with_load]
)
def test_trace_stages_are_the_objects_the_pipeline_used(build, monkeypatch):
    """Each recorded stage is the very object handed to the next function, not a
    recomputation: spy on the functions compute_eom calls and compare identities."""
    import geomech.core.math.ibp as ibp_mod
    import geomech.dynamics.eom as eom_mod

    seen = {}
    orig_ibp = eom_mod._integrate_by_parts
    orig_expand = eom_mod.expand
    orig_sep = eom_mod._separate_variations

    def spy_ibp(expr, targets, record=None):
        seen["ibp_in"] = expr
        out = orig_ibp(expr, targets, record)
        seen["ibp_out"] = out
        return out

    def spy_expand(expr):
        out = orig_expand(expr)
        seen["expand_out"] = out
        return out

    def spy_sep(dS, vecs, record=None):
        seen["sep_in"] = dS
        return orig_sep(dS, vecs, record)

    monkeypatch.setattr(eom_mod, "_integrate_by_parts", spy_ibp)
    monkeypatch.setattr(eom_mod, "expand", spy_expand)
    monkeypatch.setattr(eom_mod, "_separate_variations", spy_sep)
    assert ibp_mod  # imported for clarity of what is spied
    s = build()
    tr = EOMTrace()
    eom = compute_eom(s.L, s.dW, s.variables, trace=tr)
    assert tr.after_simplify is seen["ibp_in"]
    assert tr.after_ibp is seen["ibp_out"]
    assert tr.after_expand is seen["expand_out"] is seen["sep_in"]
    for key, (_, eq) in eom.items():
        assert tr.equations[key] is eq


def test_trace_after_failure_holds_the_stages_before_and_empty_dicts_after():
    R = SO3("R")
    Om = R.get_tangent_vector()
    w, A = Vector("w", attr=["Constant"]), Matrix("A", attr=["Constant"])
    L = Dot(w, Vee(Hat(Om) * A))
    tr = EOMTrace()
    with pytest.raises(AlgebraicError):
        compute_eom(L, None, SystemVariables(matrices=[R]), trace=tr)
    assert tr.stage_names()[-1] == "expand"
    assert tr.after_expand is not None
    assert tr.raw_equations == {} and tr.equations == {}


def test_trace_is_empty_when_validation_fails():
    tr = EOMTrace()
    x = Vector("x")
    with pytest.raises(TypeError):
        compute_eom(Dot(x, x), None, SystemVariables(vectors=[x + Vector("y")]), trace=tr)
    assert tr.stages == [] and tr.dL is None


def test_reused_eom_trace_is_reset_or_rejected():
    tr = EOMTrace()
    a, b = point_mass_r3(), rigid_pendulum_pivot()
    compute_eom(a.L, a.dW, a.variables, trace=tr)
    try:
        compute_eom(b.L, b.dW, b.variables, trace=tr)
    except ValueError:
        return  # rejecting a used trace is an acceptable fix
    assert len(tr.stages) == 13
    assert list(tr.equations) == ["\\eta_{R}"]


def test_reused_linearization_trace_is_reset_or_rejected():
    s = spherical_pendulum_s2()
    eom = compute_eom(s.L, s.dW, s.variables)
    tr = LinearizationTrace()
    linearize(eom, s.variables, s.inputs, trace=tr)
    n = len(tr.stages)
    try:
        linearize(eom, s.variables, s.inputs, trace=tr)
    except ValueError:
        return
    assert len(tr.stages) == n


# ---------------------------------------------------------------------------
# Held: zero cost when off
# ---------------------------------------------------------------------------


def test_rule_builds_no_repr_when_tracing_is_off(monkeypatch):
    monkeypatch.setattr(trace_mod, "_TRACE_RULES", False)
    built = []

    class Probe:
        def __repr__(self):
            built.append(1)
            return "p"

        __str__ = __repr__

    assert not trace_mod.trace_enabled()
    trace_mod.rule("x.y", Probe(), Probe())
    assert built == []


def test_stage_recorder_is_none_without_trace_or_debug_logger():
    log = logging.getLogger("geomech")
    old = log.level
    log.setLevel(logging.WARNING)
    try:
        assert trace_mod.stage_recorder(None, "compute_eom") is None
    finally:
        log.setLevel(old)


def test_debug_logger_uses_lazy_formatting():
    """No f-string or str.format in a logger call of the library (lazy %-formatting)."""
    for path in (ROOT / "geomech").rglob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if (
                isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr in ("debug", "info", "warning")
                and getattr(node.func.value, "id", "").endswith("logger")
                and node.args
            ):
                first = node.args[0]
                assert not isinstance(first, ast.JoinedStr), f"{path}:{node.lineno}"


# ---------------------------------------------------------------------------
# Named rules
# ---------------------------------------------------------------------------


def _load_rule_names_test():
    import tests.test_rule_names as trn

    return trn


def test_sync_test_catches_a_named_rule_missing_from_the_docstring(tmp_path, monkeypatch):
    trn = _load_rule_names_test()
    src = (ROOT / "geomech/core/transformations/vector_identities.py").read_text(encoding="utf-8")
    src = src.replace(
        "        case Hat() if expr.expr.is_zero:",
        "        case Transpose() if expr.expr.is_zero:\n"
        '            return rule("vector.transpose_zero", expr, expr.expr)\n\n'
        "        case Hat() if expr.expr.is_zero:",
    )
    monkeypatch.setattr(trn, "_source", lambda name: src)
    with pytest.raises(AssertionError, match="vector.transpose_zero"):
        trn.test_docstring_and_code_list_the_same_rules(
            "geomech.core.transformations.vector_identities", "vector"
        )


def test_sync_test_is_blind_to_a_rewrite_that_bypasses_rule():
    """Documents the limit (finding F6): an unnamed rewrite passes the sync test."""
    trn = _load_rule_names_test()
    src = (ROOT / "geomech/core/transformations/vector_identities.py").read_text(encoding="utf-8")
    src = src.replace(
        "        case Hat() if expr.expr.is_zero:",
        "        case Transpose() if expr.expr.is_zero:\n"
        "            return expr.expr\n\n"
        "        case Hat() if expr.expr.is_zero:",
    )
    old = trn._source
    trn._source = lambda name: src
    try:
        trn.test_docstring_and_code_list_the_same_rules(
            "geomech.core.transformations.vector_identities", "vector"
        )
    finally:
        trn._source = old


def test_vee_over_vadd_is_unreachable():
    """``Vee``'s constructor rejects a vector operand, so the rule's VAdd branch is dead."""
    from geomech import ExpressionMismatchError

    a, b = Vector("a"), Vector("b")
    with pytest.raises(ExpressionMismatchError):
        Vee(a + b)


def test_rule_hypotheses_sampled():
    """Ten rules checked against the hypothesis written next to them."""
    q, p = S2("q"), S2("p")
    w, wp = q.get_tangent_vector(), p.get_tangent_vector()
    a = Vector("a")
    J = Matrix("J", attr=["SymmetricMatrix"])
    A = Matrix("A")
    R = SO3("R")
    from geomech.core.transformations.simplification import _eliminate
    from geomech.core.transformations.vector_identities import vector_rules

    # vector.dot_tangent_orthogonal: only for the TS2's own point
    assert str(vector_rules(Dot(w, q))) == "0"
    assert str(vector_rules(Dot(wp, q))) != "0"
    # vector.s2_tangent_recovery: needs tangency
    assert str(vector_rules(Cross(q, Cross(w, q)))) == str(w)
    assert str(vector_rules(Cross(q, Cross(a, q)))) != str(a)
    # vector.dot_unit_norm: both unit norm and equal
    assert str(vector_rules(Dot(q, q))) == "1"
    assert str(vector_rules(Dot(q, p))) != "1"
    # simplify.transpose_symmetric: symmetric flag only
    assert str(_eliminate(Transpose(J))) == "J"
    assert str(_eliminate(Transpose(A))) != "A"
    # simplify.orthogonal_transpose: orthogonal and equal by name
    assert str(_eliminate(Transpose(R) * R)) == str(geomech.IdentityMatrix)
    assert str(_eliminate(Transpose(A) * A)) != str(geomech.IdentityMatrix)
    # simplify.hat_vee_skew: skew operand only
    assert str(_eliminate(Hat(Vee(Hat(a))))) == str(Hat(a))
    assert str(_eliminate(Hat(Vee(A)))) != "A"
    # vector.cross_self, vector.cross_scaled_self
    assert vector_rules(Cross(a, a)).is_zero
    assert vector_rules(Cross(a, a * Scalar("s"))).is_zero


def test_rule_trace_names_the_rule_that_changed_the_node(rules_log):
    a, w, v = Vector("a"), Vector("w"), Vector("v")
    from geomech.core.transformations.simplification import _eliminate

    _eliminate(Scalar("k") + Scalar("m") * Dot(a, MVMul(Hat(w), v)))
    assert any(line.startswith("simplify.mvmul_hat_to_cross:") for line in rules_log)


def test_fold_rules_are_not_logged_for_a_child_change(rules_log):
    a, w, v = Vector("a"), Vector("w"), Vector("v")
    from geomech.core.transformations.simplification import _eliminate

    _eliminate(Scalar("k") + Scalar("m") * Dot(a, MVMul(Hat(w), v)))
    names = [line.split(":")[0] for line in rules_log]
    assert names == ["simplify.mvmul_hat_to_cross"]


def test_set_trace_does_not_duplicate_lines_under_basic_config():
    code = (
        "import logging; logging.basicConfig(level=logging.DEBUG, format='ROOT %(message)s')\n"
        "import geomech; from geomech import S2, Dot\n"
        "from geomech.core.transformations.vector_identities import vector_rules\n"
        "geomech.set_trace(True); q = S2('q'); vector_rules(Dot(q.get_tangent_vector(), q))\n"
    )
    out = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, cwd=ROOT, check=True
    ).stderr
    lines = [ln for ln in out.splitlines() if "dot_tangent_orthogonal" in ln]
    assert len(lines) == 1, lines


# ---------------------------------------------------------------------------
# check_tree and the conftest fixture
# ---------------------------------------------------------------------------


def _mutated(node, index, child):
    # nodes are frozen (phase 5b, D3): bypass to build the malformed tree
    nodes = list(node.nodes)
    nodes[index] = child
    object.__setattr__(node, "nodes", tuple(nodes))
    return node


@pytest.mark.parametrize(
    "make, what, path",
    [
        (lambda: Add([Add(Scalar("m"), Scalar("g")), Scalar("k")]), "nested directly", "root.0"),
        (
            lambda: _mutated(Scalar("m") + Scalar("k"), 1, Vector("x")),
            "addend has type VECTOR",
            "root.1",
        ),
        (lambda: Cross(Vector("x"), Transpose(Vector("y"))), "row vector", "root.R"),
        (lambda: Hat(Transpose(Vector("y"))), "row vector", "root.x"),
        (lambda: Transpose(Vector("x")) + Vector("y"), "mixes row and column", "root"),
        (lambda: TS2("w"), "no S2 point", "root"),
        (lambda: Scalar("(-1)"), "printed like a number", "root"),
        (
            lambda: Scalar("m") * Dot(Vector("a"), Variation(Variation(Vector("z")).expr)),
            "",
            "",
        ),
    ],
)
def test_check_tree_names_the_invariant_and_path(make, what, path):
    try:
        tree = make()
    except Exception:
        pytest.skip("constructor already rejects this tree")
    if not what:
        check_tree(tree)  # well-formed control
        return
    with pytest.raises(InvariantError) as err:
        check_tree(tree)
    assert what in str(err.value)
    assert f" at {path}:" in str(err.value)


def test_check_tree_flags_variation_over_a_manifold_point_and_deep_path():
    q = S2("q")
    node = TimeDerivative.__new__(TimeDerivative)
    node.__dict__.update(TimeDerivative(Vector("z")).__dict__)
    object.__setattr__(node, "nodes", (q,))  # frozen node: bypass to build the bad tree
    tree = Scalar("m") * Dot(Vector("a"), Cross(Vector("b"), node))
    with pytest.raises(InvariantError, match=r"manifold point 'q' .* at root\.R\.R\.R:"):
        check_tree(tree)


def test_check_tree_has_no_false_positive_on_fuzz_trees():
    """3000 random trees of the engine fuzzer and their five rewrites (18000 trees)."""
    np = pytest.importorskip("numpy")
    sys.path.insert(0, str(ROOT / "scripts"))
    try:
        import engine_fuzz as F
    finally:
        sys.path.pop(0)
    from geomech import expand, pull, simplify
    from geomech.core.transformations.vector_identities import vector_rules

    rng = np.random.default_rng(4)
    gen = F.Gen(rng, F.Leaves())
    checked = 0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for _ in range(3000):
            expr = gen.any(int(rng.integers(1, 5)))
            for fn in (lambda e: e, expand, pull, vector_rules, simplify, full_simplify):
                try:
                    out = fn(expr)
                except Exception:
                    continue
                check_tree(out)
                checked += 1
    assert checked > 15000


def test_conftest_fixture_rejects_a_bad_tree_through_the_patched_pipeline():
    """The autouse fixture really checks: a row in a Dot of the Lagrangian passes the
    pipeline unchecked (it derives a mixed row/column sum) but the wrapper raises."""
    x, a = Vector("x"), Vector("a", attr=["Constant"])
    L = Dot(x.t_diff(), x.t_diff()) * 0.5 + Dot(Transpose(a), x)
    from geomech.dynamics.eom import compute_eom as defining

    assert defining is compute_eom  # the patched name, also in the defining module
    with pytest.raises(InvariantError, match="row vector"):
        compute_eom(L, None, SystemVariables(vectors=[x]))


# ---------------------------------------------------------------------------
# Callable modules (4.6)
# ---------------------------------------------------------------------------


def test_package_simplify_still_calls_through():
    from geomech.core.transformations import simplify

    x = Scalar("x")
    assert str(simplify(x + x)) == str(geomech.simplify(x + x))


def test_package_simplify_is_picklable():
    from geomech.core.transformations import simplify

    pickle.dumps(simplify)


def test_package_simplify_has_the_function_signature():
    import inspect

    from geomech.core.transformations import simplify

    assert isinstance(simplify, types.FunctionType)
    assert str(inspect.signature(simplify)) == "(expr)"


# ---------------------------------------------------------------------------
# Error context (4.2)
# ---------------------------------------------------------------------------


def test_reraises_keep_type_and_cause_and_do_not_repeat_the_engine_text():
    R = SO3("R")
    Om = R.get_tangent_vector()
    w, A = Vector("w", attr=["Constant"]), Matrix("A", attr=["Constant"])
    with pytest.raises(AlgebraicError) as err:
        compute_eom(Dot(w, Vee(Hat(Om) * A)), None, SystemVariables(matrices=[R]))
    e = err.value
    assert type(e.__cause__) is AlgebraicError
    assert str(e).count(str(e.__cause__)) == 1
    assert str(e).count("compute_eom:") == 1


def test_nonlinear_reraise_keeps_type():
    x = Vector("x")
    a = Vector("a", attr=["Constant"])
    with pytest.raises(NonLinearError) as err:
        compute_eom(
            Dot(x.t_diff(), x.t_diff()),
            Dot(x.delta(), a) * Dot(x.delta(), a),
            SystemVariables(vectors=[x]),
        )
    assert isinstance(err.value.__cause__, NonLinearError)


def test_huge_term_is_truncated_in_the_message():
    x = Vector("x")
    L = Dot(x.t_diff(), x.t_diff()) * 0.5
    big = Vector("a0", attr=["Constant"])
    for i in range(1, 60):
        big = big + Vector(f"a{i}", attr=["Constant"])
    L = L + Dot(big, x.t_diff().t_diff())
    with pytest.raises(AlgebraicError) as err:
        compute_eom(L, None, SystemVariables(vectors=[x]))
    assert len(str(err.value)) < 2000


def test_second_order_lagrangian_error_names_the_real_cause():
    x = Vector("x")
    L = Dot(x.t_diff().t_diff(), x.t_diff().t_diff()) * 0.5
    tr = EOMTrace()
    with pytest.raises(AlgebraicError) as err:
        compute_eom(L, None, SystemVariables(vectors=[x]), trace=tr)
    msg = str(err.value)
    stage = re.search(r"look at the stage '([^']*)'", msg)
    assert stage is None or stage.group(1) in tr.stage_names()
    assert "collect could not expose" not in msg


# ---------------------------------------------------------------------------
# Docs: every error text in troubleshooting.md matches the code
# ---------------------------------------------------------------------------


def _doc_messages():
    text = (ROOT / "docs/troubleshooting.md").read_text(encoding="utf-8")
    text = text.split("## Typed errors", 1)[1]
    rows = {}
    for line in text.splitlines():
        if not line.startswith("| ") or line.startswith("|---"):
            continue
        first, _, rest = line[2:].partition(" | ")
        cell = rest.rstrip().rstrip("|").rstrip()
        if not cell.startswith("`"):
            continue
        cell = cell.split("` — ")[0]
        msg = cell[1:].rstrip("`")
        rows[first.strip()] = msg
    return rows


DOC = _doc_messages()


def _msg(func):
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            func()
    except Exception as e:  # noqa: BLE001 - the whole point is the message
        return f"{type(e).__name__}: {e}"
    return "<no error>"


def _syms():
    m, g, k = Scalar("m"), Scalar("g"), Scalar("k")
    x, y = Vector("x"), Vector("y")
    J = Matrix("J")
    return m, g, k, x, y, J


def _q():
    q = S2("q")
    return q, q.get_tangent_vector(), q.get_variation_vector()


def _ke(x):
    return Dot(x.t_diff(), x.t_diff()) * 0.5


def _jerk():
    xs = Vector("xs")
    L = _ke(xs) + Dot(xs.t_diff().t_diff(), xs.t_diff().t_diff()) * 0  # placeholder
    del L
    eom = {
        "\\delta{xs}": (
            xs.delta(),
            TimeDerivative(xs.t_diff().t_diff()),
        )
    }
    to_standard_form(eom, SystemVariables(vectors=[xs]))


def _bilinear():
    x1, x2, m = Vector("x1"), Vector("x2"), Scalar("m")
    eom = {"\\delta{x1}": (x1.delta(), Dot(x1.t_diff().t_diff(), x2.t_diff().t_diff() * m))}
    to_standard_form(eom, SystemVariables(vectors=[x1, x2]))


def _lin_missing_input():
    x, f = Vector("x"), Vector("f")
    eom = compute_eom(_ke(x), Dot(x.delta(), f), SystemVariables(vectors=[x]))
    linearize(eom, SystemVariables(vectors=[x]))


def _matrix_factor():
    R = SO3("R")
    Om = R.get_tangent_vector()
    c1, c2, c3 = (Vector(n, attr=["Constant"]) for n in ("c1", "c2", "c3"))
    J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
    L = Dot(Om, J * Om) * 0.5 + Dot(c1, (Om * Transpose(c2)) * c3)
    compute_eom(L, None, SystemVariables(matrices=[R]))


def _falling_mass():
    """``L = ½ m ẋ·ẋ − m g e₃·x`` with ``g`` constant and the force ``f`` in δW."""
    m, g = Scalar("m", attr=["Constant"]), Scalar("g", attr=["Constant"])
    e3, x, f = Vector("e3", attr=["Constant"]), Vector("x"), Vector("f")
    L = Scalar("0.5", value=0.5) * m * Dot(x.t_diff(), x.t_diff()) - m * g * Dot(e3, x)
    variables = SystemVariables(vectors=[x])
    return compute_eom(L, Dot(x.delta(), f), variables), variables, f


def _input_named_like_a_constant_linearize():
    eom, variables, f = _falling_mass()
    linearize(eom, variables, [f, Scalar("g")])


def _input_named_like_a_constant_standard_form():
    eom, variables, _ = _falling_mass()
    to_standard_form(eom, variables, inputs=[Scalar("g")])


def _input_of_another_class():
    m, x = Scalar("m", attr=["Constant"]), Vector("x")
    variables = SystemVariables(vectors=[x])
    L = Scalar("0.5", value=0.5) * m * Dot(x.t_diff(), x.t_diff())
    eom = compute_eom(L, Dot(x.delta(), Vector("F", attr=["Constant"])), variables)
    linearize(eom, variables, [Scalar("F")])


def _constant_input_standard_form():
    eom, variables, _ = _falling_mass()
    to_standard_form(eom, variables, inputs=[Scalar("g", attr=["Constant"])])


def _duplicate_input_standard_form():
    eom, variables, f = _falling_mass()
    to_standard_form(eom, variables, inputs=[f, f])


def _set_of_inputs_standard_form():
    eom, variables, f = _falling_mass()
    to_standard_form(eom, variables, inputs={f, Vector("M")})


def _input_printed_like_the_velocity():
    eom, variables, f = _falling_mass()
    x = variables.vectors[0]
    to_standard_form(eom, variables, inputs=[f, Vector(str(x.t_diff()))])


def _mass_with_a_rate_param():
    """``L = ½ m ẋ·ẋ + (x·e₁) l̇`` with ``l`` a param: the equation holds ``d/dt(l)``."""
    m, x, f, l = Scalar("m", attr=["Constant"]), Vector("x"), Vector("f"), Scalar("l")
    e1 = Vector("e1", attr=["Constant"])
    L = Scalar("0.5", value=0.5) * m * Dot(x.t_diff(), x.t_diff()) + Dot(x, e1) * l.t_diff()
    variables = SystemVariables(vectors=[x])
    return compute_eom(L, Dot(x.delta(), f), variables, params=[l]), variables, f, l


def _input_printed_like_a_node():
    eom, variables, f, l = _mass_with_a_rate_param()
    to_standard_form(eom, variables, inputs=[f, Scalar(str(l.t_diff()))])


def _param_with_a_rate_as_input():
    eom, variables, f, l = _mass_with_a_rate_param()
    to_standard_form(eom, variables, inputs=[f, l])


def _input_of_another_size():
    eom, variables, _ = _falling_mass()
    to_standard_form(eom, variables, inputs=[Vector("f", size=(6,))])


def _input_named_like_a_numeric_leaf():
    m, x, f, k = Scalar("m", attr=["Constant"]), Vector("x"), Vector("f"), Scalar("k", value=2.0)
    half = Scalar("0.5", value=0.5)
    L = half * m * Dot(x.t_diff(), x.t_diff()) - half * k * Dot(x, x)
    variables = SystemVariables(vectors=[x])
    eom = compute_eom(L, Dot(x.delta(), f), variables)
    to_standard_form(eom, variables, inputs=[f, Scalar("k")])


def _single_input_linearize():
    eom, variables, f = _falling_mass()
    linearize(eom, variables, f)


def _rigid_body_ke():
    R = SO3("R")
    J = Matrix("J", attr=["Constant", "SymmetricMatrix"])
    return Dot(R.get_tangent_vector(), J * R.get_tangent_vector()) * 0.5


def _param_named_like_a_constant():
    m, l, x = Scalar("m", attr=["Constant"]), Scalar("l", attr=["Constant"]), Vector("x")
    L = Scalar("0.5", value=0.5) * m * l * l * Dot(x.t_diff(), x.t_diff())
    compute_eom(L, None, SystemVariables(vectors=[x]), params=[Scalar("l")])


# --- code generation and simulation: one StateLayout for both modules (final review R1-R3)


def _two_inputs():
    """A point mass on R³ with a force ``f`` and a scalar thrust ``T`` along ``e3``."""
    m, e3 = Scalar("m", attr=["Constant"]), Vector("e3", attr=["Constant"])
    x, f, T = Vector("x"), Vector("f"), Scalar("T")
    L = Dot(x.t_diff(), x.t_diff()) * m * 0.5
    dW = Dot(x.delta(), f) + Dot(x.delta(), e3) * T
    V = SystemVariables(vectors=[x])
    return L, dW, V, compute_eom(L, dW, V), f, T


def _lin_inputs_in_another_order():
    from geomech.codegen import StateLayout, generate_module

    _, _, V, eom, f, T = _two_inputs()
    generate_module(linearize(eom, V, [T, f]), StateLayout(V, [f, T], []))


def _lin_variables_in_another_order():
    from geomech.codegen import StateLayout, generate_module

    m = Scalar("m", attr=["Constant"])
    x, y = Vector("x"), Vector("y")
    L = (Dot(x.t_diff(), x.t_diff()) + Dot(y.t_diff(), y.t_diff())) * m * 0.5
    eom = compute_eom(L, None, SystemVariables(vectors=[x, y]))
    lin = linearize(eom, SystemVariables(vectors=[y, x]))
    generate_module(lin, StateLayout(SystemVariables(vectors=[x, y]), [], [m]))


def _lqr_plant_with_another_layout():
    import numpy as np

    from geomech.codegen import StateLayout, generate_module, load_module
    from geomech.simulation import lqr

    _, _, V, eom, f, T = _two_inputs()
    params = [Vector("e3", attr=["Constant"]), Scalar("m", attr=["Constant"])]
    lin = load_module(generate_module(linearize(eom, V, [f, T]), StateLayout(V, [f, T], params)))
    sf = to_standard_form(eom, V, [T, f])
    plant = load_module(generate_module(sf, StateLayout(V, [T, f], params)))
    lqr(lin, np.zeros(6), np.zeros(4), {"m": 1.0, "e3": [0, 0, 1]}, 1.0, 1.0, module=plant)


def _variation_vector_as_param():
    from geomech.codegen import StateLayout

    q = S2("q")
    StateLayout(SystemVariables(vectors=[q]), [], [q.get_variation_vector()])


def _collect_params_of_a_string():
    from geomech.codegen import collect_params

    L, _, V, _, f, T = _two_inputs()
    collect_params([L, "dW"], V, [f, T])


TRIGGERS = {
    "`Dot(m, x)`": lambda: Dot(_syms()[0], _syms()[3]),
    "`m + x`": lambda: _syms()[0] + _syms()[3],
    "`x - m`": lambda: _syms()[3] - _syms()[0],
    "`x * J`": lambda: _syms()[3] * _syms()[5],
    "`J * x.T()`": lambda: _syms()[5] * _syms()[3].T(),
    "`x.T() * y.T()`": lambda: _syms()[3].T() * _syms()[4].T(),
    "`x * y`": lambda: _syms()[3] * _syms()[4],
    "`Hat(m)`": lambda: Hat(_syms()[0]),
    "`Vee(x)`": lambda: Vee(_syms()[3]),
    '`Vector("a", size=(2,)) + x`': lambda: Vector("a", size=(2,)) + _syms()[3],
    '`m * "x"`': lambda: _syms()[0] * "x",
    "`m * True`": lambda: _syms()[0] * True,
    "`x / m`": lambda: _syms()[3] / _syms()[0],
    "`1 / m`": lambda: 1 / _syms()[0],
    "`x / 0`": lambda: _syms()[3] / 0,
    "`Number([1])`": lambda: Number([1]),
    '`Scalar("m", attr=["constant"])`': lambda: Scalar("m", attr=["constant"]),
    '`Scalar("m", attr="Constant")`': lambda: Scalar("m", attr="Constant"),
    "`Transpose()`": lambda: Transpose(),
    "`Add()`": lambda: Add(),
    "`x.delta().delta()`": lambda: _syms()[3].delta().delta(),
    "`xi.delta()`": lambda: _q()[2].delta(),
    "`TimeDerivative(q)`": lambda: TimeDerivative(S2("q")),
    "`Variation(R)`": lambda: Variation(SO3("R")),
    "`m.t_integrate()`, `m` constant": lambda: Scalar("m", attr=["Constant"]).t_integrate(),
    "`SystemVariables(scalars=[x])`, `x` a Vector": lambda: compute_eom(
        Scalar("k"), None, SystemVariables(scalars=[Vector("x")])
    ),
    "`SystemVariables(vectors=x)`": lambda: SystemVariables(vectors=Vector("x")),
    "a `Constant` variable": lambda: compute_eom(
        _ke(Vector("x")), None, SystemVariables(vectors=[Vector("e3", attr=["Constant"])])
    ),
    "`ω` listed as a variable": lambda: compute_eom(
        _ke(Vector("x")), None, SystemVariables(vectors=[_q()[1]])
    ),
    "`ẋ` listed": lambda: compute_eom(
        _ke(Vector("x")), None, SystemVariables(vectors=[Vector("x").t_diff()])
    ),
    "`x + y` listed": lambda: compute_eom(
        _ke(Vector("x")), None, SystemVariables(vectors=[Vector("x") + Vector("y")])
    ),
    "a plain `Matrix` in `matrices`": lambda: compute_eom(
        _ke(Vector("x")), None, SystemVariables(matrices=[Matrix("A")])
    ),
    "a name listed twice": lambda: compute_eom(
        _ke(Vector("x")), None, SystemVariables(vectors=[Vector("x"), Vector("x")])
    ),
    "a variable named like `ω_q`": lambda: compute_eom(
        _ke(Vector("x")),
        None,
        SystemVariables(vectors=[Vector("\\omega_{q}"), S2("q")]),
    ),
    "`SystemVariables()`": lambda: compute_eom(
        Scalar("k", attr=["Constant"]), None, SystemVariables()
    ),
    "a vector Lagrangian": lambda: compute_eom(
        Vector("x").t_diff(), None, SystemVariables(vectors=[Vector("x")])
    ),
    "`inf_work=True`": lambda: compute_eom(
        _ke(Vector("x")), True, SystemVariables(vectors=[Vector("x")])
    ),
    "`y` used but not listed (or `mm` not `Constant`)": lambda: compute_eom(
        _ke(Vector("x")) + Dot(Vector("x"), Vector("y")),
        None,
        SystemVariables(vectors=[Vector("x")]),
    ),
    "`Dot(x.t_diff().delta(), f)` in δW": lambda: compute_eom(
        _ke(Vector("x")),
        Dot(Vector("x").t_diff().delta(), Vector("f")),
        SystemVariables(vectors=[Vector("x")]),
    ),
    "a hand-built `Variation(x.t_diff())`": lambda: compute_eom(
        _ke(Vector("x")),
        Dot(Variation(Vector("x").t_diff()), Vector("f")),
        SystemVariables(vectors=[Vector("x")]),
    ),
    "`params=l`": lambda: compute_eom(
        _ke(Vector("x")), None, SystemVariables(vectors=[Vector("x")]), params=Scalar("l")
    ),
    "a `Constant` param": lambda: compute_eom(
        _ke(Vector("x")),
        None,
        SystemVariables(vectors=[Vector("x")]),
        params=[Scalar("m", attr=["Constant"])],
    ),
    "a param that is a variable": lambda: compute_eom(
        _ke(Vector("x")), None, SystemVariables(vectors=[Vector("x")]), params=[Vector("x")]
    ),
    "a velocity inside a matrix factor, `a·((Ω cᵀ) b)`": _matrix_factor,
    "a term quadratic in a variation": lambda: compute_eom(
        _ke(Vector("x")),
        VVMul_q(),
        SystemVariables(vectors=[Vector("x")]),
    ),
    "a term bilinear in two accelerations": _bilinear,
    "an input missing from `inputs` (or a symbol missing `Constant`)": _lin_missing_input,
    "an input that is a variable": lambda: linearize(
        compute_eom(_ke(Vector("x")), None, SystemVariables(vectors=[Vector("x")])),
        SystemVariables(vectors=[Vector("x")]),
        [Vector("x")],
    ),
    "`ω` as an input": lambda: linearize(
        compute_eom(Dot(_q()[1], _q()[1]) * 0.5, None, SystemVariables(vectors=[_q()[0]])),
        SystemVariables(vectors=[_q()[0]]),
        [_q()[1]],
    ),
    "not the `compute_eom` dict": lambda: linearize([], SystemVariables(vectors=[Vector("x")])),
    "other `SystemVariables` than for `compute_eom`": lambda: linearize(
        compute_eom(_ke(Vector("x")), None, SystemVariables(vectors=[Vector("x")])),
        SystemVariables(vectors=[_q()[0]]),
    ),
    "`extract_linear_coeff(a*a, a)`": lambda: extract_linear_coeff(
        Scalar("a") * Scalar("a"), Scalar("a")
    ),
    "`extract_linear_coeff(x * x.T(), x)`": lambda: extract_linear_coeff(
        Vector("x") * Vector("x").T(), Vector("x")
    ),
    "`extract_linear_coeff(Dot(ẋ, f), x)`": lambda: extract_linear_coeff(
        Dot(Vector("x").t_diff(), Vector("f")), Vector("x")
    ),
    "`extract_linear_coeff(Hat(x), x)`": lambda: extract_linear_coeff(
        Hat(Vector("x")), Vector("x")
    ),
    "`collect(x, x)`": lambda: geomech.collect(Vector("x"), Vector("x")),
    "`Add([Add(m, g), k])` (a list is not flattened)": lambda: check_tree(
        Add([Add(Scalar("m"), Scalar("g")), Scalar("k")])
    ),
    "`Dot(x.T(), y)`": lambda: check_tree(Dot(Vector("x").T(), Vector("y"))),
    "`x.T() + y`": lambda: check_tree(Vector("x").T() + Vector("y")),
    '`TS2("w")`': lambda: check_tree(TS2("w")),
    '`Scalar("(-1)")`': lambda: check_tree(Scalar("(-1)")),
    # an input or param printed like a different leaf of the equations
    'an input named like a constant, `linearize(eom, V, [f, Scalar("g")])`, `g` `Constant`'
    " in `L`": _input_named_like_a_constant_linearize,
    'an input named like a constant of `L` (`inputs=[Scalar("g")]`, `g` `Constant` in `L`)': (
        _input_named_like_a_constant_standard_form
    ),
    'an input of another class than its leaf (`Scalar("F")`, the equations hold the constant'
    ' `Vector("F")`)': _input_of_another_class,
    'a param named like a constant of `L` (`params=[Scalar("l")]`, `l` `Constant` in `L`)': (
        _param_named_like_a_constant
    ),
    # the input checks shared by to_standard_form and linearize (check_inputs)
    "a `Constant` input (`inputs=[g]`, `g` `Constant` in `L`)": _constant_input_standard_form,
    "an input listed twice (`inputs=[f, f]`)": _duplicate_input_standard_form,
    "a set (or dict) of inputs, `inputs={f, M}`": _set_of_inputs_standard_form,
    "an input printed like a velocity (`Vector(str(x.t_diff()))`)": (
        _input_printed_like_the_velocity
    ),
    "an input printed like an operation node of the equations (`Scalar(str(l.t_diff()))`,"
    " `l` a param)": _input_printed_like_a_node,
    "a param of `compute_eom` whose rate is in the equations given as an input"
    " (`inputs=[f, l]`)": _param_with_a_rate_as_input,
    'an input of another size than its leaf (`Vector("f", size=(6,))`)': _input_of_another_size,
    'an input named like a named numeric leaf (`Scalar("k")`, `L` holds'
    ' `Scalar("k", value=2.0)`)': _input_named_like_a_numeric_leaf,
    "not the `compute_eom` dict, `to_standard_form([], V)`": lambda: to_standard_form(
        [], SystemVariables(vectors=[Vector("x")])
    ),
    "`ξ`, `Ω` or `η` as an input": lambda: linearize(
        compute_eom(_rigid_body_ke(), None, SystemVariables(matrices=[SO3("R")])),
        SystemVariables(matrices=[SO3("R")]),
        [SO3("R").get_variation_vector()],
    ),
    "a plain `Vector` printed like `Ω` as an input": lambda: linearize(
        compute_eom(_rigid_body_ke(), None, SystemVariables(matrices=[SO3("R")])),
        SystemVariables(matrices=[SO3("R")]),
        [Vector(str(SO3("R").get_tangent_vector()))],
    ),
    "an input that is not a plain leaf (`ẋ`, `F + u`, a `Matrix`, an `S2`/`SO3` point)": (
        lambda: linearize(
            compute_eom(_ke(Vector("x")), None, SystemVariables(vectors=[Vector("x")])),
            SystemVariables(vectors=[Vector("x")]),
            [Vector("x").t_diff()],
        )
    ),
    "one input instead of a list, `linearize(eom, V, f)`": _single_input_linearize,
    # code generation and simulation
    "a linearization with other inputs than the layout (`linearize(eom, V, [T, f])`,"
    " `StateLayout(V, [f, T])`)": _lin_inputs_in_another_order,
    "a linearization with another variable order than the layout"
    " (`SystemVariables(vectors=[y, x])` for `linearize`, `[x, y]` for the layout)": (
        _lin_variables_in_another_order
    ),
    "`lqr(lin_module, ..., module=plant)` with modules of two layouts": (
        _lqr_plant_with_another_layout
    ),
    "a variation vector as a parameter (`collect_params` leaves them out)": (
        _variation_vector_as_param
    ),
    "something without symbols in `collect_params`": _collect_params_of_a_string,
}


def VVMul_q():  # noqa: N802 - the δx·δx term of the quadratic-variation row
    x = Vector("x")
    return x.delta().T() * x.delta()


@pytest.mark.parametrize("label", sorted(TRIGGERS))
def test_troubleshooting_message_matches_the_code(label):
    assert label in DOC, f"row {label!r} not found in docs/troubleshooting.md"
    assert _msg(TRIGGERS[label]) == DOC[label]


def test_troubleshooting_rows_without_a_trigger_are_few():
    """Rows not reproduced above: the jerk/IBP rows (a section reference, checked by
    test_error_context), the params row of linearize and the full_simplify warning (its
    message ends in a repr placeholder)."""
    missing = sorted(set(DOC) - set(TRIGGERS))
    assert len(missing) <= 4, missing


def test_troubleshooting_full_simplify_warning_prefix():
    from geomech import Cross as C

    a, b = Vector("a"), Vector("b")
    expr = Dot(a, C(b, a + b)) + Dot(a, a) * Scalar("m") * (Scalar("k") + Scalar("k"))
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        full_simplify(expr, max_iter=1)
    texts = [str(w.message) for w in caught if w.category is RuntimeWarning]
    doc = DOC["`full_simplify(expr, max_iter=1)` at the cap"]
    prefix = doc.split("RuntimeWarning: ", 1)[1].split("<repr>")[0]
    assert texts and texts[0].startswith(prefix)


def test_jerk_message_matches_troubleshooting_block():
    text = (ROOT / "docs/troubleshooting.md").read_text(encoding="utf-8")
    block = text.split("**2. A jerk**", 1)[1].split("```", 2)[1]
    expected = " ".join(block.split())
    got = " ".join(_msg(_jerk).split())
    assert got == expected


def test_ibp_survivor_message_matches_troubleshooting_block():
    text = (ROOT / "docs/troubleshooting.md").read_text(encoding="utf-8")
    block = text.split("**1. A `d/dt` of a variation vector", 1)[1].split("```", 2)[1]
    R = SO3("R")
    Om = R.get_tangent_vector()
    w, A = Vector("w", attr=["Constant"]), Matrix("A", attr=["Constant"])
    got = _msg(lambda: compute_eom(Dot(w, Vee(Hat(Om) * A)), None, SystemVariables(matrices=[R])))
    # the doc wraps lines inside identifiers; compare with all whitespace removed
    assert "".join(got.split()) == "".join(block.split())


def test_public_integrate_by_parts_logs_its_stages(caplog):
    from geomech import integrate_by_parts

    x = Vector("x")
    with caplog.at_level(logging.DEBUG, logger="geomech"):
        integrate_by_parts(Dot(x.delta().t_diff(), x.t_diff()), [x.delta().t_diff()])
    assert any("ibp" in r.getMessage() for r in caplog.records)
