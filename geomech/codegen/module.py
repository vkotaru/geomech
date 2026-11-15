"""``generate_module``: a self-contained Python module for a derived system.

Standard form (``to_standard_form``)::

    LAYOUT                                  # StateLayout.to_dict()
    mass(state, p)          -> (n, n_a)     # block matrix of M[aⱼ] per equation
    bias(state, p)          -> (n,)         # stacked f
    input_matrix(state, p)  -> (n, n_u)     # block matrix of G[uₖ]
    kinematics(state)       -> (n_config,)  # ṡ, ẋ, ω e × q, ω × q, R hat(Ω) (row-major);
                                            # kinematics(state, p) with an S¹ point

so that ``mass @ a + bias + input_matrix @ u = 0``. ``compute_eom`` output gives
``residual(state, a, u, p)`` (the stacked equations ``E``) and ``kinematics``; a
``LinearizedSystem`` gives ``M_lin``, ``A_lin``, ``B_lin`` of ``(state, a, u, p)`` with the
row and column order of ``LinearizedSystem.rows/dstates/states/inputs`` (recorded in
``LINEARIZATION``); they must follow the layout's variables and inputs, see
``generate_module``. With ``lagrangian=L`` the module also has ``energy(state, p)``.
Equations are ordered like the configuration variables of the layout.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable

from geomech.codegen.emitter import Code, Emitter
from geomech.codegen.energy import energy_source
from geomech.codegen.layout import StateLayout
from geomech.codegen.source import function_source, get_target, namespace, pyliteral
from geomech.core.base.expressions import MANIFOLD_POINTS, SO3, Expr
from geomech.core.base.types import ExprType
from geomech.dynamics.standard_form import StandardFormEquation
from geomech.dynamics.vbl import LinearizedSystem, perturbation_variables

__all__ = ["generate_module", "load_module", "module_source"]

_STATE_KINDS = frozenset({"configuration", "velocity", "param"})
_FULL = ("state", "a", "u", "p")


def _variation_key(var) -> str:
    if isinstance(var, MANIFOLD_POINTS):
        return str(var.get_variation_vector())
    return str(var.delta())


def _size(expr_type) -> int:
    return 1 if expr_type == ExprType.SCALAR else 3


def _spans(sizes):
    """``[(start, size)]`` of consecutive blocks of the given sizes."""
    out, pos = [], 0
    for n in sizes:
        out.append((pos, n))
        pos += n
    return out


def _check_block(expr, nr, nc, where):
    expected = {
        (1, 1): ExprType.SCALAR,
        (1, 3): ExprType.VECTOR,
        (3, 1): ExprType.VECTOR,
        (3, 3): ExprType.MATRIX,
    }[(nr, nc)]
    if expr.type != expected:
        raise ValueError(
            f"generate_module: block {where} is a {expr.type.name} '{expr}', expected a"
            f" {expected.name} for a {nr}×{nc} block"
        )


# ---------------------------------------------------------------------------
# Function builders
# ---------------------------------------------------------------------------


def _block_function(name, args, doc, layout, tgt, grid, row_spans, col_spans, kinds):
    """A function returning the block matrix of ``grid`` (``{(i, j): expr}``)."""
    keys = sorted(grid)
    exprs = [grid[k] for k in keys]

    def finish(em: Emitter, codes: list[Code]) -> list[str]:
        return tgt.block_matrix(row_spans, col_spans, dict(zip(keys, codes)))

    return function_source(name, args, doc, layout, tgt, exprs, finish, kinds=kinds, what=name)


def _stack_function(name, args, doc, layout, tgt, parts, n, kinds, flatten=()):
    """A function returning the stacked vector of ``parts`` (``[(start, size, expr)]``);
    the parts whose index is in *flatten* are 3×3 matrices written row-major."""

    def finish(em: Emitter, codes: list[Code]) -> list[str]:
        items = []
        for k, ((start, size, _), code) in enumerate(zip(parts, codes)):
            if k in flatten:
                code = em.assign(tgt.flatten_matrix(code))
            items.append((start, size, code))
        return tgt.stack(n, items)

    exprs = [e for _, _, e in parts]
    return function_source(name, args, doc, layout, tgt, exprs, finish, kinds=kinds, what=name)


def _kinematics(layout: StateLayout, tgt) -> str:
    """``kinematics(state)``, or ``kinematics(state, p)`` when the layout has an S¹ point:
    its rate ``ω e × q`` needs the axis ``e``, a parameter."""
    parts, flatten = [], set()
    for k, (var, slot) in enumerate(zip(layout.config_variables, layout.configuration)):
        if isinstance(var, SO3):
            flatten.add(k)
        parts.append((slot.start, slot.size, var.t_diff()))
    doc = "Configuration rates: ds/dt, dx/dt, omega x q, R hat(Omega) (row-major)."
    kinds = {"configuration", "velocity"}
    args: tuple[str, ...] = ("state",)
    if any(slot.manifold == "S1" for slot in layout.configuration):
        doc = (
            "Configuration rates: ds/dt, dx/dt, omega e x q (S1), omega x q (S2), R hat(Omega)"
            " (row-major)."
        )
        kinds.add("param")
        args = ("state", "p")
    return _stack_function(
        "kinematics", args, doc, layout, tgt, parts, layout.n_config, kinds, flatten
    )


# ---------------------------------------------------------------------------
# Systems
# ---------------------------------------------------------------------------


def _eom_rows(system: dict, layout: StateLayout) -> list[tuple[str, Any]]:
    """``[(key, value)]`` of the equations in layout order (one per configuration var)."""
    keys = [_variation_key(v) for v in layout.config_variables]
    if len(system) != len(keys) or set(system) != set(keys):
        raise ValueError(
            f"generate_module: the equations {sorted(system)} do not match the variation"
            f" vectors {sorted(keys)} of the layout's variables"
        )
    return [(k, system[k]) for k in keys]


def _row_spans(layout: StateLayout):
    """One block row per equation, of the size of its variation: 1 for a scalar and for
    an S¹ point (γ), 3 otherwise."""
    return _spans(_size(_variation_type(v)) for v in layout.config_variables)


def _variation_type(var):
    if isinstance(var, MANIFOLD_POINTS):
        return var.get_variation_vector().type
    return var.type


def _standard_form(sf: dict, layout: StateLayout, tgt) -> list[str]:
    rows = _eom_rows(sf, layout)
    row_spans = _row_spans(layout)
    acc = {s.key: j for j, s in enumerate(layout.acceleration)}
    inp = {s.key: j for j, s in enumerate(layout.input_slots)}
    acc_spans = [(s.start, s.size) for s in layout.acceleration]
    in_spans = [(s.start, s.size) for s in layout.input_slots]
    M, G, f = {}, {}, []
    for i, (key, eq) in enumerate(rows):
        M.update(_grid(eq.M, acc, i, row_spans, acc_spans, f"M[{key}]", "acceleration"))
        G.update(_grid(eq.G, inp, i, row_spans, in_spans, f"G[{key}]", "input"))
        if eq.f is not None and not eq.f.is_zero:
            _check_block(eq.f, row_spans[i][1], 1, f"f[{key}]")
            f.append((row_spans[i][0], row_spans[i][1], eq.f))
    n = sum(size for _, size in row_spans)
    args = ("state", "p")
    return [
        _block_function(
            "mass", args, _MASS_DOC, layout, tgt, M, row_spans, acc_spans, _STATE_KINDS
        ),
        _stack_function("bias", args, _BIAS_DOC, layout, tgt, f, n, _STATE_KINDS),
        _block_function(
            "input_matrix", args, _INPUT_DOC, layout, tgt, G, row_spans, in_spans, _STATE_KINDS
        ),
        _kinematics(layout, tgt),
    ]


_MASS_DOC = "Mass matrix M(state): blocks M[a_j] per equation (M = -physical inertia)."
_BIAS_DOC = "Bias f(state): the terms of each equation free of accelerations and inputs."
_INPUT_DOC = "Input matrix G(state): blocks G[u_k] per equation."


def _grid(blocks: dict, columns: dict, i, row_spans, col_spans, where, what) -> dict:
    """``{(i, j): block}`` of one block row, skipping zero blocks."""
    out = {}
    for name, block in blocks.items():
        if name not in columns:
            raise ValueError(
                f"generate_module: {where} has a block for '{name}', which is not an"
                f" {what} of the layout"
            )
        if block.is_zero:
            continue
        j = columns[name]
        _check_block(block, row_spans[i][1], col_spans[j][1], f"{where}[{name}]")
        out[(i, j)] = block
    return out


def _eom(eom: dict, layout: StateLayout, tgt) -> list[str]:
    rows = _eom_rows(eom, layout)
    parts = []
    for (start, size), (key, (_, eq)) in zip(_row_spans(layout), rows):
        _check_block(eq, size, 1, f"E[{key}]")
        parts.append((start, size, eq))
    n = sum(size for _, size, _ in parts)
    doc = "Stacked equations E(state, a, u) of compute_eom (E = 0 on a trajectory)."
    residual = _stack_function("residual", _FULL, doc, layout, tgt, parts, n, _ALL)
    return [residual, _kinematics(layout, tgt)]


_ALL = frozenset({"configuration", "velocity", "acceleration", "input", "param"})


_SAME_LAYOUT = (
    "; build the StateLayout with the SystemVariables and the inputs (same order) passed to"
    " linearize and to_standard_form, and generate both modules from that one StateLayout"
)


def _names(symbols) -> str:
    return "[" + ", ".join(f"'{s}'" for s in symbols) + "]"


def _check_linearization_order(lin: LinearizedSystem, layout: StateLayout) -> None:
    """Raise ``ValueError`` unless the rows, states and inputs of *lin* follow *layout*.

    ``B_lin`` has its columns in ``lin.inputs`` order and every block its rows and states
    in the ``SystemVariables`` order given to ``linearize``, while ``u`` and ``state``
    follow the layout: a different order would silently pair a column with the wrong
    input or state (a controller designed on it acts on the wrong channels)."""
    cfg_vars, cfg_perts, vel_perts = perturbation_variables(layout.variables)
    lin_vars = [r[len("kin:") :] for r in lin.rows if r.startswith("kin:")]
    lay_vars = [str(v) for v in cfg_vars]
    rows = ["kin:" + k for k in lay_vars] + [_variation_key(v) for v in cfg_vars]
    states = [str(s) for s in cfg_perts + vel_perts]
    if list(lin.rows) != rows or [str(s) for s in lin.states] != states:
        raise ValueError(
            f"generate_module: the linearization's configuration variables {_names(lin_vars)}"
            f" (rows {_names(lin.rows)}) are not the layout's variables {_names(lay_vars)}"
            f" in that order (rows {_names(rows)})" + _SAME_LAYOUT
        )
    lin_inputs = [str(u) for u in lin.inputs]
    expected = [str(u.delta()) for u in layout.inputs]
    if lin_inputs != expected:
        raise ValueError(
            f"generate_module: the linearization's inputs {_names(lin_inputs)} are not the"
            f" variations of the layout's inputs {_names(layout.inputs)} in that order"
            f" (expected {_names(expected)})" + _SAME_LAYOUT
        )


def _linearized(lin: LinearizedSystem, layout: StateLayout, tgt) -> list[str]:
    _check_linearization_order(lin, layout)
    row_spans = _spans(_size(lin.row_types[r]) for r in lin.rows)
    out = []
    for name, blocks, cols, doc in (
        ("M_lin", lin.M, lin.dstates, "M of M d/dt(dx) = A dx + B du (columns: dstates)."),
        ("A_lin", lin.A, lin.states, "A of M d/dt(dx) = A dx + B du (columns: states)."),
        ("B_lin", lin.B, lin.inputs, "B of M d/dt(dx) = A dx + B du (columns: inputs)."),
    ):
        col_spans = _spans(_size(lin.col_types[str(c)]) for c in cols)
        columns = {str(c): j for j, c in enumerate(cols)}
        grid = {}
        for i, row in enumerate(lin.rows):
            grid.update(
                _grid(blocks[row], columns, i, row_spans, col_spans, f"{name}[{row}]", "column")
            )
        out.append(
            _block_function(name, _FULL, doc, layout, tgt, grid, row_spans, col_spans, _ALL)
        )
    return out + [_kinematics(layout, tgt)]


def _linearization_data(lin: LinearizedSystem) -> dict:
    def entries(symbols, types):
        spans = _spans(_size(types[str(s)]) for s in symbols)
        return [{"symbol": str(s), "slice": [a, a + n]} for s, (a, n) in zip(symbols, spans)]

    rows = entries(lin.rows, lin.row_types)
    return {
        "rows": [{"label": r.pop("symbol"), **r} for r in rows],
        "dstates": entries(lin.dstates, lin.col_types),
        "states": entries(lin.states, lin.col_types),
        "inputs": entries(lin.inputs, lin.col_types),
    }


# ---------------------------------------------------------------------------
# Module
# ---------------------------------------------------------------------------


def _kind(system) -> str:
    if isinstance(system, LinearizedSystem):
        return "linearized"
    if isinstance(system, dict) and system:
        values = list(system.values())
        if all(isinstance(v, StandardFormEquation) for v in values):
            return "standard_form"
        if all(isinstance(v, tuple) and len(v) == 2 for v in values):
            return "eom"
    raise TypeError(
        "generate_module: system must be the dict of compute_eom, the dict of"
        f" to_standard_form or a LinearizedSystem, got {type(system).__name__}"
    )


_HEADER = '''"""Generated by geomech.codegen ({target} target); do not edit.

From {what}.
{description}State layout: LAYOUT (geomech.codegen.StateLayout.to_dict()).
"""'''

_WHAT = {
    "standard_form": "to_standard_form: mass @ a + bias + input_matrix @ u = 0",
    "eom": "compute_eom: residual(state, a, u, p) = 0",
    "linearized": "linearize: M_lin d/dt(dx) = A_lin dx + B_lin du",
}

_LAYOUT_FN = '''def layout():
    """A copy of LAYOUT."""
    return copy.deepcopy(LAYOUT)'''


def module_source(
    system: dict | LinearizedSystem,
    layout: StateLayout,
    target: str = "numpy",
    *,
    lagrangian: Expr | None = None,
    description: str = "",
) -> str:
    """The source of the module ``generate_module`` writes (see there)."""
    kind = _kind(system)
    tgt = get_target(target)
    builders: dict[str, Callable] = {
        "standard_form": _standard_form,
        "eom": _eom,
        "linearized": _linearized,
    }
    functions = builders[kind](system, layout, tgt)
    if lagrangian is not None:
        functions.append(energy_source(lagrangian, layout, target))
    desc = "\n" + description.strip() + "\n\n" if description.strip() else "\n"
    head = _HEADER.format(target=tgt.name, what=_WHAT[kind], description=desc)
    imports = "import copy\n\n" + "\n".join(tgt.imports)
    parts = [head + "\n\n" + imports + "\n\n" + "LAYOUT = " + pyliteral(layout.to_dict())]
    if isinstance(system, LinearizedSystem):
        parts.append("LINEARIZATION = " + pyliteral(_linearization_data(system)))
    parts.append("_PARAMS = " + pyliteral(tuple(layout.param_names)))
    parts += [h.strip() for h in tgt.helpers.strip().split("\n\n\n")]
    parts.append(_LAYOUT_FN)
    parts += [f.rstrip("\n") for f in functions]
    return "\n\n\n".join(parts) + "\n"


def generate_module(
    system: dict | LinearizedSystem,
    layout: StateLayout,
    path: str | Path | None = None,
    target: str = "numpy",
    *,
    lagrangian: Expr | None = None,
    description: str = "",
) -> str:
    """Write a Python module computing the blocks of *system* and return its source.

    Args:
        system: the dict returned by ``compute_eom``, the ``dict[str,
            StandardFormEquation]`` of ``to_standard_form``, or a ``LinearizedSystem``.
        layout: the ``StateLayout`` of the system (variables, inputs, params). Generate
            the standard-form and the linearization module of one system from the same
            layout, built with the ``SystemVariables`` and ``inputs`` (same order) passed
            to ``to_standard_form`` and ``linearize``.
        path: the ``.py`` file to write (``None``: only return the source).
        target: ``"numpy"`` (default) or ``"casadi"`` (the module imports ``casadi``).
        lagrangian: if given, the module also defines ``energy(state, p)``, the total
            energy ``Σ (∂L/∂v)·v − L`` (``energy_function``).
        description: text for the module docstring.

    The module needs only ``numpy`` (and ``casadi`` for that target), not ``geomech``.
    See the module docstring of ``geomech.codegen.module`` for the functions it defines;
    signs and shapes follow ``docs/conventions.md`` (``M`` = −physical inertia, a scalar
    equation with a vector acceleration contributes a 1×3 block). An S² equation holds
    only on the tangent plane, and its mass block may be singular (``Hat(q) Hat(q) m l²``)
    or already full rank (``−m l² I``); the module does not solve for the accelerations:
    ``geomech.simulation.accelerations`` does, with the tangency constraint ``q·ω̇ = 0``.

    Raises:
        TypeError: *system* is none of the three kinds.
        ValueError: equations that do not match the layout's variables, a block of the
            wrong shape, or a ``LinearizedSystem`` whose variables (rows and states) or
            inputs are not those of the layout in the same order (``B_lin`` would
            otherwise pair its columns with the wrong entries of ``u``).
    """
    src = module_source(system, layout, target, lagrangian=lagrangian, description=description)
    if path is not None:
        Path(path).write_text(src)
    return src


def load_module(source_or_path, name: str = "geomech_generated"):
    """Execute generated source (or a generated file) and return it as a module object."""
    import types

    text = str(source_or_path)
    if isinstance(source_or_path, Path) or "\n" not in text:
        text = Path(text).read_text()
    module = types.ModuleType(name)
    target = "casadi" if "import casadi as ca" in text else "numpy"
    if target == "casadi":
        namespace("casadi")  # raises ImportError with an install hint when absent
    exec(compile(text, f"<{name}>", "exec"), module.__dict__)
    return module
