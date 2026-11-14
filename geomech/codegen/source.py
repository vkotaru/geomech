"""Generated functions: one expression (``to_source``, ``lambdify``) or several blocks.

``function_source`` is the one code path every generated function goes through: it runs
an ``Emitter`` over the expressions, binds the slots they use and lets a *finish* callback
write the statements that assemble and return the result.
"""

from __future__ import annotations

from typing import Any, Callable, Sequence

from geomech.codegen.casadi_target import CasadiTarget, require_casadi
from geomech.codegen.emitter import ALL_KINDS, Code, Emitter, pystr
from geomech.codegen.layout import StateLayout
from geomech.codegen.numpy_target import NumpyTarget
from geomech.core.base.expressions import Expr

__all__ = ["to_source", "lambdify", "get_target", "function_source", "pyliteral"]

_INDENT = "    "


def get_target(target: str):
    """The target object for ``"numpy"`` or ``"casadi"``."""
    if target == "numpy":
        return NumpyTarget()
    if target == "casadi":
        return CasadiTarget()
    raise ValueError(f"codegen: unknown target {target!r}; use 'numpy' or 'casadi'")


def pyliteral(obj: Any, indent: str = "") -> str:
    """Python source of plain data (dict, list, tuple, str, int, float), in the layout
    ``ruff format`` keeps: collections expanded one item per line with a trailing comma,
    short lists of numbers or strings inline."""
    inner = indent + _INDENT
    if isinstance(obj, dict):
        if not obj:
            return "{}"
        items = [f"{inner}{pystr(str(k))}: {pyliteral(v, inner)}," for k, v in obj.items()]
        return "{\n" + "\n".join(items) + f"\n{indent}}}"
    if isinstance(obj, (list, tuple)):
        opener, closer = ("[", "]") if isinstance(obj, list) else ("(", ")")
        if not obj:
            return opener + closer
        flat = [pyliteral(v) for v in obj]
        if all(not isinstance(v, (dict, list, tuple)) for v in obj):
            text = ", ".join(flat) + ("," if isinstance(obj, tuple) and len(obj) == 1 else "")
            if len(text) <= 40:
                return opener + text + closer
        items = [f"{inner}{pyliteral(v, inner)}," for v in obj]
        return opener + "\n" + "\n".join(items) + f"\n{indent}{closer}"
    if isinstance(obj, str):
        return pystr(obj)
    if isinstance(obj, bool) or obj is None:
        return repr(obj)
    if isinstance(obj, (int, float)):
        return repr(obj)
    raise TypeError(f"pyliteral: cannot write {type(obj).__name__}")


def _docstring(doc: str) -> list[str]:
    lines = doc.strip().splitlines()
    if len(lines) == 1:
        return [f'"""{lines[0]}"""']
    return [f'"""{lines[0]}', *lines[1:], '"""']


def function_source(
    name: str,
    args: Sequence[str],
    doc: str,
    layout: StateLayout,
    target,
    exprs: Sequence[Expr],
    finish: Callable[[Emitter, list[Code]], list[str]],
    *,
    kinds=ALL_KINDS,
    what: str = "expression",
    params_const: str = "_PARAMS",
) -> str:
    """Source of ``def name(args)`` computing *exprs* and returning what *finish* builds.

    *finish* receives the emitter and the code of each expression and returns the final
    statements (assembly and ``return``). A slot whose array is not among *args* raises
    ``ValueError`` (an energy that depends on an input).
    """
    em = Emitter(layout, target, kinds=kinds, what=what)
    em.prepare(exprs)
    codes = [em.emit(e) for e in exprs]
    tail = finish(em, codes)
    arg_names = {a.split("=")[0] for a in args}
    for slot in em.used.values():
        if slot.array not in arg_names:
            raise ValueError(
                f"codegen: {what} depends on the {slot.kind} '{slot.key}', but {name}() has"
                f" no '{slot.array}' argument"
            )
    body = _docstring(doc)
    if "p" in arg_names:
        body.append(target.check_params(params_const))
    body += em.preamble() + em.lines + tail
    head = f"def {name}({', '.join(args)}):"
    return "\n".join([head] + [_INDENT + line if line else "" for line in body]) + "\n"


def _return(em: Emitter, codes: list[Code]) -> list[str]:
    return [f"return {codes[0].text}"]


_ARGS = ("state", "u", "p", "a=None")


def to_source(expr: Expr, layout: StateLayout, target: str = "numpy", *, name: str = "f") -> str:
    """Python source of ``def name(state, u, p, a=None)`` returning the value of *expr*.

    ``state``, ``u``, ``a`` and ``p`` follow *layout* (``StateLayout``). ``a``, the
    stacked accelerations, is required only when *expr* holds an acceleration (an equation
    of ``compute_eom``), and ``u`` only when it holds an input; a required argument given
    as ``None`` raises ``ValueError``, an unused one may be ``None``. ``None`` never means
    zero here. The source is preceded by ``_PARAMS``, the parameter names, and uses the
    helpers ``_arr``, ``_scalar``, ``_vec``, ``_mat``, ``_hat``, ``_vee``, ``_cross``
    (NumPy) and ``_check_params`` that ``generate_module`` emits (``lambdify`` provides
    them). A vector is returned as a
    ``(3,)`` array, a matrix as ``(3, 3)``.

    Raises ``ValueError`` naming the symbol for a leaf that is not in the layout, a
    ``TimeDerivative`` of a non-leaf, a variation or a time integral.
    """
    tgt = get_target(target)
    const = f"_PARAMS = {pyliteral(tuple(layout.param_names))}\n\n\n"
    doc = f"Generated by geomech.codegen: {name}."
    return const + function_source(name, _ARGS, doc, layout, tgt, [expr], _return)


def namespace(target: str) -> dict:
    """The globals a generated function runs in: the array module(s) and the helpers."""
    import numpy as np

    ns: dict[str, Any] = {"np": np}
    if target == "casadi":
        ns["ca"] = require_casadi()
    exec(get_target(target).helpers, ns)
    return ns


def lambdify(expr: Expr, layout: StateLayout, target: str = "numpy") -> Callable:
    """``to_source`` compiled: a callable ``f(state, u, p, a=None)``.

    The callable carries its source in ``f.source``. With ``target="casadi"`` the
    arguments may be CasADi expressions and the result is a CasADi object; ``casadi`` must
    be installed (``ImportError`` otherwise).
    """
    src = to_source(expr, layout, target)
    ns = namespace(target)
    exec(compile(src, "<geomech.codegen>", "exec"), ns)
    fn = ns["f"]
    fn.source = src
    return fn
