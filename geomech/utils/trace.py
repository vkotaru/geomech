"""Debugging aids: the ``geomech`` logger, named rewrite rules and pipeline stage traces.

Three switches, all off by default and free when off:

``logging.getLogger("geomech")``
    One ``DEBUG`` line per pipeline stage of ``compute_eom``, ``integrate_by_parts`` (per
    target), ``to_standard_form`` (per equation) and ``linearize`` (per row), giving the
    stage name, the node count and the ``repr`` cut at 200 characters. Enable it with
    ``logging.getLogger("geomech").setLevel(logging.DEBUG)`` plus a handler, or with the
    environment variable ``GEOMECH_DEBUG=1`` (read when ``geomech`` is imported), which
    sets the level, attaches a ``stderr`` handler if neither the logger nor any ancestor
    (the root logger included) has one, and also makes
    every pipeline run ``geomech.check_tree`` on its input and output.

``GEOMECH_TRACE=1`` or ``geomech.set_trace(True)``
    Every firing of a named rewrite rule (see ``rule``) is logged on the
    ``geomech.rules`` logger as ``name: before ==> after``. The rule names and the
    hypotheses each rule needs are listed in the docstrings of the modules
    ``vector_identities``, ``simplification``, ``pulling``, ``expansion``
    (``geomech.core.transformations``), ``collecting``, ``ibp`` (``geomech.core.math``)
    and ``manifold_rules``.

``compute_eom(..., trace=EOMTrace())`` / ``linearize(..., trace=LinearizationTrace())``
    The trace object is filled stage by stage while the pipeline runs, so after an
    exception it holds every stage up to the failing one and ``failed_at`` names the stage
    that raised. A trace is cleared when a pipeline starts, so reusing one object keeps
    only the last run. ``summary()`` prints one line
    per stage (name, node count, start of the ``repr``); ``show(stage)`` prints one stage
    with ``print_tree(style="indent")``. See ``docs/pipeline.md``.
"""

from __future__ import annotations

import logging
import os
import sys
from dataclasses import MISSING, dataclass, field
from typing import Any, Callable, TextIO, TypeVar

from geomech.utils.printing import print_tree

_T = TypeVar("_T")

logger = logging.getLogger("geomech")
rules_logger = logging.getLogger("geomech.rules")

#: Characters of ``repr`` kept in a stage summary or a debug line.
STAGE_REPR_CHARS = 200
#: Characters of ``repr`` kept for each side of a rule-trace line.
RULE_REPR_CHARS = 300

_TRACE_RULES = False


def _env_flag(name):
    return os.environ.get(name, "").strip() not in ("", "0", "false", "False", "no")


def debug_enabled() -> bool:
    """True when ``GEOMECH_DEBUG`` is set (pipelines then run ``check_tree``)."""
    return _env_flag("GEOMECH_DEBUG")


def _ensure_handler(log):
    """Attach a ``stderr`` handler to the ``geomech`` logger if no handler would print *log*.

    ``log.hasHandlers()`` looks at *log* and every ancestor it propagates to, the root
    logger included, so after ``logging.basicConfig()`` nothing is added and every line
    is printed once (by the root handler). Configure logging *before* calling
    ``set_trace(True)``/importing with ``GEOMECH_DEBUG=1``: a handler added here and a
    root handler configured later would both print each line.
    """
    if log.hasHandlers():
        return
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(logging.Formatter("%(name)s: %(message)s"))
    logger.addHandler(handler)


def set_trace(on: bool = True) -> None:
    """Log every firing of a named rewrite rule (``name: before ==> after``).

    The lines go to the ``geomech.rules`` logger at ``DEBUG`` level; ``set_trace(True)``
    sets that level and attaches a ``stderr`` handler only if no handler is configured on
    that logger or any ancestor (so ``logging.basicConfig(level=logging.DEBUG)`` called
    first prints each line once, through the root handler).
    ``set_trace(False)`` stops the tracing (the logger configuration is left as is).
    Equivalent to running with ``GEOMECH_TRACE=1``.
    """
    global _TRACE_RULES
    _TRACE_RULES = bool(on)
    if _TRACE_RULES:
        rules_logger.setLevel(logging.DEBUG)
        _ensure_handler(rules_logger)


def trace_enabled() -> bool:
    """True while rule firings are being logged (``set_trace``/``GEOMECH_TRACE``)."""
    return _TRACE_RULES


def short_repr(expr: object, limit: int = STAGE_REPR_CHARS) -> str:
    """``repr(expr)`` cut to *limit* characters (with ``...``)."""
    s = repr(expr)
    return s if len(s) <= limit else s[: limit - 3] + "..."


def rule(name: str, before: object, after: _T) -> _T:
    """Return *after*; when rule tracing is on, log ``name: before ==> after``.

    Every rewriting ``case`` of the transformation passes returns through this helper,
    so the dotted *name* (``vector.dot_tangent_orthogonal``, ``simplify.hat_vee_skew``,
    ``collect.peel_matrix`` ...) identifies the rule that produced a node. *before* is
    the node as it entered the pass (for bottom-up passes its children may already have
    been rewritten by rules that log their own lines first); the ``simplify.*_fold``
    rules pass the node rebuilt from the rewritten children instead, so a fold is logged
    only for its own work. A call whose *after* prints like *before* (same ``repr`` and
    same ``str``) is not logged; when only the ``str`` differs (the ``repr`` hides
    grouping, ``(a·b)·c`` and ``a·(b·c)``), the line shows the two ``str`` forms, and a
    side longer than ``RULE_REPR_CHARS`` is cut around the first difference, so the two
    sides of a logged line never print the same. With tracing off this is one global
    read.
    """
    if _TRACE_RULES:
        before_s, after_s = repr(before), repr(after)
        if before_s == after_s:
            before_s, after_s = str(before), str(after)
        if before_s != after_s:
            rules_logger.debug("%s: %s ==> %s", name, *_cut_pair(before_s, after_s))
    return after


def _cut(s, limit):
    return s if len(s) <= limit else s[: limit - 3] + "..."


def _cut_pair(a, b, limit=RULE_REPR_CHARS):
    """*a* and *b* cut to *limit* characters so that the cut strings still differ.

    When the two agree on more than the first ``limit - 60`` characters, both are shown
    from shortly before their first difference (``...`` marks the dropped start).
    """
    a_cut, b_cut = _cut(a, limit), _cut(b, limit)
    if a_cut != b_cut:
        return a_cut, b_cut
    common = 0
    for ca, cb in zip(a, b):
        if ca != cb:
            break
        common += 1
    start = max(0, common - 60)
    return "..." + _cut(a[start:], limit - 3), "..." + _cut(b[start:], limit - 3)


def node_count(expr: object) -> int:
    """Number of nodes (leaves included) of the tree *expr*."""
    count = 0
    stack = [expr]
    while stack:
        node = stack.pop()
        count += 1
        stack.extend(getattr(node, "nodes", None) or ())
    return count


def configure_from_env() -> None:
    """Apply ``GEOMECH_DEBUG`` and ``GEOMECH_TRACE`` (called once when ``geomech`` is imported)."""
    if debug_enabled():
        logger.setLevel(logging.DEBUG)
        _ensure_handler(logger)
    if _env_flag("GEOMECH_TRACE"):
        set_trace(True)


# ---------------------------------------------------------------------------
# Stage recording shared by the pipelines
# ---------------------------------------------------------------------------


def stage_recorder(trace: Any, where: str) -> Callable[[str, object], None] | None:
    """A ``record(stage, value)`` callable for one pipeline run, or ``None`` when idle.

    *trace* is an ``EOMTrace``/``LinearizationTrace`` (or ``None``); *where* prefixes the
    debug lines (``compute_eom``). Returns ``None`` when there is no trace and the
    ``geomech`` logger is not enabled for ``DEBUG``, so the pipelines skip every
    recording call on the happy path.
    """
    log = logger.isEnabledFor(logging.DEBUG)
    if trace is None and not log:
        return None

    def record(stage, value):
        if trace is not None:
            trace._add(stage, value)
        if log:
            if hasattr(value, "nodes") or hasattr(value, "type"):
                logger.debug(
                    "%s: %s [%d nodes] %s",
                    where,
                    stage,
                    node_count(value),
                    short_repr(value),
                )
            else:
                logger.debug("%s: %s %s", where, stage, _cut(repr(value), STAGE_REPR_CHARS))

    return record


_NOT_STAGES = ("stages", "failed_at", "error")


@dataclass
class _Trace:
    """Ordered ``(stage name, expression)`` pairs recorded by one pipeline run.

    The pipeline clears the trace when it starts (``_reset``), so a trace object passed
    to several runs holds the last run only. ``failed_at`` is ``None`` after a run that
    returned; after a run that raised it is the name of the stage that was running (that
    stage is not in ``stages``), and ``error`` is the exception's class name.
    """

    stages: list = field(default_factory=list)
    failed_at: str | None = None
    error: str | None = None

    def _reset(self):
        """Clear every recorded stage and attribute (called at pipeline entry)."""
        for name, spec in self.__dataclass_fields__.items():
            if spec.default_factory is not MISSING:
                setattr(self, name, spec.default_factory())
            else:
                setattr(self, name, spec.default)

    def _fail(self, stage, err):
        self.failed_at = stage
        self.error = type(err).__name__

    def _add(self, stage, value):
        self.stages.append((stage, value))
        self._store(stage, value)

    def _store(self, stage, value):  # pragma: no cover - overridden
        pass

    def stage_names(self) -> list[str]:
        """The recorded stage names, in pipeline order."""
        return [name for name, _ in self.stages]

    def __getitem__(self, stage):
        """The expression recorded under *stage* (the last one if recorded twice).

        An attribute name (``after_ibp``, ``equations``, ``dE`` ...) is accepted as well.
        """
        for name, value in reversed(self.stages):
            if name == stage:
                return value
        if stage not in _NOT_STAGES and stage in getattr(self, "__dataclass_fields__", {}):
            return getattr(self, stage)
        raise KeyError(
            f"{type(self).__name__}: no stage {stage!r}; recorded stages: {self.stage_names()}"
        )

    def format_summary(self, chars: int = STAGE_REPR_CHARS) -> str:
        """One line per stage: name, node count and the first *chars* characters of repr.

        After a run that raised, the first line is ``FAILED at <stage>: <exception>``.
        """
        head = []
        if self.failed_at is not None:
            head.append(
                f"FAILED at {self.failed_at}: {self.error} (the stages below ran before it)"
            )
        if not self.stages:
            return "\n".join(head + [f"{type(self).__name__}: no stage recorded"])
        width = max(len(name) for name, _ in self.stages)
        lines = []
        for name, value in self.stages:
            if hasattr(value, "nodes") or hasattr(value, "type"):
                count = f"{node_count(value):5d} nodes"
                text = short_repr(value, chars)
            else:
                count = " " * 11
                text = _cut(repr(value), chars)
            lines.append(f"{name:<{width}}  {count}  {text}")
        return "\n".join(head + lines)

    def summary(self, chars: int = STAGE_REPR_CHARS, file: TextIO | None = None) -> None:
        """Print ``format_summary()`` (one line per stage) to *file* (default stdout)."""
        print(self.format_summary(chars), file=file or sys.stdout)

    def show(self, stage: str, style: str = "indent") -> None:
        """Print the tree of one stage with ``print_tree(..., style="indent")``."""
        value = self[stage]
        if isinstance(value, dict):
            for key, item in value.items():
                print(f"--- {stage}[{key}]")
                print_tree(item, style=style)
            return
        print_tree(value, style=style)

    def __str__(self):
        return self.format_summary()


@dataclass(repr=False)
class EOMTrace(_Trace):
    """Intermediate results of one ``compute_eom`` run (``compute_eom(..., trace=tr)``).

    Attributes follow the pipeline order (``docs/pipeline.md``); each is ``None`` until
    its stage has run, so after an exception the trace holds every stage before it.

    dL, dW, dS
        The variation of the Lagrangian, the virtual work (``Zero`` for ``None``) and
        their sum δS = δL + δW.
    after_manifold_rules
        δS with ``δω``/``δΩ`` replaced by their kinematic expansions.
    after_simplify
        ``full_simplify`` of the previous stage.
    variation_vectors, ibp_targets
        The variation vectors (``δx``, ``ξ_q``, ``η_R``) and their time derivatives.
    after_ibp_prepare
        The ``expand`` + ``full_simplify`` that ``integrate_by_parts`` does first.
    ibp
        ``{printed target: {"collect": ..., "integrate": ..., "simplify": ...}}``, the
        expression after ``collect``, after the integration-by-parts pass and after
        ``full_simplify``, one entry per target in order. Not to be confused with the
        *stage* ``"ibp"``: ``tr["ibp"]`` and ``tr.show("ibp")`` give ``after_ibp``; one
        per-target step is the stage ``"ibp.collect[<target>]"`` (and ``integrate``,
        ``simplify``).
    after_ibp, after_expand
        δS after integration by parts, and after the final ``expand``.
    raw_equations, equations
        ``{printed variation vector: expression}``: the extracted coefficient before and
        after its final ``full_simplify``.
    """

    dL: object = None
    dW: object = None
    dS: object = None
    after_manifold_rules: object = None
    after_simplify: object = None
    variation_vectors: list | None = None
    ibp_targets: list | None = None
    after_ibp_prepare: object = None
    ibp: dict = field(default_factory=dict)
    after_ibp: object = None
    after_expand: object = None
    raw_equations: dict = field(default_factory=dict)
    equations: dict = field(default_factory=dict)

    _SIMPLE = {
        "dL": "dL",
        "dW": "dW",
        "dS": "dS",
        "manifold_rules": "after_manifold_rules",
        "simplify": "after_simplify",
        "ibp.prepare": "after_ibp_prepare",
        "ibp": "after_ibp",
        "expand": "after_expand",
    }

    def _store(self, stage, value):
        attr = self._SIMPLE.get(stage)
        if attr is not None:
            setattr(self, attr, value)
        elif stage.startswith("ibp.") and "[" in stage:
            step, target = stage[4:].split("[", 1)
            self.ibp.setdefault(target[:-1], {})[step] = value
        elif stage.startswith("raw["):
            self.raw_equations[stage[4:-1]] = value
        elif stage.startswith("equation["):
            self.equations[stage[9:-1]] = value

    def __getitem__(self, stage):
        attr = {v: k for k, v in self._SIMPLE.items()}.get(stage)
        return super().__getitem__(attr or stage)


@dataclass(repr=False)
class LinearizationTrace(_Trace):
    """Intermediate results of one ``linearize`` run (``linearize(..., trace=tr)``).

    Per dynamic row ``label`` (the printed variation vector of the equation):

    dE
        ``{label: δE}``, the expanded, simplified variation of the equation.
    buckets
        ``{label: {printed perturbation symbol: expression}}``, the partition of δE's
        terms by the perturbation symbol each carries (``_partition_terms``); each bucket
        is the sum of its terms.
    raw_blocks
        ``{label: {printed symbol: coefficient}}``, the coefficient straight out of
        ``extract_linear_coeff``, before ``full_simplify`` and before the mass blocks are
        negated.
    blocks
        ``{label: {printed symbol: block}}``, the final entry of ``M``, ``A`` or ``B``.

    Stage names in ``stages`` are ``dE[label]``, ``bucket[label][sym]`` (the combined
    terms), ``raw[label][sym]`` and ``block[label][sym]``.
    """

    dE: dict = field(default_factory=dict)
    buckets: dict = field(default_factory=dict)
    raw_blocks: dict = field(default_factory=dict)
    blocks: dict = field(default_factory=dict)

    def _store(self, stage, value):
        kind, _, rest = stage.partition("[")
        if kind == "dE":
            self.dE[rest[:-1]] = value
            return
        label, _, sym = rest.partition("][")
        sym = sym[:-1]
        target = {"bucket": self.buckets, "raw": self.raw_blocks, "block": self.blocks}.get(kind)
        if target is not None:
            target.setdefault(label, {})[sym] = value
