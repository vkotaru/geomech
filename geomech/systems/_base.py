"""``System``, ``Reference``, the registry and the numeric helpers shared by the builders."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Callable, NamedTuple

import numpy as np

from geomech.core.base import S2, Expr, Scalar
from geomech.dynamics.variables import SystemVariables

if TYPE_CHECKING:  # pragma: no cover
    from geomech.codegen import StateLayout

__all__ = ["Reference", "System", "SystemInfo", "get", "list_systems", "register"]

#: ½ as a constant scalar leaf, written the way the rest of the package writes it.
HALF = Scalar("0.5", value=0.5, attr=["Constant"])

#: The largest size a parametric builder accepts (``n_links``, ``n_quads``).
MAX_SIZE = 20


def frozen_array(values) -> np.ndarray:
    """A read-only float copy of *values*."""
    out = np.array(values, dtype=float)
    out.setflags(write=False)
    return out


E3 = frozen_array([0.0, 0.0, 1.0])
#: Two inertias about the centre of mass (symmetric positive definite, principal moments
#: satisfying the triangle inequality) that are not diagonal, so that no product of
#: inertia is silently assumed zero.
J_A = frozen_array([[0.02, 0.001, 0.0], [0.001, 0.025, 0.002], [0.0, 0.002, 0.04]])
J_B = frozen_array([[0.5, 0.0, 0.05], [0.0, 0.7, 0.0], [0.05, 0.0, 0.4]])


class FrozenDict(dict):  # type: ignore[type-arg]
    """A ``dict`` that refuses writes, and still copies and pickles.

    ``System`` fields hold their parameters in one so that changing a default in place
    raises instead of changing every system built afterwards. Copy it to change a value:
    ``p = dict(sys.params); p["m"] = 2.0``.
    """

    def _readonly(self, *args: Any, **kwargs: Any) -> None:
        raise TypeError(
            "System mappings are read-only; copy first: p = dict(sys.params); p['m'] = 2.0"
        )

    __setitem__ = __delitem__ = _readonly  # type: ignore[assignment]
    clear = pop = popitem = setdefault = update = _readonly  # type: ignore[assignment]
    __ior__ = _readonly  # type: ignore[assignment]

    def __reduce__(self) -> tuple[Any, ...]:
        return (type(self), (dict(self),))

    def __copy__(self) -> FrozenDict:
        return self

    def __hash__(self) -> int:  # type: ignore[override]
        raise TypeError("unhashable type: 'FrozenDict'")


@dataclass(frozen=True)
class Reference:
    """Where a system comes from.

    Attributes:
        source: a citation key of ``docs/references.bib`` (``"lee2018global"``), or a short
            text when the system is not taken from the literature (``"elementary"``,
            ``"check_eom"``: verified numerically only).
        equation: the closed-form equation(s) of *source* that the test suite compares the
            derived equations with (``"eq. 5.38"``), or ``""`` when there is none. Only
            equations with such a test are named here.
        section: where in *source* the model is treated (background reading, not a
            comparison): ``"§5.6.1; general form eq. 5.14"``, or ``""``.
        note: one or two sentences on the model and on how it relates to the source.
    """

    source: str
    equation: str = ""
    section: str = ""
    note: str = ""

    def __str__(self) -> str:
        return ", ".join(x for x in (self.source, self.section, self.equation) if x)


@dataclass(frozen=True, eq=False)
class System:
    """A mechanical system ready for ``compute_eom``, ``to_standard_form``, code generation,
    simulation and visualization.

    Attributes:
        name: the registry name (``get(name)`` rebuilds it).
        configuration: the configuration manifold, written like ``"R³ × SO(3) × S²"``.
        L: the Lagrangian.
        dW: the virtual work of the inputs.
        variables: the configuration variables.
        inputs: the input leaves, in the order of the columns of ``G`` and of ``u``.
        params: default numeric value of every parameter of the system, keyed by printed
            symbol (``float`` for a scalar, ``(3,)`` for a vector, ``(3, 3)`` for a matrix);
            the ``p`` of ``simulate``. A read-only mapping of read-only arrays, owned by
            this system: change a value on a copy, ``p = dict(sys.params); p["m"] = 2.0``.
        x0: a consistent initial state (``|q| = 1``, ``q·ω = 0``, ``R ∈ SO(3)``), packed with
            ``StateLayout.pack``.
        reference: the source of the model (``Reference``).
        symbols: the leaves and point expressions a user may want (``sys.symbols["xQ"]``), a
            read-only mapping (lists of symbols are tuples).
        scene_builder: ``scene_builder(layout, symbols)`` returns a ``geomech.viz.Scene``;
            ``None`` when the system has no scene. Call ``scene(layout)`` instead.
        options: the keyword arguments the builder was called with (read-only mapping).
    """

    name: str
    configuration: str
    L: Expr
    dW: Expr
    variables: SystemVariables
    inputs: tuple
    params: Mapping[str, Any]
    x0: np.ndarray
    reference: Reference
    symbols: Mapping[str, Any] = field(default_factory=lambda: FrozenDict({}))
    scene_builder: Callable | None = field(default=None, repr=False)
    options: Mapping[str, Any] = field(default_factory=lambda: FrozenDict({}))

    # ------------------------------------------------------------------ derived data

    @property
    def dof(self) -> int:
        """Degrees of freedom: 1 per scalar, 3 per vector, 2 per S² point, 3 per SO(3)."""
        v = self.variables
        n = len(v.scalars) + 3 * len(v.matrices)
        n += sum(2 if isinstance(q, S2) else 3 for q in v.vectors)
        return n

    @property
    def has_scene(self) -> bool:
        return self.scene_builder is not None

    def derive(self) -> dict:
        """``compute_eom(L, dW, variables)`` (not cached: each call derives again)."""
        from geomech.dynamics.eom import compute_eom

        return compute_eom(self.L, self.dW, self.variables)

    def layout(self) -> StateLayout:
        """The ``StateLayout`` of the system, with every parameter of ``L`` and ``dW``."""
        from geomech.codegen import StateLayout, collect_params

        exprs = [self.L, self.dW]
        params = collect_params(exprs, self.variables, list(self.inputs))
        return StateLayout(self.variables, list(self.inputs), params)

    def scene(self, layout: StateLayout | None = None):
        """A ``geomech.viz.Scene`` of the bodies on *layout* (``self.layout()`` by
        default), or ``None`` for a system without one."""
        if self.scene_builder is None:
            return None
        return self.scene_builder(self.layout() if layout is None else layout, self.symbols)

    def as_tuple(self) -> tuple:
        """``(L, dW, variables, inputs)``, the ``system=`` argument of ``Report``."""
        return (self.L, self.dW, self.variables, list(self.inputs))


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------


class SystemInfo(NamedTuple):
    """One row of ``list_systems()``."""

    name: str
    configuration: str
    dof: int
    reference: Reference
    parameters: tuple  # keyword arguments of the builder, e.g. ("n_links",)
    summary: str  # first sentence of the builder's docstring


#: The order of ``list_systems()`` (the table of ``docs/systems.md``).
ORDER = (
    "point_mass",
    "spring_mass",
    "spherical_pendulum",
    "particle_on_rotating_sphere",
    "double_spherical_pendulum",
    "n_link_spherical_chain",
    "elastic_spherical_pendulum",
    "cart_spherical_pendulum",
    "free_rigid_body",
    "heavy_top",
    "rigid_pendulum",
    "double_rigid_pendulum",
    "two_rigid_bodies_spherical_joint",
    "rigid_body",
    "quadrotor",
    "quadrotor_offset_com",
    "quadrotor_suspended_load",
    "quadrotor_elastic_load",
    "quadrotor_flexible_cable",
    "multi_quadrotor_point_mass",
    "multi_quadrotor_rigid_load",
)

_REGISTRY: dict[str, Callable[..., System]] = {}
_CONFIG: dict[str, str] = {}


def register(configuration: str):
    """Register a builder under its function name (which must be listed in ``ORDER``)."""

    def deco(fn):
        if fn.__name__ not in ORDER:
            raise ValueError(f"register: add {fn.__name__!r} to geomech.systems._base.ORDER")
        _REGISTRY[fn.__name__] = fn
        _CONFIG[fn.__name__] = configuration
        return fn

    return deco


def get(name: str, **kwargs) -> System:
    """The system *name* built with *kwargs* (``get("n_link_spherical_chain", n_links=4)``).

    Raises ``KeyError`` listing the registered names for an unknown *name*.
    """
    try:
        build = _REGISTRY[name]
    except KeyError:
        raise KeyError(
            f"geomech.systems.get: no system named {name!r}; the systems are {names()}"
        ) from None
    return build(**kwargs)


def _summary(fn) -> str:
    doc = " ".join((fn.__doc__ or "").strip().split("\n\n")[0].split())
    return doc.split(". ")[0].rstrip(".") + "." if doc else ""


def list_systems() -> list[SystemInfo]:
    """Every registered system built with its default arguments: name, configuration,
    degrees of freedom, reference, the builder's keyword arguments and a one-line summary,
    in a fixed order (the order of the table in the docs)."""
    import inspect

    out = []
    for name in names():
        build = _REGISTRY[name]
        s = build()
        kw = tuple(inspect.signature(build).parameters)
        out.append(SystemInfo(name, _CONFIG[name], s.dof, s.reference, kw, _summary(build)))
    return out


def names() -> list[str]:
    """The registered names, in ``list_systems()`` order."""
    return [n for n in ORDER if n in _REGISTRY]


# ---------------------------------------------------------------------------
# Helpers for the builders
# ---------------------------------------------------------------------------


def make_system(
    name,
    configuration,
    L,
    dW,
    variables,
    inputs,
    params,
    reference,
    *,
    x0: dict,
    symbols=None,
    scene=None,
    options=None,
) -> System:
    """Build a ``System``, packing *x0* (slot name → value) with the system's layout.

    Every array parameter is copied and made read-only, and ``params``, ``symbols`` and
    ``options`` are wrapped in read-only mappings, so a system never shares mutable state
    with another one or with this module.
    """
    from geomech.codegen import StateLayout

    state = StateLayout(variables, list(inputs)).pack(**x0)
    state.setflags(write=False)
    params = {
        k: frozen_array(v) if np.ndim(v) else float(np.asarray(v, dtype=float))
        for k, v in params.items()
    }
    symbols = {k: tuple(v) if isinstance(v, list) else v for k, v in (symbols or {}).items()}
    return System(
        name=name,
        configuration=configuration,
        L=L,
        dW=dW,
        variables=variables,
        inputs=tuple(inputs),
        params=FrozenDict(params),
        x0=state,
        reference=reference,
        symbols=FrozenDict(symbols),
        scene_builder=scene,
        options=FrozenDict(dict(options or {})),
    )


def check_count(value: Any, what: str, minimum: int = 1, maximum: int = MAX_SIZE) -> int:
    """*value* as an ``int`` in ``[minimum, maximum]``, else ``ValueError`` naming *what*."""
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, np.integer))
        or not minimum <= value <= maximum
    ):
        raise ValueError(f"{what} must be an int from {minimum} to {maximum}, got {value!r}")
    return int(value)


def unit(v) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v)


def tangent(q, v) -> np.ndarray:
    """*v* projected onto the tangent plane at the unit vector *q*."""
    q = np.asarray(q, dtype=float)
    v = np.asarray(v, dtype=float)
    return v - (v @ q) * q


def rot(v) -> np.ndarray:
    """The rotation ``expm(hat(v))``."""
    from geomech.numeric import expm_so3

    return expm_so3(np.asarray(v, dtype=float))


#: Short text of the citation keys of ``docs/references.bib`` and of the two plain sources.
CITATIONS = {
    "lee2018global": "Lee et al.",
    "lee2010geometric": "Lee, Leok, McClamroch (2010)",
    "sreenath2013geometric": "Sreenath, Lee, Kumar (2013)",
    "elementary": "elementary",
    "check_eom": "—",
}


def source_text(ref: Reference) -> str:
    """``"Lee et al. §5.6.1; general form eq. 5.14"``: the source of *ref* as the docs
    table shows it (not the compared equation; see ``compared_text``)."""
    text = CITATIONS.get(ref.source, ref.source)
    if ref.source == "check_eom":
        return "—"
    return f"{text} {ref.section}" if ref.section else text


def compared_text(ref: Reference) -> str:
    """What the test suite compares the system's equations with, as the docs table shows
    it: ``"Lee et al. eq. 5.38"``, or ``"—"`` when no closed form is compared (every
    system is still checked by ``check_eom`` and against a model of its bodies)."""
    if ref.equation:
        return f"{CITATIONS.get(ref.source, ref.source)} {ref.equation}"
    return "—"


def systems_table() -> str:
    """The Markdown table of ``list_systems()`` (name, configuration, dof, source, what the
    equations are compared with), as checked in to ``docs/systems.md``."""
    rows = [
        "| Name | Configuration | dof | Source | Closed form compared |",
        "|---|---|---|---|---|",
    ]
    for info in list_systems():
        call = f"{info.name}({', '.join(f'{p}=…' for p in info.parameters)})"
        call = call if info.parameters else info.name
        rows.append(
            f"| `{call}` | {info.configuration} | {info.dof} | {source_text(info.reference)}"
            f" | {compared_text(info.reference)} |"
        )
    return "\n".join(rows) + "\n"
