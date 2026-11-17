"""Generate the math of the model pages in ``docs/models/`` from the library.

Each page has hand-written prose, a diagram and a code block, plus generated blocks
between ``<!-- model:<block> -->`` and ``<!-- /model:<block> -->`` markers:

* ``energy``: the kinetic and potential terms of ``L`` (addends with a velocity are
  kinetic, the others are ``−V``), the Lagrangian and the virtual work;
* ``eom``: the derived equations, one per configuration variable (long ones folded);
* ``checks``: ``check_eom``, the energy drift of an unforced run from the page's initial
  state, and the closed form of the literature the test suite compares, if any.

``docs/models/index.md`` is generated whole from :data:`MODELS`.
``tests/test_model_pages.py`` fails when a page no longer matches the library.

    python scripts/model_pages.py                     # rewrite every page that exists
    python scripts/model_pages.py --only quadrotor    # one page
    python scripts/model_pages.py --reports DIR       # write the interactive reports

``docs/conf.py`` calls :func:`write_reports` at build time, so the animations are built
from the library each time and are not stored in the repository.
"""

from __future__ import annotations

import argparse
import html
import re
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PAGES = ROOT / "docs" / "models"
BLOCK = re.compile(r"(<!-- model:(\w+) -->\n)(.*?)(<!-- /model:\2 -->)", re.S)
#: An equation block longer than this many lines in all is shown folded, one per variable.
FOLD_LINES = 40


@dataclass(frozen=True)
class Page:
    """One model page.

    Attributes:
        name: the library system (``geomech.systems.get(name)``), also the page file name.
        title: the page and report title.
        group: the section of the models index.
        t_end: simulated time of the animation and of the energy check.
        x0: ``x0(system, layout)``, the initial state when the library's default start
            makes a dull animation; ``None`` uses ``system.x0``.
        u: ``u(system, layout)``, constant inputs of the animation (a quadrotor needs
            thrust); ``None`` runs unforced. The energy check always runs unforced.
        scene: a scene for a system the library draws none for: ``(params, build)``,
            where *params* maps extra drawing-only parameters (symbol → value) added to
            the layout and *build(layout, symbols)* returns the ``Scene``.
        options: keyword arguments of the library builder (``n_links=5``): the page,
            its animation and its checks use ``get(name, **options)``.
        indexed: the model of ``geomech.indexed.systems`` with a symbolic number of
            bodies n that this system is an instance of; the page shows its generic
            equations, derived once for every n.
    """

    name: str
    title: str
    group: str
    t_end: float = 10.0
    x0: Callable | None = None
    u: Callable | None = None
    scene: tuple | None = None
    options: dict = field(default_factory=dict)
    indexed: str | None = None

    def system(self):
        from geomech.systems import get

        return get(self.name, **self.options)

    @property
    def size(self) -> int | None:
        """The number of bodies n of an indexed page."""
        return next(iter(self.options.values())) if self.indexed else None

    def initial_state(self, system, layout) -> np.ndarray:
        return system.x0 if self.x0 is None else np.asarray(self.x0(system, layout))

    def inputs(self, system, layout) -> np.ndarray | None:
        return None if self.u is None else np.asarray(self.u(system, layout), dtype=float)


# ------------------------------------------------------------ page settings


def _spring_scene():
    """The block slides along ``e₁`` at ``x e₁``; a spring joins it to a wall at
    ``−e₁``. ``e₁`` is a drawing-only parameter: ``L`` is one-dimensional."""
    from geomech import Vector

    e1 = Vector("e1", attr=["Constant"])

    def build(layout, s):
        from geomech.viz import PointMass, Rod, Scene, Trace

        scene = Scene(layout)
        block = e1 * s["x"]
        scene.add(Rod([-1.0, 0.0, -0.3], [-1.0, 0.0, 0.3], radius=0.03, name="wall"))
        scene.add(Rod([-1.0, 0.0, 0.0], block, radius=0.015, color="cable", name="spring"))
        scene.add(PointMass(block, radius=0.12, name="m"))
        scene.add(Trace(block))
        return scene

    return ({e1: np.array([1.0, 0.0, 0.0])}, build)


def _swirl(system, layout):
    """1 rad off the vertical, spinning about e₂ and about q × e₂: a 3-D rosette, where
    the library's default start (spin about e₂ only) stays in one vertical plane."""
    from geomech.systems.spheres import Q0

    n = np.cross(Q0, [0.0, 1.0, 0.0])
    n /= np.linalg.norm(n)
    return layout.pack(q=Q0, omega_q=1.2 * np.array([0.0, 1.0, 0.0]) + 2.0 * n)


def _elastic_swirl(system, layout):
    """As :func:`_swirl`, with the spring stretched to l = 1.1 and lengthening at 0.2 m/s:
    the library's default start stays nearly in one vertical plane."""
    from geomech.systems.spheres import Q0

    n = np.cross(Q0, [0.0, 1.0, 0.0])
    n /= np.linalg.norm(n)
    w = 1.2 * np.array([0.0, 1.0, 0.0]) + 2.0 * n
    return layout.pack(l=1.1, l_dot=0.2, q=Q0, omega_q=w)


def weight(system, layout):
    """``F = m g e₃`` at the centre of mass cancels gravity; no moment."""
    p = system.params
    u = np.zeros(layout.n_inputs)
    u[layout.input_slices["F"]] = p["m"] * p["g"] * np.asarray(p["e3"])
    return u


def _drift(system, layout):
    """The library's start, drifting at 0.1 m/s along e₁ instead of thrown upward."""
    state = layout.unpack(system.x0)
    state["x_dot"] = np.array([0.1, 0.0, 0.0])
    return layout.pack(**state)


def _toss(system, layout):
    """Body 1 thrown up at 4.9 m/s, both bodies tumbling at three times the library's
    rates: a one-second flight instead of a 470 m fall."""
    state = layout.unpack(system.x0)
    state["x_dot"] = np.array([0.2, 0.0, 4.9])
    state["Omega_R1"] = 3 * np.asarray(state["Omega_R1"])
    state["Omega_R2"] = 3 * np.asarray(state["Omega_R2"])
    return layout.pack(**state)


def _step(system, layout):
    """One planar step: stance leg 0.25 rad behind the vertical vaulting forward, swing
    foot on the ground behind, swinging forward (touch-down at t ≈ 0.775 s)."""
    from geomech.systems._base import rot

    return layout.pack(
        R1=rot([0.0, -0.25, 0.0]),
        R2=rot([0.0, 0.25, 0.0]),
        Omega_R1=[0.0, 1.0, 0.0],
        Omega_R2=[0.0, -1.0, 0.0],
    )


def wheel_torque(system, layout):
    """0.02 N·m on wheel 1 for the whole run; the bus takes the reaction."""
    u = np.zeros(layout.n_inputs)
    u[layout.input_slices["\\tau_1"]] = 0.02
    return u


def hover(system, layout):
    """Each thrust input ``f…`` carries an equal share of the total weight; no moments."""
    p = system.params
    weight = p["g"] * sum(float(v) for k, v in p.items() if k.startswith("m") and np.ndim(v) == 0)
    thrusts = [k for k in layout.input_slices if k.startswith("f")]
    u = np.zeros(layout.n_inputs)
    for k in thrusts:
        u[layout.input_slices[k]] = weight / len(thrusts)
    return u


def level(system, layout):
    """The library's start with every quadrotor level and not rotating (attitudes
    ``R``, ``R1``, …, but not a load's ``R_L``): with :func:`hover` thrust and no moment
    the attitude equation is that of a free rigid body at rest, so the vehicles stay
    level, and the motion is the loads and cables swinging, not an uncontrolled tumble."""
    state = layout.unpack(system.x0)
    for key in [k for k in state if re.fullmatch(r"R\d*", k)]:
        state[key] = np.eye(3)
        state[f"Omega_{key}"] = np.zeros(3)
    return layout.pack(**state)


def _tethered(system, layout):
    """A static equilibrium of the vehicles at the library's cable directions: equal
    cable tensions carry the load (``R_L`` level), and each vehicle is tilted so its
    thrust balances its weight plus its cable's pull. Returns ``(state, thrusts)``; the
    load and cable velocities keep the library's values, the vehicles do not rotate."""
    p, state = system.params, layout.unpack(system.x0)
    n = sum(1 for k in state if re.fullmatch(r"q\d+", k))
    if "R_L" in state:
        state["R_L"] = np.eye(3)
    e3 = np.array([0.0, 0.0, 1.0])
    tension = p["m_L"] * p["g"] / sum(-state[f"q{i}"][2] for i in range(1, n + 1))
    thrusts = {}
    for i in range(1, n + 1):
        F = p[f"m{i}"] * p["g"] * e3 - tension * state[f"q{i}"]
        b3 = F / np.linalg.norm(F)
        k = np.cross(e3, b3)
        s, c = np.linalg.norm(k), b3[2]
        K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]]) / s
        state[f"R{i}"] = np.eye(3) + s * K + (1 - c) * K @ K
        state[f"Omega_R{i}"] = np.zeros(3)
        thrusts[f"f{i}"] = float(np.linalg.norm(F))
    return state, thrusts


def tethered(system, layout):
    return layout.pack(**_tethered(system, layout)[0])


def tethered_thrust(system, layout):
    u = np.zeros(layout.n_inputs)
    for key, f in _tethered(system, layout)[1].items():
        u[layout.input_slices[key]] = f
    return u


PARTICLES = "Particles and pendulums"
RIGID = "Rigid bodies"
SPACE = "Spacecraft"
AERIAL = "Aerial vehicles"

MODELS = [
    Page("point_mass", "Point mass", PARTICLES),
    Page("spring_mass", "Spring–mass", PARTICLES, scene=_spring_scene()),
    Page("planar_pendulum", "Planar pendulum", PARTICLES),
    Page("double_planar_pendulum", "Double planar pendulum", PARTICLES),
    Page("cart_planar_pendulum", "Cart with a planar pendulum", PARTICLES),
    Page("spherical_pendulum", "Spherical pendulum", PARTICLES, x0=_swirl),
    Page("particle_on_rotating_sphere", "Particle on a rotating sphere", PARTICLES),
    Page("double_spherical_pendulum", "Double spherical pendulum", PARTICLES),
    Page(
        "n_link_spherical_chain",
        "Chain of spherical pendulums",
        PARTICLES,
        options={"n_links": 5},
        indexed="spherical_chain",
    ),
    Page("elastic_spherical_pendulum", "Elastic spherical pendulum", PARTICLES, x0=_elastic_swirl),
    Page("cart_spherical_pendulum", "Cart with a spherical pendulum", PARTICLES),
    Page("free_rigid_body", "Free rigid body", RIGID),
    Page("heavy_top", "Heavy top", RIGID),
    Page("rigid_pendulum", "Rigid pendulum", RIGID),
    Page("double_rigid_pendulum", "Double rigid pendulum", RIGID),
    Page(
        "n_link_rigid_chain",
        "Chain of rigid bodies",
        RIGID,
        options={"n_links": 5},
        indexed="rigid_chain",
    ),
    Page(
        "two_rigid_bodies_spherical_joint",
        "Two rigid bodies, spherical joint",
        RIGID,
        1.0,
        x0=_toss,
    ),
    Page("rigid_body", "Rigid body in space", RIGID, x0=_drift, u=weight),
    Page("biped_3d", "3-D biped", RIGID, 0.775, x0=_step),
    Page("spacecraft_reaction_wheels", "Spacecraft with reaction wheels", SPACE, u=wheel_torque),
    Page(
        "spacecraft_reaction_wheels_s1",
        "Spacecraft with reaction wheels on S¹",
        SPACE,
        u=wheel_torque,
    ),
    Page("spacecraft_solar_arrays", "Spacecraft with solar arrays", SPACE),
    Page("quadrotor", "Quadrotor", AERIAL, 6.0, x0=level, u=hover),
    Page(
        "quadrotor_offset_com", "Quadrotor, offset centre of mass", AERIAL, 2.0, x0=level, u=hover
    ),
    Page(
        "quadrotor_suspended_load",
        "Quadrotor with a suspended load",
        AERIAL,
        6.0,
        x0=level,
        u=hover,
    ),
    Page(
        "quadrotor_elastic_load", "Quadrotor with an elastic load", AERIAL, 6.0, x0=level, u=hover
    ),
    Page(
        "quadrotor_flexible_cable",
        "Quadrotor with a flexible cable",
        AERIAL,
        6.0,
        x0=level,
        u=hover,
        options={"n_links": 5},
        indexed="quadrotor_flexible_cable",
    ),
    Page(
        "multi_quadrotor_point_mass",
        "Quadrotors carrying a point mass",
        AERIAL,
        6.0,
        x0=tethered,
        u=tethered_thrust,
        options={"n_quads": 5},
        indexed="multi_quadrotor_point_mass",
    ),
    Page(
        "multi_quadrotor_rigid_load",
        "Quadrotors carrying a rigid body",
        AERIAL,
        6.0,
        x0=tethered,
        u=tethered_thrust,
        options={"n_quads": 5},
        indexed="multi_quadrotor_rigid_load",
    ),
    Page("multi_quadrotor_hose", "Two quadrotors carrying a hose", AERIAL, 6.0, x0=level, u=hover),
]
BY_NAME = {p.name: p for p in MODELS}


def existing() -> list[Page]:
    """The pages that have been written (the others are listed but not generated)."""
    return [p for p in MODELS if (PAGES / f"{p.name}.md").exists()]


# ---------------------------------------------------------------------- LaTeX


#: Library parameter names that ``latex_names`` would print as a product or with a broken
#: subscript: ``m_Q0`` → ``m_{Q0}``, ``Ia_1`` → ``I_{a,1}``.
_PARAM_NAMES = [
    (re.compile(r"^([A-Za-z])_([A-Za-z]+)(\d+)$"), r"\1_{\2\3}"),
    (re.compile(r"^([A-Z])([a-z])_(\d+)$"), r"\1_{\2,\3}"),
]


def _names(system, eom):
    from geomech.viz.latex import latex_names

    v = system.variables
    names = _library_names(system, eom, v, latex_names)
    for key in system.params:
        for pat, rep in _PARAM_NAMES:
            if pat.match(key):
                names[key] = pat.sub(rep, key)
    return names


def _library_names(system, eom, v, latex_names):
    return latex_names(
        system.L,
        system.dW,
        eom,
        list(v.scalars),
        list(v.vectors),
        list(v.matrices),
        list(system.inputs),
    )


def _config(system):
    v = system.variables
    return [*v.scalars, *v.vectors, *v.matrices]


def _velocities(system):
    from geomech.codegen.layout import _kinematic_symbols

    return [_kinematic_symbols(var)[1] for var in _config(system)]


def _split_energy(system):
    """``(kinetic addends, potential addends of L)``; the second list is ``−V``."""
    from geomech.viz.report import _terms

    vel = _velocities(system)
    T, minus_V = [], []
    for term in _terms(system.L):
        (T if any(term.has(w) for w in vel) else minus_V).append(term)
    return T, minus_V


def _lines(terms, names):
    from geomech.utils.printing import to_latex

    out = []
    for i, term in enumerate(terms):
        s = to_latex(term, names=names)
        if i > 0 and not s.startswith("-"):
            s = "+ " + s
        out.append(s)
    return out or ["0"]


def _negated(lines):
    """``−(a + b − c)`` written term by term."""
    out = []
    for i, s in enumerate(lines):
        s = s.removeprefix("+ ")
        s = s[1:] if s.startswith("-") else ("-" + s)
        if i > 0 and not s.startswith("-"):
            s = "+ " + s
        out.append(s)
    return out


def _term_lines(expr, names):
    from geomech.viz.report import _term_lines

    return _term_lines(expr, names)


def _math(lhs, lines):
    from geomech.viz.report import _aligned

    return "$$\n" + _aligned(lhs, lines) + "\n$$\n"


def energy_block(system, eom) -> str:
    names = _names(system, eom)
    T, minus_V = _split_energy(system)
    parts = [
        "Kinetic energy (the terms of $L$ with a velocity):\n",
        _math("T", _lines(T, names)),
        "Potential energy (the other terms, with the sign changed):\n",
        _math("V", _negated(_lines(minus_V, names)) if minus_V else ["0"]),
        "Lagrangian, $L = T - V$, as the library holds it:\n",
        _math("L", _term_lines(system.L, names)),
        "Virtual work of the inputs:\n",
        _math("\\delta W", _term_lines(system.dW, names)),
    ]
    return "\n".join(parts)


def eom_block(system, eom) -> str:
    from geomech.utils.printing import to_latex

    names = _names(system, eom)
    eqs = []
    for var in _config(system):
        _, eq = eom[str(var.get_variation_vector())]
        label = (
            f"Coefficient of ${to_latex(var.get_variation_vector(), names=names)}$, the "
            f"equation of ${to_latex(var, names=names)}$"
        )
        eqs.append((label, _term_lines(eq, names)))
    return _equations(eqs)


def _equations(eqs) -> str:
    """``0 = …`` for each ``(label, lines)``, folded one per equation when long."""
    fold = sum(len(lines) for _, lines in eqs) > FOLD_LINES
    parts = []
    for label, lines in eqs:
        if fold:
            # raw HTML: MyST leaves $...$ alone there, and Sphinx's MathJax only typesets
            # elements of class "math", so write the inline math as Sphinx itself does
            summary = re.sub(
                r"\$(.+?)\$",
                lambda m: (
                    '<span class="math notranslate nohighlight">'
                    f"\\({html.escape(m.group(1))}\\)</span>"
                ),
                label,
            )
            parts.append(
                f"<details>\n<summary>{summary} ({len(lines)} terms)</summary>\n\n"
                + _math("0", lines)
                + "\n</details>\n"
            )
        else:
            parts.append(f"{label}:\n\n" + _math("0", lines))
    return "\n".join(parts)


def _bound(x: float) -> str:
    """A power of ten two decades above *x*, so that the page stays true across
    platforms and NumPy builds."""
    return f"10^{{{int(np.floor(np.log10(max(x, 1e-16)))) + 2}}}"


def energy_drift(system, eom, page: Page) -> tuple[float, float]:
    """``(E(0), max |E(t) − E(0)|)`` of an unforced run from the page's initial state."""
    E = np.asarray(simulate_page(system, eom, page, forced=False).energy(), dtype=float)
    return float(E[0]), float(np.max(np.abs(E - E[0])))


def checks_block(system, eom, page: Page) -> str:
    from geomech.numeric import check_eom
    from geomech.systems._base import compared_text

    err = check_eom(system.L, system.dW, system.variables, eom, inputs=system.inputs)
    E0, drift = energy_drift(system, eom, page)
    ref = system.reference
    rows = [
        "| Check | Result |",
        "|---|---|",
        "| `check_eom`: the equations against a numeric Euler–Lagrange evaluation at random "
        f"states | below ${_bound(err)}$ |",
        f"| Energy $H = T + V$, {page.t_end:g} s unforced from the animation's initial "
        "state, RKMK4 with the default step: "
        f"$\\max_t \\lvert H(t) - H(0) \\rvert$ | below ${_bound(drift)}$, "
        f"with $H(0) = {E0:.3g}$ |",
    ]
    if page.indexed:
        gen_err, l_err = indexed_checks(page)
        rows += [
            f"| Equations for every $n$ above, at $n = {page.size}$: `check_eom` against the "
            f"Euler–Lagrange evaluation of the generic $L$, $\\delta W$ at $n = {page.size}$ "
            f"| below ${_bound(gen_err)}$ |",
            f"| Generic $L$ at $n = {page.size}$ against the library system's $L$, at random "
            f"states | relative difference below ${_bound(l_err)}$ |",
        ]
    if ref is not None and ref.equation:
        rows.append(
            f"| Closed form, {compared_text(ref)} | equal numerically at random states "
            "(`tests/test_systems_book.py`) |"
        )
    else:
        rows.append(
            "| Closed form | none compared; the equations are verified numerically (`check_eom`) |"
        )
    return "\n".join(rows) + "\n"


# ------------------------------------------------------------- any n (indexed)


def indexed_model(page: Page):
    from geomech.indexed.systems import MODELS as INDEXED_MODELS

    return INDEXED_MODELS[page.indexed]()


#: ``geomech.indexed`` names → the names of the rest of the page (``\Omega_{R_{i}}`` → Ωᵢ).
_TIDY = [
    (re.compile(r"\\dot\{\\(Omega|omega)_\{[Rq]_\{(\w+)\}\}\}"), r"\\dot{\\\1_{\2}}"),
    (re.compile(r"\\(Omega|omega|xi|eta)_\{[Rq]_\{(\w+)\}\}"), r"\\\1_{\2}"),
    (re.compile(r"(?<![\\\w])rho_\{"), r"\\rho_{"),
    (re.compile(r"(?<![\\\w])e3(?!\w)"), r"e_3"),
]


def _tidy(s: str) -> str:
    for pat, rep in _TIDY:
        s = pat.sub(rep, s)
    return s


def _addends(expr):
    from geomech.core.transformations.utils import addends

    return list(addends(expr))


def _summand(term, free: str = "i") -> str:
    """``Σ_j … body`` of a summand whose body is one product: no brackets, and a leading
    minus sign of the body written before the sums."""
    from geomech.indexed import FREE, rename_indices
    from geomech.indexed.printing import _layout
    from geomech.utils.printing import to_latex

    _, names, bounds, guards = _layout(term, {FREE: free})
    prefix = ""
    rel = {"≤": r"\le", "<": "<"}
    for a, op, b in guards:
        prefix += rf"[{a} {rel[op]} {b}]\, "
    for letter, lows, highs in bounds:
        lo = lows[0] if len(lows) == 1 else rf"\max({', '.join(lows)})"
        hi = highs[0] if len(highs) == 1 else rf"\min({', '.join(highs)})"
        prefix += rf"\sum_{{{letter}={lo}}}^{{{hi}}} "
    body = _tidy(to_latex(rename_indices(term.body, names)))
    sign = ""
    if body.startswith("-"):
        sign, body = "-", body[1:].lstrip()
    return sign + prefix + body


def _indexed_lines(terms) -> list[str]:
    """One line per addend of each summand, keeping the summand's sums:
    ``Σ_j (a + b)`` is written ``Σ_j a + Σ_j b``."""
    from geomech.indexed import Term

    out = []
    for term in terms:
        for body in _addends(term.body):
            s = _summand(Term(term.bound, term.cons, body))
            if out and not s.startswith("-"):
                s = "+ " + s
            out.append(s)
    return out or ["0"]


def _expr_lines(expr) -> list[str]:
    from geomech.indexed import IndexedExpr

    if isinstance(expr, IndexedExpr):
        return _indexed_lines(expr.terms)
    from geomech.utils.printing import to_latex

    lines = []
    for a in _addends(expr):
        s = _tidy(to_latex(a))
        lines.append(s if not lines or s.startswith("-") else "+ " + s)
    return lines


def indexed_block(page: Page) -> str:
    from geomech.indexed.eom import _variation
    from geomech.utils.printing import to_latex

    model = indexed_model(page)
    eom = model.eom()
    parts = [
        f"Lagrangian for every $n$, as `geomech.indexed.systems.{page.indexed}()` holds it "
        "(the velocities expanded, so the inner sums are the inertia couplings):\n",
        _math("L", _expr_lines(model.lagrangian)),
        "Virtual work:\n",
        _math("\\delta W", _expr_lines(model.inf_work)),
        "Equations from `compute_indexed_eom`, derived once for every $n$; $i$ is the "
        "body, $1 \\le i \\le n$:\n",
    ]
    eqs = []
    for eq in eom.values():
        if eq.family is not None:
            sym = _tidy(to_latex(eq.family[1])).replace("_{1}", "_i")
            var = _tidy(to_latex(eq.variation)).replace("ι0", "i")
            label = f"Coefficient of ${var}$, the equation of ${sym}$"
        else:
            owner = next(
                v for v in eom.variables._globals() if str(_variation(v)) == str(eq.variation)
            )
            var, sym = _tidy(to_latex(eq.variation)), _tidy(to_latex(owner))
            label = f"Coefficient of ${var}$, the equation of ${sym}$"
        eqs.append((label, _indexed_lines(eq.terms)))
    parts.append(_equations(eqs))
    return "\n".join(parts)


def _norm(name: str) -> str:
    """A leaf name with the decoration removed: ``\\rho_1``, ``rho_{1}`` → ``rho1``."""
    return re.sub(r"[\\{}_ ]", "", name)


def _leaves(expr, out: dict) -> dict:
    """``{str(leaf): leaf}`` of *expr*, a velocity ``TimeDerivative(x)`` being a leaf."""
    from geomech.codegen.layout import is_literal_leaf
    from geomech.core.operations.calculus import TimeDerivative

    if isinstance(expr, TimeDerivative) or not getattr(expr, "nodes", None):
        if not is_literal_leaf(expr):
            out[str(expr)] = expr
    else:
        for c in expr.nodes:
            _leaves(c, out)
    return out


def _random_value(leaf, rng):
    from geomech.core.base import SO3, Matrix, Scalar

    if isinstance(leaf, SO3):
        Q, R = np.linalg.qr(rng.normal(size=(3, 3)))
        return Q * np.sign(np.diag(R))
    if isinstance(leaf, Matrix):
        A = rng.normal(size=(3, 3))
        return A + A.T
    if isinstance(leaf, Scalar):
        return float(rng.uniform(0.5, 1.5))
    v = rng.normal(size=(3, 1))
    return v / np.linalg.norm(v)


def indexed_checks(page: Page, seeds=(0, 1, 2)) -> tuple[float, float]:
    """``(check_eom of the generic equations at n = size, max relative difference of the
    generic L at n = size and the library system's L at random states)``."""
    from geomech.numeric import check_eom, evaluate

    model, system = indexed_model(page), page.system()
    L, dW, variables, inputs = model.unroll(page.size)
    gen_err = check_eom(L, dW, variables, model.eom().instantiate(page.size), inputs=inputs)
    lib = _leaves(system.L, {})
    gen = _leaves(L, {})
    by_norm = {_norm(k): k for k in lib}
    missing = [k for k in gen if _norm(k) not in by_norm]
    if missing or len(gen) != len(lib):
        raise ValueError(f"{page.name}: generic and library leaves differ: {missing}")
    diff = 0.0
    for seed in seeds:
        rng = np.random.default_rng(seed)
        env_lib = {k: _random_value(v, rng) for k, v in lib.items()}
        env_gen = {k: env_lib[by_norm[_norm(k)]] for k in gen}
        a, b = float(evaluate(system.L, env_lib)), float(evaluate(L, env_gen))
        diff = max(diff, abs(a - b) / max(abs(a), 1e-300))
    return float(gen_err), diff


# ---------------------------------------------------------------- simulation


def _plant(system, eom, page: Page | None = None):
    from geomech import to_standard_form
    from geomech.codegen import StateLayout, load_module, module_source

    sf = to_standard_form(eom, system.variables, system.inputs)
    layout = system.layout()
    if page is not None and page.scene is not None:
        extra = list(page.scene[0])
        layout = StateLayout(system.variables, list(system.inputs), [*layout.params, *extra])
    return sf, layout, load_module(module_source(sf, layout, lagrangian=system.L), system.name)


def simulate_page(system, eom, page: Page, forced: bool = True):
    """The page's run: its initial state, and its inputs when *forced*."""
    from geomech.simulation import simulate

    _, layout, plant = _plant(system, eom, page)
    u = page.inputs(system, layout) if forced else None
    x0 = page.initial_state(system, layout)
    p = dict(system.params)
    if page.scene is not None:
        p.update({str(k): v for k, v in page.scene[0].items()})
    return simulate(plant, x0, page.t_end, p, u=u)


def write_reports(outdir, pages=None) -> list[Path]:
    """Write one interactive report per written page to *outdir*; return the paths. A
    system without a scene gets no report."""
    from geomech.viz import Report

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    paths = []
    for page in existing() if pages is None else pages:
        system = page.system()
        eom = system.derive()
        sf, layout, _ = _plant(system, eom, page)
        if page.scene is not None:
            scene = page.scene[1](layout, system.symbols)
        else:
            scene = system.scene(layout)
        if scene is None:
            continue
        report = Report(
            page.title,
            system.as_tuple(),
            eom=eom,
            standard_form=sf,
            scene=scene,
            trajectory=simulate_page(system, eom, page),
            units={"g": "m/s²"} if "g" in system.params else {},
        )
        paths.append(report.write(outdir / f"{page.name}.html"))
    return paths


# ---------------------------------------------------------------- the pages


def blocks(page: Page) -> dict[str, str]:
    system = page.system()
    eom = system.derive()
    out = {
        "energy": energy_block(system, eom),
        "eom": eom_block(system, eom),
        "checks": checks_block(system, eom, page),
    }
    if page.indexed:
        out["indexed"] = indexed_block(page)
    return out


def render(text: str, new: dict[str, str]) -> str:
    """*text* with every generated block replaced by ``new[block]``."""
    seen = set()

    def sub(m):
        seen.add(m.group(2))
        return m.group(1) + new[m.group(2)] + m.group(4)

    out = BLOCK.sub(sub, text)
    missing = set(new) - seen
    if missing:
        raise ValueError(f"page has no <!-- model:… --> block for {sorted(missing)}")
    return out


def index_page(pages=None) -> str:
    """``docs/models/index.md``: one toctree per group, in :data:`MODELS` order."""
    pages = existing() if pages is None else pages
    out = [
        "# Models\n",
        "Each page animates a system, draws it and names its variables, and shows the",
        "kinetic and potential energy, the Lagrangian and virtual work, the equations",
        "`compute_eom` derives from them, the code, and the checks the equations pass.",
        "Five pages of systems with $n$ bodies also derive the equations once for every $n$",
        "([sums over n bodies](../indexed.md)) and then substitute $n = 5$. The math on",
        "every page is generated from the library and tested against it",
        "(`scripts/model_pages.py`, `tests/test_model_pages.py`).\n",
    ]
    for group in dict.fromkeys(p.group for p in MODELS):
        names = [p.name for p in pages if p.group == group]
        if names:
            out += ["```{toctree}", f":caption: {group}", ":maxdepth: 1", "", *names, "```\n"]
    return "\n".join(out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--reports", metavar="DIR", help="write the interactive reports to DIR")
    ap.add_argument("--only", nargs="+", metavar="NAME", help="these pages only")
    args = ap.parse_args(argv)
    pages = [BY_NAME[n] for n in args.only] if args.only else existing()
    if args.reports:
        for path in write_reports(args.reports, pages):
            print(f"wrote {path} ({path.stat().st_size / 1e3:.0f} kB)")
        return
    for page in pages:
        path = PAGES / f"{page.name}.md"
        path.write_text(render(path.read_text(), blocks(page)))
        print(f"updated {path.relative_to(ROOT)}")
    (PAGES / "index.md").write_text(index_page())
    print("updated docs/models/index.md")


if __name__ == "__main__":
    main()
