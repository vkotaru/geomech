"""Closed-loop LQR examples on the systems of ``geomech.systems``.

Each example stabilises an equilibrium with an LQR gain designed from the variation-based
linearization (``linearize``) and simulates the nonlinear generated model on the manifold
(RKMK4, dt = 1 ms). See ``control_helper.py`` for the chain and ``README.md`` for the
results.

    python examples/control/run.py spherical_pendulum_inverted
    python examples/control/run.py --all          # every example, then the README table
    python examples/control/run.py --list

Needs scipy (``pip install "geomech[sim]"``).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from control_helper import (  # noqa: E402
    DOWN,
    FLIP,
    RESULTS,
    UP,
    Example,
    rotate,
    run,
    tangent,
    tilt,
    weights,
)

I3 = np.eye(3)


def rest(**slots):
    """Equilibrium function returning fixed slot values."""
    return lambda sys_, layout: slots


def perturbed(**slots):
    """Perturbation function: the equilibrium with the given slots replaced (a callable
    value is applied to the equilibrium value of its slot)."""

    def f(sys_, layout, x_eq):
        vals = layout.unpack(x_eq)
        for k, v in slots.items():
            vals[k] = v(vals[k]) if callable(v) else v
        return layout.pack(**vals)

    return f


def rot_to(d) -> np.ndarray:
    """The rotation of least angle taking e₃ to the unit vector *d*."""
    from geomech.numeric import expm_so3

    d = np.asarray(d, float) / np.linalg.norm(d)
    axis = np.cross(UP, d)
    s = np.linalg.norm(axis)
    if s < 1e-15:
        return I3.copy()
    return expm_so3(axis / s * np.arctan2(s, d @ UP))


def team_equilibrium(n, tilt_rad, rigid):
    """Symmetric team equilibrium: cable i leaves the load at the angle *tilt_rad* from the
    vertical towards the azimuth 2π(i−1)/n; each quadrotor's thrust balances its weight and
    the cable tension ``T = m_L g / (n cos α)``."""

    def f(sys_, layout):
        p = sys_.params
        g, mL = p["g"], p["m_L"]
        T = mL * g / (n * np.cos(tilt_rad))
        vals = {"x_L": np.zeros(3)}
        if rigid:
            vals["R_L"] = I3
        for i in range(1, n + 1):
            th = 2 * np.pi * (i - 1) / n
            q = np.array(
                [-np.sin(tilt_rad) * np.cos(th), -np.sin(tilt_rad) * np.sin(th), -np.cos(tilt_rad)]
            )
            vals[f"q{i}"] = q
            vals[f"R{i}"] = rot_to(p[f"m{i}"] * g * UP - T * q)
        return vals

    return f


#: The commanded sun-pointing angle of the solar-array panels (deg): both panel normals
#: turned from b₃ about b₁ by this angle, towards −b₂.
SUN_ANGLE = 40.0


def wheels_at_x0(**slots):
    """Equilibrium function: the S¹ wheel directions of the system's ``x0`` (any direction
    on each wheel's circle is an equilibrium) and the given slots."""

    def f(sys_, layout):
        vals = layout.unpack(sys_.x0)
        return {**{k: vals[k] for k in ("q_1", "q_2", "q_3")}, **slots}

    return f


def sun_pointing(angle_deg):
    """Equilibrium of ``spacecraft_solar_arrays``: ``R = I``, the wheels at ``x0`` and both
    panel normals at *angle_deg* about b₁ from b₃."""

    def f(sys_, layout):
        n = tilt(UP, angle_deg, sys_.params["b1"])
        return wheels_at_x0(R=I3, n_1=n, n_2=n)(sys_, layout)

    return f


def _axis(*v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v)


def units(system: str, layout_inputs) -> dict:
    """The unit of each input of *system*: forces in N, moments and torques in N·m."""
    force = {
        "cart_spherical_pendulum": {"f"},
        "cart_planar_pendulum": {"f"},
        "rigid_body": {"F"},
        "two_rigid_bodies_spherical_joint": {"F"},
        "elastic_spherical_pendulum": {"f_l"},
    }
    out = {}
    for key in layout_inputs:
        if system.startswith(("quadrotor", "multi_quadrotor")):
            out[key] = "N" if key.startswith("f") else "N·m"
        else:
            out[key] = "N" if key in force.get(system, set()) else "N·m"
    return out


EXAMPLES = [
    Example(
        "spherical_pendulum_inverted",
        "spherical_pendulum",
        rest(q=UP),
        perturbed(
            q=tilt(UP, 30, (1, 1, 0)),
            omega_q=lambda w: tangent(tilt(UP, 30, (1, 1, 0)), [0.5, -0.3, 0]),
        ),
        weights(10, 1),
        1.0,
        8.0,
        description="Inverted spherical pendulum (q = e₃), moment f about the pivot; 30° tilt.",
    ),
    Example(
        "double_spherical_pendulum_inverted",
        "double_spherical_pendulum",
        rest(q1=UP, q2=UP),
        perturbed(q1=tilt(UP, 10, (1, 0, 0)), q2=tilt(UP, 15, (0, 1, 0))),
        weights(10, 1),
        1.0,
        10.0,
        description="Both links upright; moments f₁, f₂ on the links; 10° and 15° tilts.",
    ),
    Example(
        "n_link_spherical_chain_inverted",
        "n_link_spherical_chain",
        rest(q1=UP, q2=UP, q3=UP),
        perturbed(
            q1=tilt(UP, 8, (1, 0, 0)), q2=tilt(UP, 10, (0, 1, 0)), q3=tilt(UP, 12, (1, 1, 0))
        ),
        weights(10, 1),
        1.0,
        10.0,
        kwargs={"n_links": 3},
        description="Three-link spherical chain upright; 8°, 10°, 12° tilts.",
    ),
    Example(
        "elastic_spherical_pendulum_inverted",
        "elastic_spherical_pendulum",
        # upright: the spring carries the weight, k (l − l₀) = −m g → l = l₀ − m g / k
        lambda s, lay: {
            "q": UP,
            "l": s.params["l0"] - s.params["m"] * s.params["g"] / s.params["k"],
        },
        perturbed(l=lambda l: l + 0.1, q=tilt(UP, 20, (0, 1, 0))),
        weights(10, 1),
        1.0,
        8.0,
        description="Upright elastic pendulum at the compressed spring length"
        " l = l₀ − m g / k; inputs f_l (along the rod) and M_q (moment); 20° tilt, +0.1 m.",
    ),
    Example(
        "cart_planar_pendulum_inverted",
        "cart_planar_pendulum",
        rest(s=0.0, q=UP),
        perturbed(s=0.3, q=tilt(UP, 10, (0, 1, 0))),
        weights(10, 1),
        0.1,
        6.0,
        description="Planar cart-pole upright with the rail force only: the start, weights"
        " and horizon of cart_spherical_pendulum_inverted_in_plane. Fully controllable"
        " (rank 4 of 4), so the library's LQR applies.",
    ),
    Example(
        "cart_spherical_pendulum_inverted_in_plane",
        "cart_spherical_pendulum",
        rest(s=0.0, q=UP),
        perturbed(s=0.3, q=tilt(UP, 10, (0, 1, 0))),
        weights(10, 1),
        0.1,
        6.0,
        description="3-D cart-pole upright with the rail force only. NOT stabilizable: the"
        " pendulum mode out of the rail's vertical plane is uncontrollable and unstable."
        " LQR on the controllable (in-plane) subspace; the perturbation is in-plane.",
        notes="Converges only because the out-of-plane perturbation is exactly zero; a"
        " horizon of a few seconds keeps round-off from growing visibly (rate √(g/l)).",
    ),
    Example(
        "cart_spherical_pendulum_out_of_plane",
        "cart_spherical_pendulum",
        rest(s=0.0, q=UP),
        perturbed(s=0.3, q=tilt(tilt(UP, 10, (0, 1, 0)), 1, (1, 0, 0))),
        weights(10, 1),
        0.1,
        4.0,
        description="Negative example: the same controller with a 1° out-of-plane tilt"
        " diverges (the uncontrollable mode is unstable).",
        notes="Expected to fail; kept to show that the stabilizability check is real.",
    ),
    Example(
        "rigid_pendulum_inverted",
        "rigid_pendulum",
        rest(R=FLIP),
        perturbed(R=rotate(FLIP, 25, (1, 2, 0.5)), Omega_R=np.array([0.2, -0.1, 1.0])),
        weights(10, 1),
        1.0,
        8.0,
        description="3-D rigid pendulum inverted (R = half turn about e₁); body moment M;"
        " 25° rotation and a spin about the link.",
    ),
    Example(
        "double_rigid_pendulum_inverted",
        "double_rigid_pendulum",
        rest(R1=FLIP, R2=FLIP),
        perturbed(R1=rotate(FLIP, 10, (1, 0, 0)), R2=rotate(FLIP, 15, (0, 1, 0.3))),
        weights(10, 1),
        1.0,
        8.0,
        description="Both links upright; body moments M₁, M₂; 10° and 15° rotations.",
    ),
    Example(
        "n_link_rigid_chain_inverted",
        "n_link_rigid_chain",
        rest(R1=FLIP, R2=FLIP, R3=FLIP),
        perturbed(
            R1=rotate(FLIP, 8, (1, 0, 0)),
            R2=rotate(FLIP, 10, (0, 1, 0)),
            R3=rotate(FLIP, 12, (1, 1, 0.5)),
        ),
        weights(10, 1),
        100.0,
        10.0,
        kwargs={"n_links": 3},
        description="Three rigid links upright, joint torques τ₁ (ground–link 1), τ₂, τ₃;"
        " 8°, 10°, 12° rotations.",
    ),
    Example(
        "free_rigid_body_attitude",
        "free_rigid_body",
        rest(R=I3),
        perturbed(R=rotate(I3, 110, (1, -2, 1))),
        weights(10, 1),
        1.0,
        10.0,
        description="Attitude stabilisation to R = I from a 110° rotation at rest; body moment M.",
    ),
    Example(
        "spacecraft_reaction_wheels_slew",
        "spacecraft_reaction_wheels",
        rest(R=I3, theta_1=0.0, theta_2=0.0, theta_3=0.0),
        perturbed(R=rotate(I3, 60, (1, 1, 1))),
        weights(1000, 1, delta_theta=0.1),
        1e4,
        25.0,
        drop=("delta_theta_1", "delta_theta_2", "delta_theta_3"),
        description="Rest-to-rest 60° slew with the three wheel torques only (R = 1e4:"
        " peak wheel torque about 0.19 N·m).",
        notes="Wheel angles are cyclic (zero columns of A and B) and dropped from the design;"
        " the three components of the conserved angular momentum are uncontrollable"
        " (eigenvalue 0) and LQR is designed on the controllable subspace. The rest-to-rest"
        " (zero total momentum) start is essential: internal torques cannot change the"
        " total angular momentum, so a start with momentum cannot come to rest at R = I.",
    ),
    Example(
        "spacecraft_reaction_wheels_s1_slew",
        "spacecraft_reaction_wheels_s1",
        wheels_at_x0(R=I3),
        perturbed(R=rotate(I3, 60, (1, 1, 1))),
        weights(1000, 1, delta_omega_q=0.1),
        1e4,
        25.0,
        drop=("gamma_q_1", "gamma_q_2", "gamma_q_3"),
        description="The rest-to-rest 60° slew of spacecraft_reaction_wheels_slew on the S¹"
        " model (wheel directions qᵢ on circles): same weights, same input weight, same"
        " horizon.",
        notes="The S¹ wheel angles are cyclic like θᵢ and dropped; the design and the closed"
        " loop are those of the R³ model up to the order of the coordinates, so t_s and ν"
        " match it.",
    ),
    Example(
        "spacecraft_solar_arrays_slew",
        "spacecraft_solar_arrays",
        sun_pointing(SUN_ANGLE),
        perturbed(R=rotate(I3, 60, (1, 1, 1)), n_1=UP, n_2=UP),
        weights(1000, 1, delta_omega_q=0.1, gamma_n=10, delta_omega_n=100),
        [1e4, 1e4, 1e4, 1e5, 1e5],
        30.0,
        drop=("gamma_q_1", "gamma_q_2", "gamma_q_3"),
        description="Rest-to-rest 60° slew of the spacecraft with two solar-array panels,"
        f" wheel torques τ₁–τ₃ and array-drive torques μ₁, μ₂, ending with both panel"
        f" normals at {SUN_ANGLE:g}° about b₁ from b₃ (the commanded sun-pointing angle);"
        " the panels start at 0°. The drive torques are weighted 10× the wheel torques"
        " (R = 1e5 against 1e4) and the panel rates 100, which keeps the drive torque"
        " near 0.01 N·m.",
        notes="Wheel angles are cyclic and dropped; the panel angles are not (they are"
        " controlled through μ). The total angular momentum is conserved under the internal"
        " wheel and drive torques: three uncontrollable modes at eigenvalue 0, LQR on the"
        " 10-D controllable subspace, rest-to-rest start.",
    ),
    Example(
        "biped_3d_balance",
        "biped_3d",
        rest(R1=I3, R2=rotate(I3, 20, (0, 1, 0))),
        perturbed(R1=rotate(I3, 5, (1, 1, 0)), R2=rotate(I3, 25, (0, 1, 0.2))),
        weights(100, 1),
        1.0,
        10.0,
        description="Biped balance: stance leg vertical, swing leg held 20° forward by the"
        " hip torque, the ankle torque balancing the stance leg; the stance leg starts 5° and"
        " the swing leg 6.6° from the equilibrium.",
        notes="Contact is not modelled: the stance foot is a spherical joint that can apply"
        " any ankle torque, and the swing foot may pass through the ground.",
    ),
    Example(
        "rigid_body_pose",
        "rigid_body",
        rest(x=np.zeros(3), R=I3),
        perturbed(x=np.array([1.0, -0.5, 0.5]), R=rotate(I3, 90, (0, 0, 1))),
        weights(10, 1),
        1.0,
        10.0,
        description="6-DOF position and attitude regulation, force F and moment M; 1.2 m and 90°.",
    ),
    Example(
        "quadrotor_hover",
        "quadrotor",
        rest(x=np.zeros(3), R=I3),
        perturbed(x=np.array([1.5, -1.0, 0.5]), R=rotate(I3, 60, (0, 0, 1))),
        weights(10, 1),
        1.0,
        10.0,
        description="Hover at the origin from a 1.9 m offset and a 60° yaw error.",
    ),
    Example(
        "quadrotor_offset_com_hover",
        "quadrotor_offset_com",
        rest(x=np.zeros(3), R=I3),
        perturbed(x=np.array([1.0, 0.5, -0.5]), R=rotate(I3, 45, (0, 0, 1))),
        weights(10, 1),
        1.0,
        10.0,
        description="Hover with the centre of mass offset from the thrust point: level"
        " attitude, thrust m g, and a constant moment M = m g ρ × e₃ (from u_eq) holding it.",
    ),
    Example(
        "quadrotor_suspended_load_hover",
        "quadrotor_suspended_load",
        rest(x_L=np.zeros(3), q=DOWN, R=I3),
        perturbed(
            x_L=np.array([0.6, -0.4, 0.3]),
            q=tilt(DOWN, 28.6, (0, 1, 0)),
            R=rotate(I3, 45, (0, 0, 1)),
            Omega_R=np.array([0.0, 0.0, 1.5]),
        ),
        weights(1, 1),
        1.0,
        10.0,
        description="The documented example: load at the origin, cable hanging; the load"
        " starts 0.8 m away, swung 0.5 rad, the quadrotor yawed 45° and spinning.",
    ),
    Example(
        "quadrotor_elastic_load_hover",
        "quadrotor_elastic_load",
        lambda s, lay: {
            "x_L": np.zeros(3),
            "q": DOWN,
            "R": I3,
            "l": s.params["l0"] + s.params["m_L"] * s.params["g"] / s.params["k"],
        },
        perturbed(
            x_L=np.array([0.4, -0.2, 0.15]),
            q=tilt(DOWN, 12, (1, 0, 0)),
            l=lambda l: l + 0.03,
            R=rotate(I3, 20, (0, 0, 1)),
        ),
        weights(1, 1, delta_l_dot=10.0),
        1.0,
        15.0,
        description="Hover with an elastic cable stretched by the load weight"
        " (l = l₀ + m_L g / k); 0.47 m offset, 12° swing, +3 cm stretch, 20° yaw. The"
        " cable-spring mode near 28 rad/s is lightly damped: damping ratio 0.015 with"
        " uniform weights, 0.039 with the stretch rate weighted 10 (used here). Heavier"
        " stretch weights damp it more but need thrust deviations of 12–28 N (hover"
        " 17.7 N); these weights keep the peak deviation at 2.5 N (14 % of hover).",
    ),
    Example(
        "quadrotor_flexible_cable_hover",
        "quadrotor_flexible_cable",
        rest(x_Q=np.zeros(3), R=I3, q1=DOWN, q2=DOWN, q3=DOWN),
        perturbed(
            x_Q=np.array([0.3, -0.2, 0.1]),
            q1=tilt(DOWN, 6, (1, 0, 0)),
            q2=tilt(DOWN, 8, (0, 1, 0)),
            q3=tilt(DOWN, 6, (1, 1, 0)),
        ),
        weights(1, 1),
        [1.0, 10.0, 10.0, 10.0],
        40.0,
        kwargs={"n_links": 3},
        description="Hover with a three-link cable hanging straight down; 0.37 m offset and"
        " 6–8° link swings. Moments weighted R = 10 to keep the peak moment below 1 N·m;"
        " the two cable modes near 25 rad/s are left lightly damped (damping ratio 0.007,"
        " LQR's choice for these weights), so the horizon is 40 s.",
    ),
    Example(
        "multi_quadrotor_point_mass_hover",
        "multi_quadrotor_point_mass",
        team_equilibrium(3, 0.5, rigid=False),
        perturbed(x_L=np.array([0.3, -0.2, 0.2]), q1=lambda q: tilt(q, 8, (0, 1, 0))),
        weights(1, 1),
        1.0,
        12.0,
        kwargs={"n_quads": 3},
        description="Three quadrotors hold a point mass with cables 0.5 rad from vertical;"
        " each thrust balances its weight and the cable tension m_L g / (3 cos 0.5).",
    ),
    Example(
        "multi_quadrotor_rigid_load_hover",
        "multi_quadrotor_rigid_load",
        team_equilibrium(3, 0.3, rigid=True),
        perturbed(
            x_L=np.array([0.3, -0.2, 0.2]),
            R_L=rotate(I3, 15, (0, 0, 1)),
            q2=lambda q: tilt(q, 6, (1, 0, 0)),
        ),
        weights(1, 1),
        1.0,
        12.0,
        kwargs={"n_quads": 3},
        description="Three quadrotors hold a rigid plate level; 0.4 m offset, 15° yaw of"
        " the plate, one cable swung 6°.",
    ),
    Example(
        "two_rigid_bodies_spherical_joint_hold",
        "two_rigid_bodies_spherical_joint",
        rest(x=np.zeros(3), R1=I3, R2=rotate(I3, 90, (0, 1, 0))),
        perturbed(
            x=np.array([0.3, 0.2, -0.2]),
            R1=rotate(I3, 20, (0, 0, 1)),
            R2=rotate(I3, 70, (0.2, 1, 0)),
        ),
        weights(10, 1),
        1.0,
        8.0,
        description="Body 1 held level in mid-air with body 2 hanging below the joint;"
        " force F on body 1 and moments M₁, M₂ (fully actuated); u_eq holds both weights.",
    ),
]

BY_NAME = {ex.name: ex for ex in EXAMPLES}


def _set_units():
    from control_helper import build

    for ex in EXAMPLES:
        if not ex.units:
            layout = build(ex.system, ex.kwargs)[1]
            ex.units = units(ex.system, [slot.key for slot in layout.input_slots])


def _fmt(x, spec=".3g"):
    return "—" if x is None else format(x, spec)


def _peaks(m) -> str:
    """The largest input deviation per unit: ``"f 2.3 N, M 0.5 N·m"`` style, by unit."""
    if not m.get("peak_inputs"):
        return "—"
    by_unit: dict = {}
    for v in m["peak_inputs"].values():
        by_unit[v["unit"]] = max(by_unit.get(v["unit"], 0.0), v["peak_deviation"])
    return ", ".join(f"{val:.3g} {unit}" for unit, val in sorted(by_unit.items()))


def table(results: list[dict]) -> str:
    """The README tables: stabilizable examples, then the ones that are not."""
    head = (
        "| Example | n | rank | max Re λ | min ζ | e₀ pos. (m) | e₀ att. (rad) | t_s |"
        " peak ‖u − u_eq‖∞ | drift | converged |\n"
        "|---|---|---|---|---|---|---|---|---|---|---|\n"
    )

    def row(m):
        e0 = m["initial_error"]
        ts = "—" if m.get("settling_time") is None else f"{m['settling_time']:.2f} s"
        ok = "yes" if m["converged"] else ("diverged" if m.get("failure") else "no")
        mark = " †" if m["uncontrollable_eigenvalues"] else ""
        return (
            f"| `{m['name']}`{mark} | {m['n_reduced']} | {m['controllability_rank']}"
            f" | {m['closed_loop_max_real']:.3g} | {m['closed_loop_min_damping']:.2f}"
            f" | {_fmt(e0['position'])} | {_fmt(e0['attitude'])} | {ts} | {_peaks(m)}"
            f" | {_fmt(m.get('constraint_drift'), '.1e')} | {ok} |\n"
        )

    def unstable(m):
        return any(re > 1e-9 for re, _ in m["uncontrollable_eigenvalues"])

    good = [m for m in results if not unstable(m)]
    bad = [m for m in results if unstable(m)]
    out = head + "".join(row(m) for m in good)
    if any(m["uncontrollable_eigenvalues"] for m in good):
        out += (
            "\n† Not stabilizable in the strict sense: uncontrollable modes with eigenvalue 0"
            " (conserved quantities). LQR on the controllable subspace converges from initial"
            " states on the invariant subspace (for the spacecraft, zero total angular"
            " momentum).\n"
        )
    if bad:
        out += "\n**Not stabilizable** (uncontrollable modes listed; LQR on the controllable"
        out += " subspace):\n\n"
        out += "| Example | n | rank | uncontrollable eigenvalues | max Re λ (controlled) |"
        out += " e₀ att. (rad) | t_s | converged |\n|---|---|---|---|---|---|---|---|\n"
        for m in bad:
            unc = ", ".join(f"{re:+.3g}" for re, im in m["uncontrollable_eigenvalues"])
            ts = "—" if m.get("settling_time") is None else f"{m['settling_time']:.2f} s"
            ok = (
                "yes (in-plane start only)"
                if m["converged"]
                else ("diverged" if m.get("failure") else "no")
            )
            out += (
                f"| `{m['name']}` | {m['n_reduced']} | {m['controllability_rank']} | {unc}"
                f" | {m['closed_loop_max_real']:.3g} | {_fmt(m['initial_error']['attitude'])}"
                f" | {ts} | {ok} |\n"
            )
    return out


def write_table() -> str:
    """Regenerate ``results/table.md`` from the stored JSON files."""
    stored = [json.loads((RESULTS / f"{ex.name}.json").read_text()) for ex in EXAMPLES]
    text = table(stored)
    (RESULTS / "table.md").write_text(text)
    return text


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("names", nargs="*", help="examples to run")
    parser.add_argument("--all", action="store_true", help="run every example")
    parser.add_argument("--list", action="store_true", help="list the examples")
    args = parser.parse_args(argv)
    if args.list:
        for ex in EXAMPLES:
            print(f"{ex.name:42s} {ex.system}")
        return 0
    names = [ex.name for ex in EXAMPLES] if args.all else args.names
    if not names:
        parser.error("give example names, --all or --list")
    _set_units()
    for name in names:
        m = run(BY_NAME[name])
        ts = m["settling_time"]
        if m.get("failure"):
            print(
                f"{name:42s} rank {m['controllability_rank']}/{m['n_reduced']}  maxRe"
                f" {m['closed_loop_max_real']:.3f}  DIVERGED: {m['failure']}",
                flush=True,
            )
            continue
        print(
            f"{name:42s} res {m['equilibrium_residual']:.1e}  rank {m['controllability_rank']}"
            f"/{m['n_reduced']}  maxRe {m['closed_loop_max_real']:.3f}  zeta"
            f" {m['closed_loop_min_damping']:.2f}  nu_end {m['final_normalized_error']:.1e}"
            f"  ts {'—' if ts is None else f'{ts:.2f}'}  peak {_peaks(m)}"
            f"  drift {m['constraint_drift']:.1e}  {'OK' if m['converged'] else 'NOT CONVERGED'}",
            flush=True,
        )
    if args.all:
        write_table()
    return 0


if __name__ == "__main__":
    sys.exit(main())
