"""Export the control examples of geomech (examples/control/results) for the paper.

    python export.py <geomech>/examples/control/results

Writes, next to this script:
- nu_<example>.dat: time and normalized error nu(t), decimated from 100 Hz to 20 Hz;
- track_<example>.dat: time, position error (m), attitude error (rad), 20 Hz;
- table_rows.tex: one row per example of the control table, from the JSON metrics;
- macros.tex: \\ctrl{<example>}{<field>} for the catalog entries of the appendix
  (fields ts, peak, rank, n, maxre, zeta);
- hist_<example>.dat and lqr_<example>.tex: for each example with a time history, the
  norms of the position, attitude and velocity errors and the input deviations (20 Hz),
  and the weights and the plot shown in its catalog entry.
Only stored results are used; nothing is simulated here. The weights are read from the
example definitions (``run.EXAMPLES`` of examples/control), evaluated on the stored labels.
"""

import json
import re
import sys
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent
DECIMATE = 5  # 100 Hz -> 20 Hz

# paper names (Table I) and the order of the control table
ROWS = [
    ("point_mass_hold", "Point mass, hold"),
    ("spring_mass_rest", "Spring--mass, rest"),
    ("planar_pendulum_inverted", "Planar pendulum, inverted"),
    ("double_planar_pendulum_inverted", "Double planar pendulum, inverted"),
    ("cart_planar_pendulum_inverted", "Cart, planar pendulum, inverted"),
    ("spherical_pendulum_inverted", "Spherical pendulum, inverted"),
    ("particle_on_rotating_sphere_hold", "Particle, rotating sphere, hold"),
    ("double_spherical_pendulum_inverted", "Double spherical pendulum, inverted"),
    ("n_link_spherical_chain_inverted", "Spherical chain ($n = 3$), inverted"),
    ("elastic_spherical_pendulum_inverted", "Elastic spherical pendulum, inverted"),
    (
        "cart_spherical_pendulum_inverted_in_plane",
        "Cart, sph.\\ pendulum, in-plane start$^{a}$",
    ),
    (
        "cart_spherical_pendulum_out_of_plane",
        "Cart, sph.\\ pendulum, $1^\\circ$ out of plane$^{a}$",
    ),
    ("free_rigid_body_attitude", "Free rigid body, attitude"),
    ("heavy_top_inverted", "Heavy top, inverted"),
    ("rigid_pendulum_inverted", "Rigid pendulum, inverted"),
    ("double_rigid_pendulum_inverted", "Double rigid pendulum, inverted"),
    ("n_link_rigid_chain_inverted", "Rigid chain ($n = 3$), inverted"),
    ("two_rigid_bodies_spherical_joint_hold", "Two bodies, ball joint, hold"),
    ("rigid_body_pose", "Rigid body, pose"),
    ("spacecraft_reaction_wheels_slew", "Spacecraft, wheels in $\\RealR$, slew$^{b}$"),
    ("spacecraft_reaction_wheels_s1_slew", "Spacecraft, wheels on $S^1$, slew$^{b}$"),
    ("spacecraft_solar_arrays_slew", "Spacecraft, solar arrays, slew$^{b}$"),
    ("biped_3d_balance", "Biped, balance of a posture"),
    ("quadrotor_hover", "Quadrotor, hover"),
    ("quadrotor_offset_com_hover", "Quadrotor, offset c.o.m., hover"),
    ("quadrotor_suspended_load_hover", "Quadrotor, suspended load, hover"),
    ("quadrotor_elastic_load_hover", "Quadrotor, elastic cable, hover"),
    ("quadrotor_flexible_cable_hover", "Quadrotor, flexible cable ($n = 3$), hover"),
    ("multi_quadrotor_point_mass_hover", "3 quadrotors, point mass, hover"),
    ("multi_quadrotor_rigid_load_hover", "3 quadrotors, rigid body, hover"),
    ("multi_quadrotor_hose_hover", "2 quadrotors, hose ($n = 3$), hover"),
]


def sig(x, n=3):
    return f"{x:.{n}g}"


def zeta(z):
    return f"{z:.2f}" if z >= 0.1 else f"{z:.2g}"


def unit_tex(u):
    return {"N": "N", "N·m": "N\\,m"}.get(u, u)


def main(results):
    W = example_weights(results)
    results = Path(results)
    lines, macros = [], []
    for name, label in ROWS:
        d = json.loads((results / f"{name}.json").read_text())
        e0 = d["initial_error"]
        pos = "--" if e0.get("position") in (None, 0) else sig(e0["position"])
        att = "--" if e0.get("attitude") in (None, 0) else sig(e0["attitude"])
        ts = d.get("settling_time")
        ts = "--" if ts is None else f"{ts:.2f}"
        peaks = {}
        for v in (d.get("peak_inputs") or {}).values():
            u = v["unit"]
            peaks[u] = max(peaks.get(u, 0.0), v["peak_deviation"])
        peak = ", ".join(f"{sig(v)}\\,{unit_tex(u)}" for u, v in peaks.items()) or "--"
        if not d["converged"]:
            peak = "diverges"
        lines.append(
            f"{label} & {d['n_reduced']} & {d['controllability_rank']} & "
            f"${sig(d['closed_loop_max_real'])}$ & {zeta(d['closed_loop_min_damping'])} & "
            f"{pos} & {att} & {ts} & {peak} & {W[name][0]} \\\\"
        )
        for key, val in (
            ("ts", ts),
            ("peak", peak),
            ("rank", d["controllability_rank"]),
            ("n", d["n_reduced"]),
            ("maxre", sig(d["closed_loop_max_real"])),
            ("zeta", zeta(d["closed_loop_min_damping"])),
        ):
            macros.append(f"\\expandafter\\def\\csname ctrl@{name}@{key}\\endcsname{{{val}}}")
        npz = results / f"{name}.npz"
        if npz.exists():
            z = np.load(npz)
            t, nu = z["t"][::DECIMATE], z["normalized_error"][::DECIMATE]
            np.savetxt(
                OUT / f"nu_{name}.dat",
                np.column_stack([t, np.maximum(nu, 1e-12)]),
                header="t nu",
                comments="",
                fmt="%.6g",
            )
    (OUT / "table_rows.tex").write_text("\n".join(lines) + "\n")
    head = [
        "% written by sim/control/export.py",
        "\\makeatletter",
        "\\newcommand{\\ctrl}[2]{\\@ifundefined{ctrl@#1@#2}{\\textbf{??}}"
        "{\\csname ctrl@#1@#2\\endcsname}}",
    ]
    (OUT / "macros.tex").write_text("\n".join(head + macros + ["\\makeatother"]) + "\n")
    for name in ("tracking_quadrotor_figure_eight", "tracking_suspended_load_circle"):
        z = np.load(results / f"{name}.npz")
        k = slice(None, None, DECIMATE)
        np.savetxt(
            OUT / f"track_{name.removeprefix('tracking_')}.dat",
            np.column_stack(
                [
                    z["t"][k],
                    np.maximum(z["position_error"][k], 1e-12),
                    np.maximum(z["attitude_error"][k], 1e-12),
                ]
            ),
            header="t pos att",
            comments="",
            fmt="%.6g",
        )


GROUPS = (
    ("position", "pos", "position (m)"),
    ("attitude", "att", "attitude (rad)"),
    ("velocity", "vel", "velocity"),
)


def symbol(var):
    """TeX symbol of a reduced coordinate name (``delta_omega_q1`` -> ``\\delta\\omega_{q_1}``)."""

    def sub(name):
        name = re.sub(r"^([A-Za-z]+)_?(\d+)$", r"\1_\2", name)
        name = name.replace("R_L", "R_L")
        return name

    m = re.fullmatch(r"(gamma|xi|eta)_(.+)", var)
    if m:
        return f"\\{m.group(1)}_{{{sub(m.group(2))}}}"
    m = re.fullmatch(r"delta_(omega|Omega)_(.+)", var)
    if m:
        return f"\\delta\\{m.group(1)}_{{{sub(m.group(2))}}}"
    m = re.fullmatch(r"delta_(.+?)(_dot)?", var)
    base = m.group(1)
    head, _, idx = base.partition("_")
    head = {"theta": "\\theta"}.get(head, head)
    body = f"\\dot{{{head}}}" if m.group(2) else head
    return f"\\delta {body}" + (f"_{{{sub(idx)}}}" if idx else "")


def num(w):
    w = float(w)
    if w >= 1000:
        return f"10^{{{int(round(np.log10(w)))}}}"
    return f"{w:g}"


def weights_tex(labels, Q, R, is_velocity):
    by_var = {}
    for lab, w in zip(labels, Q):
        by_var.setdefault(re.sub(r"\[\d+\]$", "", lab), w)
    groups = {True: [], False: []}
    for var, w in by_var.items():
        groups[is_velocity(var + "[0]")].append((var, w))
    allw = [w for _, w in groups[False] + groups[True]]
    parts = []
    uniform = len(set(allw)) == 1
    for vel, name in ((False, "configuration"), (True, "velocity")):
        ws = [w for _, w in groups[vel]]
        if not ws:
            continue
        base = max(set(ws), key=ws.count)
        other = [f"${num(w)}$ on ${symbol(v)}$" for v, w in groups[vel] if w != base]
        text = f"${num(base)}$ on the {name} errors"
        if other:
            text += " except " + ", ".join(other)
        parts.append(text)
    if uniform:
        wx = "$W_x = I$" if allw[0] == 1 else f"$W_x = {num(allw[0])}\\, I$"
    else:
        wx = "$W_x$ diagonal, " + " and ".join(parts)
    Rs = np.atleast_1d(np.asarray(R, float))
    if np.all(Rs == Rs[0]):
        wu = "$W_u = I$" if Rs[0] == 1 else f"$W_u = {num(Rs[0])}\\, I$"
    else:
        wu = "$W_u = \\mathrm{diag}(" + ", ".join(num(r) for r in Rs) + ")$"
    return wx + ", " + wu


def weights_short(labels, Q, R, is_velocity):
    """Table cell of the weights: ``c, v; r`` when W_x is c on every configuration error and
    v on every velocity error and W_u = r I, with a dagger when some coordinate differs
    (the entry then lists the weights); the second value is whether it differs."""
    groups = {True: [], False: []}
    for lab, w in zip(labels, Q):
        groups[is_velocity(lab)].append(float(w))
    Rs = np.atleast_1d(np.asarray(R, float))
    base = [max(set(g), key=g.count) for g in (groups[False], groups[True]) if g]
    custom = any(len(set(g)) > 1 for g in groups.values()) or not np.all(Rs == Rs[0])
    cell = "$" + ", ".join(num(b) for b in base) + "; " + num(Rs[0])
    return cell + ("{}^\\dagger$" if custom else "$"), custom


def example_weights(results):
    """``{example name: (cell, custom, full sentence)}`` for every example of run.py; an
    example without stored histories takes the labels of another of the same system."""
    here = Path(results).resolve().parent
    sys.path.insert(0, str(here))
    import run
    from control_helper import is_velocity

    labels_of = {}
    for ex in run.EXAMPLES:
        npz = Path(results) / f"{ex.name}.npz"
        if npz.exists():
            labels_of.setdefault(ex.system, [str(x) for x in np.load(npz)["labels"]])
    out = {}
    for ex in run.EXAMPLES:
        npz = Path(results) / f"{ex.name}.npz"
        labels = [str(x) for x in np.load(npz)["labels"]] if npz.exists() else labels_of[ex.system]
        cell, custom = weights_short(labels, ex.Q(labels), ex.R, is_velocity)
        out[ex.name] = (
            cell,
            custom,
            weights_tex(labels, ex.Q(labels), ex.R, is_velocity),
        )
    return out


def histories(results):
    """Time histories and weights of each example for its catalog entry."""
    here = Path(results).resolve().parent
    sys.path.insert(0, str(here))
    import run
    from control_helper import label_group

    import geomech.systems as S

    W = example_weights(results)
    for ex in run.EXAMPLES:
        npz = Path(results) / f"{ex.name}.npz"
        if not npz.exists():
            continue
        d = json.loads((Path(results) / f"{ex.name}.json").read_text())
        z = np.load(npz)
        labels = [str(x) for x in z["labels"]]
        k = slice(None, None, DECIMATE)
        cols, names, plots = [z["t"][k]], ["t"], []
        for group, col, legend in GROUPS:
            idx = [i for i, lab in enumerate(labels) if label_group(lab) == group]
            if idx:
                cols.append(np.linalg.norm(z["z"][k][:, idx], axis=1))
                names.append(col)
                if group != "velocity":
                    plots.append((col, legend))
        du = np.abs(z["u"] - z["u_eq"])[k]
        units = {v["unit"] for v in d["peak_inputs"].values()}
        start = 0
        unit_of = []
        for u in S.get(d["system"], **d["kwargs"]).inputs:
            size = 3 if type(u).__name__ == "Vector" else 1
            unit_of += [d["peak_inputs"][str(u)]["unit"]] * size
            start += size
        assert start == du.shape[1], ex.name
        inputs = []
        for unit, col in (("N", "uf"), ("N·m", "um")):
            idx = [i for i, un in enumerate(unit_of) if un == unit]
            if idx:
                cols.append(du[:, idx].max(axis=1))
                names.append(col)
                inputs.append((col, "force (N)" if unit == "N" else "moment (N\\,m)"))
        assert {u for u in unit_of} == units
        np.savetxt(
            OUT / f"hist_{ex.name}.dat",
            np.column_stack(cols),
            header=" ".join(names),
            comments="",
            fmt="%.5g",
        )
        _, custom, w = W[ex.name]
        dat = f"sim/control/hist_{ex.name}.dat"
        style = ["blue, thick", "red, thick, dashed"]
        left = "\n".join(
            f"\\addplot[no marks, {s}] table[x=t, y={c}] {{{dat}}}; \\addlegendentry{{{lg}}}"
            for (c, lg), s in zip(plots, style)
        )
        right = "\n".join(
            f"\\addplot[no marks, {s}] table[x=t, y={c}] {{{dat}}}; \\addlegendentry{{{lg}}}"
            for (c, lg), s in zip(inputs, ["black, thick", "orange, thick, dashed"])
        )
        tex = [
            f"% written by sim/control/export.py from {ex.name}",
            (
                f"Weights: {w}; horizon {ex.horizon:g}\\,s.\\par"
                if custom
                else f"Horizon {ex.horizon:g}\\,s.\\par"
            ),
            "\\noindent\\begin{tikzpicture}",
            "\\begin{groupplot}[group style={group size=2 by 1, horizontal sep=1.1cm}, lqrplot]",
            "\\nextgroupplot[ylabel={error norm}]",
            left,
            "\\nextgroupplot[ylabel={$|u - u_\\mathrm{eq}|_\\infty$}]",
            right,
            "\\end{groupplot}",
            "\\end{tikzpicture}\\par",
        ]
        (OUT / f"lqr_{ex.name}.tex").write_text("\n".join(tex) + "\n")


if __name__ == "__main__":
    main(sys.argv[1])
    histories(sys.argv[1])
