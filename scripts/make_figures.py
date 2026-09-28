#!/usr/bin/env python3
"""
Make the README figures: neutron-star sequences for every EOS table.

Runs an enthalpy-formalism solver for each EOS, keeps the stable branch up to
M_max, and writes light and dark versions of each figure to figures/. Also
prints the summary table (values at M_max and 1.4 Msun) used in the README.
The solvers are bit-identical, so both give the same figures:
    --solver c       c/tov_h_c.x (default; build it with `make -C c`)
    --solver python  python/tovsolve.py

Usage (from the repository root):
    python scripts/make_figures.py [--solver c|python]

The functions are also used by python/tovsolve.ipynb.
"""

import argparse
import os
import subprocess
import sys

import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXE = os.path.join(ROOT, "c", "tov_h_c.x")
OUT = os.path.join(ROOT, "figures")
SOLVERS = {"c": "c/tov_h_c.x", "python": "python/tovsolve.py"}

# Central number densities (fm^-3): 281 stars per EOS
RHO_START, RHO_END, NSTEPS = 0.10, 1.50, 281
M_MIN = 0.2  # Msun; lighter configurations have R >> 20 km

# ---------------------------------------------------------------------------
# EOS tables. The MDI family is ordered by x, so it takes an ordinal blue
# ramp (lightness = x); the four unrelated EOSs take categorical hues. All
# colors are validated (dataviz validate_palette.js): MDI ramp --ordinal,
# the four categorical hues --pairs all, in light and dark mode.
# ---------------------------------------------------------------------------
EOS = [
    # key, file, label, group, light color, dark color
    ("MDI_x-2.0", "eos_MDI_x-2.0.in", "MDI  x = −2", "mdi", "#0d366b", "#9ec5f4"),
    ("MDI_x-1.0", "eos_MDI_x-1.0.in", "MDI  x = −1", "mdi", "#1c5cab", "#5598e7"),
    ("MDI_x0.0", "eos_MDI_x0.0.in", "MDI  x = 0", "mdi", "#3987e5", "#256abf"),
    ("MDI_x0.3", "eos_MDI_x0.3.in", "MDI  x = 0.3", "mdi", "#86b6ef", "#184f95"),
    ("APR", "eos_APR.in", "APR", "other", "#eda100", "#c98500"),
    ("DBHF", "eos_DBHF_BonnB.in", "DBHF (Bonn B)", "other", "#e87ba4", "#d55181"),
    ("FPS", "eos_FPS.in", "FPS", "other", "#008300", "#008300"),
    ("SLy4", "eos_SLY4.in", "SLy4", "other", "#4a3aa7", "#9085e9"),
]

THEMES = {
    "light": dict(surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e", muted="#898781",
                  grid="#e1e0d9", axis="#c3c2b7", context="#d3d2cb", band="#0b0b0b"),
    "dark": dict(surface="#1a1a19", ink="#ffffff", ink2="#c3c2b7", muted="#898781",
                 grid="#2c2c2a", axis="#383835", context="#454542", band="#ffffff"),
}


def rows_c(fname):
    """Output rows (as printed) of c/tov_h_c.x for one EOS."""
    if not os.path.exists(EXE):
        sys.exit(f"{EXE} not found: run `make -C c` first")
    res = subprocess.run(
        [EXE, fname, str(RHO_START), str(RHO_END), str(NSTEPS)],
        cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        universal_newlines=True, check=True)
    return res.stdout.splitlines()[1:]


def rows_python(fname):
    """Output rows (formatted as printed) of python/tovsolve.py for one EOS."""
    sys.path.insert(0, os.path.join(ROOT, "python"))
    import tovsolve
    seq = tovsolve.sequence(tovsolve.Eos(os.path.join(ROOT, fname)),
                            RHO_START, RHO_END, NSTEPS)
    stars = zip(*(seq[k] for k in tovsolve.Star._fields))
    return [tovsolve.format_row(s) for s in stars]


def stable_branch(rows):
    """Arrays from printed output rows, stable branch (M >= M_MIN, up to M_max).
    Both solvers go through the 6-decimal printed values, so they give
    identical figures and tables."""
    parsed = []
    for line in rows:
        if "*" in line:  # f11.6 overflow at the lowest densities
            continue
        parsed.append([float(v) for v in line.split()])
    a = np.array(parsed)
    d = dict(M=a[:, 0], R=a[:, 1], k2=a[:, 2], lam=a[:, 3], I=a[:, 4],
             beta=a[:, 5], rhoc=a[:, 6])
    d["Lam"] = (2.0 / 3.0) * d["k2"] / d["beta"] ** 5  # dimensionless Lambda
    imax = int(np.argmax(d["M"]))
    keep = np.arange(len(d["M"])) <= imax
    keep &= d["M"] >= M_MIN
    return {k: v[keep] for k, v in d.items()}


def load_all(solver="c"):
    """Stable-branch data for every EOS: {key: dict of arrays}."""
    rows = rows_c if solver == "c" else rows_python
    return {key: stable_branch(rows(fname)) for key, fname, *_ in EOS}


def at_mass(d, key, m=1.4):
    """Interpolate a quantity at mass m on the stable branch."""
    return float(np.interp(m, d["M"], d[key]))


def style_axes(ax, t):
    ax.set_facecolor(t["surface"])
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(t["axis"])
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(colors=t["muted"], labelcolor=t["ink2"], labelsize=10,
                   length=3, width=0.8)
    ax.grid(True, which="major", color=t["grid"], linewidth=0.8, linestyle="-")
    ax.set_axisbelow(True)


def place_labels(ax, pts, side, t, min_gap_px=15):
    """Direct labels in a gutter beside the data, joined to the curve ends
    by leader lines; labels repel vertically so they never collide."""
    fig = ax.figure
    fig.canvas.draw()
    trans = ax.transData
    inv = trans.inverted()
    disp = [trans.transform((x, y)) for (x, y, _, _) in pts]
    order = np.argsort([p[1] for p in disp])
    ys = [disp[i][1] for i in order]
    for k in range(1, len(ys)):  # push up
        ys[k] = max(ys[k], ys[k - 1] + min_gap_px)
    shift = (np.mean([disp[i][1] for i in order]) - np.mean(ys))
    ys = [y + shift for y in ys]
    for k in range(len(ys) - 2, -1, -1):  # keep order after the shift
        ys[k] = min(ys[k], ys[k + 1] - min_gap_px)
    bbox = ax.get_window_extent()
    xg = bbox.x1 - 4 if side == "right" else bbox.x0 + 4
    for k, i in enumerate(order):
        x0, y0 = disp[i]
        lx, ly = inv.transform((xg, ys[k]))
        cx, cy = pts[i][0], pts[i][1]
        ax.annotate(pts[i][2], xy=(cx, cy), xytext=(lx, ly),
                    textcoords="data", fontsize=9.5, color=t["ink"],
                    ha="right" if side == "right" else "left", va="center",
                    arrowprops=dict(arrowstyle="-", color=t["muted"], lw=0.7,
                                    shrinkA=2, shrinkB=5,
                                    connectionstyle="arc3,rad=0"),
                    zorder=6)


def draw_panel(ax, data, group, xkey, ykey, t, mode, logy=False):
    for key, _, label, grp, cl, cd in EOS:  # gray context: the other group
        if grp != group:
            d = data[key]
            ax.plot(d[xkey], d[ykey], color=t["context"], lw=1.1,
                    solid_capstyle="round", zorder=2)
    handles = []
    for key, _, label, grp, cl, cd in EOS:
        if grp != group:
            continue
        d = data[key]
        c = cl if mode == "light" else cd
        (h,) = ax.plot(d[xkey], d[ykey], color=c, lw=2.0, label=label,
                       solid_capstyle="round", solid_joinstyle="round", zorder=4)
        ax.plot(d[xkey][-1], d[ykey][-1], "o", ms=6.5, color=c,
                mec=t["surface"], mew=1.6, zorder=5)
        handles.append(h)
    if logy:
        ax.set_yscale("log")
    return handles


def make_figure(data, spec, mode, solver="c", save=True):
    """Draw one figure. save=True writes figures/<name>[_dark].png and returns
    the path; save=False returns the matplotlib Figure (e.g. for a notebook)."""
    t = THEMES[mode]
    plt.rcParams.update({"font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans"})
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.9), sharey=True,
                             gridspec_kw=dict(wspace=0.08))
    fig.patch.set_facecolor(t["surface"])
    fig.subplots_adjust(left=0.075, right=0.985, top=0.80, bottom=0.13)

    titles = ("MDI family, x = −2 … 0.3", "APR · DBHF (Bonn B) · FPS · SLy4")
    for ax, group, sub in zip(axes, ("mdi", "other"), titles):
        style_axes(ax, t)
        handles = draw_panel(ax, data, group, spec["x"], spec["y"], t, mode,
                             spec.get("logy", False))
        ax.set_xlim(*spec["xlim"])
        ax.set_ylim(*spec["ylim"])
        ax.set_title(sub, loc="left", fontsize=11, color=t["ink2"], pad=8)
        ax.set_xlabel(spec["xlabel"], fontsize=11, color=t["ink2"])
        if "extra" in spec:
            spec["extra"](ax, t)
        leg = ax.legend(handles=handles, loc=spec.get("legend", "best"),
                        frameon=False, fontsize=9, handlelength=1.6,
                        labelcolor=t["ink2"])
        leg.set_zorder(7)
        pts = []
        for key, _, label, grp, _, _ in EOS:
            if grp == group:
                d = data[key]
                pts.append((d[spec["x"]][-1], d[spec["y"]][-1], label, key))
        place_labels(ax, pts, spec["label_side"], t)
    axes[0].set_ylabel(spec["ylabel"], fontsize=11, color=t["ink2"])

    fig.text(0.075, 0.945, spec["title"], fontsize=15, color=t["ink"],
             fontweight="semibold", ha="left", va="center")
    fig.text(0.075, 0.885, spec["subtitle"].replace("SOLVER", SOLVERS[solver]),
             fontsize=10.5, color=t["ink2"], ha="left", va="center")
    if not save:
        return fig
    suffix = "" if mode == "light" else "_dark"
    path = os.path.join(OUT, f"{spec['name']}{suffix}.png")
    fig.savefig(path, dpi=160, facecolor=t["surface"])
    plt.close(fig)
    return path


def psr_band(ax, t, where="right"):
    """PSR J0740+6620: M = 2.08 +/- 0.07 Msun (Fonseca et al. 2021)."""
    ax.axhspan(2.01, 2.15, color=t["band"], alpha=0.06, lw=0, zorder=1)
    x0, x1 = ax.get_xlim()
    pad = 0.015 * (x1 - x0)
    ax.text(x1 - pad if where == "right" else x0 + pad, 2.08, "PSR J0740+6620",
            fontsize=8.5, color=t["muted"], ha=where, va="center", zorder=3)


def gw_bar(ax, t):
    """GW170817: Lambda(1.4) = 190 +390 -120 at 90% (Abbott et al. 2018)."""
    ax.errorbar([1.4], [190.0], yerr=[[120.0], [390.0]], fmt="s", ms=5,
                color=t["ink2"], mec=t["surface"], mew=1.2, elinewidth=1.2,
                capsize=3, zorder=6)
    ax.text(1.36, 115.0, "GW170817\n(90%)", fontsize=8.5, color=t["muted"],
            ha="right", va="center", linespacing=1.1, zorder=6)


SUB = "Stable branch up to M$_{max}$ (dot), 281 central densities per EOS · SOLVER"
FIGS = [
    dict(name="mass_radius", x="R", y="M", xlim=(7.4, 16.0), ylim=(0.2, 2.45),
         xlabel="Radius R (km)", ylabel="Mass M (M$_\\odot$)",
         title="Mass–radius relation", subtitle=SUB, label_side="left",
         legend="lower left", extra=psr_band),
    dict(name="love_number", x="M", y="k2", xlim=(0.2, 2.95), ylim=(0.0, 0.16),
         xlabel="Mass M (M$_\\odot$)", ylabel="Love number k$_2$",
         title="Tidal Love number k$_2$", subtitle=SUB, label_side="right",
         legend="upper right"),
    dict(name="tidal_deformability", x="M", y="Lam", xlim=(0.2, 2.95),
         ylim=(1.0, 3.0e5), logy=True, xlabel="Mass M (M$_\\odot$)",
         ylabel="Dimensionless Λ = (2/3) k$_2$ (c$^2$R/GM)$^5$",
         title="Tidal deformability Λ", subtitle=SUB, label_side="right",
         legend="upper right", extra=gw_bar),
    dict(name="moment_of_inertia", x="M", y="I", xlim=(0.2, 2.95), ylim=(0.0, 3.0),
         xlabel="Mass M (M$_\\odot$)", ylabel="Moment of inertia I (10$^{45}$ g cm$^2$)",
         title="Moment of inertia (slow rotation)", subtitle=SUB, label_side="right",
         legend="upper left"),
    dict(name="mass_density", x="rhoc", y="M", xlim=(0.1, 1.95), ylim=(0.2, 2.45),
         xlabel="Central baryon density n$_c$ (fm$^{-3}$)", ylabel="Mass M (M$_\\odot$)",
         title="Mass versus central density", subtitle=SUB, label_side="right",
         legend="lower right", extra=lambda ax, t: psr_band(ax, t, "left")),
]


def summary_table(data):
    """Markdown table: values at M_max and interpolated at 1.4 Msun."""
    lines = ["| EOS | M_max (M☉) | R at M_max (km) | n_c at M_max (fm⁻³) "
             "| R_1.4 (km) | k2_1.4 | Λ_1.4 | I_1.4 (10⁴⁵ g cm²) |",
             "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for key, _, label, *_ in EOS:
        d = data[key]
        lines.append(f"| {label.replace('  ', ' ')} | {d['M'][-1]:.3f} | {d['R'][-1]:.2f} "
                     f"| {d['rhoc'][-1]:.3f} | {at_mass(d, 'R'):.2f} | {at_mass(d, 'k2'):.4f} "
                     f"| {at_mass(d, 'Lam'):.0f} | {at_mass(d, 'I'):.3f} |")
    return "\n".join(lines)


def main(argv=None):
    ap = argparse.ArgumentParser(description="Make the README figures.")
    ap.add_argument("--solver", choices=sorted(SOLVERS), default="c",
                    help="enthalpy solver to run (default: %(default)s)")
    args = ap.parse_args(argv)

    import matplotlib
    matplotlib.use("Agg")
    os.makedirs(OUT, exist_ok=True)
    data = load_all(args.solver)

    for spec in FIGS:
        for mode in ("light", "dark"):
            path = make_figure(data, spec, mode, solver=args.solver)
            print("wrote", os.path.relpath(path, ROOT))

    print()
    print(summary_table(data))


if __name__ == "__main__":
    main()
