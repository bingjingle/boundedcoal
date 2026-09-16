#!/usr/bin/env python3
"""The publication 3x3 simulation panel, as a single PDF.

One figure, three rows (Ne_1..Ne_3) x three columns (RI-BM, RI-SE, Discrete).
Reads ONLY ../coords/panel_<tag>_<method>.csv, so it never touches a sampler
output and can be re-run in a second to restyle.

    python3 plot_panel3x3.py            -> figures/panels_3x3.pdf and .png

Layout notes, since alignment was the point of this file:
  * Rows share both axes (sharey/sharex per row).  Each row has its own Ne scale
    and its own bound tau, but the three panels within a row are identical axes,
    so their plot areas line up exactly.
  * Tick labels are drawn once per row, on the left column only.  subplots_adjust
    (not tight_layout) fixes every axes to the same width, so the "120" labels of
    row 3 cannot push that row's plot area rightwards relative to row 1's "5".
    This is what was misaligned before.
  * Axis labels appear once each: the trajectory (with its tau) as the y label of
    column 0, "Time" under the bottom row only.
  * Column titles on the top row, one shared legend under the figure.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import bcstyle as S

BASE = os.environ.get("BC_BASE") or os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))
CO = os.path.join(BASE, "coords")
FIG = os.path.join(BASE, "figures")
os.makedirs(FIG, exist_ok=True)

# tag, row label (trajectory and bound), tau, y limit
ROWS = [
    ("syn1", r"$N_{e,1}(t)=1,\ \tau=1$",            1.00,   5.0),
    ("syn2", r"$N_{e,2}(t)=3e^{-t},\ \tau=0.7$",    0.70,  14.0),
    ("syn3", r"$N_{e,3}(t)=25e^{-5t},\ \tau=0.71$", 0.71, 120.0),
]
COLS = [("RI_BM", "RI-BM",              False),
        ("RI_SE", "RI-SE",              False),
        ("DIS",   "Discrete (100 grids)", True)]


def draw_panel(ax, tag, key, step, tau, ylim, band=None):
    f = os.path.join(CO, f"panel_{tag}_{key}.csv")
    if not os.path.exists(f):
        ax.text(.5, .5, "not run", ha="center", va="center",
                transform=ax.transAxes, color="#999999", fontsize=9)
        return False
    d = pd.read_csv(f)
    x = d["x"].values
    for lik, col in (("sc", S.SC), ("bc", S.BC)):      # BC drawn last, on top
        if f"{lik}_med" not in d:
            continue
        which = band or S.BAND
        lo = f"{lik}_low_hpd" if (which == "hpd" and f"{lik}_low_hpd" in d) else f"{lik}_low_eq"
        hi = f"{lik}_hi_hpd" if (which == "hpd" and f"{lik}_hi_hpd" in d) else f"{lik}_hi_eq"
        S.draw(ax, x, d[f"{lik}_med"].values, d[lo].values, d[hi].values,
               col, step=step)
    if "truth" in d:
        ax.plot(x, d["truth"].values, color=S.TRUTH, lw=2.0, zorder=5)
    ev = pd.read_csv(os.path.join(CO, f"events_{tag}.csv"))["coal_time"].values
    S.rug(ax, ev, y=0.0)
    ax.set_xlim(0, tau)
    ax.set_ylim(0, ylim)
    return True


def main(band=None):
    suf = "" if (band or S.BAND) == "equal" else "_hpd"
    fig, axes = plt.subplots(3, 3, figsize=(9.6, 8.4),
                             sharex="row", sharey="row")
    # Fixed geometry: every axes gets an identical box, decorations live in the
    # margins.  tight_layout would size each axes around its own tick labels and
    # reintroduce the row-to-row misalignment.
    fig.subplots_adjust(left=0.105, right=0.985, top=0.935, bottom=0.115,
                        wspace=0.10, hspace=0.30)

    n = 0
    for r, (tag, rlab, tau, ylim) in enumerate(ROWS):
        for c, (key, clab, step) in enumerate(COLS):
            ax = axes[r, c]
            n += draw_panel(ax, tag, key, step, tau, ylim, band)
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
            ax.tick_params(labelsize=9)
            if r == 0:
                ax.set_title(clab, fontsize=11, pad=8)
            if c == 0:
                ax.set_ylabel(rlab, fontsize=10.5, labelpad=6)
            if r == len(ROWS) - 1:
                ax.set_xlabel("Time", fontsize=10.5)

    h, l = S.legend_handles()
    fig.legend(h, l, loc="lower center", ncol=4, frameon=False,
               fontsize=10, bbox_to_anchor=(0.5, 0.012))

    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"panels_3x3{suf}.{ext}"), dpi=200)
    plt.close(fig)
    print(f"panels_3x3{suf}.pdf  --  {n}/9 panels have data")


if __name__ == "__main__":
    import sys
    main("hpd" if "--hpd" in sys.argv else None)
