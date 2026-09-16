#!/usr/bin/env python3
"""The simulation figure: 4 rows (Ne_1..Ne_4) x 3 columns.

Columns, left to right, as asked:
    RI-BM   random integral, Brownian-motion kernel      (Bingjing's sampler)
    RI-SE   random integral, squared-exponential kernel  (Bingjing's sampler)
    DIS     discretised, 100 grid points                 (Julia's phylodyn)

This replaces the published 4-column layout, which interleaved the two kernels
(RI-BM, Discrete-BM, RI-SE, Discrete-SE).

Every number drawn here also exists as CSV in ../coords/, so the figure can be
rebuilt in R without touching Python.  Writes the combined grid and, because
main2.tex includes panels individually at 0.3\\textwidth, each panel on its own.
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bcstyle as S

BASE = os.environ.get("BC_BASE") or os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))
CO, FIG = os.path.join(BASE, "coords"), os.path.join(BASE, "figures")
os.makedirs(FIG, exist_ok=True)

ROWS = [("syn1", r"$N_{e,1}(t)=1$",          1.00,   5.0),
        ("syn2", r"$N_{e,2}(t)=3e^{-t}$",    0.70,  14.0),
        ("syn3", r"$N_{e,3}(t)=25e^{-5t}$",  0.71, 120.0),
        ("syn4", r"$N_{e,4}(t)=25e^{-5t}$",  0.55, 120.0)]
COLS = [("RI_BM", "RI-BM",              False),
        ("RI_SE", "RI-SE",              False),
        ("DIS",   "Discrete (phylodyn)", True)]


def panel(ax, tag, key, step, tau, ylim, band=None, show_y=True, title=None):
    f = os.path.join(CO, f"panel_{tag}_{key}.csv")
    if not os.path.exists(f):
        ax.text(.5, .5, "not yet run", ha="center", va="center",
                transform=ax.transAxes, color="#999999", fontsize=9)
        ax.set_xlim(0, tau); ax.set_ylim(0, ylim); S.tidy(ax, title=title)
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
    ax.set_xlim(0, tau); ax.set_ylim(0, ylim)
    S.tidy(ax, xlabel="Time", ylabel=(r"$N_e(t)$" if show_y else None), title=title)
    return True


def main(band=None):
    suf = "" if (band or S.BAND) == "equal" else "_hpd"
    fig, axes = plt.subplots(4, 3, figsize=(10.5, 12.0))
    n = 0
    for r, (tag, lab, tau, ylim) in enumerate(ROWS):
        for c, (key, clab, step) in enumerate(COLS):
            ok = panel(axes[r, c], tag, key, step, tau, ylim, band,
                       show_y=(c == 0), title=(clab if r == 0 else None))
            n += ok
        axes[r, 0].text(-0.30, 0.5, lab, transform=axes[r, 0].transAxes,
                        rotation=90, va="center", ha="center", fontsize=11)
    h, l = S.legend_handles()
    fig.legend(h, l, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.004))
    fig.tight_layout(rect=(0.02, 0.035, 1, 1))
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"panels_3col{suf}.{ext}"), dpi=200,
                    bbox_inches="tight")
    plt.close(fig)

    for tag, lab, tau, ylim in ROWS:            # individual panels for main2.tex
        for key, clab, step in COLS:
            f, ax = plt.subplots(figsize=(3.6, 3.0))
            if panel(ax, tag, key, step, tau, ylim, band):
                # Each panel carries its own legend: main2.tex includes these one
                # at a time at 0.3\textwidth, so the shared legend on the combined
                # panels_3col figure never reaches the reader.  Matches the
                # published panels, which are each self-contained.
                lh, ll = S.legend_handles()
                ax.legend(lh, ll, fontsize=6.0, loc="best", framealpha=0.85)
                f.tight_layout()
                f.savefig(os.path.join(FIG, f"{tag}_{key}{suf}.pdf"),
                          bbox_inches="tight")
                f.savefig(os.path.join(FIG, f"{tag}_{key}{suf}.png"), dpi=200,
                          bbox_inches="tight")
            plt.close(f)
    print(f"panels_3col{suf}.pdf  --  {n}/12 panels have data")


if __name__ == "__main__":
    main("hpd" if "--hpd" in sys.argv else None)
