#!/usr/bin/env python3
"""The three COVID real-data figures, one per method.

Axis convention is taken from Julia's covid_plot.R:
    t0 = 2020-06-08, the sampling date
    time runs backwards from the present, so the x axis is reversed (tau .. 0)
    ticks are labelled as calendar months, t0 - 365*x
There is no truth curve -- this is real data -- so each panel shows only the two
posterior summaries, BC in red and SC in blue, with the coalescent events as a rug.

Coordinates for all three panels are in ../coords/panel_covid_*.csv, including a
`date` column, so these can be redrawn in R directly.
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

T0 = np.datetime64("2020-06-08")            # covid_plot.R: date of sampling
TAU = 0.5858767153195582                    # TMRCA of the UPGMA tree
METHODS = [("RI_BM", "Random integral — BM",  False),
           ("RI_SE", "Random integral — SE",  False),
           ("DIS",   "Discrete (phylodyn)",        True)]


def date_axis(ax, tau):
    """Reversed time axis labelled with calendar months, as in covid_plot.R."""
    ticks = np.arange(0.0, tau + 1e-9, 0.1)
    labs = []
    for t in ticks:
        d = (T0 - np.timedelta64(int(round(t * 365)), "D")).astype("datetime64[D]")
        labs.append(pd.Timestamp(d).strftime("%b'%y"))
    ax.set_xticks(ticks); ax.set_xticklabels(labs, fontsize=8)
    ax.set_xlim(tau, 0.0)                    # present on the right


def one(key, label, step, ylim=None, band=None):
    f = os.path.join(CO, f"panel_covid_{key}.csv")
    if not os.path.exists(f):
        print(f"  skip {key}: no results yet"); return False
    d = pd.read_csv(f); x = d["x"].values
    fig, ax = plt.subplots(figsize=(4.4, 3.4))
    for lik, col in (("sc", S.SC), ("bc", S.BC)):
        if f"{lik}_med" not in d:
            continue
        which = band or S.BAND
        lo = f"{lik}_low_hpd" if (which == "hpd" and f"{lik}_low_hpd" in d) else f"{lik}_low_eq"
        hi = f"{lik}_hi_hpd" if (which == "hpd" and f"{lik}_hi_hpd" in d) else f"{lik}_hi_eq"
        S.draw(ax, x, d[f"{lik}_med"].values, d[lo].values, d[hi].values, col, step=step)
    ev = pd.read_csv(os.path.join(CO, "events_covid.csv"))["coal_time"].values
    S.rug(ax, ev, y=0.0)
    if ylim is None:
        # Fixed y range so the three COVID panels are directly comparable, matching
        # ylim=c(0,20) in Julia's covid_plot.R.  Set BC_COVID_YLIM=auto to go back to
        # the data-driven 97th-percentile scaling.
        env = os.environ.get("BC_COVID_YLIM", "20")
        if env.lower() == "auto":
            hi_all = [d[c].values for c in d.columns if c.endswith(("_hi_eq", "_hi_hpd"))]
            ylim = float(np.nanpercentile(np.concatenate(hi_all), 97)) if hi_all else 10.0
        else:
            ylim = float(env)
    ax.set_ylim(0, ylim)
    date_axis(ax, TAU)
    S.tidy(ax, xlabel="Time", ylabel=r"$N_e(t)$", title=label)
    h, l = S.legend_handles()
    ax.legend(h[1:], l[1:], fontsize=7.5, loc="upper left")
    fig.tight_layout()
    suf = "" if (band or S.BAND) == "equal" else "_hpd"
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"covid_{key}{suf}.{ext}"), dpi=200,
                    bbox_inches="tight")
    plt.close(fig)
    print(f"  covid_{key}{suf}.pdf")
    return True


def main(band=None):
    n = sum(one(k, l, st, band=band) for k, l, st in METHODS)
    print(f"{n}/3 COVID figures written")


if __name__ == "__main__":
    main("hpd" if "--hpd" in sys.argv else None)
