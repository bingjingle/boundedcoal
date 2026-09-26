"""Shared plotting style for the historical bounded-coalescent figures.

BC is red, SC blue, truth black, and coalescent events use a black rug.
The default bands are pointwise equal-tailed 2.5/97.5% posterior quantiles,
as used in the submitted paper. Minimal-width 95% intervals are also retained
in the coordinate files for exploratory redraws.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BC = "#FF0000"       # bounded coalescent  -- base-R "red"
SC = "#0000FF"       # standard coalescent -- base-R "blue"
TRUTH = "#000000"
FILL_ALPHA = 0.25    # adjustcolor(..., alpha.f = 0.25)
BAND = "equal"       # "equal" (submitted paper) or "hpd" (exploratory)

plt.rcParams.update({
    "font.size": 10,
    "axes.edgecolor": "#333333",
    "axes.linewidth": 0.9,
    "axes.labelcolor": "#111111",
    "xtick.color": "#333333",
    "ytick.color": "#333333",
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "legend.frameon": False,
    "pdf.fonttype": 42,
})


def hpd(x, alpha=0.05):
    """Minimal-width 95% interval retained for exploratory comparison."""
    x = np.sort(np.asarray(x, float))
    n = len(x)
    k = int(np.floor((1 - alpha) * n))
    j = np.argmin(x[k:] - x[: n - k])
    return x[j], x[j + k]


def bands(d, which=None):
    """Pick the band convention out of a curve dict."""
    which = which or BAND
    if which == "hpd" and d.get("hpd_lo") is not None:
        return d["hpd_lo"], d["hpd_hi"]
    return d["low"], d["hi"]


def draw(ax, x, med, low, high, colour, label=None, step=False):
    """One method's posterior summary: heavy median, light band edges, faint fill."""
    if step:
        g = np.r_[x[0] - (x[1] - x[0]), x]
        ax.step(g, np.r_[med[0], med], color=colour, lw=2.0, label=label, where="post")
        ax.step(g, np.r_[low[0], low], color=colour, lw=0.9, ls="--", where="post")
        ax.step(g, np.r_[high[0], high], color=colour, lw=0.9, ls="--", where="post")
        ax.fill_between(g, np.r_[low[0], low], np.r_[high[0], high], color=colour,
                        alpha=FILL_ALPHA, lw=0, step="post")
    else:
        ax.plot(x, med, color=colour, lw=2.0, label=label)
        ax.plot(x, low, color=colour, lw=0.9, ls="--")
        ax.plot(x, high, color=colour, lw=0.9, ls="--")
        ax.fill_between(x, low, high, color=colour, alpha=FILL_ALPHA, lw=0)


def rug(ax, events, y=0.0):
    ax.plot(events, np.full(len(np.asarray(events)), y), marker="|", ls="none",
            color="black", ms=5, mew=0.6, clip_on=False, zorder=6)


def tidy(ax, xlabel="Time", ylabel=None, title=None):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, fontsize=10)
    ax.tick_params(labelsize=9)


def legend_handles():
    h = [plt.Line2D([], [], color=TRUTH, lw=2.0),
         plt.Line2D([], [], color=BC, lw=2.0),
         plt.Line2D([], [], color=SC, lw=2.0),
         plt.Line2D([], [], color="black", marker="|", ls="none", ms=8)]
    return h, ["Truth", "BC (bounded)", "SC (standard)", "Coalescent events"]
