#!/usr/bin/env python3
"""Shrink finished RI sampler outputs in place, keeping everything a figure or a
table entry could need.

The samplers save the entire posterior chain: for RI-BM that is 200,000 thinned
draws x 100 grid points, for RI-SE 200,000 unthinned draws, plus a second array of
equal size holding Ne at the coalescent event times.  Across the panel and the
table run that is ~45 GB, and this machine has 50 GB free.

What is kept:
  * both credible-band conventions, computed from the post-burn-in draws --
    equal-tailed 2.5/97.5% (what Bingjing's scripts and phylodyn compute, hence
    what the published figures/table show) and minimal-width HPD (what the
    manuscript captions describe);
  * the posterior median;
  * a 5,000-draw thinned copy of the grid chain, so the bands can be recomputed
    or a trace inspected later;
  * every scalar setting and the timing.
What is dropped: the full chain and the event-time array (`aaa`), neither of which
enters any figure or any table column.

Safe to run repeatedly and while jobs are in flight: a file is touched only if it
is older than 90 s and loads cleanly, and the rewrite is atomic.
"""
import os, sys, time
import numpy as np

KEEP_DRAWS = 5000


def hpd_cols(S, alpha=0.05):
    S = np.sort(S, axis=0)
    n = S.shape[0]
    k = int(np.floor((1 - alpha) * n))
    w = S[k:] - S[: n - k]
    j = np.argmin(w, axis=0)
    c = np.arange(S.shape[1])
    return S[j, c], S[j + k, c]


def reduce_one(path):
    if time.time() - os.path.getmtime(path) < 90:
        return "young"
    try:
        z = np.load(path, allow_pickle=True)
        keys = set(z.files)
    except Exception as e:
        return f"unreadable ({e})"
    if "slim" in keys:
        return "already"
    chain = z["ne_grid"] if "ne_grid" in keys else z["aa"]
    chain = np.asarray(chain, float)
    nsim1, nsim2 = int(z["nsim1"]), int(z["nsim2"])
    rem = int(z["rem"]) if "rem" in keys else 10
    total, n = nsim1 + nsim2, chain.shape[0]
    if n == total:                      # RI-SE: appended every iteration
        burn = nsim1
    elif abs(n - total / rem) <= 2:     # RI-BM: appended every rem-th iteration
        burn = int(nsim1 / rem)
    else:
        burn = int(round(n * nsim1 / total))
    S = chain[burn:]
    if S.shape[0] < 50:
        return f"too few draws ({S.shape[0]})"
    med = np.quantile(S, 0.5, axis=0)
    lo_eq = np.quantile(S, 0.025, axis=0)
    hi_eq = np.quantile(S, 0.975, axis=0)
    lo_hp, hi_hp = hpd_cols(S)
    step = max(1, S.shape[0] // KEEP_DRAWS)

    out = dict(slim=True, grid=np.asarray(z["grid"] if "grid" in keys else z["d"]),
               med=med, low_eq=lo_eq, hi_eq=hi_eq, low_hpd=lo_hp, hi_hpd=hi_hp,
               chain_thin=S[::step].astype(np.float32), chain_thin_step=step,
               n_draws=S.shape[0], burn_rows=burn)
    for k in ("coal", "c", "nsim1", "nsim2", "rem", "secs", "per10k", "p",
              "bounded", "jitter", "const", "ntip", "T", "proposal_sigma_l",
              "t", "theta1_init", "floor"):
        if k in keys:
            out[k] = z[k]
    if "s" in keys:                      # theta chain
        s = np.asarray(z["s"], float); out["theta"] = s[:: max(1, len(s) // 20000)]
    if "theta" in keys:
        s = np.asarray(z["theta"], float); out["theta"] = s[:: max(1, len(s) // 20000)]
    z.close()

    before = os.path.getsize(path)
    tmp = path + ".tmp.npz"
    np.savez_compressed(tmp, **out)
    os.replace(tmp, path)
    return f"{before/1e6:6.0f} MB -> {os.path.getsize(path)/1e6:5.1f} MB"


def sweep(dirs):
    for d in dirs:
        if not os.path.isdir(d):
            continue
        for f in sorted(os.listdir(d)):
            if f.endswith(".npz") and not f.endswith(".tmp.npz"):
                p = os.path.join(d, f)
                r = reduce_one(p)
                if r not in ("already", "young"):
                    print(f"  {f:28s} {r}", flush=True)


if __name__ == "__main__":
    BASE = os.environ.get("BC_BASE") or os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))
    dirs = [os.path.join(BASE, "out"), os.path.join(BASE, "table", "out")]
    if "--loop" in sys.argv:
        while True:
            sweep(dirs)
            # Sweep interval is env-tunable: a tighter interval lowers PEAK disk,
            # which matters when the machine is shared with another run.
            time.sleep(int(os.environ.get("BC_REDUCE_EVERY", "300")))
    else:
        sweep(dirs)
