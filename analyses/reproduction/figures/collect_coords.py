#!/usr/bin/env python3
"""Export sampler summaries as reusable CSV coordinates.

Each panel file contains x, truth where known, posterior medians, and both
credible-band conventions. Also exports coalescent-event times, a combined
long table, and a run index. Figures can then be redrawn without sampling.
"""
import glob, json, os
import os
import numpy as np
import pandas as pd

BASE = os.environ.get("BC_BASE") or os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))
OUT, CO = os.path.join(BASE, "out"), os.path.join(BASE, "coords")
os.makedirs(CO, exist_ok=True)

TRUTH = {"syn1": lambda t: np.ones_like(t),
         "syn2": lambda t: 3 * np.exp(-t),
         "syn3": lambda t: 25 * np.exp(-5 * t),
         "syn4": lambda t: 25 * np.exp(-5 * t),
         "covid": None}
TAU = {"syn1": 1.0, "syn2": 0.7, "syn3": 0.71, "syn4": 0.55,
       "covid": 0.5858767153195582}
METHODS = [("RI_BM", "Random integral — BM", False),
           ("RI_SE", "Random integral — SE", False),
           ("DIS",   "Discrete (phylodyn)",       True)]
ROWS = ["syn1", "syn2", "syn3", "syn4"]
DATA_CSV = {"syn1": "syn1_data24.csv", "syn2": "syn2_data3.csv",
            "syn3": "syn3_data12.csv", "syn4": "syn4_data20.csv",
            "covid": "covid_coal_times.csv"}


def load_npz(p):
    z = np.load(p, allow_pickle=True)
    k = set(z.files)
    g = np.asarray(z["grid"] if "grid" in k else z["d"], float)
    if "low_eq" in k:                                   # slimmed
        d = dict(x=g, med=z["med"], low=z["low_eq"], hi=z["hi_eq"],
                 hpd_lo=z["low_hpd"], hpd_hi=z["hi_hpd"])
    else:                                               # not yet slimmed
        ch = np.asarray(z["ne_grid"] if "ne_grid" in k else z["aa"], float)
        n1, n2 = int(z["nsim1"]), int(z["nsim2"])
        rem = int(z["rem"]) if "rem" in k else 10
        tot, n = n1 + n2, ch.shape[0]
        burn = n1 if n == tot else (int(n1 / rem) if abs(n - tot / rem) <= 2
                                    else int(round(n * n1 / tot)))
        S = ch[burn:]
        Ss = np.sort(S, axis=0); m = Ss.shape[0]
        kk = int(np.floor(0.95 * m)); w = Ss[kk:] - Ss[: m - kk]
        j = np.argmin(w, axis=0); c = np.arange(Ss.shape[1])
        d = dict(x=g, med=np.quantile(S, .5, axis=0),
                 low=np.quantile(S, .025, axis=0), hi=np.quantile(S, .975, axis=0),
                 hpd_lo=Ss[j, c], hpd_hi=Ss[j + kk, c])
    meta = {kk2: (z[kk2].item() if z[kk2].shape == () else z[kk2].tolist())
            for kk2 in ("nsim1", "nsim2", "rem", "secs", "per10k", "p", "jitter",
                        "const", "ntip", "T", "proposal_sigma_l", "bounded", "t")
            if kk2 in k}
    coal = np.asarray(z["coal"], float) if "coal" in k else (
        np.asarray(z["c"], float) if "c" in k else None)
    z.close()
    return d, meta, coal


def load_csv(p):
    t = pd.read_csv(p)
    return dict(x=t["x"].values, med=t["med"].values, low=t["low"].values,
                hi=t["hi"].values,
                hpd_lo=t["hpd_lo"].values if "hpd_lo" in t else None,
                hpd_hi=t["hpd_hi"].values if "hpd_hi" in t else None), {}, None


# The exploratory driver runs both priors. Export 0.01 by default, matching
# the submitted paper. Set BC_DIS_SUFFIX="" to inspect the historical 0.1 runs.
DIS_SUFFIX = os.environ.get("BC_DIS_SUFFIX", "_p001")


def find(tag, key, lik):
    sfx = DIS_SUFFIX if key == "DIS" else ""
    p = os.path.join(OUT, f"{tag}_{key}_{lik}{sfx}.npz")
    if os.path.exists(p):
        return load_npz(p)
    p = os.path.join(OUT, f"{tag}_{key}_{lik}{sfx}.csv")
    if os.path.exists(p):
        return load_csv(p)
    return None, None, None


def main():
    long_rows, index_rows, made = [], [], []
    for tag in ROWS + ["covid"]:
        ev = None
        for key, label, _ in METHODS:
            cur = {}
            for lik in ("bc", "sc"):
                d, meta, coal = find(tag, key, lik)
                if d is None:
                    continue
                cur[lik] = d
                if coal is not None and ev is None:
                    ev = coal
                index_rows.append(dict(tag=tag, method=key, likelihood=lik,
                                       **{k: v for k, v in (meta or {}).items()}))
                for i in range(len(d["x"])):
                    long_rows.append(dict(
                        dataset=tag, method=key, likelihood=lik, x=d["x"][i],
                        med=d["med"][i], low=d["low"][i], hi=d["hi"][i],
                        hpd_lo=(d["hpd_lo"][i] if d.get("hpd_lo") is not None else np.nan),
                        hpd_hi=(d["hpd_hi"][i] if d.get("hpd_hi") is not None else np.nan)))
            if not cur:
                continue
            x = cur[next(iter(cur))]["x"]
            tb = pd.DataFrame({"x": x})
            if TRUTH[tag] is not None:
                tb["truth"] = TRUTH[tag](x)
            for lik in ("bc", "sc"):
                if lik in cur:
                    d = cur[lik]
                    tb[f"{lik}_med"] = d["med"]
                    tb[f"{lik}_low_eq"] = d["low"]
                    tb[f"{lik}_hi_eq"] = d["hi"]
                    if d.get("hpd_lo") is not None:
                        tb[f"{lik}_low_hpd"] = d["hpd_lo"]
                        tb[f"{lik}_hi_hpd"] = d["hpd_hi"]
            f = os.path.join(CO, f"panel_{tag}_{key}.csv")
            tb.to_csv(f, index=False); made.append(os.path.basename(f))
        if ev is None:
            c = os.path.join(BASE, "data", DATA_CSV[tag])
            ev = pd.read_csv(c)["coal_times"].values
        pd.DataFrame({"coal_time": np.sort(ev)}).to_csv(
            os.path.join(CO, f"events_{tag}.csv"), index=False)

    if long_rows:
        pd.DataFrame(long_rows).to_csv(os.path.join(CO, "all_curves_long.csv"), index=False)
    if index_rows:
        pd.DataFrame(index_rows).to_csv(os.path.join(CO, "run_index.csv"), index=False)
    json.dump({"tau": TAU, "methods": {k: l for k, l, _ in METHODS}},
              open(os.path.join(CO, "meta.json"), "w"), indent=1)
    print(f"wrote {len(made)} panel files, {len(long_rows)} long rows")
    for m in made:
        print("  ", m)


if __name__ == "__main__":
    main()
