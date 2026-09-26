#!/usr/bin/env python3
"""Summarize the historical Ne2-only posterior spot-check.

Reference numbers below come from an earlier draft, not submitted Table 2.
Metrics use 100 evaluation points and pointwise equal-tailed 95% intervals:
SSE = sum((median - truth)**2), coverage = mean(low <= truth <= high), and
width = mean(high - low). Across datasets report median (q25, q75); runtime
is reported as mean +/- standard deviation per 10,000 iterations.
"""
import glob, os, re, sys
import numpy as np
import pandas as pd

BASE = os.environ.get("BC_BASE") or os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(BASE, "table", "out")
TAU, NG = 0.7, 100
X = np.linspace(0.0, TAU, NG + 1)[1:]
TRUTH = 3 * np.exp(-X)

# Earlier-draft posterior table, Ne_2(t), 100 tips.  (sse_med, q25, q75), (cov...), (width...), (time, sd)
PUBLISHED = {
    ("RI_BM", "bc"): ((21.93, 10.93, 41.43), (100, 100, 100), (5.01, 4.71, 5.75), (2.47, 0.08)),
    ("RI_BM", "sc"): ((46.76, 34.54, 64.46), (85, 60, 100),   (1.82, 1.62, 2.21), (3.93, 0.73)),
    ("RI_SE", "bc"): ((25.29, 17.98, 40.82), (100, 100, 100), (6.97, 6.38, 7.99), (1421.92, 249.37)),
    ("RI_SE", "sc"): ((82.71, 65.31, 93.72), (44, 39, 53),    (1.35, 1.26, 1.51), (1561.98, 415.60)),
    ("DIS",   "bc"): ((87.59, 68.26, 114.94), (85, 47, 100),  (2.71, 2.13, 3.41), (95.38, 18.03)),
    ("DIS",   "sc"): ((50.39, 35.80, 101.51), (82, 51, 100),  (1.72, 1.54, 1.97), (11.01, 3.36)),
}
LABEL = {("RI_BM", "bc"): "BC-RI-BM", ("RI_BM", "sc"): "SC-RI-BM",
         ("RI_SE", "bc"): "BC-RI-SE", ("RI_SE", "sc"): "SC-RI-SE",
         ("DIS", "bc"): "BC-Discrete-BM", ("DIS", "sc"): "SC-Discrete-BM"}


def metrics(med, low, hi, per10k):
    return dict(sse=float(np.sum((med - TRUTH) ** 2)),
                cov=float(np.mean((TRUTH >= low) & (TRUTH <= hi)) * 100),
                width=float(np.sum(hi - low) / NG), per10k=per10k)


def read_runs():
    rows = []
    for f in sorted(glob.glob(os.path.join(OUT, "*.npz"))):
        m = re.match(r"d(\d+)_(RI_BM|RI_SE)_(bc|sc)\.npz$", os.path.basename(f))
        if not m:
            continue
        z = np.load(f, allow_pickle=True); k = set(z.files)
        if "low_eq" in k:
            med, lo, hi = z["med"], z["low_eq"], z["hi_eq"]
        else:                                   # not yet slimmed
            ch = np.asarray(z["ne_grid"] if "ne_grid" in k else z["aa"], float)
            n1, n2 = int(z["nsim1"]), int(z["nsim2"])
            rem = int(z["rem"]) if "rem" in k else 10
            tot, n = n1 + n2, ch.shape[0]
            burn = n1 if n == tot else (int(n1 / rem) if abs(n - tot / rem) <= 2
                                        else int(round(n * n1 / tot)))
            S = ch[burn:]
            med = np.quantile(S, .5, axis=0); lo = np.quantile(S, .025, axis=0)
            hi = np.quantile(S, .975, axis=0)
        per = float(z["per10k"]) if "per10k" in k else (
            float(z["p"]) if "p" in k else np.nan)
        z.close()
        rows.append(dict(data=int(m.group(1)), method=m.group(2), lik=m.group(3),
                         **metrics(np.asarray(med), np.asarray(lo), np.asarray(hi), per)))
    for f in sorted(glob.glob(os.path.join(OUT, "*_DIS_*.csv"))):
        m = re.match(r"d(\d+)_DIS_(bc|sc)\.csv$", os.path.basename(f))
        if not m:
            continue
        d = pd.read_csv(f)
        log = os.path.join(BASE, "table", "logs", os.path.basename(f)[:-4] + ".log")
        per = np.nan
        if os.path.exists(log):
            mm = re.search(r"per10k\s+([0-9.]+)", open(log).read())
            if mm:
                per = float(mm.group(1))
        rows.append(dict(data=int(m.group(1)), method="DIS", lik=m.group(2),
                         **metrics(d["med"].values, d["low"].values, d["hi"].values, per)))
    return pd.DataFrame(rows)


def q(v):
    v = np.asarray(v, float); v = v[~np.isnan(v)]
    return (np.median(v), np.quantile(v, .25), np.quantile(v, .75)) if len(v) else (np.nan,)*3


def main():
    df = read_runs()
    if df.empty:
        print("no table runs finished yet"); return
    df.to_csv(os.path.join(BASE, "table", "historical_ne2_per_dataset.csv"), index=False)
    out = []
    print(f"{'Method':16s} {'n':>3s}  {'SSE  rerun (earlier draft)':>44s}"
          f"  {'Coverage % rerun (earlier draft)':>34s}  {'Width rerun (earlier draft)':>30s}"
          f"  {'s/10k rerun (earlier draft)':>26s}")
    for (meth, lik), pub in PUBLISHED.items():
        g = df[(df["method"] == meth) & (df["lik"] == lik)]
        if g.empty:
            print(f"{LABEL[(meth,lik)]:16s}   -   (not run yet)"); continue
        s, c, w = q(g["sse"]), q(g["cov"]), q(g["width"])
        t = np.asarray(g["per10k"], float); t = t[~np.isnan(t)]
        tm, ts = (t.mean(), t.std(ddof=1) if len(t) > 1 else 0.0) if len(t) else (np.nan, np.nan)
        ps, pc, pw, pt = pub
        print(f"{LABEL[(meth,lik)]:16s} {len(g):3d}  "
              f"{s[0]:8.2f} ({s[1]:7.2f},{s[2]:8.2f}) vs {ps[0]:7.2f} ({ps[1]:6.2f},{ps[2]:7.2f})"
              f"  {c[0]:5.0f} ({c[1]:3.0f},{c[2]:4.0f}) vs {pc[0]:4.0f} ({pc[1]:3.0f},{pc[2]:4.0f})"
              f"  {w[0]:6.2f} ({w[1]:5.2f},{w[2]:6.2f}) vs {pw[0]:5.2f}"
              f"  {tm:8.1f}+-{ts:7.1f} vs {pt[0]:7.1f}")
        out.append(dict(method=LABEL[(meth, lik)], n=len(g),
                        sse=s[0], sse_q25=s[1], sse_q75=s[2], sse_pub=ps[0],
                        sse_pub_q25=ps[1], sse_pub_q75=ps[2],
                        cov=c[0], cov_q25=c[1], cov_q75=c[2], cov_pub=pc[0],
                        width=w[0], width_q25=w[1], width_q75=w[2], width_pub=pw[0],
                        per10k=tm, per10k_sd=ts, per10k_pub=pt[0], per10k_pub_sd=pt[1],
                        sse_in_pub_iqr=bool(ps[1] <= s[0] <= ps[2])))
    if out:
        pd.DataFrame(out).to_csv(os.path.join(BASE, "table", "historical_ne2_check.csv"), index=False)
        print("\nwrote table/historical_ne2_check.csv and table/historical_ne2_per_dataset.csv")


if __name__ == "__main__":
    main()
