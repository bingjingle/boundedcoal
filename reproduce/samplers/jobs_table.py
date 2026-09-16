#!/usr/bin/env python3
"""Spot-check of one block of Table 1: Ne_2(t) = 3exp(-t), tau = 0.7, 100 tips.

Purpose is narrow: decide whether OUR runs reproduce Bingjing's published numbers
closely enough that the whole table does not need re-running.  So every setting is
hers, not the manuscript's, wherever the two disagree:

  RI-BM     jitter 1e-8,  const 1,          1,000,000 + 1,000,000   (her syn2 script)
  RI-SE     step 0.000005, jitter 1e-11,      100,000 +   100,000   (her syn2 script)
  Discrete  prec_alpha = prec_beta = 0.01,    100,000 +   200,000
            ^ phylodyn's default, which is what her dis_syn2_ntip100.R actually used.
              The manuscript says 0.1.  Using 0.1 here would compare our run against
              a different prior and tell us nothing about reproducibility.

RI-SE runs on datasets 0-9 only; the other rows use all 30.
"""
import os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor

BASE = os.environ.get("BC_BASE") or os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(BASE, "data", "syn2_ntip100")
OUT  = os.path.join(BASE, "table", "out")
LOG  = os.path.join(BASE, "table", "logs")
os.makedirs(OUT, exist_ok=True); os.makedirs(LOG, exist_ok=True)

WORKERS = int(os.environ.get("BC_WORKERS") or os.cpu_count() or 4)

ENV = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
           MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1", NUMEXPR_NUM_THREADS="1")

NTIP, TAU = 100, 0.7
N_ALL, N_SE = 30, 10
# Subset knobs for a partial spot-check: BC_TABLE_N datasets overall, of which the
# first BC_TABLE_NSE also get the (very expensive) squared-exponential sampler.
# Defaults reproduce the full 140-run table exactly.  Iteration counts are NEVER
# reduced here -- a cheaper table means fewer datasets, not a different sampler,
# so the numbers stay comparable to the published ones.
N_ALL = int(os.environ.get("BC_TABLE_N", N_ALL))
N_SE  = int(os.environ.get("BC_TABLE_NSE", N_SE))
SMOKE = bool(os.environ.get("BC_SMOKE"))
BM_IT = ("200", "200") if SMOKE else ("1000000", "1000000")
SE_IT = ("100", "100") if SMOKE else ("100000", "100000")
if SMOKE:
    N_ALL, N_SE = 2, 1

def jobs():
    J = []
    for i in range(N_ALL):
        csv = os.path.join(DATA, f"syn2_data{i}.csv")
        for lik, flag in (("bc", []), ("sc", ["--standard"])):
            if i < N_SE:
                J.append((f"d{i:02d}_RI_SE_{lik}", 5, [
                    sys.executable, os.path.join(BASE, "code", "ri_se.py"),
                    "--proposal-sigma-l", "5e-06", "--nsim1", SE_IT[0],
                    "--nsim2", SE_IT[1], "--jitter", "1e-11", "--input-file", csv,
                    "--ntip", str(NTIP), "--bound", str(TAU), "--output-dir", OUT,
                    "--tag", f"d{i:02d}_RI_SE_{lik}"] + flag))
            J.append((f"d{i:02d}_RI_BM_{lik}", 3, [
                sys.executable, os.path.join(BASE, "code", "ri_bm.py"),
                "--const", "1", "--nsim1", BM_IT[0], "--nsim2", BM_IT[1],
                "--jitter", "1e-8", "--data", csv, "--ntip", str(NTIP),
                "--bound", str(TAU), "--output-dir", OUT,
                "--tag", f"d{i:02d}_RI_BM_{lik}"] + flag))
            J.append((f"d{i:02d}_DIS_{lik}", 1, [
                "Rscript", os.path.join(BASE, "code", "discrete.R"), lik, "0.01",
                f"d{i:02d}_DIS_{lik}", csv, str(NTIP), str(TAU), OUT]))
    # Method filter for a cheap partial spot-check, e.g. BC_TABLE_METHODS=RI_BM.
    # Default (unset) keeps every method, i.e. the published table.
    want = os.environ.get("BC_TABLE_METHODS", "").strip()
    if want:
        keys = [w.strip() for w in want.split(",") if w.strip()]
        J = [j for j in J if any(f"_{k}_" in j[0] for k in keys)]
    J.sort(key=lambda j: -j[1])
    return J

def run(job):
    tag, _, cmd = job
    done = os.path.join(OUT, tag + (".npz" if "RI" in tag else ".csv"))
    if os.path.exists(done):
        return tag, 0, 0.0
    t0 = time.time()
    with open(os.path.join(LOG, tag + ".log"), "w") as fh:
        rc = subprocess.call(cmd, stdout=fh, stderr=subprocess.STDOUT, env=ENV)
    el = time.time() - t0
    print(f"{'OK  ' if rc==0 else 'FAIL'} {tag:20s} {el/60:7.1f} min", flush=True)
    return tag, rc, el

if __name__ == "__main__":
    J = jobs()
    print(f"{len(J)} jobs", flush=True)
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        res = list(ex.map(run, J))
    bad = [t for t, rc, _ in res if rc != 0]
    print("FAILED:", bad if bad else "none", flush=True)
