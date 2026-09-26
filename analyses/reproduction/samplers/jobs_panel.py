#!/usr/bin/env python3
"""Build the historical exploratory simulation/COVID job list.

Includes RI-BM, RI-SE, and Discrete at two precision priors, together with a
fourth synthetic scenario. The submitted paper uses three synthetic scenarios,
RI-BM, and Discrete with precision prior 0.01. This workflow retains earlier
settings for provenance; it is not an exact regeneration of every paper result.
"""
import json, os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor

BASE = os.environ.get("BC_BASE") or os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))
OUT  = os.path.join(BASE, "out")
LOG  = os.path.join(BASE, "logs")
os.makedirs(OUT, exist_ok=True); os.makedirs(LOG, exist_ok=True)

WORKERS = int(os.environ.get("BC_WORKERS") or os.cpu_count() or 4)

ENV = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
           MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1", NUMEXPR_NUM_THREADS="1")

# tag, csv, ntip, tau, BM jitter, SE step, SE jitter
ROWS = [
    ("syn1", "syn1_data24.csv", 100, 1.00, 1e-8,  0.05,     1e-11),
    ("syn2", "syn2_data3.csv",  100, 0.70, 1e-8,  0.000005, 1e-11),
    ("syn3", "syn3_data12.csv", 100, 0.71, 1e-7,  0.003,    1e-14),
    ("syn4", "syn4_data20.csv", 100, 0.55, 1e-7,  0.003,    1e-14),   # syn3 settings
]
# COVID: 103 tips, tau = TMRCA.  Settings retained from the real-data workflow.
COVID = ("covid", "covid_coal_times.csv", 103, 0.5858767153195582, 1e-9, 0.001, 1e-13)

RI_BM_IT = (1_000_000, 1_000_000)   # burn-in, sample
RI_SE_IT = (100_000, 100_000)

# BC_SMOKE=1 runs the identical pipeline at tiny iteration counts (~2 min) so the
# environment can be proved end-to-end before committing many hours.  Results are
# statistically meaningless -- it only checks that everything runs and plots.
SMOKE = bool(os.environ.get("BC_SMOKE"))
if SMOKE:
    RI_BM_IT, RI_SE_IT = (200, 200), (100, 100)

def jobs():
    J = []
    todo = [ROWS[0], COVID] if SMOKE else ROWS + [COVID]
    for tag, csv, ntip, tau, jbm, step, jse in todo:
        data = os.path.join(BASE, "data", csv)
        for lik, flag in (("bc", []), ("sc", ["--standard"])):
            J.append((f"{tag}_RI_SE_{lik}", 5, [
                sys.executable, os.path.join(BASE, "code", "ri_se.py"),
                "--proposal-sigma-l", repr(step), "--nsim1", str(RI_SE_IT[0]),
                "--nsim2", str(RI_SE_IT[1]), "--jitter", repr(jse),
                "--input-file", data, "--ntip", str(ntip), "--bound", repr(tau),
                "--output-dir", OUT, "--tag", f"{tag}_RI_SE_{lik}"] + flag))
            J.append((f"{tag}_RI_BM_{lik}", 3, [
                sys.executable, os.path.join(BASE, "code", "ri_bm.py"),
                "--const", "1", "--nsim1", str(RI_BM_IT[0]), "--nsim2", str(RI_BM_IT[1]),
                "--jitter", repr(jbm), "--data", data, "--ntip", str(ntip),
                "--bound", repr(tau), "--output-dir", OUT,
                "--tag", f"{tag}_RI_BM_{lik}"] + flag))
            # Retain both historical prior variants. Coordinate export uses
            # the 0.01 variant, matching the submitted paper's specification.
            for prec, suff in ((0.1, ""), (0.01, "_p001")):
                J.append((f"{tag}_DIS_{lik}{suff}", 1, [
                    "Rscript", os.path.join(BASE, "code", "discrete.R"), lik, str(prec),
                    f"{tag}_DIS_{lik}{suff}", data, str(ntip), repr(tau), OUT]))
    J.sort(key=lambda j: -j[1])          # slowest first
    return J

def run(job):
    tag, _, cmd = job
    done = os.path.join(OUT, tag + (".npz" if "RI" in tag else ".csv"))
    if os.path.exists(done):
        print(f"SKIP {tag}", flush=True); return tag, 0, 0.0
    t0 = time.time()
    with open(os.path.join(LOG, tag + ".log"), "w") as fh:
        rc = subprocess.call(cmd, stdout=fh, stderr=subprocess.STDOUT, env=ENV)
    el = time.time() - t0
    print(f"{'OK  ' if rc==0 else 'FAIL'} {tag:24s} {el/60:7.1f} min", flush=True)
    return tag, rc, el

if __name__ == "__main__":
    J = jobs()
    print(f"{len(J)} jobs, {os.cpu_count()} cpus", flush=True)
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        res = list(ex.map(run, J))
    bad = [t for t, rc, _ in res if rc != 0]
    print("FAILED:", bad if bad else "none", flush=True)
    sys.exit(1 if bad else 0)
