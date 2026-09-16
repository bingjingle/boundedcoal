#!/usr/bin/env bash
# One entry point for the whole rerun.  Resumable: every stage skips work that is
# already finished, so if the machine sleeps or you Ctrl-C, just run it again.
#
#   ./run_all.sh              everything: panel -> figures -> table -> summary
#   ./run_all.sh panel        just the simulation + COVID runs and their figures
#   ./run_all.sh table        just the Table 1 spot-check
#   ./run_all.sh figures      redraw from whatever results already exist
#   ./run_all.sh smoke        ~2 min end-to-end test at tiny iteration counts.
#                             RUN THIS FIRST -- it proves the environment works
#                             before you commit many hours.  Results are junk by
#                             design; delete $BC_WORK afterwards.
#
# Chains are written to a LOCAL work directory, never into iCloud: the raw
# posterior samples are tens of GB, syncing them would be slow and costly, and
# two machines writing the same folder would produce conflict files.  Only the
# small deliverables (figures, coordinates, summaries) are copied back next to
# this script.
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# FIX: this Mac's default python3 is 3.14 (Homebrew, PEP 668 externally-managed) with
# no scientific stack, and scipy<1.16 -- a hard requirement of ri_se.py, which reaches
# the Genz MVNDST routine through statsmodels -- publishes no 3.14 wheels.  The pinned
# stack therefore lives in a dedicated Python 3.11 venv.  Put it on PATH here so that
# `./run_all.sh` and `caffeinate -i ./run_all.sh` work as documented, and so a resume
# after sleep/Ctrl-C picks up the same interpreter without anything to remember.
BC_VENV="${BC_VENV:-$HOME/boundedcoal_venv}"
if [ -x "$BC_VENV/bin/python3" ]; then
  export VIRTUAL_ENV="$BC_VENV"
  export PATH="$BC_VENV/bin:$PATH"
else
  echo "!! expected Python venv not found at $BC_VENV" >&2
  echo "   recreate with: python3.11 -m venv $BC_VENV && $BC_VENV/bin/pip install 'numpy>=1.26' 'pandas>=2.0' 'matplotlib>=3.7' 'scipy<1.16' 'statsmodels>=0.14' packaging" >&2
  exit 1
fi
# FIX: the work dir must be chosen AFTER the stage is known.  This assignment
# used to sit above the smoke branch, which left BC_WORK already set, so the
# branch's "${BC_WORK:-$HOME/boundedcoal_smoke}" could never substitute and a
# smoke run wrote its 200-iteration output into the REAL work dir.  The resume
# logic in jobs_panel.py / jobs_table.py keys only on the tag, so the full run
# would then SKIP those tags and quietly publish smoke junk (syn1 row, all three
# COVID figures, table datasets d00/d01).  An explicitly exported BC_WORK still wins.
BC_WORK_PRESET="${BC_WORK:-}"
STAGE="${1:-all}"
if [ "$STAGE" = "smoke" ]; then
  export BC_SMOKE=1
  export BC_WORK="${BC_WORK_PRESET:-$HOME/boundedcoal_smoke}"
else
  export BC_WORK="${BC_WORK_PRESET:-$HOME/boundedcoal_work}"
fi
export BC_BASE="$BC_WORK"
export BC_WORKERS="${BC_WORKERS:-$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4)}"

mkdir -p "$BC_WORK/code" "$BC_WORK/data" "$BC_WORK/out" "$BC_WORK/logs" \
         "$BC_WORK/coords" "$BC_WORK/figures" "$BC_WORK/table/out" "$BC_WORK/table/logs"
# Some of Bingjing's original files are mode 444; cp onto an existing read-only
# file fails, which would break every re-run.  Make the working copies writable.
cp -Rf "$HERE/code/." "$BC_WORK/code/" 2>/dev/null || true
cp -Rf "$HERE/data/." "$BC_WORK/data/" 2>/dev/null || true
chmod -R u+w "$BC_WORK/code" "$BC_WORK/data" 2>/dev/null || true

echo "work dir : $BC_WORK   (local, not synced)"
echo "workers  : $BC_WORKERS"
echo "stage    : $STAGE"
echo

reducer_start() {
  python3 "$BC_WORK/code/reduce_npz.py" --loop > "$BC_WORK/logs/REDUCE.log" 2>&1 &
  echo $! > "$BC_WORK/.reducer.pid"
  echo "reducer running (pid $(cat "$BC_WORK/.reducer.pid")) -- slims finished chains every 5 min"
}
reducer_stop() {
  if [ -f "$BC_WORK/.reducer.pid" ]; then kill "$(cat "$BC_WORK/.reducer.pid")" 2>/dev/null; fi
  rm -f "$BC_WORK/.reducer.pid"
  python3 "$BC_WORK/code/reduce_npz.py" >/dev/null 2>&1   # final sweep
}
figures() {
  python3 "$BC_WORK/code/collect_coords.py"
  python3 "$BC_WORK/code/plot_panels3.py"
  python3 "$BC_WORK/code/plot_covid3.py"
  ( cd "$BC_WORK/code" && Rscript replot.R ) || echo "  (replot.R skipped)"
}
deliver() {
  D="$HERE/deliverables"; mkdir -p "$D/table"
  for d in figures figures_R coords; do
    if [ -d "$BC_WORK/$d" ]; then mkdir -p "$D/$d"; cp -f "$BC_WORK/$d/"* "$D/$d/" 2>/dev/null; fi
  done
  cp -f "$BC_WORK/table/"*.csv "$D/table/" 2>/dev/null
  cp -f "$BC_WORK/table/table1_check.txt" "$D/table/" 2>/dev/null
  cp -f "$BC_WORK/logs/QUEUE.log" "$D/" 2>/dev/null
  echo "deliverables copied to: $D"
  du -sh "$D" 2>/dev/null
}

if [ "$STAGE" = "panel" ] || [ "$STAGE" = "all" ] || [ "$STAGE" = "smoke" ]; then
  reducer_start
  echo "--- PANEL: 40 runs (4 simulations + COVID) x 3 methods x 2 likelihoods ---"
  python3 "$BC_WORK/code/jobs_panel.py" 2>&1 | tee "$BC_WORK/logs/QUEUE.log"
  reducer_stop
  figures
  deliver
fi

if [ "$STAGE" = "table" ] || [ "$STAGE" = "all" ] || [ "$STAGE" = "smoke" ]; then
  reducer_start
  echo "--- TABLE 1 spot-check: Ne_2(t), 100 tips, 140 runs ---"
  python3 "$BC_WORK/code/jobs_table.py" 2>&1 | tee "$BC_WORK/table/QUEUE.log"
  reducer_stop
  python3 "$BC_WORK/code/summarize_table.py" | tee "$BC_WORK/table/table1_check.txt"
  deliver
fi

if [ "$STAGE" = "figures" ]; then
  figures
  deliver
fi

echo
echo "DONE ($STAGE)"
