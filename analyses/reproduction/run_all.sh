#!/usr/bin/env bash
# Redraw archived coordinates, or run the historical exploratory workflows.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAGE="${1:-redraw}"

case "$STAGE" in
  redraw|figures)
    # Base R is sufficient; no sampler, phylodyn, or Python environment is needed.
    export BC_REPRODUCE="$HERE"
    export BC_FIGURES_OUT="${BC_FIGURES_OUT:-$HERE/deliverables/figures_R}"
    exec Rscript "$HERE/figures/replot.R"
    ;;
  panel|table|all|smoke) ;;
  *) echo "Usage: $0 {redraw|panel|table|all|smoke}" >&2; exit 2 ;;
esac

# These stages preserve the earlier exploratory settings, including RI-SE and
# the Ne2-only table check. They do not regenerate the full submitted Table 2.
if [ -n "${BC_VENV:-}" ]; then
  [ -x "$BC_VENV/bin/python3" ] || { echo "No Python in BC_VENV=$BC_VENV" >&2; exit 1; }
  export PATH="$BC_VENV/bin:$PATH"
elif [ -x "$HOME/boundedcoal_venv/bin/python3" ]; then
  export PATH="$HOME/boundedcoal_venv/bin:$PATH"
fi
command -v python3 >/dev/null || { echo "Python 3 is required." >&2; exit 1; }
command -v Rscript >/dev/null || { echo "Rscript is required." >&2; exit 1; }

if [ "$STAGE" = smoke ]; then
  # Always use a new child directory, even when BC_WORK is explicitly set.
  # Smoke runs can never satisfy the resume checks for a full run.
  SMOKE_PARENT="${BC_WORK:-${TMPDIR:-/tmp}}"
  mkdir -p "$SMOKE_PARENT"
  export BC_WORK="$(mktemp -d "$SMOKE_PARENT/boundedcoal-smoke.XXXXXX")"
  export BC_SMOKE=1
else
  export BC_WORK="${BC_WORK:-$HOME/boundedcoal_work}"
  if [ -n "${BC_SMOKE:-}" ]; then
    echo "Unset BC_SMOKE for a full run, or choose the smoke stage." >&2; exit 2
  fi
fi
export BC_BASE="$BC_WORK"
export BC_WORKERS="${BC_WORKERS:-$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4)}"
mkdir -p "$BC_WORK"/{code,data,out,logs,coords,figures,table/out,table/logs}

# Stage the two source folders in the flat layout used by the job drivers.
cp -f "$HERE/samplers/"*.py "$HERE/samplers/"*.R "$BC_WORK/code/"
cp -f "$HERE/figures/"*.py "$HERE/figures/"*.R "$BC_WORK/code/"
cp -R "$HERE/data/." "$BC_WORK/data/"
chmod -R u+w "$BC_WORK/code" "$BC_WORK/data"
if [ "$STAGE" = smoke ]; then
  DELIVERY="$BC_WORK/deliverables"
else
  DELIVERY="${BC_DELIVERABLES:-$HERE/deliverables}"
fi
printf 'Work directory: %s\nWorkers: %s\nStage: %s\n' "$BC_WORK" "$BC_WORKERS" "$STAGE"

REDUCER_PID=""
reducer_stop() {
  if [ -n "$REDUCER_PID" ]; then
    kill "$REDUCER_PID" 2>/dev/null || true
    wait "$REDUCER_PID" 2>/dev/null || true
    REDUCER_PID=""
  fi
}
trap reducer_stop EXIT
reducer_start() {
  python3 "$BC_WORK/code/reduce_npz.py" --loop > "$BC_WORK/logs/REDUCE.log" 2>&1 &
  REDUCER_PID=$!
}
figures() {
  python3 "$BC_WORK/code/collect_coords.py"
  python3 "$BC_WORK/code/plot_panels3.py"
  python3 "$BC_WORK/code/plot_covid3.py"
  # A smoke panel contains only syn1 and COVID, so a full 3-row redraw is skipped.
  if [ "$STAGE" != smoke ]; then
    BC_REPRODUCE="$BC_WORK" BC_FIGURES_OUT="$BC_WORK/figures_R" \
      Rscript "$BC_WORK/code/replot.R"
  fi
}
deliver() {
  mkdir -p "$DELIVERY"
  for d in figures figures_R coords table; do
    if [ -d "$BC_WORK/$d" ]; then
      mkdir -p "$DELIVERY/$d"
      for f in "$BC_WORK/$d/"*.pdf "$BC_WORK/$d/"*.png "$BC_WORK/$d/"*.csv "$BC_WORK/$d/"*.txt; do
        [ ! -f "$f" ] || cp -f "$f" "$DELIVERY/$d/"
      done
    fi
  done
  printf 'Outputs: %s\n' "$DELIVERY"
}
if [ "$STAGE" = panel ] || [ "$STAGE" = all ] || [ "$STAGE" = smoke ]; then
  reducer_start
  python3 "$BC_WORK/code/jobs_panel.py" 2>&1 | tee "$BC_WORK/logs/QUEUE.log"
  reducer_stop
  python3 "$BC_WORK/code/reduce_npz.py"
  figures
  deliver
fi
if [ "$STAGE" = table ] || [ "$STAGE" = all ] || [ "$STAGE" = smoke ]; then
  reducer_start
  python3 "$BC_WORK/code/jobs_table.py" 2>&1 | tee "$BC_WORK/table/QUEUE.log"
  reducer_stop
  python3 "$BC_WORK/code/reduce_npz.py"
  python3 "$BC_WORK/code/summarize_table.py" | tee "$BC_WORK/table/historical_ne2_check.txt"
  deliver
fi
printf 'DONE (%s)\n' "$STAGE"
