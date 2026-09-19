#!/usr/bin/env bash
# Check (and where possible install) everything the reruns need.  Safe to re-run.
set -u
ok=1
say() { printf "%-42s %s\n" "$1" "$2"; }

echo "=== Python ==="
PY=$(command -v python3 || true)
[ -z "$PY" ] && { say "python3" "MISSING -- install Python 3.10+"; exit 1; }
say "python3" "$($PY -V 2>&1)"
$PY - <<'PYC'
import importlib, sys
need = {"numpy":"numpy>=1.26","pandas":"pandas>=2.0","matplotlib":"matplotlib>=3.7",
        "scipy":"scipy<1.16","statsmodels":"statsmodels>=0.14"}
miss=[]
for m,spec in need.items():
    try:
        mod=importlib.import_module(m); print(f"  {m:14s} {getattr(mod,'__version__','?')}")
    except ImportError:
        miss.append(spec); print(f"  {m:14s} MISSING")
open("/tmp/bc_missing_py.txt","w").write(" ".join(miss))
PYC
MISS=$(cat /tmp/bc_missing_py.txt 2>/dev/null)
if [ -n "$MISS" ]; then
  echo "installing: $MISS"
  $PY -m pip install --quiet $MISS || { say "pip install" "FAILED"; ok=0; }
fi
# scipy < 1.16 is a hard requirement: SciPy removed the Genz MVNDST routine that the
# squared-exponential lengthscale update reaches through statsmodels.
$PY - <<'PYC'
import sys, scipy
from packaging.version import Version
if Version(scipy.__version__) >= Version("1.16"):
    print(f"  !! scipy {scipy.__version__} is too new; ri_se.py needs scipy<1.16")
    print("     fix:  python3 -m pip install 'scipy<1.16'"); sys.exit(3)
try:
    from statsmodels.sandbox.distributions.extras import mvstdnormcdf
    print(f"  mvstdnormcdf OK (scipy {scipy.__version__})")
except Exception as e:
    print("  !! mvstdnormcdf unavailable:", e); sys.exit(3)
PYC
[ $? -ne 0 ] && ok=0

echo; echo "=== R ==="
if ! command -v Rscript >/dev/null; then
  say "Rscript" "MISSING -- install R 4.x from https://cran.r-project.org"; exit 1
fi
say "Rscript" "$(Rscript -e 'cat(R.version.string)')"
Rscript - <<'RC'
ip <- rownames(installed.packages())
need <- setdiff(c("ape","spam","expm","coda","Matrix","truncnorm","Rcpp","remotes"), ip)
if (length(need)) {
  cat("  installing:", paste(need, collapse=", "), "\n")
  install.packages(need, repos="https://cloud.r-project.org", quiet=TRUE)
}
for (p in c("ape","phylodyn")) cat(sprintf("  %-12s %s\n", p,
  if (p %in% rownames(installed.packages())) as.character(packageVersion(p)) else "MISSING"))
RC

# Pinned to ba7b607 (2026-09-19): switches the standard-coalescent path from ESS
# to ESS2 and doubles the random-walk step variance, so SC now agrees with BNPR.
# The bounded path is untouched -- BC results are unchanged (verified: medians
# agree to 3e-3).  Previously d3a6e5f (2026-07-24).  The previous pin 4a3c160 (2026-03-21) predated
# 84 commits including a rewrite of R/bounded_coal_functions.R and a likelihood
# fix; under it the bounded discrete runs gave 39-43% coverage where 95% is
# nominal, and crashed outright on syn1 at prec 0.1 and syn4 at 0.01.  If you
# bump this pin, re-run every discrete result -- the RI columns are unaffected.
# phylodyn: MUST be JuliaPalacios/phylodyn.  mdkarcher/phylodyn is the upstream
# package and contains no bounded-coalescent code at all -- no bound_ESS, no
# bounded_skyline_ascent -- so the discrete bounded runs cannot work against it.
Rscript - <<'RC'
has <- "phylodyn" %in% rownames(installed.packages())
bounded <- FALSE
if (has) {
  suppressMessages(library(phylodyn))
  bounded <- exists("bounded_skyline_ascent", where=asNamespace("phylodyn"))
  d <- packageDescription("phylodyn")
  cat("  phylodyn from:", if (!is.null(d$RemoteUsername)) d$RemoteUsername else "(unknown)", "\n")
}
if (!has || !bounded) {
  cat("  installing JuliaPalacios/phylodyn (has the bounded-coalescent code)\n")
  if (!requireNamespace("remotes", quietly=TRUE))
    install.packages("remotes", repos="https://cloud.r-project.org", quiet=TRUE)
  remotes::install_github("JuliaPalacios/phylodyn", ref="ba7b607",
                          upgrade="never", quiet=TRUE)
  suppressMessages(library(phylodyn))
  bounded <- exists("bounded_skyline_ascent", where=asNamespace("phylodyn"))
}
cat("  bounded_skyline_ascent present:", bounded, "\n")
if (!bounded) { cat("  !! phylodyn lacks the bounded code -- discrete runs will fail\n"); quit(status=3) }
RC
[ $? -ne 0 ] && ok=0

echo; echo "=== disk ==="
W="${BC_WORK:-$HOME/boundedcoal_work}"
while [ ! -d "$W" ] && [ "$W" != "/" ]; do W=$(dirname "$W"); done
df -h "$W" | tail -1
echo "  need ~10 GB free for transient chains (the reducer keeps it bounded)"
echo
[ "$ok" = 1 ] && echo "PREFLIGHT OK" || { echo "PREFLIGHT FAILED -- see above"; exit 1; }
