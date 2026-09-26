#!/usr/bin/env bash
# Read-only by default. Pass --install to install missing/mismatched dependencies.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
case "${1:-}" in
  "") INSTALL=0 ;;
  --install) INSTALL=1 ;;
  *) echo "Usage: $0 [--install]" >&2; exit 2 ;;
esac
if [ -n "${BC_VENV:-}" ]; then
  [ -x "$BC_VENV/bin/python3" ] || { echo "No Python in BC_VENV=$BC_VENV" >&2; exit 1; }
  export PATH="$BC_VENV/bin:$PATH"
elif [ -x "$HOME/boundedcoal_venv/bin/python3" ]; then
  export PATH="$HOME/boundedcoal_venv/bin:$PATH"
fi
command -v python3 >/dev/null || { echo "Python 3 is missing." >&2; exit 1; }
command -v Rscript >/dev/null || { echo "Rscript is missing." >&2; exit 1; }
if [ "$INSTALL" = 1 ]; then
  python3 -m pip install -r "$HERE/requirements.txt"
fi
export BC_DEPENDENCIES_INSTALL="$INSTALL"
export BC_REQUIREMENTS="$HERE/requirements.txt"
status=0
python3 - <<'PYTHON' || status=1
import importlib, os, sys
from pathlib import Path
try:
    from packaging.requirements import Requirement
except ImportError:
    sys.exit("Missing packaging; install requirements.txt in a virtual environment.")
print("Python:", sys.version.split()[0])
failed = []
for line in Path(os.environ["BC_REQUIREMENTS"]).read_text().splitlines():
    if not line or line.startswith("#"):
        continue
    req = Requirement(line)
    try:
        module = importlib.import_module(req.name)
        version = module.__version__
        valid = version in req.specifier
        print(f"  {req.name}: {version} ({'OK' if valid else 'MISMATCH: ' + str(req.specifier)})")
        if not valid: failed.append(req.name)
    except ImportError:
        print(f"  {req.name}: MISSING")
        failed.append(req.name)
if failed:
    sys.exit("Python dependencies require attention: " + ", ".join(failed))
from statsmodels.sandbox.distributions.extras import mvstdnormcdf
print("  RI-SE multivariate-normal interface: available")
PYTHON
Rscript - <<'RSCRIPT' || status=1
wanted <- "ba7b607"
install <- identical(Sys.getenv("BC_DEPENDENCIES_INSTALL"), "1")
required <- c("ape", "spam", "expm", "coda", "Matrix", "truncnorm", "Rcpp", "remotes")
missing <- setdiff(required, rownames(installed.packages()))
if (length(missing) && install) {
  install.packages(missing, repos = "https://cloud.r-project.org")
  missing <- setdiff(required, rownames(installed.packages()))
}
check_pin <- function() {
  if (!"phylodyn" %in% rownames(installed.packages())) return(FALSE)
  d <- packageDescription("phylodyn")
  isTRUE(identical(d$RemoteUsername, "JuliaPalacios")) &&
    isTRUE(identical(d$RemoteRepo, "phylodyn")) &&
    isTRUE(startsWith(d$RemoteSha, wanted))
}
if (!check_pin() && install) {
  if (!requireNamespace("remotes", quietly = TRUE)) stop("remotes is required to install phylodyn")
  remotes::install_github("JuliaPalacios/phylodyn", ref = wanted, upgrade = "never")
}
pinned <- check_pin()
has <- requireNamespace("phylodyn", quietly = TRUE)
# Read namespace exports after checking package metadata; installation is opt-in.
bounded <- has && exists("bounded_skyline_ascent", where = asNamespace("phylodyn"), inherits = FALSE)
cat(R.version.string, "\n")
if (has) cat("phylodyn installed commit:", packageDescription("phylodyn")$RemoteSha, "\n")
cat("Required phylodyn: JuliaPalacios/phylodyn@", wanted, "\n", sep = "")
cat("Pin matches:", pinned, "; bounded functions present:", bounded, "\n")
if (length(missing)) cat("Missing R packages:", paste(missing, collapse = ", "), "\n")
if (length(missing) || !pinned || !bounded) quit(status = 1)
RSCRIPT
if [ "$status" = 0 ]; then
  echo "PREFLIGHT OK"
else
  echo "PREFLIGHT FAILED. Use a virtual environment and see README.md; --install opts into installation." >&2
fi
exit "$status"
