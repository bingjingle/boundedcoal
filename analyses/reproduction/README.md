# Posterior inference: coordinates and reproduction tools

Redraw the synthetic posterior panels and Washington State COVID-19 analysis
from the included coordinates, without running MCMC. These tools also preserve
historical exploratory workflows and the final synthetic table artifact.
Contribution information is recorded in the [analysis overview](../README.md).

## Quick start: redraw saved results

From the repository root, with base R installed:

```bash
bash analyses/reproduction/run_all.sh redraw
```

This writes the three-by-two synthetic panel (`panels.pdf`), six individual
synthetic panels, and the three-panel COVID figure (`covid_2_DIS.pdf`) to
`analyses/reproduction/deliverables/figures_R/`. The archived inputs and PDFs
are unchanged. `redraw` is also the default stage; `figures` is an alias.
No Python packages, phylodyn, network access, or posterior sampling are needed.

To choose an output folder:

```bash
BC_FIGURES_OUT=/tmp/boundedcoal-redraw bash analyses/reproduction/run_all.sh redraw
```

The saved coordinates reproduce the numerical curves. Small rendering differences
from the submitted PDFs can arise from the graphics device and fonts. To change
colours or limits, edit the plotting constants in [`figures/replot.R`](figures/replot.R).

## Find the paper results

| Paper item | Files | Scope |
|---|---|---|
| Figure 5: synthetic posterior inference | [`coords/panel_syn1_RI_BM.csv`](coords/panel_syn1_RI_BM.csv), corresponding `syn2`/`syn3` and `DIS` CSVs; [`figures_R/panels.pdf`](figures_R/panels.pdf) | Three trajectories, RI with the Brownian-motion kernel, and Discrete |
| Figure 6: Washington State COVID-19 | `coords/panel_covid_RI_BM.csv`, `coords/panel_covid_DIS.csv`; [`figures_R/covid_2_DIS.pdf`](figures_R/covid_2_DIS.pdf) | Posterior curves and the case-count panel |
| Table 2: posterior inference over 30 datasets per trajectory | [`results/table2_posterior_100tips.tex`](results/table2_posterior_100tips.tex), [`results/table2_run_manifest.txt`](results/table2_run_manifest.txt), [`results/runs.csv`](results/runs.csv) | Final numerical table and retained run records |
| Earlier posterior spot-check | `results/historical_ne2_check.*`, `results/historical_ne2_per_dataset.csv` | An incomplete Ne2-only check against an earlier draft; **not Table 1 or final Table 2** |
| Exploratory variants | `RI_SE` coordinates, `syn4` coordinates, [`output_figures/`](output_figures/) | Additional kernels/settings/layouts retained for provenance |

**Table 1 in the submitted paper compares simulation algorithms and wall-clock
costs. It is a different experiment.** The former `table1_*` filenames in this
folder referred to an earlier numbering of the posterior table. The final table
artifact is now named `table2_posterior_100tips.tex`.

Figure 4 (maximum likelihood estimation) and the original application/synthetic
source scripts are organized elsewhere under [`analyses/`](../README.md).

## Folder guide

| Folder | Contents |
|---|---|
| [`samplers/`](samplers/) | Parameterized RI and Discrete drivers, historical job lists, and a historical spot-check summarizer |
| [`figures/`](figures/) | Base-R redraw, Python exploratory layouts, coordinate export, and chain reduction |
| [`coords/`](coords/) | Saved curve coordinates, coalescent-event times, and run metadata |
| [`data/`](data/) | Coalescent-time CSVs, the CCD0 tree, and 30 Ne2 replicate CSVs |
| [`figures_R/`](figures_R/) | Archived two-column synthetic and combined COVID figures |
| [`output_figures/`](output_figures/) | Archived exploratory figures |
| [`results/`](results/) | Final Table 2 artifact, historical manifest, and earlier spot-check outputs |
| `deliverables/` | Newly generated output; ignored by Git |

Each `coords/panel_<tag>_<method>.csv` stores `x`, `truth` where applicable,
and posterior curves prefixed `bc_` or `sc_`. Suffixes `_med`, `_low_eq`, and
`_hi_eq` identify medians and equal-tailed 2.5/97.5% quantiles. The `_hpd`
columns store minimal-width intervals for exploratory use. The submitted paper
uses pointwise equal-tailed intervals. `events_<tag>.csv` supplies the event rug.

## Optional: run historical samplers

These are computationally expensive, resumable exploratory workflows. Their
settings are preserved; they do **not** regenerate the whole submitted Table 2.
In particular:

- `panel` includes RI-SE, a fourth synthetic scenario, and two Discrete precision
  priors. The submitted Figure 5 uses three scenarios and RI-BM/Discrete only.
- `table` checks only Ne2, with RI-SE on the first 10 datasets and the other methods
  on all 30. It uses shorter SC-Discrete chains than the final Table 2.
- Only the Ne2 replicate CSVs are bundled here. A complete final-table rerun would
  need the remaining replicate inputs and final run schedule.

Use Python 3.10–3.13 and a virtual environment. Python 3.11 was used for the
historical runs. The RI-SE path requires SciPy below 1.16. Direct `.rda` inputs
also need optional `pyreadr`; the bundled workflows use CSVs.

```bash
python3.11 -m venv /tmp/boundedcoal-venv
export BC_VENV=/tmp/boundedcoal-venv
bash analyses/reproduction/preflight.sh --install
bash analyses/reproduction/preflight.sh
```

`preflight.sh` is read-only unless `--install` is supplied. It validates Python
requirements, R dependencies, and the phylodyn repository/commit. Installation
uses the active R library and the selected Python environment.

The historical workflow's configured Discrete dependency is
[`JuliaPalacios/phylodyn@ba7b607`](https://github.com/JuliaPalacios/phylodyn/commit/ba7b607).
The preflight checks the recorded GitHub commit prefix, not merely the presence
of bounded-coalescent functions. The archived Table 2 run manifest records
`d3a6e5f39a6622eb252914636ef6e64d011d8568` for those earlier runs; this provenance
is retained. A newer dependency pin is not a claim that archived results were
produced with it. The earlier README's `4a3c160` pin was stale.

After dependency checks:

```bash
# Small, statistically meaningless environment check; isolated in a new temp folder.
bash analyses/reproduction/run_all.sh smoke

# Historical exploratory workflows; potentially many hours and large files.
BC_WORK=/tmp/boundedcoal-work BC_WORKERS=4 bash analyses/reproduction/run_all.sh panel
BC_WORK=/tmp/boundedcoal-work BC_WORKERS=4 bash analyses/reproduction/run_all.sh table
```

`BC_WORK` stores staged code, chains, logs, and new coordinates. It defaults to
`~/boundedcoal_work`. Normal runs copy small outputs into `deliverables/` (or
`BC_DELIVERABLES`). A smoke run always uses a fresh child of `BC_WORK`, if set,
or the system temporary directory; its outputs remain there. Smoke files cannot
cause a full run to skip sampling. For full runs, resume checks use output tags:
use a fresh `BC_WORK` after changing settings, data, or dependencies.

The background reducer replaces completed raw RI chains with summary arrays and
a thinned diagnostic chain. Preserve full chains separately if needed for further
analysis. Failed sampler jobs now produce a nonzero exit status.

## Recorded numerical considerations

Earlier phylodyn versions encountered bounded-likelihood underflow for some
settings; the historical dependency pins and comments should be read in that
version context. Seed changes do not establish convergence. The RI workflow
uses a tiny relative increase in the COVID bound to avoid a zero-width final
interval. These numerical conventions and all archived result values are
unchanged by the repository reorganization.
