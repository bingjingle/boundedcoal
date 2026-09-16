# Reproducing the figures

Everything needed to regenerate the simulation panel and the COVID-19 figures,
and to redraw either of them without re-running any sampler.

Rerun date: 2026-09-16. Nothing here modifies the original repository contents;
this is a new top-level `reproduce/` directory.

---

## Quick start

```bash
./preflight.sh        # checks/installs Python + R dependencies, ~5 min
./run_all.sh smoke    # whole pipeline at tiny iteration counts, ~2 min
./run_all.sh panel    # the real run, writes figures and coordinates
```

`run_all.sh` is resumable: every finished run is skipped, so if it is
interrupted, run it again.

**Chains are written to `$BC_WORK` (default `~/boundedcoal_work`), not into this
repository.** The raw posterior samples are tens of gigabytes. A background
reducer slims each finished chain from ~320 MB to ~5 MB, keeping the median,
both credible-band conventions, a 5,000-draw thinned chain and every setting.
Only the small outputs are copied back.

---

## What is here

| directory | contents |
|---|---|
| `samplers/` | `ri_bm.py`, `ri_se.py` (Bingjing's random-integral samplers, parameterised), `discrete.R` (phylodyn), `bc_lik.py`, and the job drivers `jobs_panel.py` / `jobs_table.py` / `summarize_table.py` |
| `figures/` | `plot_panel3x3.py` (the 3x3 simulation figure), `plot_covid3.py`, `collect_coords.py`, `replot.R` (base-R redraw), `bcstyle.py`, `reduce_npz.py` |
| `coords/` | **one CSV per panel** — the coordinates of every curve, so no figure ever has to be re-run to be redrawn |
| `data/` | the six input CSVs, `median_ccd0.tree`, and `syn2_ntip100/` (30 replicate datasets) |
| `output_figures/` | the PDFs produced by this rerun |
| `results/` | `table1_check.txt` and the per-dataset CSVs from the Table 1 spot-check |

---

## The coordinates

`coords/panel_<tag>_<method>.csv` has one row per grid point:

```
x, truth, bc_med, bc_low_eq, bc_hi_eq, bc_low_hpd, bc_hi_hpd,
          sc_med, sc_low_eq, sc_hi_eq, sc_low_hpd, sc_hi_hpd
```

`bc_` is the bounded coalescent, `sc_` the standard coalescent. Both credible
band conventions are stored: `_eq` is equal-tailed (2.5% / 97.5% posterior
quantiles) and `_hpd` is the minimal-width interval.

Also present: `events_<tag>.csv` (coalescent times, the tick marks),
`all_curves_long.csv` (everything in one long table) and `run_index.csv`
(settings and timing for every run).

To restyle any figure — colours, band convention, axis limits — edit the
constants at the top of `figures/replot.R` and run `Rscript replot.R`. No
sampler runs.

---

## Settings

| method | prior | jitter | MH step | iterations |
|---|---|---|---|---|
| RI-BM | Γ(0.1, 0.1) | 1e-8 (syn1/2), 1e-7 (syn3) | — | 1,000,000 + 1,000,000 |
| RI-SE | Γ(0.1, 0.1) | 1e-11 (syn1/2), 1e-14 (syn3) | 0.05 / 0.000005 / 0.003 | 100,000 + 100,000 |
| Discrete | **Γ(0.01, 0.01)** | n/a | — | 100,000 burn-in + 200,000 |

The discrete runs use `prec_alpha = prec_beta = 0.01`, phylodyn's default. This
is what the published discrete results used; `discrete.R` passes it explicitly
rather than relying on the default.

COVID: 103 tips, bound τ = the TMRCA of the CCD0 genealogy, following
`covid_plot.R`.

---

## phylodyn

The discretized method requires **`JuliaPalacios/phylodyn`**, pinned here to
commit `4a3c160500ccf3470c31fadf5da9e4ef99cf4bb8`. `preflight.sh` installs it.

The widely-known `mdkarcher/phylodyn` contains **no bounded-coalescent code** —
no `bound_ESS`, no `bounded_skyline_ascent` — so the discrete column cannot be
reproduced against it.

## Python

`ri_se.py` reaches the Genz multivariate-normal routine through
`statsmodels.sandbox.distributions.extras.mvstdnormcdf`, which requires
**`scipy < 1.16`**: SciPy removed the underlying MVNDST routine in 1.16.
`preflight.sh` enforces this. This rerun used Python 3.11 with scipy 1.15.3.

---

## Known numerical issues

Two are worth recording, since both are properties of the methods rather than of
this setup.

**1. phylodyn's bounded likelihood can break down.** `mcmc_sampling(alg =
"bound_ESS")` can abort with

```
Error in if (bound_prob < 0) { : missing value where TRUE/FALSE needed
Calls: ... -> loglik -> coal_loglik_bounded
```

`bound_prob` underflows to `NA`. Which datasets this hits depends on the
precision prior: at seed 123, syn1 fails at 0.1 and syn4 fails at 0.01. Neither
precision completes all four simulation scenarios. `discrete.R` accepts
`BC_SEED` so the chain can be re-drawn, but that sidesteps the problem rather
than fixing it — the underlying computation should be done in log space.

**2. The random-integral samplers cannot take τ = max(coalescent time).** The
bounded branch builds intervals `[0, t_1, ..., t_n, τ]`, so a bound sitting
exactly on the last coalescent event gives a final interval of width zero.
`ri_bm.py` raises on this; `ri_se.py` has no guard and will sample it silently.
Grid-based methods are unaffected, since they build `seq(0, tau, length.out =
ngrid + 1)` and never form per-interval widths. The COVID bound is exactly this
case, so the RI runs there use τ raised by a relative 1e-9 — about 15,000 times
smaller than the smallest gap between consecutive coalescent times in that
dataset.
