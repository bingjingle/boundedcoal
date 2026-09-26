# Maximum likelihood · Figure 4

Figure 4 compares standard-coalescent and bounded-coalescent skyline estimates
for three simulated genealogies with 100 tips. This analysis was conducted by
Julia Palacios.

| Scenario | Effective population size | Bound | Input genealogy |
|---|---|---|---|
| 1 | `1` | `1` | `../synthetic/data/syn1_bounded_ntip100.rda`, row 25 |
| 2 | `3 exp(-t)` | `0.7` | `../synthetic/data/syn2_bounded_ntip100.rda`, row 6 |
| 3 | `25 exp(-5t)` | `0.71` | `../synthetic/data/syn3_bounded_ntip100.rda`, row 13 |

Rows are one-based. The portable entrypoint rounds coalescent times to eight
decimal places, exactly reproducing the vectors embedded in the original
`legacy/fig4_syn*_ntip100.R` scripts. It retains their `ape::skyline` and
`phylodyn:::bounded_skyline_ascent` estimation routines.

## Run

Install `ape` and the bounded-coalescent version of `phylodyn`:

```r
install.packages(c("ape", "remotes"))
remotes::install_github("JuliaPalacios/phylodyn@4a3c160500ccf3470c31fadf5da9e4ef99cf4bb8")
```

From the repository root:

```bash
Rscript analyses/maximum_likelihood/figure4_mle.R
# Optional: supply an output directory as the first argument.
```

The command writes `generated/figure4_mle.pdf`, separate CSVs for the fitted
curves and truth, optimizer logs, and `sessionInfo.txt`. It does not run
posterior sampling.

The PDF redraws the legacy example estimates with explicit labels and axis
limits. Its plotting style is not a pixel-exact reproduction of the submitted
manuscript PDF.

## Provenance and limits

`legacy/` preserves the original combined MLE/posterior exploration scripts.
The scenario 1 and 2 scripts reference an undefined `bnpr` object in their
posterior sections. `plot_MLE_DIS_bounddata.R` has an incomplete `points()` call.
Use the portable entrypoint for the MLE plots.

**Dataset alignment:** the manuscript describes Figures 4 and 5 as using the
same genealogies. The retained Figure 4 sources use scenario 2 row 6, while the
current [posterior reproduction pipeline](../reproduction/) uses row 4
(`syn2_data3.csv`, with a zero-based filename). Scenarios 1 and 3 agree to the
original rounding. This reorganization preserves both input selections; it
does not reconcile that pre-existing discrepancy or change scientific results.
