# Reproducibility and provenance

This index was checked against the Overleaf manuscript on 26 September 2026. Table and figure numbers refer to that manuscript, rather than the names used in earlier code drafts.

## Simulation: Figures 2–3 and Table 1

[`simulation/`](../simulation/) contains Algorithm 1, validation and comparison plotting scripts, and partial code for the timing benchmark. Table 1 compares 3,000 simulated genealogies per setting, with 50 or 100 tips. It is distinct from the posterior-inference comparison in Table 2. The original benchmark code does not provide a complete automated reconstruction of every submitted timing entry; no missing benchmark results have been invented.

## Maximum likelihood: Figure 4

[`analyses/maximum_likelihood/`](../analyses/maximum_likelihood/) contains the original scripts. The three scenarios have population sizes `1`, `3 exp(-t)`, and `25 exp(-5t)`, with bounds `1`, `0.7`, and `0.71`. The original scripts embed their input genealogies. See the [MLE run notes](maximum-likelihood.md) for the standalone entry point and the input difference between the second MLE example and the later posterior redraw bundle.

## Synthetic posterior inference: Figure 5 and Table 2

[`analyses/synthetic/data/`](../analyses/synthetic/data/) retains the original R datasets. The submitted Table 2 uses 30 bounded-coalescent datasets with 100 tips for each of the three scenarios above, and four methods: BC-RI, SC-RI, BC-Discrete, and SC-Discrete.

| Method | Precision prior, shape/rate | Burn-in | Subsequent iterations | Thinning |
| :--- | :--- | ---: | ---: | ---: |
| RI, boundary-corrected Brownian motion | Gamma(0.1, 0.1) | 1,000,000 | 1,000,000 | 10 |
| BC-Discrete | Gamma(0.01, 0.01) | 100,000 | 200,000 | 1 |
| SC-Discrete | Gamma(0.01, 0.01) | 100,000 | 2,000,000 | 20 |

SSE, pointwise interval coverage, and mean interval width are evaluated on 100 grid points and summarized across datasets. Coverage is the fraction of evaluation points covered by 95% equal-tailed intervals, not simultaneous coverage of an entire trajectory. Because the precision priors differ, the comparison does not isolate discretization alone.

[`analyses/reproduction/`](../analyses/reproduction/) includes the recorded Table 2 LaTeX summary, saved plotting coordinates, and historical rerun tools. Earlier files named `table1_check` refer to an older posterior-inference spot check, **not** the current Table 1 simulation benchmark. The historical job drivers include additional squared-exponential and fourth-scenario experiments and do not implement the complete submitted Table 2 schedule. Consult the [reproduction run notes](reproduction.md) before launching them.

## COVID-19: Figure 6

[`analyses/covid/`](../analyses/covid/) contains the original analysis source and CCD0 genealogy for the 103 Washington State sequences collected on 8 June 2020. Saved coordinates and redraw tools are in [`analyses/reproduction/`](../analyses/reproduction/). The manuscript identifies the sequence dataset as [GISAID EPI_SET_260825mx](https://doi.org/10.55876/gis8.260825mx).

## Preservation and scope of this cleanup

The reorganization starts from commit `b89e8177953177179fe3a69b164df97dc5769aae`. The [file-move index](file-moves.json) records the original paths and current locations of retained files. Datasets and numerical results for the paper are retained. The changes organize files, clarify attribution, repair execution paths, and document existing limitations; they do not revise the manuscript's scientific conclusions or rerun its full posterior analysis.

Original exploratory scripts may still depend on intermediate objects or historical environments. The run notes distinguish them from supported entry points. A successful redraw verifies that stored coordinates can be plotted; it does not independently validate the posterior samples or reproduce a complete MCMC experiment.
