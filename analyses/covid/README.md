# COVID-19 example · Figure 6

This folder contains the genealogy and original analysis source for the
Washington State COVID-19 example. The COVID-19 analysis was conducted by
Julia Palacios.

| Item | Contents |
|---|---|
| [`data/median_ccd0.tree`](data/median_ccd0.tree) | CCD0 genealogy with 103 tips sampled on June 8, 2020 |
| [`legacy/covid_plot.R`](legacy/covid_plot.R) | Original skyline, reported-case, and discretized posterior exploration |
| [`legacy/covid_analysis.rtf`](legacy/covid_analysis.rtf) | Historical BEAST analysis notes |
| [`../reproduction/`](../reproduction/) | Portable posterior pipeline, saved coordinates, and plotting entrypoints |

The original analysis scales all tree branches by `0.001 / 0.0012`; the scaled
tree height is approximately `0.5858767` years. It sets the coalescent bound to
that height. The saved posterior coordinates and redraw workflow are documented
in the [reproduction guide](../reproduction/README.md).

The legacy R script is an exploratory source record: it has an unmatched final
brace, depends on an interactive R session, and contains plotting fragments
using undefined objects. It is not a standalone reproduction command.
The `BC_COVID_TREE` and `BC_COVID_COORDS` environment variables can override its
input tree and output CSV paths; defaults are relative to the repository root.
The historical RTF notes include an inconsistent mutation-rate exponent in a draft
paragraph; consult the manuscript for the final scientific description.

The original BEAST XML, raw sequence alignment, and full BEAST posterior tree
sample are not included here. This folder supports analysis of the supplied
CCD0 genealogy; it does not independently reproduce the upstream sequence-to-tree
analysis.
