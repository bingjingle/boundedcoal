# Bounded-coalescent simulation

This folder contains Algorithm 1, the simulation validation associated with
Figures 2–3, and partial timing code for Table 1. The synthetic **posterior
inference** experiments are documented separately in [`../analyses/`](../analyses).

| File | Purpose |
|---|---|
| [`our algorithm.R`](../simulation/our%20algorithm.R) | Time-transformation and thinning implementation; includes a 3,000-genealogy example |
| [`simul_bounded.R`](../simulation/simul_bounded.R) | `coalsim_bounded()` and experimental bounded skyline functions |
| [`plots.R`](../simulation/plots.R) | Standard versus bounded simulation comparisons and rejection-sampling validation, Figures 2–3 |
| [`comparison of three methods.R`](../simulation/comparison%20of%20three%20methods.R) | Partial Table 1 timing comparisons: Carson's implementation, thinning, and naive rejection |
| [`data/`](../simulation/data) | Two saved 3,000-replicate simulations, with 100 tips and bound `0.5` |

## Scope of the available code

Table 1 code is present, but its original header explicitly marks it as partial:
only selected scenarios are implemented. There is no complete table-generation
driver or stored timing result for every row. No missing benchmark result has
been reconstructed during the repository reorganization.

The two saved `.rda` files contain validation simulations; older notes attached
draft figure/table numbers to them. Use the paper map above for the submitted
manuscript instead of those historical labels.

## Dependencies and use

The simulator uses `ape` and `phylodyn`. Plotting additionally uses `ggplot2`,
`cowplot`, and `covalchemy`; the timing comparison uses `BoundedCoalescent`.
The bounded `phylodyn` routines require the version linked in the
[reproduction guide](../analyses/reproduction/README.md).

For a small simulator example, run from the repository root:

```r
source("simulation/simul_bounded.R")
set.seed(1)
coalsim_bounded(0, 10, function(t) rep(1, length(t)), bound = 0.5)
```

The original simulator explicitly notes that its heterochronous sampling-time
case is not functional; the example above uses contemporaneous sampling.
`plots.R` assumes `coalsim_bounded()` has already been loaded, contains an
interactive comparison using undefined `resulta` and `resultb` objects, and
does not save a final manuscript figure automatically. Treat it as a source
record requiring section-by-section execution. The timing comparison requires
its external packages and can be expensive; its output directory is controlled
by `BC_SIMULATION_OUTPUT` (default: `generated` under the working directory).
