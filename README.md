# Bounded coalescent

Code and data accompanying *Phylodynamic inference with the bounded coalescent: a point process perspective*. The posterior inference analysis and interpretation were conducted by Shuangping Li and Julia Palacios.  

The bounded coalescent conditions a genealogy on its time to the most recent common ancestor being less than a specified bound. This repository contains the simulation algorithms and statistical analyses described in the paper: maximum-likelihood estimation, posterior inference using simulated genealogies, and the Washington State COVID-19 example.

The inference experiments compare the bounded coalescent (BC) and standard coalescent (SC) likelihoods using random-integral (RI) and discretized methods.

## Contents

- [Figures and tables](#figures-and-tables)
- [Quick start](#quick-start)
- [Dependencies and running the analyses](#dependencies-and-running-the-analyses)
- [Repository layout](#repository-layout)
- [Data and reproducibility](#data-and-reproducibility)

## Figures and tables

Numbering follows the submitted manuscript. Each directory README links to the corresponding code and data.

| Paper result | Directory |
| --- | --- |
| Figures 2–3: simulation validation and comparison | [simulation](simulation/README.md) |
| Table 1: simulation timing | [simulation](simulation/README.md) |
| Figure 4: maximum-likelihood estimation | [analyses/maximum_likelihood](analyses/maximum_likelihood/README.md) |
| Figure 5 and Table 2: synthetic posterior inference | [analyses/synthetic](analyses/synthetic/README.md) |
| Figure 6: COVID-19 analysis | [analyses/covid](analyses/covid/README.md) |

The [analysis index](analyses/README.md) collects these analyses. Shared samplers, plotting code, and recorded results are in [analyses/reproduction](analyses/reproduction/README.md).

## Quick start

With Git and R installed, clone the repository and redraw the saved posterior results:

```bash
git clone https://github.com/bingjingle/boundedcoal.git
cd boundedcoal
bash analyses/reproduction/run_all.sh redraw
```

This uses base R to plot the stored coordinates without running MCMC. Output is written to `analyses/reproduction/deliverables/figures_R/`, including `panels.pdf` for the synthetic examples and `covid_2_DIS.pdf` for the COVID-19 example. The original result files are preserved.

An alternative output directory can be specified with `BC_FIGURES_OUT`:

```bash
BC_FIGURES_OUT=/path/to/output bash analyses/reproduction/run_all.sh redraw
```

## Dependencies and running the analyses

Posterior inference uses Python and R. The discretized method uses [JuliaPalacios/phylodyn](https://github.com/JuliaPalacios/phylodyn). Dependency versions and installation instructions are given in the [reproduction notes](docs/reproduction.md); Python requirements are listed in [requirements.txt](analyses/reproduction/requirements.txt). The squared-exponential RI implementation requires SciPy below 1.16.

The reproduction dependency check can be run from the repository root:

```bash
bash analyses/reproduction/preflight.sh
```

This checks the installed packages and phylodyn revision without installing or changing them. The figure redraw above does not require these additional inference dependencies.

For individual analyses, see the instructions for [simulation](docs/simulation.md), [maximum likelihood](docs/maximum-likelihood.md), [synthetic posterior inference](docs/synthetic.md), and [COVID-19](docs/covid.md). Full posterior runs can require substantial time and memory; the [reproduction notes](docs/reproduction.md) describe the run stages, settings, and output locations.

## Repository layout

```text
boundedcoal/
├── simulation/                 # Figures 2–3 and Table 1
├── analyses/
│   ├── maximum_likelihood/     # Figure 4
│   ├── synthetic/              # Figure 5 and Table 2
│   │   ├── data/               # Simulated genealogies
│   │   ├── bounded/            # Fits to bounded-coalescent data
│   │   └── standard/           # Fits to standard-coalescent data
│   ├── covid/                  # Figure 6
│   └── reproduction/
│       ├── samplers/           # Shared inference programs
│       ├── data/               # Inputs for the recorded runs
│       ├── coords/             # Saved posterior curves
│       ├── figures/            # Plotting programs
│       ├── figures_R/          # Saved manuscript-style figures
│       └── results/            # Table 2 and run records
└── docs/                       # Run instructions and reproducibility notes
```

## Data and reproducibility

The submitted synthetic posterior comparison uses 30 genealogies with 100 tips for each of three population-size trajectories: `1`, `3 exp(-t)`, and `25 exp(-5t)`, with respective bounds `1`, `0.7`, and `0.71`. The datasets are stored in [analyses/synthetic/data](analyses/synthetic/data/). The COVID-19 analysis uses the supplied [CCD0 genealogy](analyses/covid/data/median_ccd0.tree) for 103 sequences sampled in Washington State.

The saved coordinates support figure redraws. The original scripts and historical rerun drivers are also retained, but the latter do not implement the complete submitted Table 2 run schedule. Table 1 has partial benchmark code. These distinctions, dependency revisions, and the recorded Figure 4/5 dataset difference are described in the [reproducibility notes](docs/reproducibility.md).
