# Synthetic posterior inference

Data and historical experiment scripts for comparing bounded and standard coalescent posterior inference. The submitted manuscript's synthetic posterior comparison is **Table 2**; the corresponding illustration is **Figure 5**. Table 1 concerns simulation timing and is outside this directory.

> The posterior inference analysis and interpretation were conducted by Shuangping Li and Julia Palacios.

Shuangping Li led the synthetic posterior-inference experiments. This statement describes responsibility for the analysis; it does not assign authorship of every historical source file.

## Contents

| Path | Contents |
| --- | --- |
| [`data/`](data/) | Bundled coalescent times, with 30 datasets per file |
| [`generate_datasets.R`](generate_datasets.R) | Historical generator for scenarios 1–3, with 50 and 100 tips |
| [`bounded/`](bounded/) | Bounded-data experiments, including bounded and standard likelihood fits |
| [`standard/`](standard/) | Experiments on standard-coalescent data, including RI and INLA comparisons |

**Reproduction status.** These files preserve the available implementation and exploratory runs. Their run settings are mixed: the BM scripts use 1,000,000 burn-in and 1,000,000 further iterations with thinning 10, while SE variants use shorter runs. The discrete scripts currently request 100,000 burn-in and 200,000 samples for both likelihoods. The submitted manuscript instead specifies 2,000,000 further iterations with thinning 20 for the standard discrete fit. These historical scripts do not by themselves provide a final table-aggregation pipeline or a verified end-to-end reproduction of Table 2/Figure 5. See the repository guide for the separate reproduction bundle. No numerical settings or algorithms were changed during reorganization.

## Data

Files follow `syn{scenario}_{bounded|std}_ntip{tips}.rda`. Each contains `coal_time`, a matrix with 30 rows and `tips − 1` columns; a row gives the ordered coalescent times for one simulated tree.

| Scenario | Effective population size | Bound for bounded data |
| --- | --- | --- |
| `syn1` | $N_e(t)=1$ | $\tau=1$ |
| `syn2` | $N_e(t)=3e^{-t}$ | $\tau=0.7$ |
| `syn3` | $N_e(t)=25e^{-5t}$ | $\tau=0.71$ |

The submitted comparison uses 100 tips and 30 datasets per scenario. The 50-tip datasets and `syn4_bounded_*` files are retained as additional historical material. The generator here does not define scenario 4, so those files are not assigned a manuscript result.

Python scripts select a dataset with a **zero-based index, 0–29**. The discrete R scripts use a **one-based index, 1–30**. Thus Python index `0` and R index `1` select the same row.

## Running an individual experiment

From the repository root, after installing the relevant scientific Python/R dependencies:

```bash
# Bounded data, bounded likelihood, Brownian-motion prior; first dataset.
python analyses/synthetic/bounded/RI_BM/boundedcoal_BM_bound_syn1_100_jitter8.py 0

# Same bounded dataset, discrete bounded and standard likelihood fits.
Rscript analyses/synthetic/bounded/DISCRETE/dis_syn1_ntip100.R 1

# Standard-coalescent data, Brownian-motion RI fit; first dataset.
python analyses/synthetic/standard/RI_BM/stdcoal_BMsyn1_100.py 0
```

These are full experiments and can take substantial time and memory. Python needs NumPy, SciPy, Matplotlib, pyreadr, and (for SE variants) statsmodels. R scripts use `ape`, `expm`, and a version of `phylodyn` supporting `mcmc_sampling(..., alg="bound_ESS")`; `INLA.R` additionally requires the INLA dependencies used by `phylodyn::BNPR`. The original batch examples used Python 3.9 and R 4.2. No complete version-locked environment accompanied these files.

Inputs default to this directory's `data/`; generated results default to `outputs/synthetic/` at the repository root. Paths no longer depend on an individual user's cluster account. Optional environment overrides:

```bash
export BOUNDEDCOAL_SYNTHETIC_DATA_DIR=/absolute/path/to/data
export BOUNDEDCOAL_SYNTHETIC_OUTPUT_DIR=/absolute/path/to/results
```

`generate_datasets.R` writes to `outputs/synthetic/generated_data/` by default, keeping the bundled datasets intact. Override `BOUNDEDCOAL_SYNTHETIC_GENERATED_DATA_DIR` to choose another location. Invoke the R files with `Rscript` so their locations can be resolved reliably.

The `.slurm` files are **site-specific batch templates**. Submit them from the directory containing the named Python/R script, and adapt the partition, module names, resource limits, and array range to your cluster. Output paths use the same environment overrides as local runs.

## Saved results and provenance

Python runs save `.npz` files with the original field names, including `c` (observed coalescent times), `d` (evaluation grid), `e` (truth), `f` (coverage), `j` (squared-error summary), `p` (runtime per 10,000 iterations), `q` (interval-width summary), and chain fields such as `aa`, `aaa`, and `s`. Consult each producing script for exact array shapes and the additional SE fields.

The discrete R scripts save both fitted objects and their timing, squared-error, coverage, and width summaries in `mcmc_results_{index}.rda`. Existing filenames, source citations, dataset contents, random seeds, and numerical routines have been retained. The former `synthetic examples/` material is organized here by data-generating model and inference method.
