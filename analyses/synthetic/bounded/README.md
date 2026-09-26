# Inference on bounded-coalescent data

These scripts load `../data/syn*_bounded_ntip*.rda` through the shared path configuration described in the [synthetic guide](../README.md).

| Folder | Purpose |
| --- | --- |
| `RI_BM/` | Random-integral inference with a Brownian-motion prior |
| `SE/` | Random-integral inference with a squared-exponential prior |
| `DISCRETE/` | Discrete fits of both bounded and standard likelihoods |
| `discrete_test/` | Earlier standalone checks using selected datasets |

In `RI_BM/`, filenames containing `_bound_` fit the bounded likelihood and those containing `_std_` fit the standard likelihood to the same bounded-data scenarios. In `SE/`, filenames containing `_std_` indicate the standard-likelihood fit; the other scripts fit the bounded likelihood. The `50`/`100` suffix records the tip count.

The scripts include exploratory settings as well as settings used in the manuscript. See the parent guide's reproduction-status note before using results as Table 2 or Figure 5 reproductions.
