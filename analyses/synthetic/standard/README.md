# Inference on standard-coalescent data

`RI_BM/` contains Brownian-motion random-integral fits and `INLA.R`, an INLA comparison across scenarios and tip counts. `SE/` contains squared-exponential random-integral fits. Inputs come from `../data/syn*_std_ntip*.rda` through the shared path configuration.

[`Figure2.py`](Figure2.py) is an older comparison-plot script. Its filename is historical and does **not** identify Figure 2 in the submitted manuscript. It reads the corresponding standard-data BM and SE `.npz` outputs, using the dataset indices embedded in the script, and writes PDFs beneath `outputs/synthetic/figures/`. The required posterior output files are not bundled.

See the [synthetic guide](../README.md) for dependencies, index conventions, output locations, and reproduction limits.
