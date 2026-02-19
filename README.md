# Bounded Coalescent: Synthetic Experiments

This repository contains all datasets and code used for the synthetic experiments in Sections 5.2 and 5.3.

---

## 📁 Directory Structure

- boundedcoal/
  - synthetic examples/
    - data/
    - codes/
      - standardcoaldata/
      - boundedcoaldata/


📊 Synthetic Datasets

All synthetic datasets are stored in:

boundedcoal/synthetic examples/data

File naming convention:

- syn_i_std_ntip_j.rda  
  Contains 30 datasets simulated from the standard coalescent likelihood  
  under effective population size trajectory Ne_i(t)  
  with j tips.

- syn_i_bounded_ntip_j.rda  
  Contains 30 datasets simulated from the bounded coalescent likelihood  
  under the same trajectory Ne_i(t)  
  with j tips.

---

## Section 5.2 – Standard Coalescent Experiments

Code location:

boundedcoal/synthetic examples/codes/standardcoaldata

All datasets here are simulated from the standard coalescent model.

Subfolders:

- ./RI_BM  
  Brownian motion kernels under:
  - Ne = 1
  - Ne(t) = 3 exp(-t)
  - Ne(t) = 25 exp(-5t)

- ./SE  
  Squared exponential kernels under the same three trajectories.

INLA implementation:

boundedcoal/synthetic examples/codes/standardcoaldata/RI_BM/INLA.R




---

## Section 5.3 – Bounded Coalescent Experiments

Code location:

boundedcoal/synthetic examples/codes/boundedcoaldata

All datasets here are simulated from the bounded coalescent model.

Subfolders:

- ./RI_BM  
  Brownian motion kernels under:
  - Ne = 1, tau = 1
  - Ne(t) = 3 exp(-t), tau = 0.7
  - Ne(t) = 25 exp(-5t), tau = 0.71

- ./SE  
  Squared exponential kernels under the same three scenarios.



### Python File Naming

- ****_**_bound_syn*_**_**.py  
  Inference under bounded coalescent likelihood

- ****_**_std_syn*_**_**.py  
  Inference under standard coalescent likelihood

---

## Reproducibility

Each dataset file contains 30 simulated genealogies.

Scripts are organized by:
- Model type (standard vs bounded)
- Kernel type (RI_BM vs SE)
- Effective population size trajectory
