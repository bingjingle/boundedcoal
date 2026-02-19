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


🔬 Section 5.2: Standard Coalescent Experiments
All code is located at: boundedcoal/synthetic examples/codes/standardcoaldata
All datasets in this directory are simulated from the standard coalescent model.

*******Kernel Types
1️⃣ Brownian Motion Kernels
./RI_BM
Effective population size trajectories:
𝑁𝑒=1; 𝑁𝑒(𝑡)=3 𝑒xp(-t); 𝑁𝑒(𝑡)=25 𝑒xp(-5t).

2️⃣ Squared Exponential Kernels
./SE
Uses the same three trajectories as above.

*******INLA Implementation
The INLA method is implemented at: boundedcoal/synthetic examples/codes/standardcoaldata/RI_BM/INLA.R.






🔬 Section 5.3: Bounded Coalescent Experiments
All code is located at: boundedcoal/synthetic examples/codes/boundedcoaldata
All datasets in this directory are simulated from the bounded coalescent model.


###########Kernel Types
1️⃣ Brownian Motion Kernels
./RI_BM
Scenarios:
𝑁𝑒=1 𝜏=1; 𝑁𝑒(𝑡)=3 𝑒xp(-t) 𝜏=0.7; 𝑁𝑒(𝑡)=25 𝑒xp(-5t) 𝜏=0.71.

2️⃣ Squared Exponential Kernels
./SE
Uses the same three scenarios listed above.


##########Python File Naming Convention
****_**_bound_syn*_**_**.py
→ Inference under the bounded coalescent likelihood

****_**_std_syn*_**_**.py
→ Inference under the standard coalescent likelihood







🧪 Reproducibility

Each dataset file contains 30 simulated genealogies.

All scripts are organized by:

Model type (standard vs bounded)

Kernel type (RI_BM vs SE)

Effective population size trajectory
