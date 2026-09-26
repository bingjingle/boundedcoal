# Historical discrete-method checks

These standalone scripts are retained for exploratory numerical checks; they are not the submitted Table 2 aggregation pipeline.

| Script | Likelihood | Bound/data |
| --- | --- | --- |
| `test_BC_07.R` | Bounded | $\tau=0.7$; times embedded in the script |
| `test_SC_07.R` | Standard | Same embedded-data setting |
| `test_BC_15.R` | Bounded | $\tau=1.5$; bundled `test100_samplesize2_tau1.5.rda` |
| `test3_SC_15.R` | Standard | Same bundled-data setting |

The trajectory is scenario 2. `bound`, `Ngrid`, and `noise_var` are set in the scripts. Run using `Rscript`; the two scripts with external data resolve the bundled file relative to their own location. Full chains have not been rerun during reorganization.
