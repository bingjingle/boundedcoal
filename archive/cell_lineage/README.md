# Archived cell-lineage analysis

Historical UPGMA cell-lineage material retained from the former `realexample/` directory. This is a separate analysis from the submitted manuscript's COVID-19 example and is not presented as a reproduction of a current manuscript figure or table.

| Path | Contents |
| --- | --- |
| `codes/` | Parameterized BM/SE, random-integral/discrete inference scripts and cluster job templates |
| `arxiv/` | Earlier Python scripts and exploratory notebooks; the original folder name is retained |
| `UGPMA100new.rda` | `bt_adj`: 99 adjusted coalescent times for the 100-tip tree |
| `UGPMA3257new.rda` | `bt_adj`: 3,256 adjusted coalescent times for the 3,257-tip tree |
| `UPGMAtree_100.txt` | Tree input for the 100-tip analysis |
| `UPGMAtree3257.R` | Historical tree reconstruction/preparation script |
| `INLA_realexp.R` | INLA comparison for the cell-lineage data |
| `real100_all_settings_posterior_HPD_midpoint_traceplots_with_events (1).ipynb` | Posterior/trace plotting notebook |

## Usage and limitations

The Python scripts use bundled inputs by default and write to `outputs/cell_lineage/` at the repository root. Override `BOUNDEDCOAL_CELL_LINEAGE_DATA_DIR` or `BOUNDEDCOAL_CELL_LINEAGE_OUTPUT_DIR` to use other locations. The parameterized programs in `codes/` expose their run settings through `--help`; SE variants also accept `--input-file`, and all four accept `--output-dir`.

Run Python/R files by path, and use `Rscript` for the R files. Launch notebooks with the **repository root as the working directory**; their relative data/result paths assume that location. The Slurm files are site-specific templates: submit from `codes/` and adjust partitions, modules, resources, and settings before use.

`UPGMAtree3257.R` references `Supplementary_File_2_DataTableMOI19.csv` and `shared_edit_matrix_3257.csv`, which are not included in this repository. Rebuilding that tree therefore requires the original external inputs. The supplied adjusted coalescent-time files remain available for the inference programs. Posterior `.npz` outputs expected by the plotting notebook are not bundled.

These files retain historical numerical routines, settings, notebook outputs, and source citations. Machine-specific input/output paths were replaced with repository-relative defaults. Printed personal-directory paths in the plotting notebook were anonymized; its numerical outputs remain historical. No full inference run or scientific validation of this archived analysis was performed during reorganization.
