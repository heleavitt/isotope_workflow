# Pt. Fourchon Isotope Workflow

## Overview
This repository contains the batch scripts, input data, and R utilities that support the 2016-2022 Pt. Fourchon isotope comparison project. The workflow centers on running MixSIAR mixing models on the Georgia Advanced Computing Resource Center (GACRC) cluster, followed by local post-processing to summarise posterior draws, compare model structures, and evaluate isotopic niche hypervolumes.

## Repository Layout
- `5_data_analyses/gacrc_packets/fa22_extra_error_BA/data/` - Input files grouped into `mix_files`, `source_files`, and `tef_files`.
- `5_data_analyses/gacrc_packets/fa22_extra_error_BA/scripts/` - Slurm submission helpers (`*.sh`) and analysis scripts (`*.R`).
- `README.md` - Project overview and operating instructions (this document).

Outputs produced by the cluster (`outputs/`, `results/`, and log files) are intentionally left out of version control because they are large and run-specific.

## Software Requirements
### Core runtimes
- Bash 4 or newer together with Slurm (tested on GACRC)
- R 4.3.x
- JAGS 4.3.1

### R packages
Install these packages before running the scripts: `MixSIAR`, `tidyverse` (dplyr, ggplot2, purrr, readr, stringr, tidyr), `hypervolume`, `fs`, `R2jags`, `scales`, `reshape2`, `dunn.test`, and `loo`.

On GACRC the batch script `fourchon_mixing_model_hypervolume.sh` loads `foss/2022a`, `R/4.3.1-foss-2022a`, and `JAGS/4.3.1-foss-2022a`. Update the module list if the environment changes.

## Input Data Organisation
Each model run combines one file from each data subfolder:
- `mix_files/edge_<edge>_buf_<buffer>_<species>_mix.csv` - Animal isotope measurements and covariates for a specific edge distance and buffer radius.
- `source_files/<species>_mix.csv` - Source isotope measurements with concentration dependencies.
- `tef_files/tef_<species>_step_<value>.csv` - Trophic enrichment factors evaluated at multiple step sizes.

Before launching jobs, copy the `data/` and `scripts/` directories to `/work/jnlab/fourchon_isotope_runs/fa22_extra_error_BA/` (or the appropriate run folder) so the absolute paths in the scripts resolve correctly.

## Running the MixSIAR Models on GACRC
1. Stage the run directory under `/work/jnlab/fourchon_isotope_runs/` and confirm that `data/` and `scripts/` contain the desired files.
2. Generate the array definition file:
   ```bash
   bash scripts/combinations.sh
   ```
   This scans `data/tef_files` and `data/mix_files`, writes `scripts/combinations.txt` (one line per species/TEF/buffer combination), and logs any missing files to `scripts/logs/missing_files.log`.
3. Inspect `scripts/combinations.txt` to confirm that every intended pairing exists. Delete and rerun the script if you need to regenerate the list.
4. Submit the Slurm array job:
   ```bash
   bash scripts/submit_array_jobs.sh
   ```
   The wrapper calculates the array length from `combinations.txt` and launches `fourchon_mixing_model_hypervolume.sh`. Each array task runs `Rscript scripts/fourchon_mixing_hv_scales.r <species> <tef> <edge> <buffer>`.
5. Monitor progress with `squeue -u <user>` or `sacct -j <jobid>` and review per-task logs in `scripts/logs/`.
6. Collect the outputs once the jobs finish. The R script saves `.rds`, `.csv`, and figure files to `outputs/<species>/tef_step_<value>/<buffer>/` within the run directory.

## Post-processing Scripts
All scripts assume the working directory is `5_data_analyses/gacrc_packets/fa22_extra_error_BA` unless otherwise noted.
- `model_outputs.R` converts each `results_*.rds` (and related models) into posterior draw tables, source contribution summaries, and trophic level estimates per buffer.
- `confidence_intervals.R` rebuilds MixSIAR objects and produces custom continuous covariate plots (resource contribution vs. percent mangrove cover) for every model, saving PNGs under `results/<species>/<tef>/<buffer>/`.
- `compare_models.R` pairs `results_*.rds` with `null_results_*.rds`, computes Akaike weights, and plots how model support changes with buffer size.
- `hypervolumes.R` and `hypervolumes_with _boxplot.R` offer prototype code for isotopic niche hypervolume plots and contribution boxplots using the CSVs created by `model_outputs.R`. Review and clean these scripts before reuse.
- `hypervolume_bootstrap.r` bootstraps hypervolume size and Sorensen overlap between low- and high-mangrove habitats, writing confidence intervals to `results/hv_bootstrap/`.
- `_on.R` and `_on-MSMARS-60014H.R` are ad-hoc DIC comparison utilities retained for historical reference.

Before running any R script, double-check the hard-coded `setwd()` or absolute paths so that they point to your local clone or the cluster workspace.

## Customising For New Runs
- Update the `run_name` constant and absolute paths inside `fourchon_mixing_hv_scales.r` (and other scripts) when creating a new packet.
- Edit `SPECIES_CODES` in `combinations.sh` if the focal species list changes.
- The path near the top of `combinations.sh` currently references `/work/jnlab/fourchon_isotope_runs/fa22_isotope_scales_extra_error/`; adjust it to match the folder that holds the data for the new run.
- Modify the Slurm directives in `fourchon_mixing_model_hypervolume.sh` (memory, walltime, CPU count, email) to fit future workloads.

## Troubleshooting
- An empty `scripts/combinations.txt` usually indicates a filename mismatch in `data/`.
- If `Rscript` fails inside the array job, check the corresponding `scripts/logs/mix_hypervolume_<job>-<task>.err` file for the error message.
- MixSIAR requires JAGS on the compute node. The job will terminate immediately if the module is missing or incompatible.
- The hypervolume scripts expect balanced sample sizes between low and high mangrove groups. If they exit early, verify that the `mixing_model_df_*.csv` files contain enough observations per habitat.

## Contacts And Attribution
Maintained by the Pt. Fourchon Food Webs research group at the University of Georgia. Cite MixSIAR (Stock et al. 2018) and the hypervolume package (Blonder 2018) when publishing results derived from this workflow.
