#!/usr/bin/env Rscript

# Load required libraries
library(MixSIAR)
library(dplyr)
library(tidyverse)
library(hypervolume)

# Get input arguments: species code and TEF step
args <- commandArgs(trailingOnly = TRUE)
species_code <- tolower(args[1])  # Ensure species code is lowercase
tef_step <- args[2]  # TEF step value (e.g., "1.0", "2.5")
edge_val <- args[3]
buf_val <- args[4]
run_name <- "fa22_extra_error_BA"


# Define output directory (include TEF step)
output_dir <- file.path("/work/jnlab/fourchon_isotope_runs", run_name,"outputs", species_code, paste0("tef_step_", tef_step), buf_val)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cat("Final output directory:", output_dir, "\n")

# Construct paths for mix, source, and TEF files
mix_file <- file.path("/work/jnlab/fourchon_isotope_runs", run_name,"data/mix_files", 
                      paste0("edge_", edge_val, "_buf_",buf_val, "_", species_code, "_mix.csv"))

source_file <- file.path("/work/jnlab/fourchon_isotope_runs", run_name,"data/source_files", 
                         paste0(species_code, "_mix.csv"))

tef_file <- file.path("/work/jnlab/fourchon_isotope_runs", run_name,"data/tef_files", 
                      paste0("tef_", species_code, "_step_", tef_step, ".csv"))

# Ensure that the paired files exist
if (!file.exists(mix_file)) stop(paste("ERROR: Mix file not found:", mix_file))
if (!file.exists(source_file)) stop(paste("ERROR: Source file not found:", source_file))
if (!file.exists(tef_file)) stop(paste("ERROR: TEF file not found:", tef_file))

# Log paths for debugging
cat("Mix file path:", mix_file, "\n")
cat("Source file path:", source_file, "\n")
cat("TEF file path:", tef_file, "\n")

# Function to safely read CSV files
read_safe <- function(file) {
  tryCatch(read.csv(file), error = function(e) {
    cat("Error reading", file, ":", conditionMessage(e), "\n")
    NULL
  })
}

# Load data
animals <- read_safe(mix_file)
producers <- read_safe(source_file)
tef <- read_safe(tef_file)

# Check if any data is missing
if (is.null(animals) || is.null(producers) || is.null(tef)) {
  stop("Missing or invalid input files.")
}

cat("N total:", nrow(animals), "\n")
cat("N complete for d13C & d34S:", sum(complete.cases(animals[, c("d13C", "d34S")])), "\n")
cat("N complete with edge_l.mangrove:", sum(complete.cases(animals[, c("d13C", "d34S", "edge_l.mangrove")])), "\n")


# Prepare MixSIAR data
mix <- load_mix_data(
  filename = mix_file, 
  iso_names = c("d13C", "d34S"), 
  factors = NULL,
  fac_random = NULL,
  fac_nested = NULL,
  cont_effects = "edge_l.mangrove"
)

sources <- load_source_data(
  filename = source_file, 
  conc_dep = TRUE, 
  data_type = "raw", 
  mix = mix
)

discr <- load_discr_data(
  filename = tef_file, 
  mix = mix
)

######## N-Mixes ############
N_mix <- load_mix_data(
  filename = mix_file, 
  iso_names = c("d13C", "d15N"), 
  factors = NULL,
  fac_random = NULL,
  fac_nested = NULL,
  cont_effects = "edge_l.mangrove"
)

N_sources <- load_source_data(
  filename = source_file, 
  conc_dep = TRUE, 
  data_type = "raw", 
  mix = N_mix
)

N_discr <- load_discr_data(
  filename = tef_file, 
  mix = N_mix
)

##### NULL Mixes #########
null_mix <- load_mix_data(
  filename = mix_file, 
  iso_names = c("d13C", "d34S"), 
  factors = NULL,
  fac_random = NULL,
  fac_nested = NULL,
  cont_effects = NULL
)

null_sources <- load_source_data(
  filename = source_file, 
  conc_dep = TRUE, 
  data_type = "raw", 
  mix = null_mix
)

null_discr <- load_discr_data(
  filename = tef_file, 
  mix = null_mix
)

# IND Mixes
ind_mix <- load_mix_data(
  filename = mix_file, 
  iso_names = c("d13C", "d34S"), 
  factors = "wrapping_inventory.vial_code",
  fac_random = FALSE,
  fac_nested = FALSE,
  cont_effects = NULL
)

ind_sources <- load_source_data(
  filename = source_file, 
  conc_dep = TRUE, 
  data_type = "raw", 
  mix = ind_mix
)

ind_discr <- load_discr_data(
  filename = tef_file, 
  mix = ind_mix
)


# Write and run the JAGS model
model_file <- file.path(output_dir, paste0("model_", species_code, "_tef_", tef_step, ".txt"))
null_model_file <- file.path(output_dir, paste0("null_model_", species_code, "_tef_", tef_step, ".txt"))
ind_model_file <- file.path(output_dir, paste0("ind_model_", species_code, "_tef_", tef_step, ".txt"))

# Write the JAGS model with corrected argument usage
write_JAGS_model(
  filename = model_file, 
  resid_err = TRUE,   # Residual error disabled
  process_err = TRUE,  # Process error enabled
  mix = mix, 
  source = sources
)

write_JAGS_model(
  filename = null_model_file, 
  resid_err = TRUE,   # Residual error disabled
  process_err = TRUE,  # Process error enabled
  mix = null_mix, 
  source = null_sources
)

write_JAGS_model(
  filename = ind_model_file, 
  resid_err = FALSE,   # Residual error disabled
  process_err = TRUE,  # Process error enabled
  mix = ind_mix, 
  source = ind_sources
)




# Run the JAGS model with corrected argument usage
results <- run_model(
  run = "very long", 
  mix = mix, 
  source = sources, 
  discr = discr, 
  model_filename = model_file,
  alpha.prior = 1,
  resid_err = TRUE,   # Include resid_err
  process_err = TRUE   # Include process_err
)

# Run the JAGS model with corrected argument usage
null_results <- run_model(
  run = "very long", 
  mix = null_mix, 
  source = null_sources, 
  discr = null_discr, 
  model_filename = null_model_file,
  alpha.prior = 1,
  resid_err = TRUE,   # Include resid_err
  process_err = TRUE   # Include process_err
)

# Run the JAGS model with corrected argument usage
ind_results <- run_model(
  run = "very long", 
  mix = ind_mix, 
  source = ind_sources, 
  discr = ind_discr, 
  model_filename = ind_model_file,
  alpha.prior = 1,
  resid_err = TRUE,   # Include resid_err
  process_err = TRUE   # Include process_err
)


saveRDS(null_results, file.path(output_dir, sprintf("null_results_%s_tef_%s_buf_%s.rds", species_code, tef_step, buf_val)))
saveRDS(results, file.path(output_dir, sprintf("results_%s_tef_%s_buf_%s.rds", species_code, tef_step, buf_val)))
saveRDS(ind_results, file.path(output_dir, sprintf("ind_results_%s_tef_%s_buf_%s.rds", species_code, tef_step, buf_val)))

# Save the current working directory
old_wd <- getwd()

# Set the working directory to the output directory
setwd(output_dir)
# Generate isospace plot and save to output directory
plot_data_filename <- file.path(output_dir, sprintf("isospace_plot_%s_tef_%s.png", species_code, tef_step))

png(filename = plot_data_filename, width = 800, height = 600)
plot_data(
  mix = N_mix,
  source = N_sources,
  discr = N_discr,
  plot_save_pdf = FALSE,
  plot_save_png = FALSE  # Don't let it save automatically
)
dev.off()

cat("Isospace plot saved to:", plot_data_filename, "\n")


# Output diagnostics and results
output_JAGS(
  results, 
  mix, 
  sources, 
  output_options = list(
    summary_save = TRUE,
    summary_name = sprintf("summary_statistics_%s_tef_%s", species_code, tef_step),
    sup_post = TRUE,
    plot_post_save_png = FALSE,
    plot_post_name = sprintf("posterior_density_%s_tef_%s", species_code, tef_step),
    sup_pairs = TRUE,
    plot_pairs_save_png = FALSE,
    plot_pairs_name = sprintf("pairs_plot_%s_tef_%s", species_code, tef_step),
    sup_xy = FALSE,
    plot_xy_save_png = FALSE,
    plot_xy_name = sprintf("xy_plot_%s_tef_%s", species_code, tef_step),
    gelman = TRUE,
    heidel = FALSE,
    geweke = TRUE,
    diag_save = TRUE,
    diag_name = sprintf("diagnostics_%s_tef_%s", species_code, tef_step),
    indiv_effect = FALSE,
    plot_post_save_pdf = FALSE,
    plot_pairs_save_pdf = FALSE,
    plot_xy_save_pdf = FALSE
  )
)

# Output diagnostics and results
output_JAGS(
  null_results, 
  null_mix, 
  null_sources, 
  output_options = list(
    summary_save = TRUE,
    summary_name = sprintf("null_summary_statistics_%s_tef_%s", species_code, tef_step),
    sup_post = TRUE,
    plot_post_save_png = FALSE,
    plot_post_name = sprintf("null_posterior_density_%s_tef_%s", species_code, tef_step),
    sup_pairs = TRUE,
    plot_pairs_save_png = FALSE,
    plot_pairs_name = sprintf("null_pairs_plot_%s_tef_%s", species_code, tef_step),
    sup_xy = FALSE,
    plot_xy_save_png = FALSE,
    plot_xy_name = sprintf("null_xy_plot_%s_tef_%s", species_code, tef_step),
    gelman = TRUE,
    heidel = FALSE,
    geweke = TRUE,
    diag_save = TRUE,
    diag_name = sprintf("null_diagnostics_%s_tef_%s", species_code, tef_step),
    indiv_effect = FALSE,
    plot_post_save_pdf = FALSE,
    plot_pairs_save_pdf = FALSE,
    plot_xy_save_pdf = FALSE
  )
)

output_JAGS(
  ind_results, 
  ind_mix, 
  ind_sources, 
  output_options = list(
    summary_save = TRUE,
    summary_name = sprintf("null_summary_statistics_%s_tef_%s", species_code, tef_step),
    sup_post = TRUE,
    plot_post_save_png = FALSE,
    plot_post_name = sprintf("null_posterior_density_%s_tef_%s", species_code, tef_step),
    sup_pairs = TRUE,
    plot_pairs_save_png = FALSE,
    plot_pairs_name = sprintf("null_pairs_plot_%s_tef_%s", species_code, tef_step),
    sup_xy = FALSE,
    plot_xy_save_png = FALSE,
    plot_xy_name = sprintf("null_xy_plot_%s_tef_%s", species_code, tef_step),
    gelman = TRUE,
    heidel = FALSE,
    geweke = TRUE,
    diag_save = TRUE,
    diag_name = sprintf("null_diagnostics_%s_tef_%s", species_code, tef_step),
    indiv_effect = FALSE,
    plot_post_save_pdf = FALSE,
    plot_pairs_save_pdf = FALSE,
    plot_xy_save_pdf = FALSE
  )
)

