library(tidyverse)
library(MixSIAR)

# 1. List all RDS model paths
rds_paths <- list.files("outputs", pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)

rds_tbl <- tibble(filepath = rds_paths) %>%
  mutate(
    filename = basename(filepath),
    type     = str_extract(filename, "^(null_results|results|ind)"),
    species_code = str_extract(filename, "(?<=^(null_results|results|ind_results)_)[^_]+"),
    tef_step     = str_extract(filename, "(?<=tef_)[^_]+"),
    buf_val      = str_extract(filename, "(?<=buf_)[0-9]+")
  )

# 3. Process each file
for (i in seq_len(nrow(rds_tbl))) {
  row <- rds_tbl[i, ]
  rds_path <- row$filepath
  species_code <- row$species_code
  tef_step <- row$tef_step
  buf_val <- row$buf_val
  model_type <- row$type
  
  output_dir <- file.path("outputs", species_code, paste0("tef_step_", tef_step), buf_val)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Build input file paths
  mix_path <- file.path("data", "mix_files", sprintf("edge_3_buf_%s_%s_mix.csv", buf_val, species_code))
  source_path <- file.path("data", "source_files", sprintf("%s_mix.csv", species_code))
  
  if (!file.exists(mix_path) || !file.exists(source_path)) {
    message("Missing input files for ", species_code, " buf=", buf_val)
    next
  }
  
  # Load model + inputs
  results <- readRDS(rds_path)
  mix_df <- read_csv(mix_path, show_col_types = FALSE)
  source_df <- read_csv(source_path, show_col_types = FALSE)
  

  if(model_type == "ind"){
    
    mix <- load_mix_data(
      filename = mix_path, 
      iso_names = c("d13C", "d34S"), 
      factors = "wrapping_inventory.vial_code",
      fac_random = FALSE,
      fac_nested = FALSE,
      cont_effects = NULL
    )
    
    source <- load_source_data(
      filename = source_path,
      conc_dep = TRUE,
      data_type = "raw",
      mix = mix
      )
    
    # Extract p.ind
    p_array <- results$BUGSoutput$sims.list$p.fac1
    

    # Optional: calculate TL if d15N present
    if ("d15N" %in% names(mix_df)) {
      animal_d15N <- mix_df %>% mutate(ind = row_number()) %>% select(ind, c_d15N = d15N)
      
      if ("d15N" %in% names(source_df)) {
        s_d15N_df <- source_df %>%
          group_by(Source) %>%
          summarize(s_d15N = mean(d15N, na.rm = TRUE), .groups = "drop")
      } else if ("Meand15N" %in% names(source_df)) {
        s_d15N_df <- source_df %>%
          select(Source, s_d15N = Meand15N)
      } else {
        stop("No d15N values found for TL calculation.")
      }

      
      posterior_samples <- as.data.frame.table(results$BUGSoutput$sims.list$p.fac1) %>%
        mutate(
          Source = source$source_names[as.integer(Var3)],
          ind = as.integer(Var2),
          run = as.integer(Var1),
          Contr = Freq,
          species_code = species_code
        ) %>%
        select(species_code, Source, Contr, ind, run)
      
      df_mm <- posterior_samples %>%
        left_join(s_d15N_df, by = "Source") %>%
        left_join(animal_d15N, by = "ind") %>%
        group_by(species_code, ind, run) %>%
        mutate(TL = (c_d15N - sum(Contr * s_d15N)) / 3.4 + 1) %>%
        ungroup() %>%
        distinct()
      covariate_df <- mix_df %>%
        mutate(ind = row_number()) %>%
        select(ind, edge_l.mangrove)
      
      df_mm <- df_mm %>%
        left_join(covariate_df, by = "ind")
      # Summaries
      summary_contr <- posterior_samples %>%
        group_by(species_code, Source) %>%
        summarize(
          mean_Contr = mean(Contr),
          sd_Contr = sd(Contr),
          p5_Contr = quantile(Contr, 0.05),
          p95_Contr = quantile(Contr, 0.95),
          .groups = 'drop'
        )
      
      summary_TL <- df_mm %>%
        group_by(species_code) %>%
        summarize(
          mean_TL = mean(TL),
          sd_TL = sd(TL),
          p5_TL = quantile(TL, 0.05),
          p95_TL = quantile(TL, 0.95),
          .groups = 'drop'
        )
      file.exists(output_dir)
      dir.exists(output_dir)
      file.access(output_dir, mode = 2)  # 0 = ok, -1 = no write access
      write_csv(posterior_samples, file.path(output_dir, sprintf("posterior_samples_preTL_%s_tef_%s_buf_%s.csv", species_code, tef_step, buf_val)))
      write_csv(df_mm, file.path(output_dir, sprintf("mixing_model_df_%s_tef_%s_buf_%s.csv", species_code, tef_step, buf_val)))
      write_csv(summary_contr, file.path(output_dir, sprintf("src_cont_summary_%s_tef_%s_buf_%s.csv", species_code, tef_step, buf_val)))
      write_csv(summary_TL, file.path(output_dir, sprintf("trophic_levels_summary_%s_tef_%s_buf_%s.csv", species_code, tef_step, buf_val)))
    }
  }
}
  