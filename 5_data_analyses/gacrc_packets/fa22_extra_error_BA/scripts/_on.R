library(tidyverse)

setwd("C:/Users/hl51981/OneDrive - University of Georgia/Leavitt_Herbert/PFFW/Pt Fourchon Food Webs R Directory/isotope_workflow/5_data_analyses/gacrc_packets/fa22_isotope_scales")

# 1) get all the .rds file paths
rds_paths <- list.files(
  path       = "outputs",
  pattern    = "\\.rds$", 
  recursive  = TRUE, 
  full.names = TRUE
)

# 2) build a tibble of metadata + filepaths
rds_tbl <- tibble(filepath = rds_paths) %>%
  mutate(
    filename = basename(filepath),
    # type is either "null_results" or "results"
    type     = str_extract(filename, "^(null_results|results)"),
    # extract species (the part before first "_")
    species = str_extract(
      filename,
      "(?<=^(null_results|results)_)[^_]+"
    ),
    # tef_step is what follows "tef_" up to the next "_"
    tef_step = str_extract(filename, "(?<=tef_)[^_]+"),
    # buf_val is the number following "buf_"
    buf_val  = str_extract(filename, "(?<=buf_)[0-9]+")
  )

# 3) pivot wider so each row has one null‐ and one real‐results path
paired_tbl <- rds_tbl %>%
  select(-filename) %>%
  pivot_wider(
    names_from  = type,
    values_from = filepath
  )
paired_tbl2 <- paired_tbl %>% 
  # rename the columns that currently hold file‐paths
  rename(
    null_path    = null_results,
    results_path = results
  )

final_tbl <- paired_tbl2 %>%
  mutate(
    null_model = map(null_path,    readRDS),
    stan_model = map(results_path, readRDS)
  )

# 1) Extract DIC from each model
dic_tbl <- final_tbl %>%
  mutate(
    DIC_null  = map_dbl(null_model,  ~ .x$BUGSoutput$DIC),   # or .x$BUGSoutput$DIC[1] if it's a vector
    DIC_model = map_dbl(stan_model,      ~ .x$BUGSoutput$DIC)
  )

# 2) Compute difference or ratio
dic_compare <- dic_tbl %>%
  transmute(
    species,
    tef_step,
    buf_val,
    DIC_null,
    DIC_model,
    delta_DIC = DIC_null - DIC_model,
    better     = if_else(DIC_model < DIC_null, "standard", "null")
  )

print(dic_compare)


# Prepare data for TEF step 1
plot_df <- dic_compare %>%
  # TEF step might be a character like "1.0"; convert or match accordingly
  filter(tef_step == "1.0") %>%
  # Create a signed delta for plotting:
  # if the null model was “better”, flip the sign
  mutate(
    buf_val    = as.numeric(buf_val)
  )

# Quick check
print(plot_df)

# Plot
ggplot(plot_df, aes(x = buf_val, y = delta_DIC)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_line() +
  geom_point(size = 2) +
  scale_x_continuous(breaks = plot_df$buf_val) +
  labs(
    title    = "Signed ΔDIC vs. Buffer Radius (TEF step 1.0)",
    x        = "Buffer radius (m)",
    y        = expression(Delta~DIC~"(null - standard)"),
    caption  = "Negative values indicate null model better"
  ) +
  theme_minimal(base_size = 14)
