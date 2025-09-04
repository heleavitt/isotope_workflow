library(tidyverse)
library(fs)
library(MixSIAR)

plot_continuous_var_custom <- function(jags.1, mix, source, output_options, alphaCI = 0.05, 
                                       exclude_sources_below = 0.1) {
  R2jags::attach.jags(jags.1)
  n.sources <- source$n.sources
  source_names <- source$source_names
  
  # Optional: modify source names here
  source_names <- recode(source_names,
                         "CRAVIR" = "POM",
                         "live spartina leaves and stems" = "Spartina",
                         "live mangrove leaves" = "Mangrove",
                         "bostrychia" = "Bostrychia",
                         "benthic algae" = "Benthic algae"
                         
)
  
  for (ce in 1:mix$n.ce) {
    label <- mix$cont_effects[ce]
    cont <- mix$CE[[ce]]
    ilr.cont <- get(paste("ilr.cont", ce, sep = ""))
    n.plot <- 200
    chain.len <- dim(p.global)[1]
    Cont1.plot <- seq(from = round(min(cont), 1), 
                      to = round(max(cont), 1), length.out = n.plot)
    
    ilr.plot <- array(NA, dim = c(n.plot, n.sources - 1, chain.len))
    for (src in 1:(n.sources - 1)) {
      for (i in 1:n.plot) {
        ilr.plot[i, src, ] <- ilr.global[, src] + ilr.cont[, src] * Cont1.plot[i]
      }
    }
    
    e <- matrix(rep(0, n.sources * (n.sources - 1)), 
                nrow = n.sources, ncol = (n.sources - 1))
    for (i in 1:(n.sources - 1)) {
      e[, i] <- exp(c(rep(sqrt(1/(i * (i + 1))), i), 
                      -sqrt(i/(i + 1)), rep(0, n.sources - i - 1)))
      e[, i] <- e[, i]/sum(e[, i])
    }
    
    cross <- array(data = NA, dim = c(n.plot, chain.len, n.sources, n.sources - 1))
    tmp <- array(data = NA, dim = c(n.plot, chain.len, n.sources))
    p.plot <- array(data = NA, dim = c(n.plot, chain.len, n.sources))
    for (i in 1:n.plot) {
      for (d in 1:chain.len) {
        for (j in 1:(n.sources - 1)) {
          cross[i, d, , j] <- (e[, j]^ilr.plot[i, j, d]) / sum(e[, j]^ilr.plot[i, j, d])
        }
        for (src in 1:n.sources) {
          tmp[i, d, src] <- prod(cross[i, d, src, ])
        }
        for (src in 1:n.sources) {
          p.plot[i, d, src] <- tmp[i, d, src] / sum(tmp[i, d, ])
        }
      }
    }
    
    get_high <- function(x) quantile(x, 1 - alphaCI / 2)
    get_low <- function(x) quantile(x, alphaCI / 2)
    
    p.low <- apply(p.plot, c(1, 3), get_low)
    p.high <- apply(p.plot, c(1, 3), get_high)
    p.median <- apply(p.plot, c(1, 3), median)
    colnames(p.median) <- source_names
    
    Cont1.plot <- Cont1.plot * mix$CE_scale + mix$CE_center
    
    df <- data.frame(
      reshape2::melt(p.median)[, 2:3],
      rep(Cont1.plot, n.sources),
      reshape2::melt(p.low)[, 3],
      reshape2::melt(p.high)[, 3]
    )
    colnames(df) <- c("source", "median", "x", "low", "high")
    df$source <- factor(df$source, levels = source_names)
    
    rm.srcs <- apply(p.median, 2, function(x) all(x < exclude_sources_below))
    df <- subset(df, source %in% source_names[!rm.srcs])
    
 
    return(
      ggplot(df, aes(x = x, y = median, color = source, fill = source)) +
        geom_line(size = 1.2) +
        geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.3, color = NA) +
        labs(
          x = "% Mangrove Cover",
          y = "Proportional Contribution",
          color = "Resource",
          fill = "Resource"
        ) +
        scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
        scale_y_continuous(limits = c(0, 0.9)) +
        theme_classic(base_size = 14) +
        scale_color_manual(values = c(
          "POM" = "#D81B60",
          "Spartina" = "#1E88E5",
          "Mangrove" = "#DC9203",
          "Benthic algae" = "#004D40"
        )) +
        scale_fill_manual(values = c(
          "POM" = "#D81B60",
          "Spartina" = "#1E88E5",
          "Mangrove" = "#DC9203",
          "Benthic algae" = "#004D40"
        )) +
        theme(legend.position = "right")
    )
    
  }
} 
#############
library(tidyverse)
library(readr)

# 1. Get all results RDS files
rds_files <- list.files("outputs", pattern = "^results_.*\\.rds$", recursive = TRUE, full.names = TRUE)

# 2. Loop through and generate + save plots
for (rds_path in rds_files) {
  # Parse metadata from filename
  filename <- basename(rds_path)
  matches <- str_match(filename, "^results_([^_]+)_tef_([^_]+)_buf_([0-9]+)\\.rds$")
  species <- matches[2]
  tef     <- matches[3]
  buf     <- matches[4]
  
  if(species %in% c("palsp")){
    next
  }else{
  # Build file paths
  mix_path    <- file.path("data", "mix_files", sprintf("edge_3_buf_%s_%s_mix.csv", buf, species))
  source_path <- file.path("data", "source_files", sprintf("%s_mix.csv", species))
  
  # Skip if mix or source files are missing
  if (!file.exists(mix_path)) {
    message("Missing mix file: ", mix_path)
    next
  }
  if (!file.exists(source_path)) {
    message("Missing source file: ", source_path)
    next
  }
  
  # Load model and input files
  model  <- readRDS(rds_path)
  mix    <- read_csv(mix_path, show_col_types = FALSE)
  source <- read_csv(source_path, show_col_types = FALSE)
  # Reconstruct MixSIAR-style objects
  mix_obj <- MixSIAR::load_mix_data(
    filename = mix_path, 
    iso_names = c("d13C", "d34S"), 
    factors = NULL,
    fac_random = NULL,
    fac_nested = NULL,
    cont_effects = "edge_l.mangrove")
  
  source_obj <- MixSIAR::load_source_data( filename = source_path, 
                                           conc_dep = TRUE, 
                                           data_type = "raw", 
                                           mix = mix_obj)
  
  # Generate the plot
  plt <- plot_continuous_var_custom(jags.1 = model,
                                    mix    = mix_obj,
                                    source = source_obj,
                                    output_options = list())
  
  # Save to PNG in same directory as the RDS file
  plot_dir <- file.path("results", species, tef, buf)
  dir.create(plot_dir)
  out_file <- file.path(plot_dir, sprintf("cont_plot_%s_tef_%s_buf_%s.png", species, tef, buf))
  ggsave(out_file, plt, width = 8, height = 4, units = 'in')
  
  message("Saved: ", out_file)
  }
}

