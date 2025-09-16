# Specify multiple species to combine in each hypervolume
species_vec <- c("pensets")  # example species codes
buffers_vec <- c(150)            # buffer radii (m)
setwd("C:/Users/hl51981/OneDrive - University of Georgia/Leavitt_Herbert/PFFW/Pt Fourchon Food Webs R Directory/isotope_workflow/5_data_analyses/gacrc_packets/fa22_extra_error_BA")
# Base directories
base_dir   <- "outputs"                             # root folder for mix and posterior data
default_output <- file.path(base_dir, "hv_by_mangrove_quartile")
dir.create(default_output, recursive = TRUE, showWarnings = FALSE)

# Load libraries
library(tidyverse)
library(stringr)
library(hypervolume)
# Helper: make numeric, scaled matrix(s) per habitat and pooled ---------
.prepare_mats <- function(df_both, species) {
  # recode to consistent axis names
  df_w <- df_both %>%
    mutate(src = recode(src,
                        "CRAVIR" = "POM",
                        "live spartina leaves and stems" = "Spartina",
                        "live mangrove leaves"         = "Mangrove",
                        "benthic algae"                = "Benthic algae")) %>%
    select(habitat, ind, species_code, src, p, TL) %>%
    filter(species_code == species)

  # wide -> drop zero-variance -> scale
  wide <- df_w %>%
    pivot_wider(id_cols = c(ind, habitat, species_code),
                names_from = src, values_from = p, values_fill = 0) %>%
    left_join(df_w[, c("ind","species_code","TL")] %>% distinct(), by = c("ind","species_code")) %>%
    distinct()

  # remove TL and non-numeric cols later; figure kept axes by variance
  numdf <- wide %>% select(-habitat, -species_code, -ind)
  vars  <- sapply(numdf, stats::sd, na.rm = TRUE)
  keep  <- names(vars)[vars > 0 & !names(vars) %in% c("TL","n_min")]

  # species-specific axis removal to mirror your rules
  if (species %in% c("pensets","palsp")) {
    keep <- setdiff(keep, "Benthic algae")
  } else { # calsap
    keep <- setdiff(keep, "Mangrove")
  }
  stopifnot(length(keep) >= 1)

  scaled <- wide %>%
    mutate(across(all_of(keep), ~ as.numeric(scale(.x)[,1])))

  list(
    low   = scaled %>% filter(habitat == "low")  %>% select(all_of(keep)),
    high  = scaled %>% filter(habitat == "high") %>% select(all_of(keep)),
    pooled= scaled %>% select(habitat, all_of(keep))
  )
}

# Helper: safe HV builder ------------------------------------------------
.make_hv <- function(X, hv_name, spp = 1000L) {
  # guard small matrices
  if (nrow(X) < 3) return(NULL)
  hypervolume_gaussian(
    data = as.data.frame(X),
    name = hv_name,
    samples.per.point = spp,
    sd.count = 3,
    quantile.requested = 0.95,
    quantile.requested.type = "probability",
    chunk.size = 1000,
    verbose = FALSE
  )
}

# Core bootstraps --------------------------------------------------------
bootstrap_hv <- function(mats, n_iter = 100L, frac_size = 2/3, frac_rand = 1/2, seed = 123) {
  set.seed(seed)

  # 1) Niche size per habitat (2/3 per species×habitat)
  size_boot <- map_dfr(seq_len(n_iter), function(i) {
    # sample per habitat independently
    take_n <- function(df, frac) {
      n <- nrow(df); k <- max(3, floor(frac * n))
      if (n < 3) return(NULL)
      df[sample.int(n, k, replace = n < k), , drop = FALSE]
    }
    low_s  <- take_n(mats$low,  frac_size)
    high_s <- take_n(mats$high, frac_size)
    if (is.null(low_s) || is.null(high_s)) return(tibble())

    hvL <- .make_hv(low_s,  "low")
    hvH <- .make_hv(high_s, "high")
    tibble(iter = i,
           habitat = c("low","high"),
           volume  = c(ifelse(is.null(hvL), NA_real_, get_volume(hvL)),
                       ifelse(is.null(hvH), NA_real_, get_volume(hvH))))
  })

  # 2) “Random individual” niche size (50% pooled)
  rand_boot <- map_dfr(seq_len(n_iter), function(i) {
    # pooled sample irrespective of habitat
    pool <- mats$pooled %>% select(-habitat)
    n <- nrow(pool); k <- max(3, floor(frac_rand * n))
    if (n < 3) return(tibble())
    idx <- sample.int(n, k, replace = n < k)
    hvR <- .make_hv(pool[idx, , drop = FALSE], "random")
    tibble(iter = i, habitat = "random",
           volume = ifelse(is.null(hvR), NA_real_, get_volume(hvR)))
  })

  # 3) Overlap bootstraps (low vs high; 2/3 each)
  overlap_boot <- map_dfr(seq_len(n_iter), function(i) {
    take_n <- function(df, frac) {
      n <- nrow(df); k <- max(3, floor(frac * n))
      if (n < 3) return(NULL)
      df[sample.int(n, k, replace = n < k), , drop = FALSE]
    }
    low_s  <- take_n(mats$low,  frac_size)
    high_s <- take_n(mats$high, frac_size)
    if (is.null(low_s) || is.null(high_s)) return(tibble())

    hvL <- .make_hv(low_s,  "low")
    hvH <- .make_hv(high_s, "high")
    if (is.null(hvL) || is.null(hvH)) return(tibble())

    hs   <- hypervolume_set(hvL, hvH, check.memory = FALSE)
    ost  <- hypervolume_overlap_statistics(hs)
    tibble(iter = i, sorensen = as.numeric(ost["sorensen"]))
  })

  list(size = bind_rows(size_boot, rand_boot), overlap = overlap_boot)
}

# Summaries + tests ------------------------------------------------------
summarize_boot <- function(size_df, overlap_df, species, tef_step, buf) {
  # 95% CIs
  ci <- function(x) quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  size_ci <- size_df %>%
    group_by(habitat) %>%
    summarise(n = sum(!is.na(volume)),
              lo = ci(volume)[1], md = ci(volume)[2], hi = ci(volume)[3],
              .groups = "drop") %>%
    mutate(species = species, tef = tef_step, buf = buf)

  olap_ci <- overlap_df %>%
    summarise(n = sum(!is.na(sorensen)),
              lo = ci(sorensen)[1], md = ci(sorensen)[2], hi = ci(sorensen)[3],
              .groups = "drop") %>%
    mutate(species = species, tef = tef_step, buf = buf)

  # Kruskal–Wallis (sizes among habitats; keep only low vs high)
  size_kw <- try({
    dd <- size_df %>% filter(habitat %in% c("low","high")) %>% drop_na(volume)
    if (n_distinct(dd$habitat) == 2) kruskal.test(volume ~ habitat, data = dd) else NULL
  }, silent = TRUE)

  size_dunn <- try({
    dd <- size_df %>% filter(habitat %in% c("low","high")) %>% drop_na(volume)
    if (n_distinct(dd$habitat) == 2) {
      dunn.test(x = dd$volume, g = droplevels(dd$habitat), method = "bonferroni", kw = FALSE)
    } else NULL
  }, silent = TRUE)

  # Kruskal–Wallis for overlaps vs a null? (optional)
  list(size_ci = size_ci, overlap_ci = olap_ci,
       size_kw = size_kw, size_dunn = size_dunn)
}

# Loop over TEF steps and buffer sizes
for (buf in buffers_vec) {
    # Collect data across all species
    df_list <- list()
    for (species in species_vec) {
      # Construct paths
      tef_step<-ifelse(species %in% c("palsp", "pensets"), "1.5", "1.0")
      sdir <- file.path(base_dir, species, paste0("tef_step_", tef_step), as.character(buf))
      mix_csv <- file.path(sdir, sprintf("mixing_model_df_%s_tef_%s_buf_%s.csv", species, tef_step, buf))
      
      out_png <-file.path("results","hv_by_mangrove_quartile", sprintf("comp_hv_%s_tef_%s_buf_%s.png", species, tef_step, buf))
      out_png2 <-file.path("results","hv_by_mangrove_quartile", sprintf("comp_box_%s_tef_%s_buf_%s.png", species, tef_step, buf))
      
      if (!file.exists(mix_csv)) next
      # Read data
      df_mm <- read_csv(mix_csv) %>% mutate(covariate = edge_l.mangrove,
                                            src = Source) %>% filter(species_code == species)
      

        qs <- quantile(df_mm$covariate, c(0.5), na.rm=TRUE)
        df_low  <- df_mm %>% filter(covariate < qs[1]) %>% mutate(habitat='low')
        df_high <- df_mm %>% filter(covariate > qs[1]) %>% mutate(habitat='high')
        if (nrow(df_low)==0 || nrow(df_high)==0) next
        # Summarize per individual, source across species
        df_both <- bind_rows(df_low, df_high) %>%
        group_by(habitat, ind, src, species_code) %>%
        summarise(p=mean(Contr,na.rm=TRUE),
                    TL = mean(TL), .groups='drop')
        
        
        df_both_plot <- df_both %>%
        mutate(src = recode(src,
                            "CRAVIR" = "POM",
                            "live spartina leaves and stems" = "Spartina",
                            "live mangrove leaves" = "Mangrove",
                            "benthic algae" = "Benthic algae"  # just in case
        ))

        # ----------------- RUN for each species within your existing loop ----------------
        # After you computed df_both for this species/buf:
        mats <- .prepare_mats(df_both, species = species)

        boots <- bootstrap_hv(
        mats,
        n_iter   = 100L,
        frac_size= 2/3,
        frac_rand= 1/2,
        seed     = 42
        )

        out <- summarize_boot(boots$size, boots$overlap, species, tef_step, buf)

        # Save tabular outputs (append per species)
        dir.create(file.path("results","hv_bootstrap"), showWarnings = FALSE, recursive = TRUE)
        readr::write_csv(out$size_ci,
                        file.path("results","hv_bootstrap",
                                sprintf("size_CI_%s_tef_%s_buf_%s.csv", species, tef_step, buf)),
                        append = FALSE)
        readr::write_csv(out$overlap_ci,
                        file.path("results","hv_bootstrap",
                                sprintf("overlap_CI_%s_tef_%s_buf_%s.csv", species, tef_step, buf)),
                        append = FALSE)

        # Optional: print quick sanity summary
        print(out$size_ci)
        print(out$overlap_ci)
        if (!inherits(out$size_kw, "try-error") && !is.null(out$size_kw)) print(out$size_kw)
        if (!inherits(out$size_dunn, "try-error") && !is.null(out$size_dunn)) print(out$size_dunn)
    }

  }

x <- na.omit(boots$overlap$sorensen)
ci <- quantile(x, c(0.025, 0.5, 0.975))
hist(x, breaks = seq(0, 1, by = 0.05),
     xlab = "Sørensen overlap", main = "Bootstrap overlap")
abline(v = ci, lty = c(2,1,2))

