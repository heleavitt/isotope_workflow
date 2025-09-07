# hypervolume_multispecies_by_habitat_quartile.R
# ----------------------------------
# User-defined settings
# Specify multiple species to combine in each hypervolume
species_vec <- c("pensets","calsap","palsp","minlon","ctebol")  # example species codes
buffers_vec <- c(20, 30, 50, 70, 100, 150)            # buffer radii (m)

# Base directories
base_dir   <- "outputs"                             # root folder for mix and posterior data
default_output <- file.path(base_dir, "hv_by_mangrove_quartile")
dir.create(default_output, recursive = TRUE, showWarnings = FALSE)

# Load libraries
library(tidyverse)
library(stringr)
library(hypervolume)

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
    png(out_png2,width = 500, height = 600)
   print( ggplot(df_both_plot, aes(x = src, y = p, fill = src, alpha = 0.85)) +
      geom_boxplot() +
      scale_x_discrete(limits = c(
        "POM",
        "Spartina",
        "Mangrove",
        "Benthic algae"
      )) +
      scale_y_continuous(limits = c(0, 0.9)) +
      theme_classic(base_size = 30) +
      scale_fill_manual(values = c(
        "POM" = "#D81B60",
        "Spartina" = "#1E88E5",
        "Mangrove" = "#DC9203",
        "Benthic algae" = "#004D40"
      )) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))+
      labs(x = "Resource", y = "Proportional contribution", fill = "Resource")
   )
    dev.off()
    
    # Pivot and drop zero-variance dimensions
    wide <- df_both %>% pivot_wider(id_cols=c("ind", "habitat", "species_code"), names_from=src, values_from=p) %>% 
      left_join(df_both[,c("ind", "species_code", "TL")], by = c("ind", "species_code")) %>% unique()
    numeric_df <- wide %>% ungroup() %>% dplyr::select(-habitat, -species_code, -ind)
    vars <- sapply(numeric_df, sd, na.rm=TRUE)
    keep_cols <- names(vars)[vars>0]
    dat_keep <- wide %>% dplyr::select(habitat, species_code, all_of(keep_cols))
    # Scale combined and split
    scaled <- dat_keep %>% mutate(across(all_of(keep_cols), ~ scale(.x)[,1]))
    lo_mat <- scaled %>% filter(habitat=='low') %>% select(all_of(keep_cols), species_code)
    hi_mat <- scaled %>% filter(habitat=='high') %>% select(all_of(keep_cols), species_code)
    # Count points per species per group
    count_lo <- lo_mat %>% count(species_code, name = "n_lo")
    count_hi <- hi_mat %>% count(species_code, name = "n_hi")
    
    # Join and compute min n
    min_n_per_species <- inner_join(count_lo, count_hi, by = "species_code") %>%
      mutate(n_min = pmin(n_lo, n_hi)) %>%
      select(species_code, n_min)
    
    abs_min <- min(min_n_per_species$n_min)
    
    # Slice both datasets to that min n per species
    lo_balanced <- lo_mat %>%
      inner_join(min_n_per_species, by = "species_code") %>%
      group_by(species_code) %>%
      slice_sample(n = abs_min) %>%
      ungroup()
    
    hi_balanced <- hi_mat %>%
      inner_join(min_n_per_species, by = "species_code") %>%
      group_by(species_code) %>%
      slice_sample(n = abs_min) %>%
      ungroup()
  if (species %in% c("pensets", "palsp")){
      sp_lo_mat <- lo_balanced %>% filter(lo_balanced$species_code == species) %>% dplyr::select(-species_code, -TL,- n_min,-`benthic algae`)
      sp_hi_mat <- hi_balanced %>% filter(hi_balanced$species_code == species) %>% dplyr::select(-species_code, -TL,- n_min,-`benthic algae`)
  }
    else{
      sp_lo_mat <- lo_balanced %>% filter(lo_balanced$species_code == species) %>% dplyr::select(-species_code, -TL,- n_min, -`live mangrove leaves`)
      sp_hi_mat <- hi_balanced %>% filter(hi_balanced$species_code == species) %>% dplyr::select(-species_code, -TL,- n_min, -`live mangrove leaves`)
    }
    
      
      hv_low <- hypervolume_gaussian(as.data.frame(sp_lo_mat), name= "low", samples.per.point=1000, sd.count=3,
                                     quantile.requested=0.95, quantile.requested.type='probability', chunk.size=1000)
      hv_high <- hypervolume_gaussian(as.data.frame(sp_hi_mat), name= "high", samples.per.point=1000, sd.count=3,
                                      quantile.requested=0.95, quantile.requested.type='probability', chunk.size=1000)
      low_vol<-get_volume(hv_low)
      high_vol<-get_volume(hv_high)
      
      hv_set<- hypervolume_join(hv_low, hv_high)
     
       png(out_png, width = 800, height = 600)
      
        hv_names <- if (species == "calsap") {
         c('POM', "Benthic\nalgae", "Spartina")
       } else {
         c('POM', "Mangrove", "Spartina")
       }
       
       plot(hv_set, show.3d = FALSE,
            show.axes = TRUE, show.frame = TRUE,
            show.random = TRUE, show.density = TRUE, show.data = TRUE,
            names = hv_names,
            show.legend = FALSE, limits = c(-5, 5),
            show.contour = TRUE, contour.lwd = 3.5,
            contour.type = "alphahull", contour.alphahull.alpha = 6,
            contour.ball.radius.factor = 1, contour.kde.level = 3,
            contour.raster.resolution = 100, show.centroid = TRUE,
            cex.centroid = 2, colors = c("deepskyblue4", "orange3"),
            cex.random = 3, cex.data = 0, cex.axis = 1.3, cex.names = 2.8,
            cex.legend = 3, num.points.max.data = 1000, num.points.max.random = 2000,
            reshuffle = TRUE, verbose = FALSE,
            main = sprintf("%s TEF %s — buffer %d", species, tef_step, buf))
        
        legend(x = 0.1, y = 0.2,
               legend = c(
                 bquote(bold("Marsh") * " (volume: " * .(round(low_vol, 1)) * ")"),
                 bquote(bold("Mangrove") * " (volume: " * .(round(high_vol, 1)) * ")")
               ),
               fill = c("deepskyblue4", "orange3"),
               cex = 1.8)
      dev.off()
      
      hv_low_list[[species]] <- hv_low
      hv_high_list[[species]]<- hv_high
      
    }
    }

  # hv_low_list is a named list of Hypervolume objects
  hv_low_joined <- Reduce(
    function(acc, next_hv) hypervolume_join(acc, next_hv, check.memory = FALSE),
    hv_low_list
  )
  hv_high_joined <- Reduce(
    function(acc, next_hv) hypervolume_join(acc, next_hv, check.memory = FALSE),
    hv_high_list
  )
  
  
  hv_all_joined <- hypervolume_join(hv_high_joined, hv_low_joined)
      png(out_png, width = 800, height = 600)
      plot(hv_high_joined, show.3d=F,plot.3d.axes.id=NULL,
           show.axes=TRUE, show.frame=TRUE,
           show.random=TRUE, show.density=TRUE,show.data=TRUE,
           show.legend=TRUE,limits = c(-5,5),
           show.contour=F, contour.lwd=5, 
           contour.type="alphahull",
           contour.alphahull.alpha=3,
           contour.ball.radius.factor=1,
           contour.kde.level=3,
           contour.raster.resolution=100,
           show.centroid=TRUE, cex.centroid=2, colors=c("deepskyblue4","orange3"),
           cex.random= 1,cex.data= 0,cex.axis= 2,cex.names=2.8, cex.legend= 3,
           num.points.max.data = 10000, num.points.max.random = 2000, reshuffle=T,
           verbose=FALSE,
           main = sprintf("%s TEF %s — buffer %d", species_code, tef_step, buf))
      legend("topright", c("low","high"), fill = c("deepskyblue4","orange3"))
      dev.off()
      
      message("Done: ", species_code, " TEF=", tef_step, " buffer=", buf)
    
  


