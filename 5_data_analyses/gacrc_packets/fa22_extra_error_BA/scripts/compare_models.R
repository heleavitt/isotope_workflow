library(tidyverse)
library(MixSIAR)

setwd("C:/Users/hl51981/OneDrive - University of Georgia/Leavitt_Herbert/PFFW/Pt Fourchon Food Webs R Directory/isotope_workflow/5_data_analyses/gacrc_packets/fa22_extra_error_BA")

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
    type     = str_extract(filename, "^(null_results|results|ind)"),
    species_code = str_extract(filename, "(?<=^(null_results|results|ind_results)_)[^_]+"),
    tef_step     = str_extract(filename, "(?<=tef_)[^_]+"),
    buf_val      = str_extract(filename, "(?<=buf_)[0-9]+")
  )




for (sp in unique(rds_tbl$species_code)) {
  for (tv in unique(rds_tbl$tef_step)) {

    # 3) pivot wider so each row has one null‐ and one real‐results path
  rds_subset <- rds_tbl %>%
      filter(species_code == sp, tef_step == tv, type != "ind") %>% 
      mutate(buf_val = as.numeric(buf_val))

  comp_scores <- data.frame() 
   
  for (buf in unique(rds_subset$buf_val)) {
  
     rds_buf_set <- rds_subset %>%
       filter(buf_val == buf, ) 
          
      paired_tbl <- rds_buf_set %>%
        select(-filename) %>%
        pivot_wider(
          names_from  = type,
          values_from = filepath
        )
      
      paired_tbl2 <- paired_tbl %>% 
        # rename the columns that currently hold file‐paths
        rename(
          null_path    = null_results,
          results_path = results)
          
          
    final_tbl <- paired_tbl2 %>%
      mutate(null_model = map(null_path,readRDS),
             stan_model = map(results_path, readRDS))
             
    model_comp <- list()
    model_comp[[1]] <- readRDS(paired_tbl2$null_path)
    model_comp[[2]] <- readRDS(paired_tbl2$results_path)
    

    comp<-compare_models(model_comp, loo = TRUE)
    
    comp_buf<-data.frame(buffer = buf, weight =  comp$weight[comp$Model == "model.2"])
    
   
     comp_scores <- rbind(comp_scores, comp_buf)
  }
     
    out_dir <- file.path("results", sp, tv)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    out_png <- file.path(out_dir, paste0("Weight_plot_tef_", tv, "_", sp, ".png"))
    png(out_png, width = 8, height = 4, units = 'in', res = 300)
    
  # build the plot object
    p <- ggplot(comp_scores, aes(x = buffer, y = weight)) +
      geom_hline(yintercept = 0.5, linetype = "dashed") +
      geom_line() +
      geom_point(size = 2) +
      scale_y_continuous(limits = c(0,1))+
      scale_x_continuous(breaks = c(20, 70, 100, 150, 200, 300, 400, 500, 600)) +
      labs(
     #   title   = paste("Signed ΔDIC vs Buffer (TEF", tv, ")"),
        x       = "Buffer radius (m)",
        y       = "Akaike weight of % mangrove model",
        #caption = "Negative = null better"
      ) +
      theme_classic(base_size = 14)
    
    # **explicitly print** it
    print(p)
    
    # close device
    dev.off()
    
    message("Wrote ", out_png, " (", nrow(comp_scores), " rows)")
  } 
  }


    
  
