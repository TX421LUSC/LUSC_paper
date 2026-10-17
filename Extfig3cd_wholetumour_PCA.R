#Script for Tx421 LUSC paper
#Ext fig 3c,d.
# PCA for pathologic features of tumour regions graded by entire tumour in LUSC

library(dplyr)
library(ggplot2)
library(reshape2)
library(ggbiplot)
library(tidyverse)

idlist <- read.csv("../df-alvin/id_conversion.csv") 
vsdfile <- readRDS("../data/vsd_median_alltumour_df.RDS") 
#this file vsd_median_alltumour_df.RDS was obtained using the exact same method which 
#resulted in the 20221108_vsd_median_df_luad.RDS df, already deposited in Zenodo. In fact,
#20221108_vsd_median_df_luad.RDS is a subset of this file, and this file contains all tumours LUAD+LUSC. 


# Load entire tumour data
squamous_entire_data <- read_csv("../df-alvin/500_entire_tumour_with_grading_with_wGIITMB_info.csv")

tumour_keratinizationdata <- readRDS("../df-alvin/tumour_regionalkeratinization.RDS")
tumour_keratinization <- tumour_keratinizationdata %>% 
  mutate(subtype_by_k = case_when(keratinization == 'have_keratin' & proportion_regions == 1 ~ "Fully_keratinized",
                                  keratinization == 'no_keratin' & proportion_regions == 1 ~ "Fully_nonkeratinized",
                                  keratinization == 'have_keratin' & proportion_regions < 1 ~ "Mixed_keratinized")) %>% 
  mutate(subtype_by_k_detail = case_when(keratinization == 'have_keratin' & proportion_regions == 1 ~ "Fully_keratinized",
                                         keratinization == 'no_keratin' & proportion_regions == 1 ~ "Fully_nonkeratinized",
                                         keratinization == 'have_keratin' & proportion_regions <0.5 ~ "Low_keratinized",
                                         keratinization == 'have_keratin' & proportion_regions >= 0.5 & proportion_regions != 1 ~ "High_keratinized")) %>% 
  filter(!is.na(subtype_by_k)) %>% 
  dplyr::select(-keratinization)

squamous_entire_data <- squamous_entire_data %>% left_join(tumour_keratinization, by = c("cruk_tumourid" = "tumour_id_muttable_cruk"))

tumour_id_region <- squamous_entire_data %>% dplyr::select(cruk_id_all, Histology, differentiation, keratinization, tumour_budding, single_cell_invasion, atypia, subtype_by_k, subtype_by_k_detail, wGII, TMB, FLOH)
# ##############################################################################################

  #Perform PCAs in regions separately--------------------------------------
  
  fullyk_regions <- tumour_id_region$cruk_id_all[tumour_id_region$subtype_by_k_detail == "Fully_keratinized" & !is.na(tumour_id_region$subtype_by_k_detail)]
  lowk_regions <- tumour_id_region$cruk_id_all[tumour_id_region$subtype_by_k_detail == "Low_keratinized" & !is.na(tumour_id_region$subtype_by_k_detail)]
  fullynonk_regions<- tumour_id_region$cruk_id_all[tumour_id_region$subtype_by_k_detail == "Fully_nonkeratinized" & !is.na(tumour_id_region$subtype_by_k_detail)]
  highk_regions <- tumour_id_region$cruk_id_all[tumour_id_region$subtype_by_k_detail == "High_keratinized" & !is.na(tumour_id_region$subtype_by_k_detail)]
  
  all_regions <- c(fullyk_regions, 
                   lowk_regions, 
                   fullynonk_regions,
                   highk_regions)
  
  region_list <- all_regions
  vst_mat <- vsdfile
  
  #Keep only regions of interest and calculate PCs
  to_pca_vst <- vst_mat[, colnames(vst_mat) %in% region_list]
  to_pca_vst <- as.data.frame(t(to_pca_vst))
  to_pca_vst$casenumber <- rownames(to_pca_vst)
  to_pca_vst1 <- to_pca_vst %>% mutate(grouping = case_when(casenumber %in% fullyk_regions ~ "FullyK", 
                                                            casenumber %in% lowk_regions ~ "LowK",
                                                            casenumber %in% fullynonk_regions ~ "FullynonK",
                                                            casenumber %in% highk_regions ~ "HighK"))
  
  pcs_region <- prcomp(dplyr::select(to_pca_vst1, -c(grouping, casenumber)))
  
  
  color_list <-  c("FullyK" = "#1C75BC",
                   "LowK" = "#7F3F98",
                   "FullynonK"= "#FB9A99",
                   "HighK" = "#8B5E3C")
  
  pcaplot1 <- ggbiplot(pcs_region, obs.scale = 1, var.scale = 1, choices = c(1,2),
                       groups = to_pca_vst1$grouping, ellipse = TRUE, circle = TRUE,
                       var.axes = FALSE) +
    scale_color_discrete(name = '') +
    theme(legend.direction = 'horizontal', legend.position = 'top') +
    scale_color_manual(name = "subtype_by_k_detail",
                       values = color_list) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.text.y = element_text(size = 13, colour = "black"), 
          axis.text.x = element_text(size = 13, colour = "black"),
          legend.title = element_text(size = 14), 
          legend.text = element_text(size = 13),
          axis.line = element_line(colour = "black"),
          panel.border = element_blank(),
          panel.background = element_blank())
  print(pcaplot1)
  #Ext fig 3d
  
  #make sure the order of casenumbers in pcs_region is the same as to_pca_vst1
  identical(rownames(pcs_region$x), to_pca_vst1$casenumber)
  #check whether the gene in first PC is irregularly high 
  head(sort(abs(pcs_region$rotation[, "PC1"]), decreasing = TRUE))

  
  
#Ext fig 3c Entire tumour PCA
  # Define parameters
  param <- "keratinization"
  param_poor_group <- "no_keratin"
  param_well_group <- "have_keratin"
  
  # Assign parameter to dataframe
  tumour_id_region$param <- tumour_id_region[[param]]
  
  # Get regions for each group
  poor_regions <- tumour_id_region$cruk_id_all[tumour_id_region$param == param_poor_group & !is.na(tumour_id_region$param)]
  well_regions <- tumour_id_region$cruk_id_all[tumour_id_region$param == param_well_group & !is.na(tumour_id_region$param)]
  
  # Combine regions for analysis
  all_regions <- c(well_regions, poor_regions)
  region_list <- all_regions
  vst_mat <- vsdfile
  
  # Keep only regions of interest and prepare data for PCA
  to_pca_vst <- vst_mat[, colnames(vst_mat) %in% region_list]
  to_pca_vst <- as.data.frame(t(to_pca_vst))
  to_pca_vst$casenumber <- rownames(to_pca_vst)
  
  # Annotate each case with its group
  to_pca_vst1 <- to_pca_vst %>% mutate(
    grouping = case_when(
      casenumber %in% well_regions ~ param_well_group,
      casenumber %in% poor_regions ~ param_poor_group,
      TRUE ~ "allothers"
    )
  )
  
  # Perform PCA
  pcs_region <- prcomp(dplyr::select(to_pca_vst1, -c(grouping, casenumber)))
  
  # Define color schemes
  color_scheme <- c( "#A6CEE3", "#FB9A99")
  # Plot with original color scheme
  pcaplot1 <- ggbiplot(pcs_region, obs.scale = 1, var.scale = 1, choices = c(1,2),
                       groups = to_pca_vst1$grouping, ellipse = FALSE, circle = TRUE,
                       var.axes = FALSE) +
    scale_color_discrete(name = '') +
    theme(legend.direction = 'horizontal', legend.position = 'top') +
    scale_color_manual(name = param, values = color_scheme) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 13, colour = "black"),
      axis.text.x = element_text(size = 13, colour = "black"),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13),
      axis.line = element_line(colour = "black"),
      panel.border = element_blank(),
      panel.background = element_blank()
    )
  
  # Ensure the order of case numbers matches
  print(identical(rownames(pcs_region$x), to_pca_vst1$casenumber))  # Should be TRUE
  
  # Check which genes have the highest loadings on PC1
  print(head(sort(abs(pcs_region$rotation[, "PC1"]), decreasing = TRUE)))
  
  # Return the PCA plot
  p2 <- pcaplot1

  print(p2)  

