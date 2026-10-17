##Script for Tx421 LUSC paper
##Ext fig 3h and 3i

#################################
### libraries and directories ###
#################################

library(tidyverse)
library(data.table)
library(fst)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(ggfittext)
library(cowplot)
library(ggpubr)
library(ggsci)
library(gtable)


########################
###  load input data ###
########################

# main data
input_plot_df_all <- readRDS("../data/20221109_LUAD_input_plot_df_all.248.RDS")
LUAD_region_df_full.ori <- readRDS("../data/20221109_LUAD_input_plot_df_all_LUAD_region_df_full.RDS")
entire_tumour_LUSC <- read_csv("../df-alvin/500_entire_tumour_with_grading_with_wGIITMB_info.csv")

patient_df <- readRDS("../data/20221109_TRACERx421_all_patient_df.rds") 

# clonal diversity and recent subcl expansion score
max_cp_df <- fread("../data/20221102_TRACERx421_tumourDivCorrectedtree.csv")
entire_tumour_LUSC -> entire_tumour_LUSC_0

#LVI data
LVI_data <- read_csv("../df-alvin/LVI_data")



########################
### helper functions ###
########################

## functions to align legends with ggplot ###
 ## (Thanks to Michelle Dietzen and Bobby Bentham) ##

# library(gtable)

#Extract the legend from a ggplot object with this function:
g_legend <- function(a.gplot){
  tmp    <- ggplot_gtable(ggplot_build(a.gplot))
  leg    <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  if(length(leg) == 0){
    return(NULL)
  }
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#Then use the following function on a list of legends to return a single grob that will contain all the legends aligned:
align.legends.fun  <- function(legend.list){
  aligned.leg <- legend.list[[1]]
  for(i in 2:length(legend.list)){
    leg1        <- legend.list[[i]]$grobs[[1]]
    leg_a       <- gtable_add_rows(aligned.leg, pos = nrow(aligned.leg) - 1, heights = sum(leg1$heights))
    leg_final   <- gtable_add_grob(leg_a, leg1, t = nrow(leg_a) - 1, l = 3)
    aligned.leg <- leg_final
  }
  return(aligned.leg)
}



## set gab_size
gab_size    <- 0.2
text_size   <- 18
legend_text_size  <- 13
legend_title_size <- 15

## finalise input data
subtype_list_full <- c("have_keratin", "no_keratin")
Sub_list <- c("have_keratin", "no_keratin")

input_plot_df_all_ready <-  entire_tumour_LUSC %>%
  mutate(cruk_id = substr(cruk_tumourid, 1, 8)) %>% 
  mutate(tumour_id_muttable_cruk = cruk_tumourid) %>% 
  filter(keratinization %in% subtype_list_full) %>%
  left_join(patient_df[,c("cruk_id", "smoking_status_merged", "pathologyTNM", "sex", "age")], by = "cruk_id") %>% 
  left_join(max_cp_df[,c("tumour_id", "maxCCF_terminalonly_tum_meanccf", "min_simpson_tum_meanccf")], by = c("tumour_id_muttable_cruk" = "tumour_id")) %>%
  mutate(keratinization = factor(keratinization, levels = c("have_keratin", "no_keratin"))) %>%
  mutate(clonal_TMB_log = log10(clonal_tmb),
         subclonal_TMB_log = log10(subclonal_tmb),
#         clonal_illusion_count_log = log10(clonal_illusion_count + 1),
#         ki67 = KI67_score/100
  ) %>%
  mutate(smoking = factor(smoking_status_merged , levels = c("Never Smoked", "Ex-Smoker", "Smoker"))) %>%
  mutate(num_clonal_gds = factor(num_clonal_gds, levels = c("0", "1", "2"))) %>% 
  mutate(num_subclonal_gds = factor(num_subclonal_gds, levels = c("0", "1", "2"))) %>%   #: Set as factors to avoid error of continuous variable
  mutate(TumourID = tumour_id_muttable_cruk) %>% 
  mutate(extensive_necrosis1 = case_when(extensive_necrosis == 1 ~ "have_necrosis", 
                                         extensive_necrosis == 0 ~ "no_necrosis"))
  
input_plot_df_all_ready <- left_join(input_plot_df_all_ready, LVI_data, by = c("cruk_tumourid" = "tumour_id_per_patient"))

gd_table <- read_tsv("../data/20221109_TRACERx421_gd_table_per_tumour.tsv")
gd_table_full <- gd_table %>% left_join(input_plot_df_all[, c("tumour_id_muttable_cruk", "histology_3")], by = c("tumour_id" = "tumour_id_muttable_cruk")) %>% 
  left_join(input_plot_df_all_ready[, c("tumour_id_muttable_cruk", "Histology", "keratinization")],by = c("tumour_id" = "tumour_id_muttable_cruk")) %>% 
  mutate(Histology = ifelse(!is.na(histology_3), histology_3, Histology))

#Ext fig 3h
x_var_list <- c("Histology")
y_var_list <- c("num_clonal_gds", "num_subclonal_gds")
in_data <- gd_table_full
# in_data$Relapse_cat_new <- factor(in_data$Relapse_cat_new, levels = c("Unknown Site", "Intra & Extra", "Extrathoracic", "Intrathoracic", "No rec"))
out_f <- paste0("../results/genomedoubling_LUADvsLUSC.pdf")
pdf(out_f, width = 8, height = 6)

for(x_var in x_var_list){
  for(y_var in y_var_list){
    
    df <- in_data
    
    df$x_var <- df[[x_var]]
    df$y_var <- as.factor(df[[y_var]])
    df <- df %>%
      filter(!is.na(x_var),!is.na(y_var))
    
    df_count <- df %>%
      group_by(y_var, x_var) %>%
      summarise(count = n())
    df_count$x_var[df_count$x_var == "Squamous cell carcinoma"] <- "LUSC"
    
    
    g4 <-ggplot(data=df_count, aes(x=x_var, y=count, fill=y_var, label = count)) +
      geom_bar(stat="identity", position = "fill", color = "black") +
      geom_bar_text( position = "fill", reflow = TRUE, outside = TRUE) +
      # scale_fill_manual(values = colorlist) +
      theme_cowplot() + 
      labs(x = x_var, y = "Proportion of patients (%)", fill = y_var)   + 
      theme(axis.title.x = element_text(size=24), 
            axis.title.y = element_text(size=24),
            axis.text.x = element_text(size=20, angle = 45,  hjust=1), 
            axis.text.y = element_text(size=16), 
            legend.position = "right",
            legend.title = element_text(size=22), #legend size
            legend.text = element_text(size=18),
            strip.text = element_text(size=16)) 
    print(g4)
  }
}
dev.off()
# Ext fig 3h


#Ext fig 3i
x_var_list <- c("keratinization")
y_var_list <- c("num_clonal_gds", "num_subclonal_gds")
in_data <- gd_table_full
out_f <- paste0("../results/genomedoubling_kvsnonKLUSC.pdf")
pdf(out_f, width = 8, height = 6)

for(x_var in x_var_list){
  for(y_var in y_var_list){
    
    df <- in_data
    
    df$x_var <- df[[x_var]]
    df$y_var <- as.factor(df[[y_var]])
    df <- df %>%
      filter(!is.na(x_var),!is.na(y_var))
    
    df_count <- df %>%
      group_by(y_var, x_var) %>%
      summarise(count = n())
    
    g4 <-ggplot(data=df_count, aes(x=x_var, y=count, fill=y_var, label = count)) +
      geom_bar(stat="identity", position = "fill", color = "black") +
      geom_bar_text( position = "fill", reflow = TRUE, outside = TRUE) +
      # scale_fill_manual(values = colorlist) +
      theme_cowplot() + 
      labs(x = x_var, y = "Proportion of patients (%)", fill = y_var)   + 
      theme(axis.title.x = element_text(size=24), 
            axis.title.y = element_text(size=24),
            axis.text.x = element_text(size=20, angle = 45,  hjust=1), 
            axis.text.y = element_text(size=16), 
            legend.position = "right",
            legend.title = element_text(size=22), #legend size
            legend.text = element_text(size=18),
            strip.text = element_text(size=16)) 
    print(g4)
  }
}
dev.off()
#Ext fig 3i
