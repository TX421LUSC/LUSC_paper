##Script for Tx421 LUSC paper
##Ext fig 7a,b,c,d

library(tidyverse)
library(ggpubr)
library(cowplot)
library(data.table)


entire_tumour_data <- readRDS("../data/all_tumour_and_regions_df/full_histology_df.20220818.RDS") 
entire_info <- entire_tumour_data[, c("Histology_per_tumour_id_muttable", "patient", "tumour_id_by.label", "wGII", "TMB", "clonal_tmb", "clonal_TMB", "subclonal_tmb", "subclonal_TMB", "FLOH", "num_clonal_gds", "num_subclonal_gds", "purity", "ploidy", "smoking_status", "pack_years")] 
neoag_info <- entire_tumour_data %>% dplyr::select(patient, tumour_id_by.label, contains("neoag"))
squamous_data <- read_csv("../df-alvin/500_entire_tumour_with_grading.csv")
all_regional_full <- readRDS("../data/all_tumour_and_regions_df/all_regional_df_20221109.RDS")
years_smoking_data <- readRDS("../data/20221109_TRACERx421_all_patient_df.rds") 
smoking_cessation_data <- read_csv("../data/all_patient_df_20221109.csv")
smoking_cessation_data <- smoking_cessation_data[, c("tumour_id_muttable_cruk", "years_after_stop")]
clonal_signature <- read_tsv("../data/signature/clonal_signature_tab.txt")
clonal_signature <- clonal_signature %>% select(ltxid, contains("SBS"))
clonal_sbscount <- read_tsv("../data/signature/clonal_signaturecount_tab.txt")
clonal_sbscount <-clonal_sbscount %>% select(ltxid, contains("SBS"))
#evometrics: ITH
wgiidata <- read_csv("../data/wgii_20221110_patient_Evo_metrics_tx421.csv")
wgiitrim <- wgiidata[, c("ltx_tumourid", "frac_abberant_genom_subcl", "perc_subclonal")]
#SBS
sbs_data <- readRDS("../data/20221110_TRACERx421_mutationSignature_weights.rds")
sbs_entiretumour <- sbs_data$signature_weights_perTumour
#subclonal expansion score  Recent_subclonal_expansion_score` = maxCCF_terminalonly_tum_meanccf
#Subclonal_diversity_within_region_(min)` = min_simpson_tum_meanccf,
max_cp_df <- fread("../data/20221102_TRACERx421_tumourDivCorrectedtree.csv")
expansion_score <- max_cp_df[, c("tumour_id", "maxCCF_terminalonly_tum_meanccf", "min_simpson_tum_meanccf")]

max_regional_data <- read_csv("../df-alvin/All_regions_500_regions-SCC_reviewed_With_max_budding_data.csv")
max_regional_data %>% dplyr::select(tumour_id_muttable_cruk, contains("max"), -max_tumour_budding) -> max_regional_data1
max_regional_data1[!duplicated(max_regional_data1), ]  -> max_regional_data2


squamous_only <- squamous_data %>% filter(!is.na(keratinization))

squamous_plot <- left_join(squamous_only, entire_info, by = c("ltx_tumour_id" = "tumour_id_by.label"))

squamous_plot1 <- left_join(squamous_plot, wgiitrim, by = c("ltx_tumour_id" = "ltx_tumourid"))
squamous_plot2 <- left_join(squamous_plot1, sbs_entiretumour, by = c("cruk_tumourid" = "tumour_id"))
squamous_plot3 <- left_join(squamous_plot2, expansion_score, by = c("cruk_tumourid" = "tumour_id"))
squamous_plot4 <- left_join(squamous_plot3, neoag_info, by = c("ltx_tumour_id" = "tumour_id_by.label")) %>% 
                  left_join(years_smoking_data, by = c("cruk_tumourid" = "tumour_id_muttable_cruk")) %>% 
                  left_join(smoking_cessation_data, by = c("cruk_tumourid" = "tumour_id_muttable_cruk")) %>% 
                  left_join(clonal_signature, by = c("ltx_clusterid" = "ltxid")) %>% 
                  left_join(clonal_sbscount, by = c("ltx_clusterid" = "ltxid"))
squamous_plot <- squamous_plot4

squamous_plot$keratinization <- factor( squamous_plot$keratinization , levels=c("have_keratin", "no_keratin"))
squamous_plot$single_cell_invasion <- factor( squamous_plot$single_cell_invasion , levels=c("no_sci", "have_sci"))

##################plot 

boxplot_function <- function(param, pval_label = 0.7){
  color_list <-  c("no_sci" = "#A6CEE3",
                   "greencolour" = "#B2DF8A",
                   "have_sci"= "#FB9A99")
  
  squamous_plot$param <- squamous_plot[[param]]
  stat_n_number <- squamous_plot %>% filter(!is.na(param))
  #compare_means(param ~ single_cell_invasion, data = squamous_plot, method = "wilcox.test")
  g1 <- ggboxplot(stat_n_number, width=0.5,lwd=1, x = "single_cell_invasion", y = param, fill = "single_cell_invasion", add = "jitter") + 
    stat_compare_means(comparisons = list( c("have_sci", "no_sci")), label.y = c(pval_label)) +
    theme_cowplot() +
    scale_fill_manual(values = color_list, na.value = '#bababa', drop = F) +
    #  scale_x_discrete(labels = c("Well-differentiated", "Poorly-differentiated")) +
    theme(  axis.title.y = element_text(size=26),
            axis.text.y = element_text(size=26),
            axis.title.x = element_text(size=26),
            axis.text.x = element_text(size=26),
            legend.position = "none",
            legend.title = element_blank()) 
  print(paste0("n number of no_sci:", nrow(stat_n_number[stat_n_number$single_cell_invasion == "no_sci",]) ))
  print(paste0("n number of have_sci:", nrow(stat_n_number[stat_n_number$single_cell_invasion == "have_sci",]) ))
  print(g1)
}


#plot these new graphs

boxplot_function("pack_years", 110)  #Ext 7a

boxplot_function("cigs_perday", 70) #Ext 7c

boxplot_function("clonalSBS92count", 400) #Ext 7d

boxplot_function("clonalSBS4count", 750)  #Ext 7d

squamous_plot1 <- squamous_plot %>% 
  left_join(years_smoking_data, by = c("cruk_tumourid" = "tumour_id_muttable_cruk")) %>% 
  left_join(smoking_cessation_data, by = c("cruk_tumourid" = "tumour_id_muttable_cruk")) %>% 
  left_join(clonal_signature, by = c("ltx_clusterid" = "ltxid")) %>% 
  left_join(clonal_sbscount, by = c("ltx_clusterid" = "ltxid")) %>% 
  left_join(max_regional_data2, by = c("cruk_tumourid" = "tumour_id_muttable_cruk"))

squamous_plot <- squamous_plot1

squamous_plot$sci_high.low <- factor( squamous_plot$sci_high.low , levels=c("no_sci", "low_sci", "high_sci"))


boxplot_function <- function(param, pval_label = 0.7){
  color_list <-  c("no_sci" = "#A6CEE3",
                   "low_sci" = "#B2DF8A",
                   "high_sci"= "#FB9A99")
  g1 <- ggboxplot(squamous_plot, width=0.5,lwd=1, x = "sci_high.low", y = param, fill = "sci_high.low", add = "jitter") + 
    stat_compare_means(comparisons = list( c("no_sci", "high_sci") , c("no_sci", "low_sci"), c("low_sci", "high_sci")), label.y = c(pval_label+10, pval_label, pval_label)) +
    theme_cowplot() +
    scale_fill_manual(values = color_list, na.value = '#bababa', drop = F) +
    #  scale_x_discrete(labels = c("Well-differentiated", "Poorly-differentiated")) +
    theme(  axis.title.y = element_text(size=26),
            axis.text.y = element_text(size=26),
            axis.title.x = element_text(size=26),
            axis.text.x = element_text(size=26),
            legend.position = "none",
            legend.title = element_blank()) 
  print(paste0("n number of have_keratin:", nrow(squamous_plot[squamous_plot$keratinization == "have_keratin",]) ))
  print(paste0("n number of no_keratin:", nrow(squamous_plot[squamous_plot$keratinization == "no_keratin",]) ))
  print(g1)
}



boxplot_function("pack_years", 110)  #Ext 7c

ggplot(squamous_plot[squamous_plot$max_sci != 0, ], aes(x = pack_years, y = max_sci)) + 
  geom_point(aes(colour = '#67001f'), alpha = 0.7) +
  geom_smooth(method = 'lm', colour = 'black', linetype = 'dashed') +
  stat_cor() +
  #scale_y_continuous( breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20 )) +
  #scale_colour_manual(name = 'Cancer Type', values = c('LUAD' = '#67001f', 'LUSC' = '#053061')) +
  #facet_grid(histology ~ variable, scales = 'free_x') + 
  xlab('Pack Years') + ylab('Regional maximal single cell invasion') +
  theme_bw() + theme(legend.position = 'top', text = element_text(size = 20))


squamous_plot$tumour_budding <- factor( squamous_plot$tumour_budding , levels=c("no_budding", "have_budding"))

boxplot_function <- function(param, pval_label = 0.7){
  color_list <-  c("no_budding" = "#A6CEE3",
                   "homo_have_keratin" = "#B2DF8A",
                   "have_budding"= "#FB9A99")
  #compare_means(param ~ keratinization, data = squamous_plot, method = "wilcox.test")
  g1 <- ggboxplot(squamous_plot, width=0.5,lwd=1, x = "tumour_budding", y = param, fill = "tumour_budding", add = "jitter") + 
    stat_compare_means(comparisons = list( c("no_budding", "have_budding")), label.y = c(pval_label)) +
    theme_cowplot() +
    scale_fill_manual(values = color_list, na.value = '#bababa', drop = F) +
    #  scale_x_discrete(labels = c("Well-differentiated", "Poorly-differentiated")) +
    theme(  axis.title.y = element_text(size=26),
            axis.text.y = element_text(size=26),
            axis.title.x = element_text(size=26),
            axis.text.x = element_text(size=26),
            legend.position = "none",
            legend.title = element_blank()) 
  print(paste0("n number of have_keratin:", nrow(squamous_plot[squamous_plot$keratinization == "have_keratin",]) ))
  print(paste0("n number of no_keratin:", nrow(squamous_plot[squamous_plot$keratinization == "no_keratin",]) ))
  print(g1)
}

boxplot_function("pack_years", 110)


