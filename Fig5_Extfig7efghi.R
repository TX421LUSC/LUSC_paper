##Script for Tx421 LUSC paper
##Fig 5a to 5d
##Ext fig 7e, 7f, 7g, 7h, 7i

### Survival only includes 132 patients
#included 131 patients with 1 LUSC lesion
#included 1 patient with 1 LUSC lesion and 1 carcinoid lesion
#Excluded 8 patients with both LUSC and LUAD lesions (double primary)


#################################
### libraries and directories ###
#################################

library(tidyverse)
library(data.table)
library(fst)
library(readxl)
library(openxlsx)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(ggfittext)
library(cowplot)
library(ggpubr)
library(survival)
library(survminer)
library(cmprsk)
library(car)
library(MASS)
library(oddsratio)
library(ggeasy)
library(nlme)
library(fgsea)
library(GSEABase)

#######################
### load input data ###
#######################

# load per patient data and per tumour data - anonymised
all_patient_df <- readRDS("../data/20221109_TRACERx421_all_patient_df.rds") 
all_tumour_df <- readRDS("../data/20221109_TRACERx421_all_tumour_df.rds") 

# load ctDNA data (cruk_id) 
Tx421_ctDNA_data_LUAD_single <- fread("../data/20221108_Tx421_ctDNA_data_LUAD_single.csv")

## time to relapse data (prepared for the primary paper) 
rec_time_df <- data.table(all_patient_df)[,.(Recurrence_time_use,cruk_id)]

# muttable 
muttable_path <- "../data/20221109_TRACERx421_mutation_table.fst"
muttable_df <- read.fst(muttable_path,
                        as.data.table = T)

#Genome doubling table
gd_df_tumour_path <- "../data/20221109_TRACERx421_gd_table_per_tumour.tsv"
gd_df_tumour_ori <- fread(gd_df_tumour_path)

#include only DriverMut == TRUE    
muttable_driver <- muttable_df %>%
  filter(DriverMut == TRUE)  

# load median vst per tumour 
vsd_median <- readRDS("../data/20221108_vsd_median_df_luad.RDS")

# load S-V gene list
SV_list.ori <- fread("../data/SV.gene.list_from_Cell2018.TableS3.csv")
# set pathway 
pathway_file <- "../data/h.all.v7.2.symbols.gmt"
pathways_hm <- gmtPathways(pathway_file)

squamous_pathdata <- read_csv("../df-alvin/500_entire_tumour_with_grading_with_wGIITMB_info.csv")

max_regional_data <- read_csv("../df-alvin/All_regions_500_regions-SCC_reviewed_With_max_budding_data.csv")
max_regional_data %>% dplyr::select(tumour_id_muttable_cruk, contains("max"), -max_tumour_budding) -> max_regional_data1
max_regional_data1[!duplicated(max_regional_data1), ]  -> max_regional_data2
  
left_join(squamous_pathdata, max_regional_data2 , by = c("cruk_tumourid" = "tumour_id_muttable_cruk")) ->
squamous_pathdata1           

sq_heterohomo <- read_csv("../df-alvin/scc_heterohomo_keratin_budding_sci.csv")

squamous_pathdata1 %>% left_join(sq_heterohomo, by = c("cruk_tumourid" = "tumour_id_muttable_cruk")) ->
  squamous_pathdata1


squamous_pathdata <- squamous_pathdata1
###########################
### organise input data ###
###########################
all_patient_df1 <- all_patient_df  %>%
  filter(histology_lesion1 == "Squamous cell carcinoma" & lesion1_sampled == TRUE) %>%
  filter(histology_multi_full != "LUAD&LUSC") %>% 
  mutate(cens_rec_cum = case_when( 
    cens_lung_event == 1 ~ 1,  # cens_lung_event == 1 : had recurrence/lung cancer death
    cens_dfs == 0 ~ 0,  # cens_dfs == 0 : no recurrence and no death 
    TRUE ~ 2),  # the rest should be non- lung cancer death 
    time_rec_cum = dfs_time) 

all_patient_df_added_luad <- all_patient_df1 %>%
  mutate(Surgery = case_when(Surgery_type %in% c("Segmentectomy or wedge resection", "Lingulectomy") ~ "sublobar resection",
                             Surgery_type %in% c("Bilobectomy", "Lobectomy", "Lobectomy and wedge resection", "Pneumonectomy") ~ "lobectomy or more" )) %>%
  left_join(squamous_pathdata , by = c("tumour_id_per_patient" = "cruk_tumourid")) 
# %>%
#   left_join(subtype_df, by = "tumour_id_per_patient")  

# For survival analysis, there are only 132 patients because two patients with double primary are excluded
# "CRUK0586_Tumour1" and "CRUK0682"
# "CRUK0134" and "CRUK0223_Tumour1"

# all_patient_df  %>%
#   filter(histology_lesion1 == "Squamous cell carcinoma" | histology_lesion2 == "Squamous cell carcinoma") %>% view()
# 
# all_patient_df_added  %>%
#   filter(!is.na(keratinization)) %>% view()
# 
# all_patient_df_added  %>%
#   filter(histology_lesion1 == "Squamous cell carcinoma" & lesion1_sampled == TRUE) %>% view()

data <- all_patient_df_added_luad %>%
  mutate(Relapse_cat_new2 = as.character(Relapse_cat_new) ) %>%
  mutate(Relapse_cat_new2 = factor(Relapse_cat_new2, levels = c("No rec", "Intrathoracic", "Intra & Extra", "Extrathoracic", "Unknown Site"))) %>%
  mutate(Relapse_cat_new = factor(Relapse_cat_new, levels = c("No rec", "Intrathoracic","Intra & Extra", "Extrathoracic"))) %>%
  mutate(TNM_stage = case_when(pathologyTNM %in% c("IA","IB") ~ "I",
                               pathologyTNM %in% c("IIA","IIB") ~ "II",
                               pathologyTNM %in% c("IIIA","IIIB") ~ "III")) %>%
  mutate(pTNM_stage = case_when(pathologyTNM %in% c("IIIA","IIIB") ~ "III",
                                TRUE ~ pathologyTNM),
         pack_years = packyears,
         adjuvant = adjuvant_treatment_YN) %>%
  mutate(surgery = case_when(Surgery_type %in% c("Segmentectomy or wedge resection", "Lingulectomy") ~ "sublobar resection",
                             Surgery_type %in% c("Bilobectomy", "Lobectomy", "Lobectomy and wedge resection", "Pneumonectomy") ~ "lobectomy or more" )) %>%
  mutate(Relapse_cat_new = case_when(cruk_id %in% c("CRUK0173","CRUK0249", "CRUK0781","CRUK0129") ~ "Unknown Site", # "CRUK0173" : relapsed before death but no CT reports available. "CRUK0249", "CRUK0781","CRUK0129": lung cancer death as first event 
                                     TRUE ~ as.character(Relapse_cat_new) )) %>%
  mutate(Relapse_cat_new = factor(Relapse_cat_new, levels = c("No rec", "Intrathoracic", "Intra & Extra", "Extrathoracic", "Unknown Site"))) %>%
  mutate(is.rec = case_when(Relapse_cat_new %in% c("Intrathoracic", "Intra & Extra", "Extrathoracic", "Unknown Site") ~ TRUE, 
                            TRUE ~ FALSE)) %>%
  mutate(time.to.rec = case_when(is.rec == TRUE ~ Recurrence_time_use)) %>%
  mutate(Relapse.Site = case_when(Relapse_cat_new %in% c("Intrathoracic") ~ "Intrathoracic",
                                  Relapse_cat_new %in% c("Intra & Extra", "Extrathoracic") ~ "Extrathoracic")) %>%
  mutate(Relapse.Site = factor(Relapse.Site , levels = c("Intrathoracic", "Extrathoracic"))) %>%
  mutate(  Age = age/10,
           PackYears = packyears/10,
           Adjuvant = adjuvant_treatment_YN) %>%
  mutate(Relapse_Site = case_when(Relapse_cat_new =="No rec" ~ "No rec before death or new primary", 
                                  Relapse_cat_new %in% c("Intra & Extra", "Extrathoracic") ~ "Extrathoracic (any)",
                                  Relapse_cat_new == "Intrathoracic" ~ "Intrathoracic (only)")) %>%
  mutate(Relapse_Site = factor(Relapse_Site, levels = c("No rec before death or new primary", "Intrathoracic (only)", "Extrathoracic (any)")))
  
  
#categorize the subclonal expansion scores and ITH values
min_simpson_tum_median <- median(data$min_simpson_tum_meanccf, na.rm = T)
max_CCF_terminal_tum_median <- median(data$maxCCF_terminalonly_tum_meanccf, na.rm = T)
frac_abberant_genom_subcl_median <- median(data$frac_abberant_genom_subcl, na.rm = T)
perc_subclonal_median <- median(data$perc_subclonal, na.rm = T)

data <- data %>%   mutate(min_simpson_tum_cat = case_when(min_simpson_tum_meanccf >= min_simpson_tum_median  ~ "high",
                                                           min_simpson_tum_meanccf < min_simpson_tum_median  ~ "low"),
                           max_CCF_terminal_tum_cat = case_when(maxCCF_terminalonly_tum_meanccf >= max_CCF_terminal_tum_median  ~ "high",
                                                                maxCCF_terminalonly_tum_meanccf < max_CCF_terminal_tum_median  ~ "low"),
                           frac_abberant_genom_subcl_cat = case_when(frac_abberant_genom_subcl >= frac_abberant_genom_subcl_median  ~ "high",
                                                                     frac_abberant_genom_subcl < frac_abberant_genom_subcl_median  ~ "low"),
                           perc_subclonal_cat = case_when(perc_subclonal >= perc_subclonal_median  ~ "high",
                                                          perc_subclonal < perc_subclonal_median  ~ "low"))


#genome doubling data
gd_df_tumour <- gd_df_tumour_ori %>%
  dplyr::select(tumour_id, num_clonal_gds, num_subclonal_gds) %>%
  mutate(subclonal_gd = case_when(is.na(num_subclonal_gds) ~ NA_character_,
                                  num_subclonal_gds > 0 ~ "Yes",
                                  num_subclonal_gds == 0 ~ "No" ),
         any_gd = case_when(is.na(num_subclonal_gds) ~ NA_character_,
                            num_subclonal_gds > 0 | num_clonal_gds > 0 ~ "Yes",
                            num_subclonal_gds == 0 & num_clonal_gds == 0 ~ "No" )) %>% 
  mutate(GD_status = case_when(subclonal_gd == "Yes" ~ "Subclonal GD",
                                    any_gd == "Yes" ~ "Clonal GD only", 
                                    any_gd == "No" ~ "No GD") )

data <- data %>%
  left_join(gd_df_tumour, by = c("tumour_id_muttable_cruk" = "tumour_id"))


##################################
### prep for survival analysis ### 
##################################


time_vec <- c( "dfs_time") 
status_vec <- c("cens_dfs") 
title_vec <- c( "DFS") 
in_data_surv <- data

# Extfig 6e, Association of LVI and budding and single cell invasion
in_data_surv <- in_data_surv %>% mutate(max_regional_budding = case_when(max_10hpf_budding >=8 ~ "high", 
                                                                         max_10hpf_budding <8 ~ "low_no")) %>% 
  mutate(max_regional_sci = case_when(max_sci >=8 ~ "high", 
                                      max_sci <8 ~ "low_no"))

in_data_surv %>% dplyr::select(LVI_per_patient, tumour_budding) %>% table

in_data_surv %>% dplyr::select(LVI_per_patient, single_cell_invasion) %>% table

table(in_data_surv$LVI_per_patient,in_data_surv$tumour_budding) -> table1

df_gathered1 <- data.frame(var1 = c(rownames(table1)[1], rownames(table1)[1], rownames(table1)[2], rownames(table1)[2] ),
                           var2 = c(colnames(table1)[1],colnames(table1)[2], colnames(table1)[1], colnames(table1)[2]),
                           count = c(table1[1], table1[3], table1[2], table1[4]))
df_gathered1$var2 <- factor( df_gathered1$var2 , levels=c("have_budding", "no_budding"))


fisher_p <- fisher.test(matrix(c(table1[1], table1[2], table1[3], table1[4]), nrow = 2))

color_list <-  c("No" = "#A6CEE3",
                 "greencolour" = "#B2DF8A",
                 "Yes"= "#FB9A99")
g1 <- ggplot(data=df_gathered1, aes(x=var2, y=count, fill=var1, label = count)) +
  geom_bar(stat="identity", position = "fill", color = "black") +
  geom_bar_text( position = "fill", reflow = TRUE) + 
  scale_fill_manual(values = color_list, na.value = '#bababa', drop = F)+
  theme_cowplot() +
  ylab("proportion") +
  ggtitle(paste0( " Fisher-exact test p-value: ", format(round(as.numeric(fisher_p$p.value), 4), nsmall=4))) +
  theme(axis.title.x = element_blank(),  
        axis.title.y = element_text(size=24),
        axis.text.x = element_text(size=20, angle = 45,  hjust=1), 
        axis.text.y = element_text(size=16),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size=18))
print(g1)
# Ext fig 6e LVI vs budding

table(in_data_surv$LVI_per_patient,in_data_surv$single_cell_invasion) -> table1

df_gathered1 <- data.frame(var1 = c(rownames(table1)[1], rownames(table1)[1], rownames(table1)[2], rownames(table1)[2] ),
                           var2 = c(colnames(table1)[1],colnames(table1)[2], colnames(table1)[1], colnames(table1)[2]),
                           count = c(table1[1], table1[3], table1[2], table1[4]))
df_gathered1$var2 <- factor( df_gathered1$var2 , levels=c("have_sci", "no_sci"))


fisher_p <- fisher.test(matrix(c(table1[1], table1[2], table1[3], table1[4]), nrow = 2))

color_list <-  c("No" = "#A6CEE3",
                 "greencolour" = "#B2DF8A",
                 "Yes"= "#FB9A99")
g2 <- ggplot(data=df_gathered1, aes(x=var2, y=count, fill=var1, label = count)) +
  geom_bar(stat="identity", position = "fill", color = "black") +
  geom_bar_text( position = "fill", reflow = TRUE) + 
  scale_fill_manual(values = color_list, na.value = '#bababa', drop = F)+
  theme_cowplot() +
  ylab("proportion") +
  ggtitle(paste0( " Fisher-exact test p-value: ", format(round(as.numeric(fisher_p$p.value), 4), nsmall=4))) +
  theme(axis.title.x = element_blank(),  
        axis.title.y = element_text(size=24),
        axis.text.x = element_text(size=20, angle = 45,  hjust=1), 
        axis.text.y = element_text(size=16),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size=18))
print(g2)
# Ext fig 6e LVI vs single cell invasion


# Plot for keratinization, fig 5a
library(survival)
library(survminer)
library(dplyr)
library(RColorBrewer)

in_df_surv <- in_data_surv  # full dataset
in_df_surv <- in_df_surv %>%
  mutate(keratinization_group = case_when(
    keratinization == "no_keratin" ~ "nonkerat_group",
    keratinization == "have_keratin" ~ "kerat_group",
    TRUE ~ NA_character_
  ))

time_vec <- "dfs_time"
status_vec <- "cens_dfs"
title_vec <- "DFS"
group_vec <- "keratinization_group"

surv_df <- in_df_surv %>%
  mutate(group = keratinization_group,
         time = dfs_time,
         status = cens_dfs) %>%
  filter(!is.na(group), !is.na(status), !is.na(time)) %>%
  mutate(group = factor(group, levels = c("kerat_group", "nonkerat_group")))

fit <- survfit(Surv(time, status) ~ group, data = surv_df)
cox_fit_uni <- coxph(Surv(time, status) ~ group, data = surv_df)
CI_tab_uni <- summary(cox_fit_uni)$conf.int
logrk <- survdiff(Surv(time, status) ~ group, surv_df)

cox_fit <- coxph(Surv(time, status) ~ age + PackYears + TNM_stage + surgery + adjuvant_treatment_YN + group, data = surv_df)
CI_tab <- summary(cox_fit)$conf.int

round_cust <- function(x, digits = 2) format(round(x, digits), nsmall = digits)

label <- paste0(
  'Hazard ratio  kerat_group vs nonkerat_group\n',
  'Uni-Cox: ', round_cust(CI_tab_uni[, 1][1]), 
  ' ; 95% CI: ', round_cust(CI_tab_uni[, 3][1]), '-', round_cust(CI_tab_uni[, 4][1]), 
  ' (P = ', signif(summary(cox_fit_uni)$coefficients[, 5][1], 2), ")\n",
  'Multi-Cox: ', round_cust(CI_tab[, 1][7]), 
  ' ; 95% CI: ', round_cust(CI_tab[, 3][7]), '-', round_cust(CI_tab[, 4][7]), 
  ' (P = ', signif(summary(cox_fit)$coefficients[, 5][7], 2), ")\n",
  'log rank p = ', round(logrk$pvalue, 3)
)

colours <- brewer.pal(3, 'Set1')[1:2]

g <- ggsurvplot(
  fit = fit,
  pval = FALSE,
  pval.method = FALSE,
  risk.table = TRUE,
  legend.title = '',
  font.legend = 16,
  palette = rev(colours),
  xlab = "Days", 
  ylab = "Disease free survival",
  xlim = c(0, 2000), 
  font.tickslab = c(18),
  font.x = c(20),
  font.y = c(20),
  break.x.by = 500,
  title = paste(title_vec, "by", group_vec, "allstage")
)

g$plot <- g$plot +  
  annotate("text", hjust = 0, x = 10, y = 0.15, label = label, size = 5.5)

g$table <- g$table + labs(title = '', y = '', x = '')

print(g)


##plot immune data
TIL_entire <- read_csv("../data/20230408-TIL_STAS_ctDNA/Pathology.TIL.scores_diagnostic.csv")
TIL_region <- read_csv("../data/20230408-TIL_STAS_ctDNA/Pathology.TIL.scores_regional.csv")
immune_data <- read_csv("../data/Tx421_MCPcounter.regional.csv")

immune_data %>% 
  filter(!is.na(`MCP_T cells`)) %>% 
  group_by(ltx_tumourid) %>% 
  summarize(mean_T_cells_MCP = mean(`MCP_T cells`, na.rm = TRUE)) -> immune_data_T

in_data_surv <- in_data_surv %>% left_join(immune_data_T, by = c("ltx_region_id" = "ltx_tumourid")) %>% 
  left_join(TIL_entire, by = c("ltx_clusterid" = "SampleID") )


# Ext fig 7i
TIL_median <- median(in_data_surv$final.TIL, na.rm = T)
in_df_surv <- in_data_surv %>% mutate(TIL_cat = case_when(final.TIL >= TIL_median  ~ "high",
                                                            final.TIL < TIL_median  ~ "low"))

time_vec <- "dfs_time"
status_vec <- "cens_dfs"
title_vec <- "DFS"
group_vec <- "TIL_cat"

surv_df <- in_df_surv %>%
  mutate(group = TIL_cat,
         time = dfs_time,
         status = cens_dfs) %>%
  filter(!is.na(group), !is.na(status), !is.na(time)) %>%
  mutate(group = factor(group, levels = c("high", "low")))

fit <- survfit(Surv(time, status) ~ group, data = surv_df)
cox_fit_uni <- coxph(Surv(time, status) ~ group, data = surv_df)
CI_tab_uni <- summary(cox_fit_uni)$conf.int
logrk <- survdiff(Surv(time, status) ~ group, surv_df)

cox_fit <- coxph(Surv(time, status) ~ age + PackYears + TNM_stage + surgery + adjuvant_treatment_YN + group, data = surv_df)
CI_tab <- summary(cox_fit)$conf.int

round_cust <- function(x, digits = 2) format(round(x, digits), nsmall = digits)

label <- paste0(
  'Hazard ratio  high vs low\n',
  'Uni-Cox: ', round_cust(CI_tab_uni[, 1][1]), 
  ' ; 95% CI: ', round_cust(CI_tab_uni[, 3][1]), '-', round_cust(CI_tab_uni[, 4][1]), 
  ' (P = ', signif(summary(cox_fit_uni)$coefficients[, 5][1], 2), ")\n",
  'Multi-Cox: ', round_cust(CI_tab[, 1][7]), 
  ' ; 95% CI: ', round_cust(CI_tab[, 3][7]), '-', round_cust(CI_tab[, 4][7]), 
  ' (P = ', signif(summary(cox_fit)$coefficients[, 5][7], 2), ")\n",
  'log rank p = ', round(logrk$pvalue, 3)
)

colours <- brewer.pal(3, 'Set1')[1:2]

g <- ggsurvplot(
  fit = fit,
  pval = FALSE,
  pval.method = FALSE,
  risk.table = TRUE,
  legend.title = '',
  font.legend = 16,
  palette = rev(colours),
  xlab = "Days", 
  ylab = "Disease free survival",
  xlim = c(0, 2000), 
  font.tickslab = c(18),
  font.x = c(20),
  font.y = c(20),
  break.x.by = 500,
  title = paste(title_vec, "by", group_vec, "allstage")
)

g$plot <- g$plot +  
  annotate("text", hjust = 0, x = 10, y = 0.15, label = label, size = 5.5)

g$table <- g$table + labs(title = '', y = '', x = '')



#Ext fig 7i
TIL_RNA_median <- median(in_data_surv$mean_T_cells_MCP, na.rm = T)
in_df_surv <- in_data_surv %>% mutate(TIL_RNA_cat = case_when(mean_T_cells_MCP >= TIL_RNA_median  ~ "high",
                                                                mean_T_cells_MCP < TIL_RNA_median  ~ "low"))

time_vec <- "dfs_time"
status_vec <- "cens_dfs"
title_vec <- "DFS"
group_vec <- "TIL_RNA_cat"

surv_df <- in_df_surv %>%
  mutate(group = TIL_RNA_cat,
         time = dfs_time,
         status = cens_dfs) %>%
  filter(!is.na(group), !is.na(status), !is.na(time)) %>%
  mutate(group = factor(group, levels = c("high", "low")))

fit <- survfit(Surv(time, status) ~ group, data = surv_df)
cox_fit_uni <- coxph(Surv(time, status) ~ group, data = surv_df)
CI_tab_uni <- summary(cox_fit_uni)$conf.int
logrk <- survdiff(Surv(time, status) ~ group, surv_df)

cox_fit <- coxph(Surv(time, status) ~ age + PackYears + TNM_stage + surgery + adjuvant_treatment_YN + group, data = surv_df)
CI_tab <- summary(cox_fit)$conf.int

round_cust <- function(x, digits = 2) format(round(x, digits), nsmall = digits)

label <- paste0(
  'Hazard ratio  high vs low\n',
  'Uni-Cox: ', round_cust(CI_tab_uni[, 1][1]), 
  ' ; 95% CI: ', round_cust(CI_tab_uni[, 3][1]), '-', round_cust(CI_tab_uni[, 4][1]), 
  ' (P = ', signif(summary(cox_fit_uni)$coefficients[, 5][1], 2), ")\n",
  'Multi-Cox: ', round_cust(CI_tab[, 1][7]), 
  ' ; 95% CI: ', round_cust(CI_tab[, 3][7]), '-', round_cust(CI_tab[, 4][7]), 
  ' (P = ', signif(summary(cox_fit)$coefficients[, 5][7], 2), ")\n",
  'log rank p = ', round(logrk$pvalue, 3)
)

colours <- brewer.pal(3, 'Set1')[1:2]

g <- ggsurvplot(
  fit = fit,
  pval = FALSE,
  pval.method = FALSE,
  risk.table = TRUE,
  legend.title = '',
  font.legend = 16,
  palette = rev(colours),
  xlab = "Days", 
  ylab = "Disease free survival",
  xlim = c(0, 2000), 
  font.tickslab = c(18),
  font.x = c(20),
  font.y = c(20),
  break.x.by = 500,
  title = paste(title_vec, "by", group_vec, "allstage")
)

g$plot <- g$plot +  
  annotate("text", hjust = 0, x = 10, y = 0.15, label = label, size = 5.5)

g$table <- g$table + labs(title = '', y = '', x = '')

print(g)



### Fibroblasts-mean-RNA score
#Ext fig 7h
immune_data %>% 
  filter(!is.na(MCP_Fibroblasts)) %>% 
  group_by(ltx_tumourid) %>% 
  summarize(mean_fibroblasts_cells_MCP = mean(MCP_Fibroblasts, na.rm = TRUE)) -> immune_data_fb

in_data_surv <- in_data_surv %>% left_join(immune_data_fb, by = c("ltx_region_id" = "ltx_tumourid"))

fibroblastscell_RNA_median <- median(in_data_surv$mean_fibroblasts_cells_MCP, na.rm = T)

in_df_surv <- in_data_surv %>% mutate(fibroblastscell_RNA_cat = case_when(mean_fibroblasts_cells_MCP >= fibroblastscell_RNA_median  ~ "high",
                                                                            mean_fibroblasts_cells_MCP < fibroblastscell_RNA_median  ~ "low"))

time_vec <- "dfs_time"
status_vec <- "cens_dfs"
title_vec <- "DFS"
group_vec <- "fibroblastscell_RNA_cat"

surv_df <- in_df_surv %>%
  mutate(group = fibroblastscell_RNA_cat,
         time = dfs_time,
         status = cens_dfs) %>%
  filter(!is.na(group), !is.na(status), !is.na(time)) %>%
  mutate(group = factor(group, levels = c("low", "high")))

fit <- survfit(Surv(time, status) ~ group, data = surv_df)
cox_fit_uni <- coxph(Surv(time, status) ~ group, data = surv_df)
CI_tab_uni <- summary(cox_fit_uni)$conf.int
logrk <- survdiff(Surv(time, status) ~ group, surv_df)

cox_fit <- coxph(Surv(time, status) ~ age + PackYears + TNM_stage + surgery + adjuvant_treatment_YN + group, data = surv_df)
CI_tab <- summary(cox_fit)$conf.int

round_cust <- function(x, digits = 2) format(round(x, digits), nsmall = digits)

label <- paste0(
  'Hazard ratio  high vs low\n',
  'Uni-Cox: ', round_cust(CI_tab_uni[, 1][1]), 
  ' ; 95% CI: ', round_cust(CI_tab_uni[, 3][1]), '-', round_cust(CI_tab_uni[, 4][1]), 
  ' (P = ', signif(summary(cox_fit_uni)$coefficients[, 5][1], 2), ")\n",
  'Multi-Cox: ', round_cust(CI_tab[, 1][7]), 
  ' ; 95% CI: ', round_cust(CI_tab[, 3][7]), '-', round_cust(CI_tab[, 4][7]), 
  ' (P = ', signif(summary(cox_fit)$coefficients[, 5][7], 2), ")\n",
  'log rank p = ', round(logrk$pvalue, 3)
)

colours <- brewer.pal(3, 'Set1')[1:2]

g <- ggsurvplot(
  fit = fit,
  pval = FALSE,
  pval.method = FALSE,
  risk.table = TRUE,
  legend.title = '',
  font.legend = 16,
  palette = rev(colours),
  xlab = "Days", 
  ylab = "Disease free survival",
  xlim = c(0, 2000), 
  font.tickslab = c(18),
  font.x = c(20),
  font.y = c(20),
  break.x.by = 500,
  title = paste(title_vec, "by", group_vec, "allstage")
)

g$plot <- g$plot +  
  annotate("text", hjust = 0, x = 10, y = 0.15, label = label, size = 5.5)

g$table <- g$table + labs(title = '', y = '', x = '')

print(g)


### Fibroblasts-max-RNA score
#Fig 5c
immune_data %>% 
  filter(!is.na(`MCP_Fibroblasts`)) %>% 
  group_by(ltx_tumourid) %>% 
  summarize(max_fibroblasts_cells_MCP = max(`MCP_Fibroblasts`, na.rm = TRUE)) -> immune_data_maxfb

in_data_surv <- in_data_surv %>% left_join(immune_data_maxfb, by = c("ltx_region_id" = "ltx_tumourid"))

fibroblastscell_RNA_max <- median(in_data_surv$max_fibroblasts_cells_MCP, na.rm = T)

in_df_surv <- in_data_surv %>% mutate(fibroblastscell_RNA_max = case_when(max_fibroblasts_cells_MCP >= fibroblastscell_RNA_max  ~ "high",
                                                                          max_fibroblasts_cells_MCP < fibroblastscell_RNA_max  ~ "low"))


time_vec <- "dfs_time"
status_vec <- "cens_dfs"
title_vec <- "DFS"
group_vec <- "fibroblastscell_RNA_max"

surv_df <- in_df_surv %>%
  mutate(group = fibroblastscell_RNA_max,
         time = dfs_time,
         status = cens_dfs) %>%
  filter(!is.na(group), !is.na(status), !is.na(time)) %>%
  mutate(group = factor(group, levels = c("low", "high")))

fit <- survfit(Surv(time, status) ~ group, data = surv_df)
cox_fit_uni <- coxph(Surv(time, status) ~ group, data = surv_df)
CI_tab_uni <- summary(cox_fit_uni)$conf.int
logrk <- survdiff(Surv(time, status) ~ group, surv_df)

cox_fit <- coxph(Surv(time, status) ~ age + PackYears + TNM_stage + surgery + adjuvant_treatment_YN + group, data = surv_df)
CI_tab <- summary(cox_fit)$conf.int

round_cust <- function(x, digits = 2) format(round(x, digits), nsmall = digits)

label <- paste0(
  'Hazard ratio  high vs low\n',
  'Uni-Cox: ', round_cust(CI_tab_uni[, 1][1]), 
  ' ; 95% CI: ', round_cust(CI_tab_uni[, 3][1]), '-', round_cust(CI_tab_uni[, 4][1]), 
  ' (P = ', signif(summary(cox_fit_uni)$coefficients[, 5][1], 2), ")\n",
  'Multi-Cox: ', round_cust(CI_tab[, 1][7]), 
  ' ; 95% CI: ', round_cust(CI_tab[, 3][7]), '-', round_cust(CI_tab[, 4][7]), 
  ' (P = ', signif(summary(cox_fit)$coefficients[, 5][7], 2), ")\n",
  'log rank p = ', round(logrk$pvalue, 3)
)

colours <- brewer.pal(3, 'Set1')[1:2]

g <- ggsurvplot(
  fit = fit,
  pval = FALSE,
  pval.method = FALSE,
  risk.table = TRUE,
  legend.title = '',
  font.legend = 16,
  palette = rev(colours),
  xlab = "Days", 
  ylab = "Disease free survival",
  xlim = c(0, 2000), 
  font.tickslab = c(18),
  font.x = c(20),
  font.y = c(20),
  break.x.by = 500,
  title = paste(title_vec, "by", group_vec, "allstage")
)

g$plot <- g$plot +  
  annotate("text", hjust = 0, x = 10, y = 0.15, label = label, size = 5.5)

g$table <- g$table + labs(title = '', y = '', x = '')

print(g)


# max tumour budding Fig 5b
in_data_surv <- in_data_surv %>% mutate(max_regional_budding_high.low.no = case_when(max_10hpf_budding >=8 ~ "high", 
                                                                                     max_10hpf_budding <8 & max_10hpf_budding >0 ~ "low",
                                                                                     max_10hpf_budding == 0 ~ "no")) %>% 
  mutate(max_regional_sci_high.low.no = case_when(max_sci >=8 ~ "high", 
                                                  max_sci <8 & max_sci >0 ~ "low", 
                                                  max_sci == 0 ~ "no"))


in_df_surv <- in_data_surv
in_df_surv$group <- in_df_surv[["max_regional_budding_high.low.no"]]

surv_df <- in_df_surv %>%
  dplyr::filter(!is.na(group), !is.na(cens_dfs)) %>%
  mutate(group = factor(group,  levels = c("no", "low", "high")))

cats <- levels(surv_df$group)

fit <- survfit(Surv(dfs_time, cens_dfs) ~ group, data = surv_df)
round_cust <- function(x, digits = 2) format(round(x, digits), nsmall = digits)

cox_fit_uni <- coxph(Surv(dfs_time, cens_dfs) ~ group, data = surv_df)
CI_tab_uni <- summary(cox_fit_uni)$conf.int

cox_fit <- coxph(Surv(dfs_time, cens_dfs) ~ age + PackYears + TNM_stage + surgery + adjuvant_treatment_YN + group, data = surv_df)
CI_tab <- summary(cox_fit)$conf.int

logrk <- survdiff(Surv(dfs_time, cens_dfs) ~ group, surv_df)

label <- paste0(
  'Hazard ratio vs ', cats[1], '\n',
  cats[2], ': Uni-Cox: ', round_cust(CI_tab_uni[, 1][1]), ' ; 95% CI: ', round_cust(CI_tab_uni[, 3][1]), '-', round_cust(CI_tab_uni[, 4][1]),
  ' (P = ', signif(summary(cox_fit_uni)$coefficients[, 5][1], 2), ")\n",
  'Multi-Cox: ', round_cust(CI_tab[, 1][7]), ' ; 95% CI: ', round_cust(CI_tab[, 3][7]), '-', round_cust(CI_tab[, 4][7]),
  ' (P = ', signif(summary(cox_fit)$coefficients[, 5][7], 2), ")\n",
  cats[3], ': Uni-Cox: ', round_cust(CI_tab_uni[, 1][2]), ' ; 95% CI: ', round_cust(CI_tab_uni[, 3][2]), '-', round_cust(CI_tab_uni[, 4][2]),
  ' (P = ', signif(summary(cox_fit_uni)$coefficients[, 5][2], 2), ")\n",
  'Multi-Cox: ', round_cust(CI_tab[, 1][8]), ' ; 95% CI: ', round_cust(CI_tab[, 3][8]), '-', round_cust(CI_tab[, 4][8]),
  ' (P = ', signif(summary(cox_fit)$coefficients[, 5][8], 2), ")\n",
  'log rank p = ', round(logrk$pvalue, 3)
)

colours <- brewer.pal(3, 'Set1')[1:3]
g <- ggsurvplot(
  fit = fit,
  risk.table = TRUE,
  legend.title = '',
  font.legend = 16,
  palette = rev(colours),
  xlab = "Days",
  ylab = "Disease free survival",
  xlim = c(0, 2000),
  font.tickslab = c(18),
  font.x = c(20),
  font.y = c(20),
  break.x.by = 500,
  title = "DFS by max_regional_budding  allstage"
)

g$plot <- g$plot +
  annotate("text", hjust = 0, x = 10, y = 0.15, label = label, size = 5.5)
g$table <- g$table + labs(title = '', y = '', x = '')

print(g)



#max single cell invasion 
#Ext fig 7f

in_df_surv <- in_data_surv
in_df_surv$group <- in_df_surv[["max_regional_sci_high.low.no"]]

surv_df <- in_df_surv %>%
  dplyr::filter(!is.na(group), !is.na(cens_dfs)) %>%
  mutate(group = factor(group,  levels = c("no", "low", "high")))

cats <- levels(surv_df$group)

fit <- survfit(Surv(dfs_time, cens_dfs) ~ group, data = surv_df)
round_cust <- function(x, digits = 2) format(round(x, digits), nsmall = digits)

cox_fit_uni <- coxph(Surv(dfs_time, cens_dfs) ~ group, data = surv_df)
CI_tab_uni <- summary(cox_fit_uni)$conf.int

cox_fit <- coxph(Surv(dfs_time, cens_dfs) ~ age + PackYears + TNM_stage + surgery + adjuvant_treatment_YN + group, data = surv_df)
CI_tab <- summary(cox_fit)$conf.int

logrk <- survdiff(Surv(dfs_time, cens_dfs) ~ group, surv_df)

label <- paste0(
  'Hazard ratio vs ', cats[1], '\n',
  cats[2], ': Uni-Cox: ', round_cust(CI_tab_uni[, 1][1]), ' ; 95% CI: ', round_cust(CI_tab_uni[, 3][1]), '-', round_cust(CI_tab_uni[, 4][1]),
  ' (P = ', signif(summary(cox_fit_uni)$coefficients[, 5][1], 2), ")\n",
  'Multi-Cox: ', round_cust(CI_tab[, 1][7]), ' ; 95% CI: ', round_cust(CI_tab[, 3][7]), '-', round_cust(CI_tab[, 4][7]),
  ' (P = ', signif(summary(cox_fit)$coefficients[, 5][7], 2), ")\n",
  cats[3], ': Uni-Cox: ', round_cust(CI_tab_uni[, 1][2]), ' ; 95% CI: ', round_cust(CI_tab_uni[, 3][2]), '-', round_cust(CI_tab_uni[, 4][2]),
  ' (P = ', signif(summary(cox_fit_uni)$coefficients[, 5][2], 2), ")\n",
  'Multi-Cox: ', round_cust(CI_tab[, 1][8]), ' ; 95% CI: ', round_cust(CI_tab[, 3][8]), '-', round_cust(CI_tab[, 4][8]),
  ' (P = ', signif(summary(cox_fit)$coefficients[, 5][8], 2), ")\n",
  'log rank p = ', round(logrk$pvalue, 3)
)

colours <- brewer.pal(3, 'Set1')[1:3]
g <- ggsurvplot(
  fit = fit,
  risk.table = TRUE,
  legend.title = '',
  font.legend = 16,
  palette = rev(colours),
  xlab = "Days",
  ylab = "Disease free survival",
  xlim = c(0, 2000),
  font.tickslab = c(18),
  font.x = c(20),
  font.y = c(20),
  break.x.by = 500,
  title = "DFS by max_single_cell_invasion  allstage"
)

g$plot <- g$plot +
  annotate("text", hjust = 0, x = 10, y = 0.15, label = label, size = 5.5)
g$table <- g$table + labs(title = '', y = '', x = '')

print(g)



#tumour budding whole tumour level
#Ext fig 7g

in_df_surv <- in_data_surv
in_df_surv$group <- in_df_surv[["budding_high.low"]]

surv_df <- in_df_surv %>%
  dplyr::filter(!is.na(group), !is.na(cens_dfs)) %>%
  mutate(group = factor(group,  levels = c("no_budding", "low_budding", "high_budding")))

cats <- levels(surv_df$group)

fit <- survfit(Surv(dfs_time, cens_dfs) ~ group, data = surv_df)
round_cust <- function(x, digits = 2) format(round(x, digits), nsmall = digits)

cox_fit_uni <- coxph(Surv(dfs_time, cens_dfs) ~ group, data = surv_df)
CI_tab_uni <- summary(cox_fit_uni)$conf.int

cox_fit <- coxph(Surv(dfs_time, cens_dfs) ~ age + PackYears + TNM_stage + surgery + adjuvant_treatment_YN + group, data = surv_df)
CI_tab <- summary(cox_fit)$conf.int

logrk <- survdiff(Surv(dfs_time, cens_dfs) ~ group, surv_df)

label <- paste0(
  'Hazard ratio vs ', cats[1], '\n',
  cats[2], ': Uni-Cox: ', round_cust(CI_tab_uni[, 1][1]), ' ; 95% CI: ', round_cust(CI_tab_uni[, 3][1]), '-', round_cust(CI_tab_uni[, 4][1]),
  ' (P = ', signif(summary(cox_fit_uni)$coefficients[, 5][1], 2), ")\n",
  'Multi-Cox: ', round_cust(CI_tab[, 1][7]), ' ; 95% CI: ', round_cust(CI_tab[, 3][7]), '-', round_cust(CI_tab[, 4][7]),
  ' (P = ', signif(summary(cox_fit)$coefficients[, 5][7], 2), ")\n",
  cats[3], ': Uni-Cox: ', round_cust(CI_tab_uni[, 1][2]), ' ; 95% CI: ', round_cust(CI_tab_uni[, 3][2]), '-', round_cust(CI_tab_uni[, 4][2]),
  ' (P = ', signif(summary(cox_fit_uni)$coefficients[, 5][2], 2), ")\n",
  'Multi-Cox: ', round_cust(CI_tab[, 1][8]), ' ; 95% CI: ', round_cust(CI_tab[, 3][8]), '-', round_cust(CI_tab[, 4][8]),
  ' (P = ', signif(summary(cox_fit)$coefficients[, 5][8], 2), ")\n",
  'log rank p = ', round(logrk$pvalue, 3)
)

colours <- brewer.pal(3, 'Set1')[1:3]
g <- ggsurvplot(
  fit = fit,
  risk.table = TRUE,
  legend.title = '',
  font.legend = 16,
  palette = rev(colours),
  xlab = "Days",
  ylab = "Disease free survival",
  xlim = c(0, 2000),
  font.tickslab = c(18),
  font.x = c(20),
  font.y = c(20),
  break.x.by = 500,
  title = "DFS by tumour budding on whole tumour level  allstage"
)

g$plot <- g$plot +
  annotate("text", hjust = 0, x = 10, y = 0.15, label = label, size = 5.5)
g$table <- g$table + labs(title = '', y = '', x = '')

print(g)



#single cell invasion whole tumour level
#Ext fig 7g

in_df_surv <- in_data_surv
in_df_surv$group <- in_df_surv[["sci_high.low"]]

surv_df <- in_df_surv %>%
  dplyr::filter(!is.na(group), !is.na(cens_dfs)) %>%
  mutate(group = factor(group,  levels = c("no_sci", "low_sci", "high_sci")))

cats <- levels(surv_df$group)

fit <- survfit(Surv(dfs_time, cens_dfs) ~ group, data = surv_df)
round_cust <- function(x, digits = 2) format(round(x, digits), nsmall = digits)

cox_fit_uni <- coxph(Surv(dfs_time, cens_dfs) ~ group, data = surv_df)
CI_tab_uni <- summary(cox_fit_uni)$conf.int

cox_fit <- coxph(Surv(dfs_time, cens_dfs) ~ age + PackYears + TNM_stage + surgery + adjuvant_treatment_YN + group, data = surv_df)
CI_tab <- summary(cox_fit)$conf.int

logrk <- survdiff(Surv(dfs_time, cens_dfs) ~ group, surv_df)

label <- paste0(
  'Hazard ratio vs ', cats[1], '\n',
  cats[2], ': Uni-Cox: ', round_cust(CI_tab_uni[, 1][1]), ' ; 95% CI: ', round_cust(CI_tab_uni[, 3][1]), '-', round_cust(CI_tab_uni[, 4][1]),
  ' (P = ', signif(summary(cox_fit_uni)$coefficients[, 5][1], 2), ")\n",
  'Multi-Cox: ', round_cust(CI_tab[, 1][7]), ' ; 95% CI: ', round_cust(CI_tab[, 3][7]), '-', round_cust(CI_tab[, 4][7]),
  ' (P = ', signif(summary(cox_fit)$coefficients[, 5][7], 2), ")\n",
  cats[3], ': Uni-Cox: ', round_cust(CI_tab_uni[, 1][2]), ' ; 95% CI: ', round_cust(CI_tab_uni[, 3][2]), '-', round_cust(CI_tab_uni[, 4][2]),
  ' (P = ', signif(summary(cox_fit_uni)$coefficients[, 5][2], 2), ")\n",
  'Multi-Cox: ', round_cust(CI_tab[, 1][8]), ' ; 95% CI: ', round_cust(CI_tab[, 3][8]), '-', round_cust(CI_tab[, 4][8]),
  ' (P = ', signif(summary(cox_fit)$coefficients[, 5][8], 2), ")\n",
  'log rank p = ', round(logrk$pvalue, 3)
)

colours <- brewer.pal(3, 'Set1')[1:3]
g <- ggsurvplot(
  fit = fit,
  risk.table = TRUE,
  legend.title = '',
  font.legend = 16,
  palette = rev(colours),
  xlab = "Days",
  ylab = "Disease free survival",
  xlim = c(0, 2000),
  font.tickslab = c(18),
  font.x = c(20),
  font.y = c(20),
  break.x.by = 500,
  title = "DFS by single cell invasion on whole tumour level  allstage"
)

g$plot <- g$plot +
  annotate("text", hjust = 0, x = 10, y = 0.15, label = label, size = 5.5)
g$table <- g$table + labs(title = '', y = '', x = '')

print(g)




############################
### Multivariable : forest plot ### 
############################

#Forest plot  Fig 5d
write_rds(in_data_surv, "../data/in_data_surv_temp.RDS")

in_df_surv_pre <- in_data_surv %>%
  mutate(TNMstage = case_when(pathologyTNM %in% c("IA","IB") ~ "I",
                              pathologyTNM %in% c("IIA","IIB") ~ "II",
                              pathologyTNM %in% c("IIIA","IIIB") ~ "III")) %>%
  mutate(Surgery = case_when(Surgery_type %in% c("Segmentectomy or wedge resection", "Lingulectomy") ~ "sublobar resection",
                             Surgery_type %in% c("Bilobectomy", "Lobectomy", "Lobectomy and wedge resection", "Pneumonectomy") ~ "lobectomy or more" )) %>%
  mutate(Adjuvant = factor(adjuvant_treatment_YN, levels = c("No adjuvant", "Adjuvant"))) %>%
  mutate(budding_high.low = factor(budding_high.low, levels = c("no_budding", "low_budding", "high_budding"))) %>% 
  mutate(differentiation = factor(differentiation, levels = c("well_dx", "poor_dx"))) %>%         # set factor level where the first one will be the reference in forest plot
  mutate(max_regional_budding_high.low.no = factor(max_regional_budding_high.low.no, levels = c("no", "low", "high"))) %>% 
  mutate(max_regional_sci_high.low.no = factor(max_regional_sci_high.low.no, levels = c("no", "low", "high"))) %>% 
  mutate(GD_status = factor(GD_status, levels = c("No GD", "Clonal GD only", "Subclonal GD"))) %>% 
  mutate(max_CCF_terminal_tum_cat = factor(max_CCF_terminal_tum_cat, levels = c("low", "high"))) %>% 
  mutate(min_simpson_tum_cat = factor(min_simpson_tum_cat, levels = c("low", "high"))) %>% 
  mutate(frac_abberant_genom_subcl_cat=factor(frac_abberant_genom_subcl_cat, levels = c("low", "high"))) %>% 
  mutate(perc_subclonal_cat=factor(perc_subclonal_cat, levels = c("low", "high"))) %>% 
  mutate(max_regional_budding = factor(max_regional_budding, levels = c("low_no", "high"))) %>% 
  mutate(max_regional_budding_high.low.no = factor(max_regional_budding_high.low.no, levels = c("no", "low", "high"))) %>% 
  mutate(max_regional_sci = factor(max_regional_sci, levels = c("low_no", "high"))) %>% 
  mutate(max_regional_sci_high.low.no = factor(max_regional_sci_high.low.no, levels = c("no", "low", "high"))) %>%
  mutate(keratinization = factor(keratinization, levels = c("no_keratin", "have_keratin"))) %>% 
  mutate(extensive_necrosis1 = case_when(extensive_necrosis == 1 ~ "Present", 
                                         extensive_necrosis == 0 ~ "Absent") ) %>% 
  mutate(extensive_necrosis1 = factor(extensive_necrosis1, levels = c("Absent", "Present"))) %>% 
  mutate(PackYears = packyears/10,
    Age = age/10 )



time_vec <- c("dfs_time")
status_vec <- c("cens_dfs")
title_vec <- c("Disease-free survival")

variables_cox <- c("Age", "TNM_stage", "PackYears", 
                   "Adjuvant" , "keratinization" , 
                   "max_regional_budding")

in_df_surv <- in_df_surv_pre
for(var in c(variables_cox) ) {
  in_df_surv <- in_df_surv[!is.na(in_df_surv[[var]]),]
}

for (m in 1:length(time_vec)) {
  
  survival_data.frame <-in_df_surv
  colnames_time <-  time_vec[m]
  colnames_status <-  status_vec[m]
  
  surv_df <- survival_data.frame
  surv_df$time <- surv_df[[colnames_time]]
  surv_df$status <- surv_df[[colnames_status]]  
  
  cox_out <- coxph(Surv(time, status) ~  Age + TNMstage + PackYears + 
                     Adjuvant + keratinization + 
                     max_regional_budding,
                   data = surv_df)
  
  fp <- ggforest(cox_out, data = NULL, main = "Hazard ratio", fontsize = 0.7, refLabel = "reference", noDigits = 2)
  print(fp)
}

