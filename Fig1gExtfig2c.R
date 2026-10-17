## Script for Tx421 LUSC paper
#Fig 1g
#Ext fig 2c

#################################
### libraries and directories ###
#################################

library(tidyverse)
library(data.table)
library(readxl)
library(ggfittext)
library(cowplot)
library(ggpubr)
library(nlme)
library(ComplexHeatmap)
library(dendsort)
library(seriation)
library(car)
library(reshape2)
library(fst)

# chr arm list
arm_list <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q","7p","7q",
              "8p","8q","9p","9q","10p","10q","11p","11q","12p","12q",
              "13q","14q","15q","16p","16q","17p","17q","18p","18q",
              "19p","19q","20p","20q","21q","22q")

LUSC_region_df_full <- read_csv("../df-alvin/All_regions_500_regions-SCC.csv")
all_tumours_genes_cn <- readRDS("../data/all_tumours_adjusted_cn_all_genes.RDS")
mutable <- read_fst("../data/20221109_TRACERx421_mutation_table.fst")
idconv <- read_csv("../df-alvin/id_conversion.csv")

driver <- mutable$Hugo_Symbol[mutable$DriverMut == TRUE] %>% unique
all_tumours_genes_cn %>% filter(gene_name %in% driver) -> all_tumours_drivergenes_cn0
all_tumours_drivergenes_cn0 %>%
  left_join(idconv[, c("pathology_id","sample_name_cruk","tumour_id_muttable_cruk")],
            by = c("sample" = "pathology_id")) ->
  all_tumours_drivergenes_cn

Tx421_weigted_gene_mean.adj.cpn_df <- dcast(
  all_tumours_drivergenes_cn,
  sample_name_cruk + sample + tumour_id_muttable_cruk ~ gene_name,
  value.var = "tot_adj_raw",
  fun.aggregate = function(x) mean(x, na.rm = TRUE)
)
Tx421_weigted_gene_mean.adj.cpn_df %>%
  select(-tumour_id_muttable_cruk) ->
  Tx421_weigted_gene_mean.adj.cpn_df

max_cp_df <- fread("../data/20221102_TRACERx421_tumourDivCorrectedtree.csv")

### arm & driver gene cn ###

tumour_subtype_df <- read_csv("../df-alvin/500_entire_tumour_with_grading.csv")

tumour_subtype_df %>%
  mutate(nonkerat_vs_kerat_group = case_when(
    keratinization %in% c("have_keratin") ~ "kerat_group",
    keratinization %in% c("no_keratin")  ~ "nonkerat_group"
  )) %>%
  mutate(nonkerat_vs_kerat_group = factor(nonkerat_vs_kerat_group,
                                          levels = c("kerat_group","nonkerat_group"))) ->
  tumour_subtype_df

LUSC_region_arm_gene_cn <- LUSC_region_df_full %>%
  left_join(Tx421_weigted_gene_mean.adj.cpn_df, by = "sample_name_cruk") %>%
  left_join(tumour_subtype_df,
            by = c("tumour_id_muttable_cruk" = "cruk_tumourid"))

var_list <- c(
  colnames(Tx421_weigted_gene_mean.adj.cpn_df)[
    !grepl("_cruk", colnames(Tx421_weigted_gene_mean.adj.cpn_df)) &
      !grepl("sample", colnames(Tx421_weigted_gene_mean.adj.cpn_df))
  ],
  arm_list
)

group_list <- c("nonkerat_vs_kerat_group")

all_res_vec <- vector()
for (group in group_list) {
  for (var in var_list) {
    tmp_df <- LUSC_region_arm_gene_cn
    tmp_df$var <- tmp_df[[var]]
    tmp_df$group <- tmp_df[[group]]
    tmp_df <- tmp_df %>% filter(!is.na(var), !is.na(group), !is.na(wGII))
    
    effect_model <- lme(var ~ group, random = ~1 | tumour_id_muttable_cruk,
                        data = tmp_df, method = "ML")
    res_table <- as.data.frame(coef(summary(effect_model)), stringsAsFactors = FALSE)
    res_vec <- c(group, var,
                 res_table$Value[2],
                 res_table$`p-value`[2])
    all_res_vec <- rbind(all_res_vec, res_vec)
  }
}
all_res_dt <- as.data.table(all_res_vec)
colnames(all_res_dt) <- c("Group","variable","coefficient","pval")
all_res_dt <- all_res_dt %>%
  mutate(coefficient = as.numeric(coefficient),
         pval        = as.numeric(pval)) %>%
  mutate(p.adj = p.adjust(pval, method = "BH"))

all_res_dt_SCNA_lme <- all_res_dt

# calculate mean in kerat_group and nonkerat_group (mean of mean tumour)
region_cn <- LUSC_region_arm_gene_cn %>%
  select(tumour_id_muttable_cruk, `1p`:`22q`, ACTN2:ZNF483)
tumour_arm_gene_cn_mean <- sapply(
  unique(region_cn$tumour_id_muttable_cruk),
  function(x) {
    tmp <- region_cn[region_cn$tumour_id_muttable_cruk == x, ]
    colMeans(tmp[ , -1], na.rm = TRUE)
  }
)

tumour_subtype_df$tumour_id_muttable_cruk <- tumour_subtype_df$cruk_tumourid
kerat_tumours <- tumour_subtype_df$tumour_id_muttable_cruk[
  tumour_subtype_df$nonkerat_vs_kerat_group == "kerat_group"
]
nonkerat_tumours <- tumour_subtype_df$tumour_id_muttable_cruk[
  tumour_subtype_df$nonkerat_vs_kerat_group == "nonkerat_group"
]

mean_kerat <- as.data.frame(t(tumour_arm_gene_cn_mean)) %>%
  rownames_to_column() %>%
  rename(tumour_id_muttable_cruk = rowname) %>%
  filter(tumour_id_muttable_cruk %in% kerat_tumours) %>%
  select(-tumour_id_muttable_cruk) %>%
  colMeans(na.rm = TRUE)

mean_nonkerat <- as.data.frame(t(tumour_arm_gene_cn_mean)) %>%
  rownames_to_column() %>%
  rename(tumour_id_muttable_cruk = rowname) %>%
  filter(tumour_id_muttable_cruk %in% nonkerat_tumours) %>%
  select(-tumour_id_muttable_cruk) %>%
  colMeans(na.rm = TRUE)

mean_kerat_df <- data.frame(variable = names(mean_kerat),
                            mean_kerat = mean_kerat) %>%
  mutate(direction_kerat = ifelse(mean_kerat > 0,
                                  "kerat_group (mean) > 0",
                                  "kerat_group (mean) < 0")) %>%
  select(variable, direction_kerat)

mean_nonkerat_df <- data.frame(variable = names(mean_nonkerat),
                               mean_nonkerat = mean_nonkerat) %>%
  mutate(direction_nonkerat = ifelse(mean_nonkerat > 0,
                                     "nonkerat_group (mean) > 0",
                                     "nonkerat_group (mean) < 0")) %>%
  select(variable, direction_nonkerat)

## create volcano plot
in_df <- all_res_dt_SCNA_lme %>%
  left_join(mean_kerat_df, by = "variable") %>%
  left_join(mean_nonkerat_df, by = "variable") %>%
  mutate(id = variable,
         sig = p.adj < 0.05,
         id.sig = ifelse(sig, id, ""),
         col_direction = case_when(
           sig & direction_nonkerat == "nonkerat_group (mean) > 0" &
             direction_kerat     == "kerat_group (mean) > 0" ~ ">0 in nonkerat and kerat",
           sig & direction_nonkerat == "nonkerat_group (mean) > 0" &
             direction_kerat     == "kerat_group (mean) < 0" ~ ">0 nonkerat, <0 kerat",
           sig & direction_nonkerat == "nonkerat_group (mean) < 0" &
             direction_kerat     == "kerat_group (mean) < 0" ~ "<0 in nonkerat and kerat",
           sig & direction_nonkerat == "nonkerat_group (mean) < 0" &
             direction_kerat     == "kerat_group (mean) > 0" ~ "<0 nonkerat, >0 kerat"
         )
  )

# Fig 1e
for (group in group_list) {
  df <- in_df %>% filter(Group == group)
  g <- ggplot(df, aes(x = coefficient, y = -log10(p.adj))) +
    geom_point(position = "jitter", size = 4,
               aes(fill = col_direction),
               colour = "black", pch = 21, alpha = 0.8) +
    labs(x = "Coefficient (vs. kerat_group)",
         y = "-log10 q-value") +
    theme_cowplot() +
    theme(axis.title = element_text(size = 20),
          axis.text = element_text(size = 18),
          legend.position = "right",
          legend.title = element_text(size = 22),
          legend.text = element_text(size = 18)) +
    scale_fill_manual(name = "Mean adj CN", values = c(
      ">0 nonkerat, <0 kerat" = "#E31A1C",
      ">0 in nonkerat and kerat" = "#FB9A99",
      "<0 in nonkerat and kerat" = "#A6CEE3",
      "<0 nonkerat, >0 kerat" = "#1F78B4"
    ), drop = FALSE) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red3") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    ggrepel::geom_text_repel(aes(label = id.sig), size = 6, max.overlaps = Inf)
  print(g)
}

### add wGII as a covariate ###
all_res_vec <- vector()
for (group in group_list) {
  for (var in var_list) {
    tmp_df <- LUSC_region_arm_gene_cn
    tmp_df$var <- tmp_df[[var]]
    tmp_df$group <- tmp_df[[group]]
    tmp_df <- tmp_df %>% filter(!is.na(var), !is.na(group), !is.na(wGII))
    
    effect_model <- lme(var ~ group + wGII, random = ~1 | tumour_id_muttable_cruk,
                        data = tmp_df, method = "ML")
    res_table <- as.data.frame(coef(summary(effect_model)), stringsAsFactors = FALSE)
    res_vec <- c(group, var,
                 res_table$Value[2],  # group coefficient
                 res_table$`p-value`[2])
    all_res_vec <- rbind(all_res_vec, res_vec)
  }
}
all_res_dt <- as.data.table(all_res_vec)
colnames(all_res_dt) <- c("Group","variable","coefficient","pval")
all_res_dt <- all_res_dt %>%
  mutate(coefficient = as.numeric(coefficient),
         pval        = as.numeric(pval)) %>%
  mutate(p.adj = p.adjust(pval, method = "BH"))
all_res_dt_SCNA_lme_wGII_adjusted <- all_res_dt


in_df2 <- all_res_dt_SCNA_lme_wGII_adjusted %>%
  left_join(mean_kerat_df, by = "variable") %>%
  left_join(mean_nonkerat_df, by = "variable") %>%
  mutate(id = variable,
         sig = p.adj < 0.05,
         id.sig = ifelse(sig, id, ""),
         col_direction = case_when(
           sig & direction_nonkerat == "nonkerat_group (mean) > 0" &
             direction_kerat     == "kerat_group (mean) > 0"       ~ ">0 in nonkerat and kerat",
           sig & direction_nonkerat == "nonkerat_group (mean) > 0" &
             direction_kerat     == "kerat_group (mean) < 0"       ~ ">0 nonkerat, <0 kerat",
           sig & direction_nonkerat == "nonkerat_group (mean) < 0" &
             direction_kerat     == "kerat_group (mean) < 0"       ~ "<0 in nonkerat and kerat",
           sig & direction_nonkerat == "nonkerat_group (mean) < 0" &
             direction_kerat     == "kerat_group (mean) > 0"       ~ "<0 nonkerat, >0 kerat"
         )
  )

# Ext fig 2c
for (group in group_list) {
  df <- in_df2 %>% filter(Group == group)
  g2 <- ggplot(df, aes(x = coefficient, y = -log10(p.adj))) +
    geom_point(position = "jitter", size = 4,
               aes(fill = col_direction),
               colour = "black", pch = 21, alpha = 0.8) +
    labs(x = "Coefficient (vs. kerat_group)",
         y = "-log10 q-value") +
    theme_cowplot() +
    theme(axis.title = element_text(size = 20),
          axis.text = element_text(size = 18),
          legend.position = "right",
          legend.title = element_text(size = 22),
          legend.text = element_text(size = 18)) +
    scale_fill_manual(name = "Mean adj CN", values = c(
      ">0 nonkerat, <0 kerat" = "#E31A1C",
      ">0 in nonkerat and kerat" = "#FB9A99",
      "<0 in nonkerat and kerat" = "#A6CEE3",
      "<0 nonkerat, >0 kerat" = "#1F78B4"
    ), drop = FALSE) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red3") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    ggrepel::geom_text_repel(aes(label = id.sig), size = 6, max.overlaps = Inf)
  print(g2)
}
