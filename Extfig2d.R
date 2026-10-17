## Script for Tx421 LUSC paper
#Fig 1g
#Ext fig 2c

library(tidyverse)
library(data.table)
library(reshape2)
library(ggpubr)
library(cowplot)
library(ggfittext)
read_csv("../data/lohhla.csv") -> loh_df

lohhla_results <-  loh_df %>% filter(!is.na(loh)) %>% as.data.table()

lohhla_results[, patient_tumour := paste0(patient, "_", tumour_id)]

lohhla_results[loh == TRUE & cn1_binned < cn2_binned, allele1_lost := TRUE]
lohhla_results[is.na(allele1_lost), allele1_lost := FALSE]

lohhla_results[loh == TRUE & cn1_binned > cn2_binned, allele2_lost := TRUE]
lohhla_results[is.na(allele2_lost), allele2_lost := FALSE]

loh_tab <- lohhla_results[,c("patient_tumour", "sample_name", "gene",
                             "allele1", "allele2", 
                             "allele1_lost", "allele2_lost", "hist")]

loh_tab[,n_regions := .N, by = c("patient_tumour", "gene")]
loh_tab[,n_allele1_lost := sum(allele1_lost), by = c("patient_tumour", "gene")]
loh_tab[,n_allele2_lost := sum(allele2_lost), by = c("patient_tumour", "gene")]

loh_tab[n_allele1_lost > 0 & n_allele1_lost == n_regions, loh_cat := "clonal"]
loh_tab[n_allele2_lost > 0 & n_allele2_lost == n_regions, loh_cat := "clonal"]

loh_tab[n_allele1_lost > 0 & n_allele1_lost < n_regions, loh_cat := "subclonal"]
loh_tab[n_allele2_lost > 0 & n_allele2_lost < n_regions, loh_cat := "subclonal"]

loh_tab[n_allele1_lost == 0 & n_allele2_lost == 0, loh_cat := "none"]
loh_tab[,.N,by = "loh_cat"]

# create only three rows for each tumour
loh_tab_per_tumour <- loh_tab %>% dplyr::select(-sample_name, -allele1_lost, -allele2_lost) %>% 
  filter(!duplicated(.)) 

read_csv("../df-alvin/LOH_summary_table.csv") -> LOH_analyse

LOH_analyse %>% group_by(ltx_tumourid, hist, loh_cat) %>% summarize(number = n()) -> tmp1
tmp1 %>% dcast(., ltx_tumourid + hist ~ loh_cat) %>% filter(!is.na(ltx_tumourid)) -> tmp2

#LUSC entire tumour data
read_csv("../df-alvin/500_entire_tumour_with_grading_with_wGIITMB_info.csv") -> squamous_entire_data
left_join(tmp2, squamous_entire_data[, c("tumour_id_by.label", "differentiation", "keratinization")], by = c("ltx_tumourid" = "tumour_id_by.label")) -> tmp3

# create any_clonal LOH, any_LOH columns
tmp3 %>% mutate(any_clonalLOH = case_when(clonal >= 1 ~ "ClonalLOH", 
                                          TRUE ~ "No_ClonalLOH")) %>% 
  mutate(any_LOH = case_when(clonal >= 1 | subclonal >= 1 ~ "AnyLOH", 
                             TRUE ~ "NoLOH")) %>% 
  mutate(subclonal_LOH_only = case_when(subclonal >= 1 & (is.na(clonal)|clonal == 0) ~ "subclonalLOH", 
                                        TRUE ~ "")) %>% 
  mutate(clonal_in_all_HLA_allele = case_when(clonal >= 1 & is.na(subclonal) & is.na(none) ~ "LOH_all_alleles", 
                                              TRUE ~ "LOH_not_all_alleles")) %>% 
  mutate(clonal_LOH_vs_noLOH = case_when(any_clonalLOH == "ClonalLOH" ~ "ClonalLOH",
                                          any_LOH == "NoLOH" ~ "NoLOH", 
                                          )) %>% 
  mutate(subclonalLOH_vs_noLOH = case_when(subclonal_LOH_only == "subclonalLOH" ~ "subclonalLOH",
                                           any_LOH == "NoLOH" ~ "NoLOH", 
  )) %>% dplyr::filter(!is.na(keratinization)) ->
  LOH_plot


#plot stacked bar chart for keratinization
#clonal LOH
table(LOH_plot$keratinization,LOH_plot$clonal_LOH_vs_noLOH) -> table1 

df_gathered1 <- data.frame(var1 = c(rownames(table1)[1], rownames(table1)[1], rownames(table1)[2], rownames(table1)[2] ),
                           var2 = c(colnames(table1)[1],colnames(table1)[2], colnames(table1)[1], colnames(table1)[2]),
                           count = c(table1[1], table1[3], table1[2], table1[4]))
fisher_p <- fisher.test(matrix(c(table1[1], table1[2], table1[3], table1[4]), nrow = 2))

color_list <-  c("have_keratin" = "#A6CEE3",
                 "greencolour" = "#B2DF8A",
                 "no_keratin"= "#FB9A99")
g1 <- ggplot(data=df_gathered1, aes(x=var2, y=count, fill=var1, label = count)) +
  geom_bar(stat="identity", position = "fill", color = "black") +
  geom_bar_text( position = "fill", reflow = TRUE) + 
  scale_fill_manual(values = color_list, na.value = '#bababa', drop = F)+
  theme_cowplot() +
  ylab("proportion") +
  # labs(x = var2, y = "proportion", fill = var1)   +
  ggtitle(paste0( " Fisher-exact test p-value: ", format(round(as.numeric(fisher_p$p.value), 4), nsmall=4))) +
  theme(axis.title.x = element_blank(),  
        axis.title.y = element_text(size=24),
        axis.text.x = element_text(size=20, angle = 45,  hjust=1), 
        axis.text.y = element_text(size=16),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size=18))
print(g1)


#Any LOH
table(LOH_plot$keratinization,LOH_plot$any_LOH) -> table1

df_gathered1 <- data.frame(var1 = c(rownames(table1)[1], rownames(table1)[1], rownames(table1)[2], rownames(table1)[2] ),
                           var2 = c(colnames(table1)[1],colnames(table1)[2], colnames(table1)[1], colnames(table1)[2]),
                           count = c(table1[1], table1[3], table1[2], table1[4]))
df_gathered1$var2 <- factor( df_gathered1$var2 , levels=c("AnyLOH", "NoLOH"))


fisher_p <- fisher.test(matrix(c(table1[1], table1[2], table1[3], table1[4]), nrow = 2))

color_list <-  c("have_keratin" = "#A6CEE3",
                 "greencolour" = "#B2DF8A",
                 "no_keratin"= "#FB9A99")
g1 <- ggplot(data=df_gathered1, aes(x=var2, y=count, fill=var1, label = count)) +
  geom_bar(stat="identity", position = "fill", color = "black") +
  geom_bar_text( position = "fill", reflow = TRUE) + 
  scale_fill_manual(values = color_list, na.value = '#bababa', drop = F)+
  theme_cowplot() +
  ylab("proportion") +
  # labs(x = var2, y = "proportion", fill = var1)   +
  ggtitle(paste0( " Fisher-exact test p-value: ", format(round(as.numeric(fisher_p$p.value), 4), nsmall=4))) +
  theme(axis.title.x = element_blank(),  
        axis.title.y = element_text(size=24),
        axis.text.x = element_text(size=20, angle = 45,  hjust=1), 
        axis.text.y = element_text(size=16),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size=18))
print(g1)
