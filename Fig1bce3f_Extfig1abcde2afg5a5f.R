## Script for Tx421 LUSC paper

#Fig 1b,c,e
#Fig 3f
#Ext fig 1a,b,c,d,e
#Ext fig 2a,f,g
#Ext fig 5a, 5f


library(tidyverse)
library(ggpubr)
library(cowplot)
library(data.table)

# === Data loading ===
entire_tumour_data <- readRDS("../data/all_tumour_and_regions_df/full_histology_df.20220818.RDS") 
entire_info <- entire_tumour_data[, c("Histology_per_tumour_id_muttable", "patient", "tumour_id_by.label", "wGII", "TMB", "clonal_tmb", "clonal_TMB", "subclonal_tmb", "subclonal_TMB", "FLOH", "num_clonal_gds", "num_subclonal_gds", "purity", "ploidy", "smoking_status", "pack_years")] 
neoag_info <- entire_tumour_data %>% dplyr::select(patient, tumour_id_by.label, contains("neoag"))
squamous_data <- read_csv("../df-alvin/500_entire_tumour_with_grading.csv")
all_regional_full <- readRDS("../data/all_tumour_and_regions_df/all_regional_df_20221109.RDS")

wgiidata <- read_csv("../data/wgii_20221110_patient_Evo_metrics_tx421.csv")
wgiitrim <- wgiidata[, c("ltx_tumourid", "frac_abberant_genom_subcl", "perc_subclonal")]

sbs_data <- readRDS("../data/20221110_TRACERx421_mutationSignature_weights.rds")
sbs_entiretumour <- sbs_data$signature_weights_perTumour

max_cp_df <- fread("../data/20221102_TRACERx421_tumourDivCorrectedtree.csv")
expansion_score <- max_cp_df[, c("tumour_id", "maxCCF_terminalonly_tum_meanccf", "min_simpson_tum_meanccf")]

years_smoking_data <- readRDS("../data/20221109_TRACERx421_all_patient_df.rds") 
smoking_cessation_data <- read_csv("../data/all_patient_df_20221109.csv")
smoking_cessation_data <- smoking_cessation_data[, c("tumour_id_muttable_cruk", "years_after_stop")]
clonal_signature <- read_tsv("../data/signature/clonal_signature_tab.txt")
clonal_signature <- clonal_signature %>% dplyr::select(ltxid, contains("SBS"))
clonal_sbscount <- read_tsv("../data/signature/clonal_signaturecount_tab.txt")
clonal_sbscount <-clonal_sbscount %>% dplyr::select(ltxid, contains("SBS"))

# === Data integration ===
squamous_only <- squamous_data %>% filter(!is.na(keratinization))
squamous_plot <- squamous_only %>%
  left_join(entire_info, by = c("ltx_tumour_id" = "tumour_id_by.label")) %>%
  left_join(wgiitrim, by = c("ltx_tumour_id" = "ltx_tumourid")) %>%
  left_join(sbs_entiretumour, by = c("cruk_tumourid" = "tumour_id")) %>%
  left_join(expansion_score, by = c("cruk_tumourid" = "tumour_id")) %>%
  left_join(neoag_info, by = c("ltx_tumour_id" = "tumour_id_by.label")) %>%
  left_join(years_smoking_data, by = c("cruk_tumourid" = "tumour_id_muttable_cruk")) %>% 
  left_join(smoking_cessation_data, by = c("cruk_tumourid" = "tumour_id_muttable_cruk")) %>% 
  left_join(clonal_signature, by = c("ltx_clusterid" = "ltxid")) %>% 
  left_join(clonal_sbscount, by = c("ltx_clusterid" = "ltxid")) 

squamous_plot$keratinization <- factor(squamous_plot$keratinization, levels = c("have_keratin", "no_keratin"))

# === Define color mapping ===
color_list <- c(
  "have_keratin" = "#A6CEE3",
  "greencolour" = "#B2DF8A",
  "no_keratin" = "#FB9A99"
)

#Fig 1b
filtered_df <- squamous_plot %>% filter(!is.na(TMB))
print(table(filtered_df$keratinization))
ggboxplot(filtered_df, width = 0.5, lwd = 1, x = "keratinization", y = "TMB",
          fill = "keratinization", add = "jitter") +
  stat_compare_means(comparisons = list(c("no_keratin", "have_keratin")), label.y = 2000) +
  theme_cowplot() +
  scale_fill_manual(values = color_list, na.value = "#bababa", drop = FALSE) +
  theme(axis.title.y = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.title.x = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        legend.title = element_blank())


filtered_df <- squamous_plot %>% filter(!is.na(clonal_TMB))
print(table(filtered_df$keratinization))
ggboxplot(filtered_df, width = 0.5, lwd = 1, x = "keratinization", y = "clonal_TMB",
          fill = "keratinization", add = "jitter") +
  stat_compare_means(comparisons = list(c("no_keratin", "have_keratin")), label.y = 2000) +
  theme_cowplot() +
  scale_fill_manual(values = color_list, na.value = "#bababa", drop = FALSE) +
  theme(axis.title.y = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.title.x = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        legend.title = element_blank())

#Fig 1c
filtered_df <- squamous_plot %>% filter(!is.na(clonal_all_neoag))
print(table(filtered_df$keratinization))
ggboxplot(filtered_df, width = 0.5, lwd = 1, x = "keratinization", y = "clonal_all_neoag",
          fill = "keratinization", add = "jitter") +
  stat_compare_means(comparisons = list(c("no_keratin", "have_keratin")), label.y = 300) +
  theme_cowplot() +
  scale_fill_manual(values = color_list, na.value = "#bababa", drop = FALSE) +
  theme(axis.title.y = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.title.x = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        legend.title = element_blank())


#Fig 1e
filtered_df <- squamous_plot %>% filter(!is.na(FLOH))
print(table(filtered_df$keratinization))
ggboxplot(filtered_df, width = 0.5, lwd = 1, x = "keratinization", y = "FLOH",
          fill = "keratinization", add = "jitter") +
  stat_compare_means(comparisons = list(c("no_keratin", "have_keratin")), label.y = 0.7) +
  theme_cowplot() +
  scale_fill_manual(values = color_list, na.value = "#bababa", drop = FALSE) +
  theme(axis.title.y = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.title.x = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        legend.title = element_blank())


#Fig 1e truncal LOH count see separate script

#Ext fig 1a
filtered_df <- squamous_plot %>% filter(!is.na(subclonal_TMB))
print(table(filtered_df$keratinization))
ggboxplot(filtered_df, width = 0.5, lwd = 1, x = "keratinization", y = "subclonal_TMB",
          fill = "keratinization", add = "jitter") +
  stat_compare_means(comparisons = list(c("no_keratin", "have_keratin")), label.y = 1500) +
  theme_cowplot() +
  scale_fill_manual(values = color_list, na.value = "#bababa", drop = FALSE) +
  theme(axis.title.y = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.title.x = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        legend.title = element_blank())


#Ext fig 1b
filtered_df <- squamous_plot %>% filter(!is.na(years_smoking))
print(table(filtered_df$keratinization))
ggboxplot(filtered_df, width = 0.5, lwd = 1, x = "keratinization", y = "years_smoking",
          fill = "keratinization", add = "jitter") +
  stat_compare_means(comparisons = list(c("no_keratin", "have_keratin")), label.y = 110) +
  theme_cowplot() +
  scale_fill_manual(values = color_list, na.value = "#bababa", drop = FALSE) +
  theme(axis.title.y = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.title.x = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        legend.title = element_blank())

filtered_df <- squamous_plot %>% filter(!is.na(pack_years))
print(table(filtered_df$keratinization))
ggboxplot(filtered_df, width = 0.5, lwd = 1, x = "keratinization", y = "pack_years",
          fill = "keratinization", add = "jitter") +
  stat_compare_means(comparisons = list(c("no_keratin", "have_keratin")), label.y = 110) +
  theme_cowplot() +
  scale_fill_manual(values = color_list, na.value = "#bababa", drop = FALSE) +
  theme(axis.title.y = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.title.x = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        legend.title = element_blank())

filtered_df <- squamous_plot %>% filter(!is.na(years_after_stop))
print(table(filtered_df$keratinization))
ggboxplot(filtered_df, width = 0.5, lwd = 1, x = "keratinization", y = "years_after_stop",
          fill = "keratinization", add = "jitter") +
  stat_compare_means(comparisons = list(c("no_keratin", "have_keratin")), label.y = 110) +
  theme_cowplot() +
  scale_fill_manual(values = color_list, na.value = "#bababa", drop = FALSE) +
  theme(axis.title.y = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.title.x = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        legend.title = element_blank())

#Ext fig 1c
filtered_df <- squamous_plot %>% filter(!is.na(purity))
print(table(filtered_df$keratinization))
ggboxplot(filtered_df, width = 0.5, lwd = 1, x = "keratinization", y = "purity",
          fill = "keratinization", add = "jitter") +
  stat_compare_means(comparisons = list(c("no_keratin", "have_keratin")), label.y = 0.8) +
  theme_cowplot() +
  scale_fill_manual(values = color_list, na.value = "#bababa", drop = FALSE) +
  theme(axis.title.y = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.title.x = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        legend.title = element_blank())


#Ext fig 1d
filtered_df <- squamous_plot %>% filter(!is.na(clonalSBS4count))
print(table(filtered_df$keratinization))
ggboxplot(filtered_df, width = 0.5, lwd = 1, x = "keratinization", y = "clonalSBS4count",
          fill = "keratinization", add = "jitter") +
  stat_compare_means(comparisons = list(c("no_keratin", "have_keratin")), label.y = 500) +
  theme_cowplot() +
  scale_fill_manual(values = color_list, na.value = "#bababa", drop = FALSE) +
  theme(axis.title.y = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.title.x = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        legend.title = element_blank())

filtered_df <- squamous_plot %>% filter(!is.na(clonalSBS92count))
print(table(filtered_df$keratinization))
ggboxplot(filtered_df, width = 0.5, lwd = 1, x = "keratinization", y = "clonalSBS92count",
          fill = "keratinization", add = "jitter") +
  stat_compare_means(comparisons = list(c("no_keratin", "have_keratin")), label.y = 400) +
  theme_cowplot() +
  scale_fill_manual(values = color_list, na.value = "#bababa", drop = FALSE) +
  theme(axis.title.y = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.title.x = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        legend.title = element_blank())

filtered_df <- squamous_plot %>% filter(!is.na(clonalSBS4))
print(table(filtered_df$keratinization))
ggboxplot(filtered_df, width = 0.5, lwd = 1, x = "keratinization", y = "clonalSBS4",
          fill = "keratinization", add = "jitter") +
  stat_compare_means(comparisons = list(c("no_keratin", "have_keratin")), label.y = 0.61) +
  theme_cowplot() +
  scale_fill_manual(values = color_list, na.value = "#bababa", drop = FALSE) +
  theme(axis.title.y = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.title.x = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        legend.title = element_blank())

filtered_df <- squamous_plot %>% filter(!is.na(clonalSBS92))
print(table(filtered_df$keratinization))
ggboxplot(filtered_df, width = 0.5, lwd = 1, x = "keratinization", y = "clonalSBS92",
          fill = "keratinization", add = "jitter") +
  stat_compare_means(comparisons = list(c("no_keratin", "have_keratin")), label.y = 0.61) +
  theme_cowplot() +
  scale_fill_manual(values = color_list, na.value = "#bababa", drop = FALSE) +
  theme(axis.title.y = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.title.x = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        legend.title = element_blank())

#Ext fig 2f
filtered_df <- squamous_plot %>% filter(!is.na(maxCCF_terminalonly_tum_meanccf))
print(table(filtered_df$keratinization))
ggboxplot(filtered_df, width = 0.5, lwd = 1, x = "keratinization", y = "maxCCF_terminalonly_tum_meanccf",
          fill = "keratinization", add = "jitter") +
  stat_compare_means(comparisons = list(c("no_keratin", "have_keratin")), label.y = 1.2) +
  theme_cowplot() +
  scale_fill_manual(values = color_list, na.value = "#bababa", drop = FALSE) +
  theme(axis.title.y = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.title.x = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        legend.title = element_blank())


#Ext fig 2g
filtered_df <- squamous_plot %>% filter(!is.na(frac_abberant_genom_subcl))
print(table(filtered_df$keratinization))
ggboxplot(filtered_df, width = 0.5, lwd = 1, x = "keratinization", y = "frac_abberant_genom_subcl",
          fill = "keratinization", add = "jitter") +
  stat_compare_means(comparisons = list(c("no_keratin", "have_keratin")), label.y = 1) +
  theme_cowplot() +
  scale_fill_manual(values = color_list, na.value = "#bababa", drop = FALSE) +
  theme(axis.title.y = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.title.x = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        legend.title = element_blank())


#Ext fig 1e
df <- squamous_plot
df <- df %>%
  filter(!is.na(clonal_TMB)) %>% 
  group_by(keratinization) %>%
  arrange(clonal_TMB) %>%
  mutate(Order_within_subtype = row_number()) %>%
  ungroup()

df2 <- df %>% mutate(final_order = case_when(keratinization == "have_keratin" ~ Order_within_subtype,
                                             keratinization == "no_keratin" ~ Order_within_subtype + length(df$keratinization[df$keratinization == "have_keratin"])))

ggplot(df2) +
  geom_col(aes(x = final_order, y = clonal_TMB, fill = keratinization), color = "black", linewidth = 0.01) +
  scale_fill_manual(values = color_list, na.value = '#bababa', drop = F) +
  theme_cowplot() +
  panel_border(color = "black", size = 1) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.x = element_blank(),
    legend.position = "top"
  ) 
# ggsave("../results/new_results/final_clonal_tmb_plot.pdf", width = 20, height = 8)


ggplot(df2) +
  geom_col(aes(x = final_order, y = pack_years, fill = keratinization), color = "black", linewidth = 0.01) +
  scale_fill_manual(values = color_list, na.value = '#bababa', drop = F) +
  theme_cowplot() +
  panel_border(color = "black", size = 1) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.x = element_blank(),
    legend.position = "top"
  ) 
# ggsave("../results/new_results/final_packyears-clonal_tmb_plot.pdf", width = 20, height = 3)

ggplot(df2) +
  geom_col(aes(x = final_order, y = SBS4, fill = keratinization), color = "black", linewidth = 0.01) +
  scale_fill_manual(values = color_list, na.value = '#bababa', drop = F) +
  theme_cowplot() +
  panel_border(color = "black", size = 1) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.x = element_blank(),
    legend.position = "top"
  ) 
# ggsave("../results/new_results/final_SBS4-clonal_tmb_plot.pdf", width = 20, height = 3)

ggplot(df2) +
  geom_col(aes(x = final_order, y = SBS92, fill = keratinization), color = "black", linewidth = 0.01) +
  scale_fill_manual(values = color_list, na.value = '#bababa', drop = F) +
  theme_cowplot() +
  panel_border(color = "black", size = 1) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.x = element_blank(),
    legend.position = "top"
  ) 
# ggsave("../results/new_results/final_SBS92-clonal_tmb_plot.pdf", width = 20, height = 3)


#Ext fig 1f
#two groups of LUSC
keratinized_lusc <- squamous_plot %>% filter(keratinization == "have_keratin")
nonkeratinized_lusc <- squamous_plot %>% filter(keratinization == "no_keratin")
ggscatter(keratinized_lusc, x = "pack_years", y = "clonalSBS4", 
          add = "reg.line", conf.int = TRUE, 
          color = '#A6CEE3',
          cor.coef = TRUE, cor.method = "pearson")

ggscatter(nonkeratinized_lusc, x = "pack_years", y = "clonalSBS4", 
          add = "reg.line", conf.int = TRUE, 
          color = '#FB9A99',
          cor.coef = TRUE, cor.method = "pearson")

ggscatter(keratinized_lusc, x = "pack_years", y = "clonalSBS92", 
          add = "reg.line", conf.int = TRUE, 
          color = '#A6CEE3',
          cor.coef = TRUE, cor.method = "pearson")

ggscatter(nonkeratinized_lusc, x = "pack_years", y = "clonalSBS92", 
          add = "reg.line", conf.int = TRUE, 
          color = '#FB9A99',
          cor.coef = TRUE, cor.method = "pearson")



#Ext Fig 2a

SBS92detection_df <- read_tsv("../data/SBS92detection_classification_perTumour.txt")

frac_SBS92detection_valueGroups <- function(group_df, gap = 5) {
  group_df <- group_df[!is.na(group_df[,2]),]
  
  min_value <- ceiling(min(group_df[,2]) / gap) * gap
  max_value <- ceiling(max(group_df[,2]) / gap) * gap
  seq_vector <- seq(min_value, max_value, by = gap)
  
  frac <- lapply(seq_vector, function(x) {
    sub <- group_df[group_df[,2] <= x,]
    counts <- sub %>%
      group_by(sub[,1]) %>%
      summarise(n = n(), .groups = 'drop') %>%
      mutate(n_group = nrow(sub))
    data.frame(max_value = x, counts)
  }) %>% bind_rows()
  
  colnames(frac) <- c('max_groupValue', 'SBS92detected', 'n', 'n_group')
  frac$frac_group <- frac$n / frac$n_group
  
  temp <- group_df %>%
    group_by(group_df[,1]) %>%
    summarise(n_SBS92group = n(), .groups = 'drop')
  colnames(temp)[1] <- 'group'
  
  frac <- left_join(frac, temp, by = c('SBS92detected' = 'group'))
  frac$frac_SBS92group <- frac$n / frac$n_SBS92group
  
  return(frac)
}


# Table version for combined plot
table_SBS92_vs_smoking <- function(df) {
  smoke_vars <- df %>%
    select(cruk_tumourid, years_smoking, cigs_perday, age, histology_multi_full, smoking_status_merged)
  
  smoke_vars <- SBS92detection_df %>%
    left_join(smoke_vars, by = c('tumour' = 'cruk_tumourid')) %>%
    filter(smoking_status_merged != 'Never Smoked') %>%
    mutate(across(c(overall_SBS92detetced, clonal_SBS92detected, subclonal_SBS92detected),
                  ~factor(.x, levels = c('YES', 'MAYBE', 'NO'))))
  
  colnames(smoke_vars) <- c('tumour', 'Patient', 'clonal_SBS92detected', 'subclonal_SBS92detected', 'total_SBS92detected',
                            'overall_SBS92detetced', 'smokeYrs', 'cigsPerDay', 'age', 'histology', 'smoking_status_merged')
  
  filtered <- smoke_vars %>%
    filter(overall_SBS92detetced != 'MAYBE') %>%
    select(overall_SBS92detetced, smokeYrs)
  
  frac_data <- frac_SBS92detection_valueGroups(filtered, gap = 5)
  
  add_data <- frac_data[1,,drop = FALSE]
  add_data$SBS92detected <- 'YES'
  add_data$frac_group <- 0
  
  bind_rows(frac_data, add_data)
}

# Combined grouped plot
nonK <- table_SBS92_vs_smoking(nonkeratinized_lusc) %>% mutate(keratinization = "no_keratin")
K    <- table_SBS92_vs_smoking(keratinized_lusc) %>% mutate(keratinization = "have_keratin")
combined <- bind_rows(nonK, K) %>% filter(SBS92detected == "YES")

ggplot(combined, aes(x = factor(max_groupValue), y = frac_group, fill = keratinization)) +
  geom_bar(stat = 'identity', position = "dodge") +
  geom_smooth(method = 'loess', aes(group = keratinization, colour = keratinization), linetype = 'dashed', se = F, size = 0.8) +
  scale_fill_manual(name = 'SBS92', values = c('have_keratin' = '#A6CEE3', 'MAYBE' = '#fdae61', 'no_keratin' = '#FB9A99'),
                    labels = c('NO' = 'not detected', 'MAYBE' = 'low confidence\ndetected', 'YES' = 'high confidence\ndetected')) +
  scale_colour_manual(name = 'SBS92', values = c('have_keratin' = '#A6CEE3', 'MAYBE' = '#fdae61', 'no_keratin' = '#FB9A99'),
                      labels = c('NO' = 'not detected', 'MAYBE' = 'low confidence\ndetected', 'YES' = 'high confidence\ndetected')) +
  scale_y_continuous( breaks = seq(0, 1, 0.2), labels = seq(0,100, 20)) +
  xlab('upper limit smoking years') + ylab('% tumours') +
  ggtitle('Probability of detecting SBS92 by smoking duration') +
  theme_bw() + 
  theme(legend.position = 'none', 
        plot.title = element_text(hjust = 0.5, size = 11.5), 
        panel.grid = element_blank(),
        text = element_text( size = 16.5 ))


########################
#TIL data
#Ext fig 5f
TIL_entire <- read_csv("../data/20230408-TIL_STAS_ctDNA/Pathology.TIL.scores_diagnostic.csv")
squamous_plot <- squamous_plot %>% mutate(SampleID = substr(ltx_tumour_id, 1, 6)) 

left_join(squamous_plot, TIL_entire , by = c("SampleID")) -> squamous_TIL 

filtered_df <- squamous_TIL %>% filter(!is.na(final.TIL))
print(table(filtered_df$keratinization))
ggboxplot(filtered_df, width = 0.5, lwd = 1, x = "keratinization", y = "final.TIL",
          fill = "keratinization", add = "jitter") +
  stat_compare_means(comparisons = list(c("no_keratin", "have_keratin")), label.y = 125) +
  theme_cowplot() +
  scale_fill_manual(values = color_list, na.value = "#bababa", drop = FALSE) +
  theme(axis.title.y = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.title.x = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        legend.title = element_blank())
# ggsave("../results/TIL_keratinization_final.pdf")

#Fig 3f
filtered_df <- squamous_TIL %>% filter(!is.na(final.TIL)) 
filtered_df$tumour_budding <- factor(filtered_df$tumour_budding, levels = c("no_budding", "have_budding"))

print(table(filtered_df$tumour_budding))
color_list <- c(
  "no_budding" = "#A6CEE3",
  "greencolour" = "#B2DF8A",
  "have_budding" = "#FB9A99"
)
ggboxplot(filtered_df, width = 0.5, lwd = 1, x = "tumour_budding", y = "final.TIL",
          fill = "tumour_budding", add = "jitter") +
  stat_compare_means(comparisons = list(c("no_budding", "have_budding")), label.y = 125) +
  theme_cowplot() +
  scale_fill_manual(values = color_list, na.value = "#bababa", drop = FALSE) +
  theme(axis.title.y = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.title.x = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        legend.title = element_blank())
# ggsave("../results/TIL_budding_final.pdf")

#Ext fig 5a
filtered_df <- squamous_TIL %>% filter(!is.na(final.TIL)) 
filtered_df$single_cell_invasion <- factor(filtered_df$single_cell_invasion, levels = c("no_sci", "have_sci"))

print(table(filtered_df$single_cell_invasion))
color_list <- c(
  "no_sci" = "#A6CEE3",
  "greencolour" = "#B2DF8A",
  "have_sci" = "#FB9A99"
)
ggboxplot(filtered_df, width = 0.5, lwd = 1, x = "single_cell_invasion", y = "final.TIL",
          fill = "single_cell_invasion", add = "jitter") +
  stat_compare_means(comparisons = list(c("no_sci", "have_sci")), label.y = 125) +
  theme_cowplot() +
  scale_fill_manual(values = color_list, na.value = "#bababa", drop = FALSE) +
  theme(axis.title.y = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.title.x = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        legend.position = "none",
        legend.title = element_blank())
# ggsave("../results/TIL_sci_final.pdf")
