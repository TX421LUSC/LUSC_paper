#Script for Tx421 LUSC paper
#Fig 3e, 3g, 3h, 
#Ext 4i, 5b, 5b, 5c

library(tidyverse)
library(ggpubr)
library(reshape2)
# library(data.table)
library(ggfittext)
library(cowplot)

LUSC_region_df_full <- read_csv("../df-alvin/All_regions_500_regions-SCC-reviewed.csv")
GSVA_data <- read_csv("../df-alvin/gsva_GOBP_keratin_df.csv")
EMT_data <- readRDS("../df-alvin/gsva_GO_df.RDS") %>% dplyr::select(short_sample, contains("MESENCHYMAL_TRAN"))
colnames(EMT_data)[1] <-  "sample_name"
immune_data <- read_csv("../data/Tx421_MCPcounter.regional.csv")
danaher_data <- read_csv("../data/danaher_scores.csv")
danaher_data <- danaher_data %>% dplyr::select(ltx_region_id, score, cell_type) %>% 
                  dcast(formula = ltx_region_id ~ cell_type, value.var = "score")
colnames(danaher_data)[-1] <- danaher_data %>% dplyr::select(-ltx_region_id) %>% colnames() %>% str_c("_danaher") 
textract_data <- read_csv("../data/Tx421_TCRA_20211217.csv")
textract_data <- textract_data %>% dplyr::select(ltx_region_id, TCRA.tcell.fraction.adj)
LUSC_region_df <- LUSC_region_df_full  %>%
  filter(!is.na(Histology_per_tumour_id_muttable)) %>% 
  filter(Histology_per_tumour_id_muttable == "Squamous cell carcinoma") %>% 
  filter(included_500_regions == 1)

scc_regionalstat <- LUSC_region_df %>% dplyr::select(tumour_id_muttable_cruk, Region_ID_cruk, ltx_region_id, regional_pattern.new, wGII, all_TMB, clonal_TMB, subclonal_TMB, FLOH, purity, ploidy, keratinization_01, keratinization, contains("budding"), contains("invasion"), contains("sci")) %>% 
                    left_join(GSVA_data, by = c("Region_ID_cruk" = "sample_name")) %>% 
                    left_join(EMT_data, by = c("Region_ID_cruk" = "sample_name"))   %>% 
                    left_join(immune_data, by = c("ltx_region_id" = "sample_name")) %>% 
                    left_join(danaher_data, by = "ltx_region_id") %>% 
                    left_join(textract_data, by = "ltx_region_id") 



plot_scatter <- function(paramX, paramY){
   in_df <- scc_regionalstat
  in_df$paramX <- in_df[[paramX]]
  in_df$paramY <- in_df[[paramY]]
  in_df <- in_df %>% filter(!is.na(paramX) & !is.na(paramY))
  g1 <- ggplot(in_df, aes(x = paramX, y = paramY)) + 
  geom_point(aes(colour = '#67001f'), alpha = 0.7) +
  geom_smooth(method = 'lm', colour = 'black', linetype = 'dashed') +
  stat_cor(method = 'spearman') +
   xlab(paramX) + ylab(paramY) +
  theme_bw() + theme(legend.position = 'top', text = element_text(size = 20))
  
  print(g1)
}

#Ext fig 4i
plot_scatter("GOBP_KERATINOCYTE_DIFFERENTIATION", "GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION" )

#Fig 3e
plot_scatter("max_tumour_budding", "GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION" )

#for plotting purity
scc_sync <- scc_regionalstat %>% 
  filter(!is.na(max_tumour_budding) & !is.na(GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION))

# 2. Plot 1: Budding vs EMT (The reference plot)
p1 <- ggplot(scc_sync, aes(x = max_tumour_budding, y = GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION)) + 
  geom_point(colour = '#67001f', alpha = 0.7) +
  geom_smooth(method = 'lm', colour = 'black', linetype = 'dashed') +
  stat_cor(method = 'spearman') +
  theme_bw() + theme(legend.position = 'top', text = element_text(size = 20))

# 3. Plot 2: Budding vs Purity (Synced to Plot 1)
# We filter (!is.na(purity)) here so the trend line doesn't break if purity is missing for a few
p2 <- ggplot(scc_sync %>% filter(!is.na(purity)), aes(x = max_tumour_budding, y = purity)) + 
  geom_point(colour = '#67001f', alpha = 0.7) +
  #geom_smooth(method = 'lm', colour = 'black', linetype = 'dashed') +
  #stat_cor(method = 'spearman') +
  theme_bw() + theme(legend.position = 'top', text = element_text(size = 20))

library(patchwork)
combined_plot <- p1 / p2 + plot_layout(heights = c(6, 1))

# Display the result
print(combined_plot)
ggsave("../results/new_buddingvspurity.pdf", height = 8, width = 9)  

#Fig 3g
plot_scatter("max_tumour_budding", "MCP_Fibroblasts")
#for plotting purity
scc_sync <- scc_regionalstat %>% 
  filter(!is.na(max_tumour_budding) & !is.na(MCP_Fibroblasts))

# 2. Plot 1: Budding vs EMT (The reference plot)
p1 <- ggplot(scc_sync, aes(x = max_tumour_budding, y = MCP_Fibroblasts)) + 
  geom_point(colour = '#67001f', alpha = 0.7) +
  geom_smooth(method = 'lm', colour = 'black', linetype = 'dashed') +
  stat_cor(method = 'spearman') +
  theme_bw() + theme(legend.position = 'top', text = element_text(size = 20))
print(p1)
# 3. Plot 2: Budding vs Purity (Synced to Plot 1)
# We filter (!is.na(purity)) here so the trend line doesn't break if purity is missing for a few
p2 <- ggplot(scc_sync %>% filter(!is.na(purity)), aes(x = max_tumour_budding, y = purity)) + 
  geom_point(colour = '#67001f', alpha = 0.7) +
  #geom_smooth(method = 'lm', colour = 'black', linetype = 'dashed') +
  #stat_cor(method = 'spearman') +
  theme_bw() + theme(legend.position = 'top', text = element_text(size = 20))

library(patchwork)
combined_plot <- p1 / p2 + plot_layout(heights = c(6, 1))
print(combined_plot)
ggsave("../results/new_buddingfibroblastvspurity.pdf", height = 8, width = 9)  

#Ext Fig 5c
plot_scatter("num_single_cell_invasion", "MCP_Fibroblasts")

#Ext Fig 5b
plot_scatter("max_tumour_budding", "MCP_T cells")

#Ext Fig 5b
plot_scatter("num_single_cell_invasion", "MCP_T cells")


# identifying correlation table  - keratinization vs B Cell , T cell, NK cell Spearman
#Fig 3h
mcp_list <- scc_regionalstat %>% dplyr::select(contains("MCP")) %>% colnames
danaher_list <- scc_regionalstat %>% dplyr::select(contains("danaher")) %>% colnames
x_var_list <- c( "GOBP_KERATINOCYTE_DIFFERENTIATION")
y_var_list <- c( mcp_list, danaher_list, "TCRA.tcell.fraction.adj")

in_df <-  scc_regionalstat 

all_res_vec <- vector()
for(x_var in x_var_list){
  for(y_var in y_var_list){
    tmp_df <- in_df
    tmp_df$x_var <- tmp_df[[x_var]]
    tmp_df$y_var <- tmp_df[[y_var]]
    tmp_df <- tmp_df %>%
      filter(!is.na(x_var),!is.na(y_var)) 
    
    res_cor.test <- cor.test(tmp_df$x_var,tmp_df$y_var, method = "spearman" )
    res_coef <- res_cor.test$estimate
    res_pval <- res_cor.test$p.value
    
    res_vec <- c(x_var,y_var, res_coef, res_pval)
    all_res_vec <- rbind(all_res_vec, res_vec)
  }
}

all_res_dt <- data.table(all_res_vec) 
colnames(all_res_dt)   <- c("grade_proportion","variable", "coef", "pval")
all_res_dt <- all_res_dt %>%
  mutate(coef = as.numeric(coef),
         pval = as.numeric(pval)) %>%
  #mutate(p.adj = p.adjust(pval, method = "bonferroni")) %>%
  mutate(p.adj = p.adjust(pval, method = "BH")) %>%
  mutate(sign_pval = case_when(pval < 0.0001 ~ "****",
                               pval < 0.001 ~ "***",
                               pval < 0.01 ~ "**",
                               pval < 0.05 ~ "*"),
         sign_p.adj = case_when(p.adj < 0.0001 ~ "****",
                                p.adj < 0.001 ~ "***",
                                p.adj < 0.01 ~ "**",
                                p.adj < 0.05 ~ "*")) %>%
  mutate(variable = factor(gsub("_"," ",variable), levels = rev(gsub("_"," ",y_var_list)) ))

all_res_dt <- all_res_dt %>%  mutate(cell_group = case_when(variable == "MCP T cells" ~ "T cell",
                                                            variable == "t cells danaher" ~ "T cell", 
                                                            variable == "MCP B lineage" ~ "B cell",
                                                            variable == "b cells danaher" ~ "B cell", 
                                                            variable == "MCP CD8 T cells" ~ "CD8 cell",
                                                            variable == "cd8 cells danaher" ~ "CD8 cell",
                                                            variable == "b cells danaher" ~ "B cell",
                                                            variable == "MCP NK cells" ~ "NK cell",
                                                            variable == "nk cells danaher" ~ "NK cell",
                                                            variable == "MCP Cytotoxic lymphocytes" ~ "Cytotoxic cell",
                                                            variable == "cytotoxic cells danaher" ~ "Cytotoxic cell",
                                                            variable == "MCP Neutrophils" ~ "Neutrophil",
                                                            variable == "neutrophils danaher" ~ "Neutrophil",
                                                            variable == "MCP Monocytic lineage" ~ "Macrophage",
                                                            variable == "macrophages danaher" ~ "Macrophage",
                                                            variable == "TCRA.tcell.fraction.adj" ~ "T cell"
                                                            ))  %>% 
                              mutate(trial_group = case_when(variable == "MCP T cells" ~ "MCP",
                                                            variable == "t cells danaher" ~ "Danaher", 
                                                            variable == "MCP B lineage" ~ "MCP",
                                                            variable == "b cells danaher" ~ "Danaher",
                                                            variable == "MCP CD8 T cells" ~ "MCP",
                                                            variable == "cd8 cells danaher" ~ "Danaher",
                                                            variable == "b cells danaher" ~ "Danaher",
                                                            variable == "MCP NK cells" ~ "MCP",
                                                            variable == "nk cells danaher" ~ "Danaher",
                                                            variable == "MCP Cytotoxic lymphocytes" ~ "MCP",
                                                            variable == "cytotoxic cells danaher" ~ "Danaher",
                                                            variable == "MCP Neutrophils" ~ "MCP",
                                                            variable == "neutrophils danaher" ~ "Danaher",
                                                            variable == "MCP Monocytic lineage" ~ "MCP",
                                                            variable == "macrophages danaher" ~ "Danaher",
                                                            variable == "TCRA.tcell.fraction.adj" ~ "Tracerx")) %>% 
                              mutate(cell_group = factor(cell_group, levels = rev(c("T cell", "CD8 cell", "Cytotoxic cell", "B cell", "NK cell", "Macrophage", "Neutrophil"))))
       
# plot table 

g <- ggplot(all_res_dt %>% filter(!is.na(cell_group)), aes(x = trial_group, 
                            y = cell_group,
                            label = sign_p.adj)) +
  geom_tile(colour = 'black', aes(fill = coef)) +
  geom_fit_text(reflow = TRUE, min.size = 6 ) +
  scale_fill_gradient2(name = "Spearman's Rho", low = 'blue',mid = "white", high = 'red',limits = c(-1,1), na.value = 'gray') +
  theme_cowplot() + 
  theme(title = element_text(size=14),
        axis.title.x =element_blank(),  
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size=14), #,angle = 90,vjust = 0.5, hjust=1),  
        axis.text.y = element_text(size=14),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line=element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=10)) + 
  scale_x_discrete(position = "top")

print(g)  


#Ext fig 5e

# scc_regionalstat$GOBP_KERATINIZATION
x_var_list <- c( "keratinization_01")
y_var_list <- c( mcp_list, danaher_list, "TCRA.tcell.fraction.adj")

in_df <-  scc_regionalstat 

all_res_vec <- vector()
for(x_var in x_var_list){
  for(y_var in y_var_list){
    tmp_df <- in_df
    tmp_df$x_var <- tmp_df[[x_var]]
    tmp_df$y_var <- tmp_df[[y_var]]
    tmp_df <- tmp_df %>%
      filter(!is.na(x_var),!is.na(y_var)) 
    
    res_cor.test <- cor.test(tmp_df$x_var,tmp_df$y_var, method = "spearman" )
    res_coef <- res_cor.test$estimate
    res_pval <- res_cor.test$p.value
    
    res_vec <- c(x_var,y_var, res_coef, res_pval)
    all_res_vec <- rbind(all_res_vec, res_vec)
  }
}

all_res_dt <- data.table(all_res_vec) 
colnames(all_res_dt)   <- c("grade_proportion","variable", "coef", "pval")
all_res_dt <- all_res_dt %>%
  mutate(coef = as.numeric(coef),
         pval = as.numeric(pval)) %>%
  #mutate(p.adj = p.adjust(pval, method = "bonferroni")) %>%
  mutate(p.adj = p.adjust(pval, method = "BH")) %>%
  mutate(sign_pval = case_when(pval < 0.0001 ~ "****",
                               pval < 0.001 ~ "***",
                               pval < 0.01 ~ "**",
                               pval < 0.05 ~ "*"),
         sign_p.adj = case_when(p.adj < 0.0001 ~ "****",
                                p.adj < 0.001 ~ "***",
                                p.adj < 0.01 ~ "**",
                                p.adj < 0.05 ~ "*")) %>%
  mutate(variable = factor(gsub("_"," ",variable), levels = rev(gsub("_"," ",y_var_list)) ))

all_res_dt <- all_res_dt %>%  mutate(cell_group = case_when(variable == "MCP T cells" ~ "T cell",
                                                            variable == "t cells danaher" ~ "T cell", 
                                                            variable == "MCP B lineage" ~ "B cell",
                                                            variable == "b cells danaher" ~ "B cell", 
                                                            variable == "MCP CD8 T cells" ~ "CD8 cell",
                                                            variable == "cd8 cells danaher" ~ "CD8 cell",
                                                            variable == "b cells danaher" ~ "B cell",
                                                            variable == "MCP NK cells" ~ "NK cell",
                                                            variable == "nk cells danaher" ~ "NK cell",
                                                            variable == "MCP Cytotoxic lymphocytes" ~ "Cytotoxic cell",
                                                            variable == "cytotoxic cells danaher" ~ "Cytotoxic cell",
                                                            variable == "MCP Neutrophils" ~ "Neutrophil",
                                                            variable == "neutrophils danaher" ~ "Neutrophil",
                                                            variable == "MCP Monocytic lineage" ~ "Macrophage",
                                                            variable == "macrophages danaher" ~ "Macrophage",
                                                            variable == "TCRA.tcell.fraction.adj" ~ "T cell"
))  %>% 
  mutate(trial_group = case_when(variable == "MCP T cells" ~ "MCP",
                                 variable == "t cells danaher" ~ "Danaher", 
                                 variable == "MCP B lineage" ~ "MCP",
                                 variable == "b cells danaher" ~ "Danaher",
                                 variable == "MCP CD8 T cells" ~ "MCP",
                                 variable == "cd8 cells danaher" ~ "Danaher",
                                 variable == "b cells danaher" ~ "Danaher",
                                 variable == "MCP NK cells" ~ "MCP",
                                 variable == "nk cells danaher" ~ "Danaher",
                                 variable == "MCP Cytotoxic lymphocytes" ~ "MCP",
                                 variable == "cytotoxic cells danaher" ~ "Danaher",
                                 variable == "MCP Neutrophils" ~ "MCP",
                                 variable == "neutrophils danaher" ~ "Danaher",
                                 variable == "MCP Monocytic lineage" ~ "MCP",
                                 variable == "macrophages danaher" ~ "Danaher",
                                 variable == "TCRA.tcell.fraction.adj" ~ "Tracerx")) %>% 
  mutate(cell_group = factor(cell_group, levels = rev(c("T cell", "CD8 cell", "Cytotoxic cell", "B cell", "NK cell", "Macrophage", "Neutrophil"))))

# plot table 

g <- ggplot(all_res_dt %>% filter(!is.na(cell_group)), aes(x = trial_group, 
                                                           y = cell_group,
                                                           label = sign_p.adj)) +
  geom_tile(colour = 'black', aes(fill = coef)) +
  geom_fit_text(reflow = TRUE, min.size = 6 ) +
  scale_fill_gradient2(name = "Spearman's Rho", low = 'blue',mid = "white", high = 'red',limits = c(-1,1), na.value = 'gray') +
  theme_cowplot() + 
  theme(title = element_text(size=14),
        axis.title.x =element_blank(),  
        axis.title.y = element_blank(), 
        axis.text.x = element_text(size=14), #,angle = 90,vjust = 0.5, hjust=1),  
        axis.text.y = element_text(size=14),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line=element_blank(),
        legend.title = element_text(size=12),
        legend.text = element_text(size=10)) +  
  scale_x_discrete(position = "top")

print(g)  


