## Script for Tx421 LUSC paper
#Fig 3e
#Ext fig 3b, 4f
#Ext fig 4g, 4h, 5d
#Fig 2c,d, Fig 3i,j, Supp fig 5g,h,i
#Supp fig 8

#################################
### libraries and directories ###
#################################

library(tidyverse)
library(data.table)
library(fst)
library(readxl)
library(RColorBrewer)
#library(circlize)
library(ggfittext)
library(cowplot)
library(ggpubr)
library(car)
library(MASS)
library(rmcorr)
library(reshape2)
##################
## scatter plot ##
##################

#library(lme4)
library(rmcorr)
LUAD_region_df_full <- read_csv("../df-alvin/All_regions_500_regions-SCC-reviewed.csv")
GSVA_data <- read_csv("../df-alvin/gsva_GOBP_keratin_df.csv")
EMT_data <- readRDS("../df-alvin/gsva_GO_df.RDS") %>% dplyr::select(short_sample, contains("MESENCHYMAL_TRAN"))
colnames(EMT_data)[1] <-  "sample_name"
immune_data <- read_csv("../data/Tx421_MCPcounter.regional.csv")
danaher_data <- read_csv("../data/danaher_scores.csv")
danaher_data <- danaher_data %>% dplyr::select(ltx_region_id, score, cell_type) %>% 
  reshape2::dcast(formula = ltx_region_id ~ cell_type, value.var = "score")
colnames(danaher_data)[-1] <- danaher_data %>% dplyr::select(-ltx_region_id) %>% colnames() %>% str_c("_danaher") 
textract_data <- read_csv("../data/Tx421_TCRA_20211217.csv")
textract_data <- textract_data %>% dplyr::select(ltx_region_id, TCRA.tcell.fraction.adj)
LUAD_region_df <- LUAD_region_df_full  %>%
  filter(!is.na(Histology_per_tumour_id_muttable)) %>% 
  filter(Histology_per_tumour_id_muttable == "Squamous cell carcinoma") %>% 
  filter(included_500_regions == 1)

scc_regionalstat <- LUAD_region_df %>% dplyr::select(tumour_id_muttable_cruk, Region_ID_cruk, ltx_region_id, regional_pattern.new, wGII, all_TMB, clonal_TMB, subclonal_TMB, FLOH, purity, ploidy, keratinization_01, keratinization, contains("budding"), contains("invasion"), contains("sci")) %>% 
  left_join(GSVA_data, by = c("Region_ID_cruk" = "sample_name")) %>% 
  left_join(EMT_data, by = c("Region_ID_cruk" = "sample_name"))   %>% 
  left_join(immune_data, by = c("ltx_region_id" = "sample_name")) %>% 
  left_join(danaher_data, by = "ltx_region_id") %>% 
  left_join(textract_data, by = "ltx_region_id") 

tumor_df <- fread("../df-alvin//500_entire_tumour_with_grading.csv")
## regional TIL scores ##
all_regional_df <- readRDS("../data/all_regional_df_20221109.RDS")


regional_df <- scc_regionalstat
paramX <- "max_tumour_budding"
paramY <- "GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION"

in_df <- regional_df
in_df$paramX <- in_df[[paramX]]
in_df$paramY <- in_df[[paramY]]
in_df <- in_df %>% filter(!is.na(paramX) & !is.na(paramY))

rmcorr(tumour_id_muttable_cruk, paramX, paramY, in_df)



#########################################
## random sample one region per tumour ##
#########################################


paramX <- "max_tumour_budding"
paramY <- "GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION"

# paramX <- "num_single_cell_invasion"
# paramY <- "GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION"

n_trial = 1000 

in_df <- regional_df
in_df$paramX <- in_df[[paramX]]
in_df$paramY <- in_df[[paramY]]
in_data <- in_df %>% 
  filter(!is.na(paramX) & !is.na(paramY)) %>%
  dplyr::select(tumour_id_muttable_cruk, paramX, paramY)

set.seed(123) 

pval_vec <- vector() 
rho_vec <- vector() 

for(i in 1:n_trial){
  print(i)
  n_rows <- nrow(in_data)
  in_data$index <- sample(1:n_rows)
  in_data_shuffled <- in_data[in_data$index,] # now this is shuffled
  in_data_single <- in_data_shuffled %>%
    distinct(tumour_id_muttable_cruk, .keep_all = T)
  
  res_cor <- cor.test(in_data_single$paramX,in_data_single$paramY, method = "spearman" )
  pval <- res_cor$p.value
  rho <- res_cor$estimate
  pval_vec <- c(pval_vec, pval)
  rho_vec <- c(rho_vec, rho)
}


all_res <- data.table(pval = pval_vec , 
                      rho = rho_vec)


# plot histograms
pval_h <- hist(all_res$pval, breaks = 100)
print(pval_h)
rho_h <- hist(all_res$rho, breaks = 100) 
print(rho_h) #Ext fig 4h



# median and range 
paramX <- "max_tumour_budding"
paramY <- "GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION"

median(all_res$pval)
# [1] 0.0002419184
range(all_res$pval)
# [1] 1.378899e-08 5.266196e-02 
median(all_res$rho)
# [1] 0.35257
range(all_res$rho)
# [1] 0.1905636 0.5214720
percentile_95 <- quantile(all_res$pval, probs = 0.95)
print(percentile_95)
# 95% 
# 0.006659472 

paramX <- "num_single_cell_invasion"
paramY <- "GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION"

median(all_res$pval)
# [1] 0.002959902
range(all_res$pval)
# [1] 5.696862e-07 3.129627e-01
median(all_res$rho)
# [1] 0.2886743
range(all_res$rho)
# [1] 0.09990269 0.46732079
percentile_95 <- quantile(all_res$pval, probs = 0.95)
print(percentile_95)
# 95% 
# 0.04737252 

set.seed(48)

plot_scatter <- function(paramX, paramY){
  in_df <- regional_df
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

#Fig 3e
plot_scatter("max_tumour_budding", "GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION")

#Ext fig 4f
plot_scatter("num_single_cell_invasion", "GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION")



#Ext fig 4g, 5d

################################################################# 
## intra-tumour regional comparison: for 3q, 19q, E2F geneset 
## by regional grade  
#################################################################  

in_data <- regional_df %>%
  mutate(budding = factor(tumour_budding, levels = c("no_budding", "have_budding")),
         budding_high = case_when(budding_high.low %in% c("high_budding") ~ "high",
                                  budding_high.low %in% c("low_budding", "no_budding") ~ "no.low") ) %>%
  mutate(budding_high = factor(budding_high, levels = c("no.low", "high"))) %>%
  mutate(single_cell_invasion = factor(single_cell_invasion, levels = c("no_sci", "have_sci"))) 
colnames(in_data) <- gsub("-","_",colnames(in_data))

library(nlme)


variable_list <- c("GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION", "MCP_Fibroblasts", 
                   "GOBP_KERATINOCYTE_DIFFERENTIATION")

var_list <- gsub("-","_",variable_list)


out_path <- "../results/new_purity"
group_column_list <- c("budding", "keratinization")

for (group in group_column_list){
  out_f <- paste0(out_path,"/",group,"_all.region.pdf")
  pdf(out_f, width = 6, height = 5)
  
  for (var in  var_list) {
    
    tmp_df <- in_data
    tmp_df$group <- factor(tmp_df[[group]])
    tmp_df$var <- tmp_df[[var]]
    tmp_df <- tmp_df %>%
      filter(!(group %in% c("not.applicable", "NA")), !is.na(group)) %>%
      filter(!is.na(var))
    
    group.var_list <- sort(unique(tmp_df$group))
    
    if(length(group.var_list) > 1) {
      
      for(i in 1:(length(group.var_list)-1)) {
        #i=1
        for(j in (i+1):length(group.var_list)) {
          
          Group1 <- levels(tmp_df$group)[i]
          Group2 <- levels(tmp_df$group)[j]
          
          df <- tmp_df %>% 
            dplyr::filter(group %in% c(as.character(Group1), as.character(Group2)))
           
          count_group <- df %>%
            distinct(tumour_id_muttable_cruk, group) %>%
            group_by(tumour_id_muttable_cruk) %>%
            summarise(count =n())  
          
          hom_tum <- count_group$tumour_id_muttable_cruk[which(count_group$count == 1)]
          het_tum <- count_group$tumour_id_muttable_cruk[which(count_group$count > 1)]
          
          
          ## heterogenous tumour
          data <- df %>%
            filter(tumour_id_muttable_cruk %in% c(het_tum))
          
          if(length(unique(data$tumour_id_muttable_cruk))>0) {
            
            # 1. Calculate Purity Stats for the plotted points
            purity_summary <- data %>%
              group_by(group) %>%
              summarise(
                count = n(),
                m = mean(purity, na.rm = TRUE),
                med = median(purity, na.rm = TRUE),
                s = sd(purity, na.rm = TRUE)
              ) %>%
              mutate(label = paste0(group, ": (N=", count, ") Mean=", round(m, 2), 
                                    " Med=", round(med, 2), " SD=", round(s, 2)))
            
            # Combine labels into one text block for the caption
            purity_caption <- paste(purity_summary$label, collapse = "\n")
            
            # calculate pvalue using lme()
            if(length(unique(data$tumour_id_muttable_cruk))>2) {  
              
              data_df <- data %>%
                dplyr::select(Region_ID_cruk, tumour_id_muttable_cruk, group, var) %>%
                mutate(group = factor(group)) #, levels = rev(levels(regional_pattern))))
              names(data_df) <- c("short_sample", "tumour_id", "group", "var")
              
              effect_model <- lme(var ~ group, random = ~ 1 |tumour_id, data = data_df, method = "ML")
              null_model <- lme(var ~ 1, random = ~ 1 | tumour_id, data = data_df, method = "ML")
              anova_res <- anova(null_model, effect_model)
              # print(var)
              # print(anova_res)
              #print(summary(effect_model))
              
              anova_pval <- anova_res$`p-value`[2]
              pval <- round(anova_pval,6)
            } else {
              pval <- "NA"
            }
            
            y_range <- max(data$var) - min(data$var)
            #my_col <- createPalette(length(het_tum),  c("#ff0000", "#00ff00", "#0000ff"))
            
            num_tum <- length(unique(data$tumour_id_muttable_cruk))
            
            
            g <-  ggplot(data, aes(x = group, y = var)) +
              geom_jitter( size = 4, height = 0.0, fill = "gray40", #aes(fill = pred_subtype),  # shape = xxx   didn't work... I don't know why but shape was simply ignored
                           alpha = 0.8, colour="black", pch=21) +
              geom_boxplot(size = 1, width=0.7, alpha = 0, outlier.shape = NA) +
              labs(x =  "Regional pattern", y = var,
                   title = paste0("N= ", num_tum, " tumours (",sum(data_df$group == Group1)," & ",sum(data_df$group == Group2)," regions)"),
                   caption = purity_caption) +
              #guides(shape = guide_legend(reverse = TRUE)) +
              # facet_wrap(~ consent, 
              #            strip.position = "bottom", ncol = 5)+   #put facet label to the bottom 
              theme_cowplot()+
              #color_palette("jco")+
              theme(#strip.placement = "outside",   #move the facet label under axis
                strip.text = element_text(size = 12),
                #strip.background = element_blank(),
                axis.title.x = element_text(size = 26),
                axis.title.y = element_text(size = 22),
                axis.text.x = element_text(size = 22),
                axis.text.y = element_text(size = 18),
                legend.position = "right",
                legend.text = element_text(size = 18),
                legend.title = element_text(size = 22)) + 
              # scale_color_manual(name = 'tumour_id', values = my_col, #color25,
              #                    drop = F) + 
              # scale_fill_manual(values = c( "Muc" ="#FDBF6F" , "Lep" = "#A6CEE3", "Pap" = "#1F78B4", "Aci" = "#B2DF8A",
              #                               "Cri" = "#33A02C","Mic" = "#FB9A99", "Sol" = "#E31A1C"), drop = F) +
              # #theme(legend.position = "none") +
              annotate("text",size = 7, x = 1.5, y = (min(data$var) + y_range*1.1),
                       label = paste("p=",pval,"(lme)"))
            
            print(g)
            
          } #if(length(unique(data$tumour_id_muttable_cruk))>0)
        } #for(j in (i+1):length(group.var_list))
      } #for(i in 1:(length(group.var_list)-1)) {
    } #if(length(group.var_list) > 1) {
  } #for (var in var_list) {
  dev.off()
}
#Ext fig 3b, 4g, 5d

################################################
# comparison of non-keratinizing regions in mixed LUSC and pure non-Keratinized LUSC
################################################

regional_df$region_id_short <- gsub("_SU", "", regional_df$ltx_region_id)

regional_df <- regional_df %>%
  left_join(all_regional_df[,c("Region_ID", "TIL_not.just.stromal")], by = c("region_id_short" = "Region_ID"))

### 

tumor_with_kerat_regions <- unique(regional_df$tumour_id_muttable_cruk[which(regional_df$keratinization == "have_keratin")])
tumor_with_nonkerat_regions <- unique(regional_df$tumour_id_muttable_cruk[which(regional_df$keratinization == "no_keratin")])

kerat_tumor <- tumor_df$cruk_tumourid[which(tumor_df$keratinization == "have_keratin")]
nonkerat_tumor <- tumor_df$cruk_tumourid[which(tumor_df$keratinization == "no_keratin")]

tumor_with_both_regions <- unique(c(intersect(tumor_with_kerat_regions, tumor_with_nonkerat_regions),
                                    intersect(tumor_with_kerat_regions, nonkerat_tumor),
                                    intersect(kerat_tumor, tumor_with_nonkerat_regions)))

tumor_with_kerat_only <- intersect(kerat_tumor, setdiff(tumor_with_kerat_regions, tumor_with_nonkerat_regions))
tumor_with_nonkerat_only <- intersect(nonkerat_tumor, setdiff(tumor_with_nonkerat_regions, tumor_with_kerat_regions))


################################################################# 
## across-tumour regional comparison: 
#1 non-kerat in mixed vs pure non-kerat
#2 kerat in mixed vs pure-kerat
#################################################################  


in_data <- regional_df %>%
  mutate(budding = factor(tumour_budding, levels = c("no_budding", "have_budding")),
         budding_high = case_when(budding_high.low %in% c("high_budding") ~ "high",
                                  budding_high.low %in% c("low_budding", "no_budding") ~ "no.low") ) %>%
  mutate(budding_high = factor(budding_high, levels = c("no.low", "high"))) %>%
  mutate(single_cell_invasion = factor(single_cell_invasion, levels = c("no_sci", "have_sci"))) %>%
  mutate(kerat_in_pure.or.mix = case_when(keratinization == "have_keratin" & tumour_id_muttable_cruk %in% tumor_with_both_regions ~ "kerat_in_mixed",
                                          keratinization == "have_keratin" & tumour_id_muttable_cruk %in% tumor_with_kerat_only  ~ "kerat_in_pure"),
         nonkerat_in_pure.or.mix = case_when(keratinization == "no_keratin" & tumour_id_muttable_cruk %in% tumor_with_both_regions ~ "nonkerat_in_mixed",
                                             keratinization == "no_keratin" & tumour_id_muttable_cruk %in% tumor_with_nonkerat_only  ~ "nonkerat_in_pure"))

colnames(in_data) <- gsub("-","_",colnames(in_data))


#################################
### plot all region version  ### 
#################################

library(nlme)
# library(Polychrome)
# #num_het <- length(het_tum)
# my_col <- createPalette(50,  c("#ff0000", "#00ff00", "#0000ff"))
# my_col <- unname(my_col)
#library(pracma)  # nthroot()


#########################
## no color ## 
#########################

variable_list <- c("GOBP_KERATINOCYTE_DIFFERENTIATION", 
                   "t_cells_danaher" , "TCRA.tcell.fraction.adj",
                    "b_cells_danaher", "MCP_T cells",  "MCP_B lineage"  , 
                   "clonal_TMB","purity")

var_list <- gsub("-","_",variable_list)


out_path <- "../results/new_purity"
group_column_list <- c("nonkerat_in_pure.or.mix")

for (group in group_column_list){
  out_f <- paste0(out_path,"/",group,"_all.region.pdf")
  pdf(out_f, width = 6, height = 5)
  
  for (var in  var_list) {
    
    tmp_df <- in_data
    tmp_df$group <- factor(tmp_df[[group]])
    tmp_df$var <- tmp_df[[var]]
    tmp_df <- tmp_df %>%
      filter(!(group %in% c("not.applicable", "NA")), !is.na(group)) %>%
      filter(!is.na(var))
    
    group.var_list <- sort(unique(tmp_df$group))
    
    if(length(group.var_list) > 1) {
      
      for(i in 1:(length(group.var_list)-1)) {
        #i=1
        for(j in (i+1):length(group.var_list)) {
          
          Group1 <- levels(tmp_df$group)[i]
          Group2 <- levels(tmp_df$group)[j]
          
          df <- tmp_df %>% 
            dplyr::filter(group %in% c(as.character(Group1), as.character(Group2)))
          
          # df_s <- df[,c("tumour_id_muttable_cruk","short_sample", "regional_pattern")] %>%
          #   filter(tumour_id_muttable_cruk %in% het_tum)
          
          # count_group <- df %>%
          #   distinct(tumour_id_muttable_cruk, group) %>%
          #   group_by(tumour_id_muttable_cruk) %>%
          #   summarise(count =n())  
          # 
          # hom_tum <- count_group$tumour_id_muttable_cruk[which(count_group$count == 1)]
          # het_tum <- count_group$tumour_id_muttable_cruk[which(count_group$count > 1)]
          # 
          # 
          ## homogenous tumour 
          ## leave this alone this time..
          
          ## heterogenous tumour
          # data <- df %>%
          #   filter(tumour_id_muttable_cruk %in% c(het_tum))
          
          data <- df
          
          if(length(unique(data$tumour_id_muttable_cruk))>0) {
            
            # 1. Calculate Purity Stats for the plotted points
            purity_summary <- data %>%
              group_by(group) %>%
              summarise(
                count = n(),
                m = mean(purity, na.rm = TRUE),
                med = median(purity, na.rm = TRUE),
                s = sd(purity, na.rm = TRUE)
              ) %>%
              mutate(label = paste0(group, ": (N=", count, ") Mean=", round(m, 2), 
                                    " Med=", round(med, 2), " SD=", round(s, 2)))
            
            # Combine labels into one text block for the caption
            purity_caption <- paste(purity_summary$label, collapse = "\n")
            
          
            # calculate pvalue using lme()
            if(length(unique(data$tumour_id_muttable_cruk))>2) {  
              
              data_df <- data %>%
                dplyr::select(Region_ID_cruk, tumour_id_muttable_cruk, group, var) %>%
                mutate(group = factor(group)) #, levels = rev(levels(regional_pattern))))
              names(data_df) <- c("short_sample", "tumour_id", "group", "var")
              
              effect_model <- lme(var ~ group, random = ~ 1 |tumour_id, data = data_df, method = "ML")
              null_model <- lme(var ~ 1, random = ~ 1 | tumour_id, data = data_df, method = "ML")
              anova_res <- anova(null_model, effect_model)
              # print(var)
              # print(anova_res)
              #print(summary(effect_model))
              
              anova_pval <- anova_res$`p-value`[2]
              pval <- round(anova_pval,6)
            } else {
              pval <- "NA"
            }
            
            y_range <- max(data$var) - min(data$var)
            #my_col <- createPalette(length(het_tum),  c("#ff0000", "#00ff00", "#0000ff"))
            
            num_tum <- length(unique(data$tumour_id_muttable_cruk))
            
            
            g <-  ggplot(data, aes(x = group, y = var)) +
              geom_jitter( size = 4, height = 0.0, fill = "gray40", #aes(fill = pred_subtype),  # shape = xxx   didn't work... I don't know why but shape was simply ignored
                           alpha = 0.8, colour="black", pch=21) +
              geom_boxplot(size = 1, width=0.7, alpha = 0, outlier.shape = NA) +
              labs(x =  "Regional pattern", y = var,
                   title = paste0("N= ", num_tum, " tumours (",sum(data_df$group == Group1)," & ",sum(data_df$group == Group2)," regions)"),
                   caption = purity_caption) +
              #guides(shape = guide_legend(reverse = TRUE)) +
              # facet_wrap(~ consent, 
              #            strip.position = "bottom", ncol = 5)+   #put facet label to the bottom 
              theme_cowplot()+
              #color_palette("jco")+
              theme(#strip.placement = "outside",   #move the facet label under axis
                strip.text = element_text(size = 12),
                #strip.background = element_blank(),
                axis.title.x = element_text(size = 26),
                axis.title.y = element_text(size = 22),
                axis.text.x = element_text(size = 22),
                axis.text.y = element_text(size = 18),
                legend.position = "right",
                legend.text = element_text(size = 18),
                legend.title = element_text(size = 22)) + 
              # scale_color_manual(name = 'tumour_id', values = my_col, #color25,
              #                    drop = F) + 
              # scale_fill_manual(values = c( "Muc" ="#FDBF6F" , "Lep" = "#A6CEE3", "Pap" = "#1F78B4", "Aci" = "#B2DF8A",
              #                               "Cri" = "#33A02C","Mic" = "#FB9A99", "Sol" = "#E31A1C"), drop = F) +
              # #theme(legend.position = "none") +
              annotate("text",size = 7, x = 1.5, y = (min(data$var) + y_range*1.1),
                       label = paste("p=",pval,"(lme)"))
            
            print(g)
            
          } #if(length(unique(data$tumour_id_muttable_cruk))>0)
        } #for(j in (i+1):length(group.var_list))
      } #for(i in 1:(length(group.var_list)-1)) {
    } #if(length(group.var_list) > 1) {
  } #for (var in var_list) {
  dev.off()
}#for (group in group_column_list){
#Fig 2c,d, Fig 3i,j, Supp fig 5g,h,i


###################
### forest plot ### 
###################
# Supp fig 8

library(survival)
library(survminer)

in_data_surv <- readRDS("../data/in_data_surv_temp.RDS")  #This file is obtained from running the script "Fig5_Extfig7efghi.R".

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
         Age = age/10 ) %>%
  mutate(TNM = case_when(TNM_stage %in% c("I","II") ~ "Stage2", 
                         TNM_stage == "III" ~ "Stage3")) %>%
  mutate(TNM = factor(TNM, levels = c("Stage2","Stage3")))  %>%
  mutate(pure.mixed.K_nonK = case_when(cruk_id %in% tumor_with_both_regions ~ "mixed",
                                       cruk_id %in% tumor_with_kerat_only ~ "pure_K",
                                       cruk_id %in% tumor_with_nonkerat_only ~ "pure_NonK")) %>%
  mutate(pure.mixed.K_nonK_full = case_when(cruk_id %in% tumor_with_both_regions & keratinization == "have_keratin" ~ "mixed_K",
                                            cruk_id %in% tumor_with_both_regions & keratinization == "no_keratin" ~ "mixed_NonK",
                                            cruk_id %in% tumor_with_kerat_only ~ "pure_K",
                                            cruk_id %in% tumor_with_nonkerat_only ~ "pure_NonK")) %>%
  mutate(pure.mixed.K_nonK_partial = case_when(pure.mixed.K_nonK_full == "pure_NonK" ~ "NonK",
                                               pure.mixed.K_nonK_full == "mixed_NonK" ~ "NonK",
                                               pure.mixed.K_nonK_full == "mixed_K" ~ "mixed_K",
                                               pure.mixed.K_nonK_full == "pure_K" ~ "pure_K")) %>%
  mutate(pure.mixed.K_nonK = factor(pure.mixed.K_nonK, levels = c("pure_NonK","mixed","pure_K"  ))) %>%
  mutate(pure.mixed.K_nonK_full = factor(pure.mixed.K_nonK_full, levels = c("pure_NonK","mixed_NonK","mixed_K","pure_K"  ))) %>% 
  mutate(pure.mixed.K_nonK_partial = factor(pure.mixed.K_nonK_partial, levels = c("NonK","mixed_K","pure_K")))


in_df_surv <- in_df_surv_pre

time_vec <- c("dfs_time")
status_vec <- c("cens_dfs")
title_vec <- c("Disease-free survival")

# each <- "pure.mixed.K_nonK" # "max_regional_budding" # "keratinization"  #  , 
# variables_cox <- c("Age", "sex","TNM_stage",   # "PackYears", #Surgery + # margin_status +
#                    "Adjuvant" , each ) #, "LVI_per_patient")
# 
# out_pdf <- paste0("../results/multiCox_",each,".pdf")
# pdf(out_pdf)


# each <- "pure.mixed.K_nonK_full" # "max_regional_budding" # "keratinization"  #  ,
# variables_cox <- c("Age", "sex","TNM_stage",   # "PackYears", #Surgery + # margin_status +
#                    "Adjuvant" , each ) #, "LVI_per_patient")
# 
# 
# out_pdf <- paste0("../results/multiCox_",each,".pdf")
# pdf(out_pdf)

each <- "pure.mixed.K_nonK_partial" # "max_regional_budding" # "keratinization"  #  ,
variables_cox <- c("Age", "sex","TNM_stage",   # "PackYears", #Surgery + # margin_status +
                   "Adjuvant" , each ) #, "LVI_per_patient")

out_pdf <- paste0("../results/multiCox_",each,".pdf")
pdf(out_pdf)



for(var in c(variables_cox) ) {
  in_df_surv <- in_df_surv[!is.na(in_df_surv[[var]]),]
}

f <- as.formula(paste("Surv(time, status)",
                      paste(variables_cox, collapse = " + "), 
                      sep = " ~ "))

for (m in 1:length(time_vec)) {
  
  survival_data.frame <-in_df_surv
  colnames_time <-  time_vec[m]
  colnames_status <-  status_vec[m]
  
  surv_df <- survival_data.frame
  surv_df$time <- surv_df[[colnames_time]]
  surv_df$status <- surv_df[[colnames_status]]  

  cox_out <- coxph(f, data = surv_df)
  
  fp <- ggforest(cox_out, data = NULL, main = "Hazard ratio", fontsize = 0.7, refLabel = "reference", noDigits = 2)
  print(fp)
}

dev.off()


