##Script for Tx421 LUSC paper
##Fig 1a
# to remove rows from the plot, pay attention to the lines with #amend here
# in particular, this line should change to the number of rows c(1:16) 
# aligend.legends1 <- align.legends.fun(legends[c(1:16)]) #: add parameter amend here


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


#########################
# plot

p_combined_list <- lapply(Sub_list, function(x){
  if(x == Sub_list[1]){
    input_plot_df <- input_plot_df_all_ready[input_plot_df_all_ready$keratinization == Sub_list[1],]
  } else if (x == Sub_list[2]){
    input_plot_df <- input_plot_df_all_ready[input_plot_df_all_ready$keratinization == Sub_list[2],]
  } 

  
  order_samples <- input_plot_df$TumourID
  input_plot_df$TumourID  <- factor(input_plot_df$TumourID, levels = order_samples)
  order_tumour <- input_plot_df$TumourID
  
  p_keratinization <- ggplot(input_plot_df, aes(x = TumourID, y = 0.6, fill = keratinization)) +
    geom_tile(colour = 'white', size = gab_size) +
    scale_fill_manual(name = 'Keratinization',
                      values = c("no_keratin" = "#FB9A99", "have_keratin" = "#A6CEE3"), na.value = 'white', drop = F) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_discrete(expand = c(0,0), drop = F) +
    ylab('Keratinization') +
    theme_bw() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),
          text = element_text(size = 20),
          legend.text = element_text(size = 15), legend.title = element_text(size = 18))
  
    
  p_age <- ggplot(input_plot_df, aes(x = TumourID, y = 0.6 , fill = age)) +
    geom_tile(colour = 'white', size = gab_size) +
    scale_fill_gradient2(name = 'Age',  low = "#4575B4" , mid = "#FFFFBF", high = "#D73027" , midpoint = 60, limits = c(30,90), na.value = 'white') +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_discrete(expand = c(0,0), drop = F) +
    ylab('Age') +
    theme_bw() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),
          text = element_text(size = 20),
          legend.text = element_text(size = 15), legend.title = element_text(size = 18))

  # plot Gender

  p_sex <- ggplot(input_plot_df, aes(x = TumourID, y = 0.6, fill = sex)) +
    geom_tile(colour = 'white', size = gab_size) +
    scale_fill_manual(name = 'Gender',
                      values = c("Male" = "#2B3990", "Female" = "#EF4136"), na.value = 'white', drop = F) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_discrete(expand = c(0,0), drop = F) +
    ylab('Gender') +
    theme_bw() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),
          text = element_text(size = 20),
          legend.text = element_text(size = 15), legend.title = element_text(size = 18))

  # plot smoking status (3 levels)

  p_smoke <- ggplot(input_plot_df, aes(x = TumourID, y = 0.6, fill = smoking)) +
    geom_tile(colour = 'white', size = gab_size) +
    scale_fill_manual(name = 'Smoking status',
                      values = c("Never Smoked"="#187FC3","Ex-Smoker" = "#878787", "Smoker" = "black"), na.value = 'white', drop = F) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_discrete(expand = c(0,0), drop = F) +
    ylab('Smoking status') +
    theme_bw() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),
          text = element_text(size = 20),
          legend.text = element_text(size = 15), legend.title = element_text(size = 18))

  # plot TNM (v7)

  p_TNM <- ggplot(input_plot_df, aes(x = TumourID, y = 0.6, fill = pathologyTNM)) +
    geom_tile(colour = 'white', size = gab_size) +
    scale_fill_manual(name = 'TNM stage',
                      values = c("IA"="#EFF3FF", "IB"="#BDD7E7","IIA"="#6BAED6", "IIB"="#3182BD", "IIIA"="#08519C", "IIIB"="midnightblue"), na.value = 'white', drop = F) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_discrete(expand = c(0,0), drop = F) +
    ylab('TNM stage') +
    theme_bw() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),
          text = element_text(size = 20),
          legend.text = element_text(size = 15), legend.title = element_text(size = 18))

  p_bud <- ggplot(input_plot_df, aes(x = TumourID, y = 0.6, fill = budding_high.low)) +
    geom_tile(colour = 'white', size = gab_size) +
    scale_fill_manual(name = 'Tumour budding',
                      values = c("no_budding"="#8DC63F","low_budding" = "#006838", "high_budding" = "#7F3F98"), na.value = 'white', drop = F) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_discrete(expand = c(0,0), drop = F) +
    ylab('Tumour budding') +
    theme_bw() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),
          text = element_text(size = 20),
          legend.text = element_text(size = 15), legend.title = element_text(size = 18))

    p_invasion <- ggplot(input_plot_df, aes(x = TumourID, y = 0.6, fill = sci_high.low)) +
    geom_tile(colour = 'white', size = gab_size) +
    scale_fill_manual(name = 'Single cell invasion',
                      values = c("no_sci"="#8DC63F","low_sci" = "#006838", "high_sci" = "#7F3F98"), na.value = 'white', drop = F) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_discrete(expand = c(0,0), drop = F) +
    ylab('Single cell invasion') +
    theme_bw() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),
          text = element_text(size = 20),
          legend.text = element_text(size = 15), legend.title = element_text(size = 18))
  
    p_atypia <- ggplot(input_plot_df, aes(x = TumourID, y = 0.6, fill = atypia)) +
      geom_tile(colour = 'white', size = gab_size) +
      scale_fill_manual(name = 'Cytological atypia',
                        values = c("high_aty" = "#EC008C", "low_aty" = "#27AAE1"), na.value = 'white', drop = F) +
      scale_y_continuous(expand = c(0,0)) +
      scale_x_discrete(expand = c(0,0), drop = F) +
      ylab('Cytological atypia') +
      theme_bw() +
      theme(axis.text = element_blank(), axis.ticks = element_blank(),
            axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),
            text = element_text(size = 20),
            legend.text = element_text(size = 15), legend.title = element_text(size = 18))
    
    p_mitoses <- ggplot(input_plot_df, aes(x = TumourID, y = 0.6 , fill = mitoses)) +
      geom_tile(colour = 'white', size = gab_size) +
      scale_fill_gradient2(name = 'Mitoses',  low = "#4575B4" , mid = "#FFFFBF", high = "#D73027" , midpoint = 30, limits = c(0,60), na.value = 'white') +
      scale_y_continuous(expand = c(0,0)) +
      scale_x_discrete(expand = c(0,0), drop = F) +
      ylab('Mitoses') +
      theme_bw() +
      theme(axis.text = element_blank(), axis.ticks = element_blank(),
            axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),
            text = element_text(size = 20),
            legend.text = element_text(size = 15), legend.title = element_text(size = 18))

    p_necrosis <- ggplot(input_plot_df, aes(x = TumourID, y = 0.6, fill = extensive_necrosis1)) +
      geom_tile(colour = 'white', size = gab_size) +
      scale_fill_manual(name = 'Extensive necrosis',
                        values = c("have_necrosis" = "#f15a29", "no_necrosis" = "#808285"), na.value = 'white', drop = F) +
      scale_y_continuous(expand = c(0,0)) +
      scale_x_discrete(expand = c(0,0), drop = F) +
      ylab('Extensive necrosis') +
      theme_bw() +
      theme(axis.text = element_blank(), axis.ticks = element_blank(),
            axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),
            text = element_text(size = 20),
            legend.text = element_text(size = 15), legend.title = element_text(size = 18))

        
    #plot clonal WGD:
    #if see Error: Continuous value supplied to discrete scale, remember to set to factor at the line "#: Set as factors to avoid error of continuous variable "
    p_clonalWGD <- ggplot(input_plot_df, aes(x = TumourID, y = 0.6, fill = num_clonal_gds)) +
      geom_tile(colour = 'white', size = gab_size) +
      scale_fill_manual(name = 'Truncal GD',
                        values = c('0' = '#d1e5f04D', '1' = '#92c5de', '2' = '#2166ac'), na.value = 'white', drop = F) +
      scale_y_continuous(expand = c(0,0)) +
      scale_x_discrete(expand = c(0,0), drop = F) +
      ylab('Truncal GD') +
      theme_bw() +
      theme(axis.text = element_blank(), axis.ticks = element_blank(),
            axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),
            text = element_text(size = 20),
            legend.text = element_text(size = 15), legend.title = element_text(size = 18))
    
  
  
    p_subclonalWGD <- ggplot(input_plot_df, aes(x = TumourID, y = 0.6, fill = num_subclonal_gds)) +
      geom_tile(colour = 'white', size = gab_size) +
      scale_fill_manual(name = 'Subclonal GD',
                        values = c('0' = '#fddbc74D', '1' = '#f4a582', '2' = '#b2182b'), na.value = 'white', drop = F) +
      scale_y_continuous(expand = c(0,0)) +
      scale_x_discrete(expand = c(0,0), drop = F) +
      ylab('Subclonal GD') +
      theme_bw() +
      theme(axis.text = element_blank(), axis.ticks = element_blank(),
            axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),
            text = element_text(size = 20),
            legend.text = element_text(size = 15), legend.title = element_text(size = 18))

    p_SNVITH <- ggplot(input_plot_df, aes(x = TumourID, y = 0.6 , fill = perc_subclonal)) +
      geom_tile(colour = 'white', size = gab_size) +
      scale_fill_gradient2(name = 'SNV ITH', low = '#fde0ef', high = '#c51b7d', limits=c(0,1), na.value = 'white') +
      scale_y_continuous(expand = c(0,0)) +
      scale_x_discrete(expand = c(0,0), drop = F) +
      ylab('SNV ITH') +
      theme_bw() +
      theme(axis.text = element_blank(), axis.ticks = element_blank(),
            axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),
            text = element_text(size = 20),
            legend.text = element_text(size = 15), legend.title = element_text(size = 18))
    
    p_SCNAITH <- ggplot(input_plot_df, aes(x = TumourID, y = 0.6 , fill = frac_abberant_genom_subcl)) +
      geom_tile(colour = 'white', size = gab_size) +
      scale_fill_gradient2(name = 'SCNA ITH', low = '#e7d4e8', high = '#3f007d', limits=c(0,1), na.value = 'white') +
      scale_y_continuous(expand = c(0,0)) +
      scale_x_discrete(expand = c(0,0), drop = F) +
      ylab('SCNA ITH') +
      theme_bw() +
      theme(axis.text = element_blank(), axis.ticks = element_blank(),
            axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),
            text = element_text(size = 20),
            legend.text = element_text(size = 15), legend.title = element_text(size = 18))
    
    p_ploidy <- ggplot(input_plot_df, aes(x = TumourID, y = 0.6 , fill = ploidy)) +
      geom_tile(colour = 'white', size = gab_size) +
      scale_fill_gradient2(name = 'Ploidy', low = '#efd1bb', high = '#442610', limits=c(1,5), na.value = 'white') +
      scale_y_continuous(expand = c(0,0)) +
      scale_x_discrete(expand = c(0,0), drop = F) +
      ylab('Ploidy') +
      theme_bw() +
      theme(axis.text = element_blank(), axis.ticks = element_blank(),
            axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),
            text = element_text(size = 20),
            legend.text = element_text(size = 15), legend.title = element_text(size = 18)) 
    
    p_LVI <- ggplot(input_plot_df, aes(x = TumourID, y = 0.6, fill = LVI_per_patient)) +
      geom_tile(colour = 'white', size = gab_size) +
      scale_fill_manual(name = 'LVI_per_patient',
                        values = c("Yes" = "#EC2227", "No" = "#2075BC"), na.value = 'white', drop = F) +
      scale_y_continuous(expand = c(0,0)) +
      scale_x_discrete(expand = c(0,0), drop = F) +
      ylab('Lymphovascular invasion') +
      theme_bw() +
      theme(axis.text = element_blank(), axis.ticks = element_blank(),
            axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),
            text = element_text(size = 20),
            legend.text = element_text(size = 15), legend.title = element_text(size = 18))
    
    
  ##########################################
  ###          combine plots             ###
  ##########################################
  #align legends
  legends         <- list(g_legend(p_keratinization), #: add parameter amend here
                          g_legend(p_age),
                          g_legend(p_sex),
                          g_legend(p_smoke),
                          g_legend(p_TNM), 
                          g_legend(p_bud), 
                          g_legend(p_invasion), 
                          g_legend(p_atypia), 
                          g_legend(p_mitoses),
                          g_legend(p_necrosis),
                          g_legend(p_clonalWGD), 
                          g_legend(p_subclonalWGD),
                          g_legend(p_SNVITH),
                          g_legend(p_SCNAITH),
                          g_legend(p_ploidy), 
                          g_legend(p_LVI)
                        
  )
 
  aligend.legends1 <- align.legends.fun(legends[c(1:16)]) #: add parameter amend here
  
  if(x == Sub_list[1]){        #choosing "have_keratin" group to produce legend  
  l_keratinization     <<- p_keratinization    
  l_age   <<-   p_age  
  l_sex     <<-   p_sex    
  l_smoke     <<-   p_smoke    
  l_TNM      <<-   p_TNM     
  l_bud  <<-   p_bud 
  l_invasion      <<-   p_invasion     
  l_atypia   <<-   p_atypia  
  l_mitoses    <<-   p_mitoses 
  l_necrosis     <<- p_necrosis
  l_clonalWGD <<-   p_clonalWGD
  l_subclonalWGD  <<-   p_subclonalWGD 
  l_SNVITH  <<-   p_SNVITH 
  l_SCNAITH  <<-   p_SCNAITH 
  l_ploidy     <<- p_ploidy
  l_LVI <<- p_LVI
  }
  

  p_keratinization       <- p_keratinization + theme(legend.position = 'none')
  p_age       <- p_age + theme(legend.position = 'none')
  p_sex       <- p_sex + theme(legend.position = 'none')
  p_smoke       <- p_smoke + theme(legend.position = 'none')
  p_TNM       <- p_TNM + theme(legend.position = 'none')
  p_bud       <- p_bud + theme(legend.position = 'none')
  p_invasion       <- p_invasion + theme(legend.position = 'none')
  p_atypia       <- p_atypia + theme(legend.position = 'none')
  p_mitoses       <- p_mitoses + theme(legend.position = 'none')
  p_necrosis <- p_necrosis + theme(legend.position = 'none')
  p_clonalWGD <- p_clonalWGD + theme(legend.position = 'none')
  p_subclonalWGD <- p_subclonalWGD + theme(legend.position = 'none')
  p_SNVITH <- p_SNVITH + theme(legend.position = 'none')
  p_SCNAITH <- p_SCNAITH + theme(legend.position = 'none')
  p_ploidy <- p_ploidy + theme(legend.position = 'none')
  p_LVI <- p_LVI + theme(legend.position = 'none')
  
  #delete y-labels except the one on the very left
  if(x != Sub_list[1]){
     p_keratinization <- p_keratinization + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
     p_age       <- p_age + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
    p_sex       <- p_sex + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
    p_smoke       <- p_smoke + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
    p_TNM       <- p_TNM + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
    p_bud       <- p_bud + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
    p_invasion       <- p_invasion + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
    p_atypia       <- p_atypia + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
    p_mitoses       <- p_mitoses + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
    p_necrosis       <- p_necrosis + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
    p_clonalWGD <- p_clonalWGD +  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
    p_subclonalWGD <- p_subclonalWGD +  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
    p_SNVITH <- p_SNVITH +  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
    p_SCNAITH <- p_SCNAITH +  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
    p_ploidy       <- p_ploidy + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
    p_LVI       <- p_LVI + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
      }
  
p_keratinization <- ggplotGrob(p_keratinization + theme(plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm")))
  p_age <- ggplotGrob(p_age + theme(plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm")))
  p_sex <- ggplotGrob(p_sex + theme(plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm")))
  p_smoke <- ggplotGrob(p_smoke + theme(plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm")))
  p_TNM <- ggplotGrob(p_TNM + theme(plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm")))  # this is the bottome this time
  p_bud <- ggplotGrob(p_bud + theme(plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm")))
  p_invasion <- ggplotGrob(p_invasion + theme(plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm")))
  p_atypia <- ggplotGrob(p_atypia + theme(plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm")))
  p_mitoses <- ggplotGrob(p_mitoses + theme(plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm")))
  p_necrosis <- ggplotGrob(p_necrosis + theme(plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm")))
  p_clonalWGD <- ggplotGrob(p_clonalWGD + theme(plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm")))
  p_subclonalWGD <- ggplotGrob(p_subclonalWGD + theme(plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm")))
  p_SNVITH <- ggplotGrob(p_SNVITH + theme(plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm")))
  p_SCNAITH <- ggplotGrob(p_SCNAITH + theme(plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm")))
  p_ploidy <- ggplotGrob(p_ploidy + theme(plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm")))  
  p_LVI <- ggplotGrob(p_LVI + theme(plot.margin = unit(c(0.1, 0.1, 0, 0.1), "cm")))
   
   grid_unit_pmax <- grid::unit.pmax(p_keratinization$widths,
                                    p_age$widths,p_sex$widths,p_smoke$widths,  #
                                    p_TNM$widths, p_bud$widths, p_invasion$widths,
                                    p_atypia$widths, p_mitoses$widths, p_necrosis$widths,
                                    p_clonalWGD$widths,
                                    p_subclonalWGD$widths,
                                    p_SNVITH$widths,
                                    p_SCNAITH$widths, p_ploidy$widths, 
                                    p_LVI$widths
                                    
                                    
  )
  

  p_keratinization$widths    <-  grid_unit_pmax
  p_age$widths    <-  grid_unit_pmax
  p_sex$widths    <-  grid_unit_pmax
  p_smoke$widths    <-  grid_unit_pmax
  p_TNM$widths    <-  grid_unit_pmax
  p_bud$widths<-  grid_unit_pmax
  p_invasion$widths<-  grid_unit_pmax
  p_atypia$widths<-  grid_unit_pmax
  p_mitoses$widths<-  grid_unit_pmax
  p_necrosis$widths <-  grid_unit_pmax
  p_clonalWGD$widths<-  grid_unit_pmax
  p_subclonalWGD$widths<-  grid_unit_pmax
  p_SNVITH$widths<-  grid_unit_pmax
  p_SCNAITH$widths<-  grid_unit_pmax
  p_ploidy$widths<-  grid_unit_pmax
  p_LVI$widths<-  grid_unit_pmax
  
  g <- gtable_rbind(p_keratinization,
                   p_age,
                    p_sex,
                    p_smoke,
                    p_TNM, 
                   p_bud,
                   p_invasion,
                   p_atypia, 
                   p_mitoses,
                   p_necrosis,
                   p_clonalWGD,
                   p_subclonalWGD,
                   p_SNVITH,
                   p_SCNAITH,
                   p_ploidy, 
                   p_LVI
                  )
                  

  id_panels_h <- unique(g$layout[g$layout$name=="panel", "t"])
  g$heights[id_panels_h] <- grid::unit(c(1, 0.8, 0.8, 0.8, 0.8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), "null") #
  
  g_combined <- g
  
  return(g_combined)
  
})


#plot

pdf(paste0('../results/Fig1a.LUSC_overview_Final.pdf'), width = 20, height = 6)
grid.arrange(p_combined_list[[1]], p_combined_list[[2]],  
             ncol = 2, widths = c(30, 15)   # add extra margin to first column    #adjust the width of each group
             )
print(l_keratinization    )
print(  l_age  )
print(  l_sex    )
print(  l_smoke    )
print(  l_TNM     )
print(  l_bud )
print(  l_invasion     )
print(  l_atypia  )
print(  l_mitoses   )
print( l_necrosis )
print(l_clonalWGD)
print(  l_subclonalWGD )
print(  l_SNVITH )
print(  l_SCNAITH )
print( l_ploidy )
print( l_LVI )
dev.off() 
# Fig 1a