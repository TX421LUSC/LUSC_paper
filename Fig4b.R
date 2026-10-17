##Script for Tx421 LUSC paper
#Fig 4b
#Ext fig 5

#################################
### libraries and directories ###
#################################

library(tidyverse)
library(data.table)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(ggfittext)
library(cowplot)
library(ggpubr)
library(ggsci)
library(gtable)
library(exactRankTests)


########################
### helper functions ###
########################


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



squamous_region_data0 <- read_csv("../df-alvin/All_regions_500_regions-SCC_reviewed_With_mean_keratinMitosis.csv")
squamous_region_data <- squamous_region_data0 %>%  
  mutate(seeding_region_label = case_when(seeding_region == "TRUE" ~ "Seeding",
                                          seeding_region == "FALSE" ~ "Nonseeding", 
                                          is.na(seeding_region) ~ "Unknown" ))
seeding_data_metpaper <- read_tsv("../data/seedingRegionInfo.txt")
metastasis_info_complete <- read_csv("../df-alvin/metastasis_info_LUSC.csv")
met_cases <- metastasis_info_complete$tumour_id
#met_cases <- c("CRUK0083", "CRUK0444", "CRUK0589", "CRUK0691") #only include cases with ancestral/descendant relation
squamous_met <- squamous_region_data %>% mutate(is_metastatic_case = case_when(tumour_id_muttable_cruk %in% met_cases ~ "TRUE" 
                                                               )) %>% 
                filter(is_metastatic_case == "TRUE")
squamous_plot <- squamous_met %>% dplyr::select(tumour_id_muttable_cruk, Region_ID_cruk, region, differentiation, keratinization, max_tumour_budding, num_single_cell_invasion, mitoses, atypia, seeding_region_label) %>% 
  mutate_all(~replace(., is.na(.), "No_data")) %>% 
  filter(!grepl("CRUK0495_T2", Region_ID_cruk)) %>% 
  filter(!grepl("CRUK0495_T3", Region_ID_cruk))

met_sites0 <- read_csv("../df-alvin/metastatic_sites_LUSC.csv")
met_sites <- met_sites0 %>% dplyr::select(tumour_id, met_site)
met_sites <- met_sites[!duplicated(met_sites), ]

squamous_plot_ready <- left_join(squamous_plot, met_sites, by = c("tumour_id_muttable_cruk" = "tumour_id"))


# check the max number of mets samples per tumour 
unique(squamous_plot_ready$region) 

squamous_plot_ready$seeding_region_label <- factor( squamous_plot_ready$seeding_region_label , levels=c("Seeding", "Nonseeding", "Unknown"))

# prep for pattern plotting of mets
primary_met_full <- data.frame(tumour_id = rep(unique(squamous_plot_ready$tumour_id_muttable_cruk), 8),
                               region = unlist(lapply(1:8, function(x){
                                 rep(paste0("R",x),  length(unique(squamous_plot_ready$tumour_id_muttable_cruk)))})), 
                               stringsAsFactors = F) %>%
  left_join(squamous_plot_ready, by = c("tumour_id" = "tumour_id_muttable_cruk", "region" = "region")) %>%
  arrange(seeding_region_label) %>% 
  arrange(tumour_id) %>% 
  mutate(region.number.by.pattern.rank = rep(paste0("r",1:8),  length(unique(squamous_plot_ready$tumour_id_muttable_cruk)))) 

  
  input_plot_df_regional <- primary_met_full
  # plot regional pattern with seeding information
  # confirm the seeding data with mets paper
  # input_plot_df_regional1 <- input_plot_df_regional %>% full_join(seeding_data_metpaper[, c("Region_ID_cruk", "Metastasizing")], by = "Region_ID_cruk")
  # input_plot_df_regional1[, c("Region_ID_cruk", "seeding_region_label", "Metastasizing")] %>%  view()

seeding_color <- c("Seeding" = "#92278f", 
                   "Nonseeding" = "#8dc63f", 
                   "Unknown" = "#bababa")  


#metastatic site
met_color <-c("lymphnode" = "#c49a6c", 
              "brain_lymphnode" = "#f15a29", 
              "lymphnode_pleura" = "#fbb040",
              "lymphnode_peritoneal" = "#9e1f63", 
              "lung" = "#c2b59b",
              "lung_lymphnode" = "#594a42", 
              "stump" = "#00a79d", 
              "adrenal" = "#f9ed32", 
              "neckLN" = "#39b54a")

input_plot_df_regional_met <- input_plot_df_regional %>% dplyr::select(tumour_id, met_site) %>% filter(!is.na(met_site)) %>% distinct()


p_met <- ggplot(input_plot_df_regional_met, aes(x = tumour_id, y = 0.6, fill = met_site)) +
  geom_tile(width = 0.8, height = 0.8, size = 1) +
  scale_fill_manual(name = 'met',
                    values = met_color, na.value = 'white', drop = F) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0.02,0.02), drop = F) +
  ylab('Metastatic site') +
  theme_void() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5),
        text = element_text(size = 20),
        legend.text = element_text(size = 15), legend.title = element_text(size = 18))
print(p_met)


  #Keratinization
  p_region.keratinization <- ggplot(input_plot_df_regional, aes(x = tumour_id, y = region.number.by.pattern.rank, 
                                                      fill = keratinization, colour = seeding_region_label, label = seeding_region_label)) + 
    geom_tile(width = 0.8, height = 0.8, size = 1) +
    scale_fill_manual(name = 'keratinization', 
                      values = c("No_data" = "#bababa", "have_keratin" = "#1F78B4", "no_keratin" = "#ed1c24"), na.value = 'white', drop = F) +
    scale_color_manual(name = 'Seeding status', 
                       values = seeding_color, na.value = 'white', drop = F) +
    scale_x_discrete(expand = c(0.02,0.02), drop = F) +
    ylab('Keratinization\n&\nseeding status') +
    theme_bw() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), 
          axis.title.x = element_blank(), 
          axis.title.y = element_text(angle = 0, vjust = 0.5),
          text = element_text(size = 20),
          legend.text = element_text(size = 15), legend.title = element_text(size = 18))
  
  print(p_region.keratinization)
  
  
  # seeding status
  
  
  p_region.seed <- ggplot(input_plot_df_regional, aes(x = tumour_id, y = region.number.by.pattern.rank, 
                                                      fill = seeding_region_label, colour = seeding_region_label, label = seeding_region_label)) + 
    geom_tile(width = 0.8, height = 0.8, size = 0.2) +
    scale_fill_manual(name = 'Seeding status', 
                      values = seeding_color , na.value = 'white', drop = F) +
    scale_color_manual(name = 'Seeding status', 
                       values = seeding_color, na.value = 'white', drop = F) +
    scale_x_discrete(expand = c(0.02,0.02), drop = F) +
    ylab("Seeding status") +
    theme_bw() +
    theme(axis.ticks = element_blank(), 
          axis.text.x =  element_text(angle = 30, hjust = 1, vjust = 1),
          axis.text.y =  element_blank(),
          axis.title.x = element_blank(), 
          axis.title.y = element_text(angle = 0, vjust = 0.5),
          text = element_text(size = 20),
          legend.text = element_text(size = 15), legend.title = element_text(size = 18))
  
  print(p_region.seed)
  
  
  
  #Budding
  input_plot_df_regional %>% mutate(budding_data = case_when(max_tumour_budding == "No_data" ~ as.numeric(NA), 
                                                             TRUE ~ as.numeric(max_tumour_budding)))  ->
    input_plot_df_regional
  
  p_region.budding <- ggplot(input_plot_df_regional, aes(x = tumour_id, y = region.number.by.pattern.rank, 
                                                         fill = budding_data, colour = seeding_region_label, label = seeding_region_label)) + 
    geom_tile(width = 0.8, height = 0.8, size = 1) +
    scale_fill_gradient(name = 'Tumour\nbudding', low = "#b7c9e2", high = '#3f007d', limits=c(0,15), na.value = 'white') +
    scale_color_manual(name = 'Seeding status', 
                       values = seeding_color, na.value = 'white', drop = F) +
    scale_x_discrete(expand = c(0.02,0.02), drop = F) +
    ylab('Tumour\nbudding\n&\nseeding status') +
    theme_bw() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), 
          axis.title.x = element_blank(), 
          axis.title.y = element_text(angle = 0, vjust = 0.5),
          text = element_text(size = 20),
          legend.text = element_text(size = 15), legend.title = element_text(size = 18))
  
  print(p_region.budding)
  

  #single cell invasion
  input_plot_df_regional %>% mutate(sci_data = case_when(num_single_cell_invasion == "No_data" ~ as.numeric(NA), 
                                                             TRUE ~ as.numeric(num_single_cell_invasion)))  ->
    input_plot_df_regional
  
  p_region.sci <- ggplot(input_plot_df_regional, aes(x = tumour_id, y = region.number.by.pattern.rank, 
                                                         fill = sci_data, colour = seeding_region_label, label = seeding_region_label)) + 
    geom_tile(width = 0.8, height = 0.8, size = 1) +
    scale_fill_gradient(name = 'Single\ncell\ninvasion', low = "#FFFFBF", high = '#c51b7d',  limits=c(0,15), na.value = 'white') +
    scale_color_manual(name = 'Seeding status', 
                       values = seeding_color, na.value = 'white', drop = F) +
    scale_x_discrete(expand = c(0.02,0.02), drop = F) +
    ylab('Single cell\ninvasion\n&\nseeding status') +
    theme_bw() +
    #theme_bw() + 
    theme(axis.text = element_blank(), axis.ticks = element_blank(), 
          axis.title.x = element_blank(), 
          axis.title.y = element_text(angle = 0, vjust = 0.5),
          text = element_text(size = 20),
          legend.text = element_text(size = 15), legend.title = element_text(size = 18))
  
  print(p_region.sci)  
  

  #plot ancestor descendant relationship
ancestor_table <-   read_csv("../df-alvin/ancestor_descent_table.csv")
ancestor_table1 <-  left_join(ancestor_table, input_plot_df_regional, by = c("assigned_region" = "Region_ID_cruk"))
p_annotated <- p_region.seed + 
    geom_text(
    data = ancestor_table1,  # Specify the new data source here
    aes(x = tumour_id, y = region.number.by.pattern.rank, label = labeling),
    #vjust = 3,
    colour = "black",  
    na.rm = TRUE
  )
print(p_annotated)


p_regionname <- p_region.seed + 
  geom_text(
    data = ancestor_table1,  # Specify the new data source here
    aes(x = tumour_id, y = region.number.by.pattern.rank, label = region),
    #vjust = 3,
    colour = "black",  
    na.rm = TRUE
  )
print(p_regionname)


# add keratin status of metastasis
met_keratin <- read_csv("../df-alvin/met_keratin_status.csv") %>% 
  dplyr::select(regionID, contains("keratin"), contains("note"))
met_keratin_plot <- metastasis_info_complete %>% left_join(met_keratin, by = c("regionID")) %>% distinct()


# Load necessary libraries
library(dplyr)
library(ggplot2)
library(ggalluvial)
library(tidyr)

# 1. Load the datasets
primary_seeding <- input_plot_df_regional
met_plot <- met_keratin_plot

# 2. Summarize Primary Seeding Status per Tumour
# We categorize the primary based on the collective status of its "Seeding" regions
primary_summary <- primary_seeding %>%
  filter(seeding_region_label == "Seeding") %>%
  group_by(tumour_id) %>%
  summarise(
    # Count the specific regional statuses
    n_seeding_kerat    = sum(keratinization == "have_keratin", na.rm = TRUE),
    n_seeding_nonkerat = sum(keratinization == "no_keratin", na.rm = TRUE),
    n_seeding_nodata   = sum(is.na(keratinization) | keratinization == "", na.rm = TRUE),
    
    primary_label = case_when(
      n_seeding_kerat > 0 & n_seeding_nonkerat == 0 ~ "Primary: Keratin",
      n_seeding_nonkerat > 0 & n_seeding_kerat == 0 ~ "Primary: No Keratin",
      n_seeding_kerat > 0 & n_seeding_nonkerat > 0 ~ "Primary: Mixed",
      TRUE                                          ~ "No keratinization data"
  ))

# 3. Merge and Clean Data
plot_df <- met_plot %>%
  left_join(primary_summary, by = "tumour_id") %>%
  filter(!is.na(met_keratinization)) %>%
  mutate(
    # Standardize metastatic labels
    met_label = ifelse(met_keratinization == "met_have_keratin", "Have Keratin", "No Keratin"),
    
    # Handle tumors that exist in met_plot but were NOT found in primary_summary
    primary_label = ifelse(is.na(primary_label), "Not in Primary File", as.character(primary_label)),
    
    # Set factor levels to include all possible categories for the plot
    primary_label = factor(primary_label, 
                           levels = c("Primary: Keratin", 
                                      "Primary: Mixed", 
                                      "Primary: No Keratin", 
                                      "No keratinization data", 
                                      "Not in Primary File")),
    met_label = factor(met_label, levels = c("Have Keratin", "No Keratin"))
  )


# --- PLOT 1 - Boxplot of Percent Keratinization ---
# Added x-axis label rotation to handle the longer category names
ggplot(plot_df, aes(x = primary_label, y = met_percent_keratinisation, fill = primary_label)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, color = "black", alpha = 0.5) +
  labs(title = "Metastatic Keratinization % by Primary Seeding Status",
       x = "Primary Seeding Status", y = "Metastatic Keratinization (%)") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# --- PLOT 2 - Stacked Bar Chart (Proportions) ---
ggplot(plot_df, aes(x = primary_label, fill = met_label)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("Have Keratin" = "#e74c3c", "No Keratin" = "#3498db")) +
  labs(title = "Proportion of Metastatic Status by Primary Seeding Status",
       x = "Primary Seeding Status", y = "Percentage of Metastases", fill = "Metastatic Status") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# --- SECTION 6: PLOT 3 - Heatmap (Counts) ---
plot_df %>%
  count(primary_label, met_label) %>%
  complete(primary_label, met_label, fill = list(n = 0)) %>% # Ensure all combinations show
  ggplot(aes(x = met_label, y = primary_label, fill = n)) +
  geom_tile() +
  geom_text(aes(label = n), color = "black", size = 5) +
  scale_fill_gradient(low = "#f1f5f9", high = "#3b82f6") +
  labs(title = "Counts: Primary Seeding vs Metastatic Status",
       x = "Metastatic Status", y = "Primary Seeding Status", fill = "Tumour Count") +
  theme_minimal()

# --- SECTION 7: PLOT 4 - Sankey / Alluvial Diagram ---
# This plot now visualizes the flow even for samples with missing primary data
ggplot(plot_df, aes(axis1 = primary_label, axis2 = met_label, y = 1)) +
  geom_alluvium(aes(fill = met_label), width = 1/12, alpha = 0.4) +
  geom_stratum(width = 1/4, fill = "grey95", color = "grey") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("Primary Seeding Status", "Metastatic Status"), expand = c(.1, .1)) +
  scale_fill_manual(values = c("Have Keratin" = "#e74c3c", "No Keratin" = "#3498db")) +
  labs(title = "Lineage Flow: Primary Seeding to Metastasis",
       y = "Number of Tumours", fill = "Metastatic Phenotype") +
  theme_minimal() +
  theme(legend.position = "bottom")

library(dplyr)
library(ggplot2)
library(ggalluvial)
library(tidyr)

# --- SECTION 2: Summarize Primary Seeding Status per Tumour ---
primary_summary <- primary_seeding %>%
  filter(seeding_region_label == "Seeding") %>%
  group_by(tumour_id) %>%
  summarise(
    n_seeding_kerat    = sum(keratinization == "have_keratin", na.rm = TRUE),
    n_seeding_nonkerat = sum(keratinization == "no_keratin", na.rm = TRUE),
    n_seeding_nodata   = sum(is.na(keratinization) | keratinization == "" | keratinization == "No_data", na.rm = TRUE),
    
    primary_label = case_when(
      n_seeding_kerat > 0 & n_seeding_nonkerat == 0 ~ "Primary: Keratin",
      n_seeding_nonkerat > 0 & n_seeding_kerat == 0 ~ "Primary: No Keratin",
      n_seeding_kerat > 0 & n_seeding_nonkerat > 0 ~ "Primary: Mixed",
      TRUE                                          ~ "No keratinization data"
    )
  )

# --- SECTION 3: Merge and Clean Data (LEFT JOIN) ---
plot_df <- met_plot %>%
  left_join(primary_summary, by = "tumour_id") %>%
  filter(!is.na(met_keratinization)) %>%
  mutate(
    met_label = ifelse(met_keratinization == "met_have_keratin", "Have Keratin", "No Keratin"),
    primary_label = ifelse(is.na(primary_label), "Not in Primary File", as.character(primary_label)),
    
    # Set factor levels
    primary_label = factor(primary_label, 
                           levels = c("Primary: Keratin", "Primary: Mixed", "Primary: No Keratin", 
                                      "No keratinization data", "Not in Primary File")),
    met_label = factor(met_label, levels = c("Have Keratin", "No Keratin"))
  )

# --- PLOT 1: Boxplot ---
ggplot(plot_df, aes(x = primary_label, y = met_percent_keratinisation, fill = primary_label)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, color = "black", alpha = 0.5) +
  labs(title = "Metastatic Keratinization % by Primary Seeding Status",
       x = "Primary Seeding Status", y = "Metastatic Keratinization (%)") +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

# --- PLOT 2: Stacked Bar Chart ---
color_list <-  c("Have Keratin" = "#A6CEE3", "No Keratin" = "#FB9A99")
ggplot(plot_df, aes(x = primary_label, fill = met_label)) +
  geom_bar(position = "fill") +
  geom_bar_text( position = "fill", reflow = TRUE) + 
  scale_fill_manual(values = color_list, na.value = '#bababa', drop = F)+
  theme_cowplot() +
  ylab("proportion") +
  scale_fill_manual(values = color_list, na.value = '#bababa', drop = F) +
  labs(title = "Proportion of Metastatic Status by Primary Seeding Status",
       x = "Primary Seeding Status", y = "Percentage", fill = "Metastasis") +
  theme(axis.title.x = element_blank(),  
        axis.title.y = element_text(size=24),
        axis.text.x = element_text(size=20, angle = 45,  hjust=1), 
        axis.text.y = element_text(size=16),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size=18))

# 1. Define your color list
color_list <- c("Have Keratin" = "#A6CEE3", "No Keratin" = "#FB9A99")

# 2. Build the plot to match your preferred style
g2 <- ggplot(plot_df, aes(x = primary_label, y = count, fill = met_label, label = count)) +
  # Add stat="identity" and color="black" for the outlines
  geom_bar(stat = "identity", position = "fill", color = "black") +
  # Text labels within the bars
  geom_bar_text(position = "fill", reflow = TRUE) + 
  # Fill scale with your specific NA color and drop=F setting
  scale_fill_manual(values = color_list, na.value = '#bababa', drop = FALSE) +
  # Apply cowplot theme
  theme_cowplot() +
  # Labels and Title
  ylab("proportion") +
  ggtitle("Proportion of Metastatic Status by Primary Seeding Status") +
  # Specific theme overrides from your preferred format
  theme(
    axis.title.x = element_blank(),  
    axis.title.y = element_text(size = 24),
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1), 
    axis.text.y = element_text(size = 16),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 18)
  )

print(g2)


# --- PLOT 3: Heatmap ---
plot_df %>%
  count(primary_label, met_label) %>%
  complete(primary_label, met_label, fill = list(n = 0)) %>%
  ggplot(aes(x = met_label, y = primary_label, fill = n)) +
  geom_tile() +
  geom_text(aes(label = n), size = 5) +
  scale_fill_gradient(low = "#f1f5f9", high = "#3b82f6") +
  labs(title = "Counts: Primary Seeding vs Metastatic Status",
       x = "Metastatic Status", y = "Primary Seeding Status") +
  theme_minimal()

# --- PLOT 4: 2-Axis Sankey with Outer ID Labels ---
ggplot(plot_df, aes(axis1 = primary_label, axis2 = met_label, y = 1)) +
  # Flows
  geom_alluvium(aes(fill = met_label), width = 1/12, alpha = 0.4) +
  
  # Categorical Strata (Boxes)
  geom_stratum(width = 1/4, fill = "grey95", color = "grey70") +
  
  # Category Labels (inside the boxes)
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, fontface = "bold") +
  
  # LEFTMOST TEXT: Tumour IDs
  # x = 1 places it at the first axis, hjust = 1.1 moves it to the left side
  geom_text(stat = "alluvium", aes(label = tumour_id), x = 1, hjust = 1, size = 1) +
  
  # RIGHTMOST TEXT: Region IDs
  # x = 2 places it at the second axis, hjust = -0.1 moves it to the right side
  geom_text(stat = "alluvium", aes(label = regionID), x = 2, hjust = -0, size = 1) +
  
  # Expansion settings ensure labels aren't cut off at the edges
  # scale_x_discrete(limits = c("Primary Seeding Status", "Metastatic Status"), 
  #                  expand = c(.4, .4)) +
  # 
  scale_fill_manual(values = c("Have Keratin" = "#e74c3c", "No Keratin" = "#3498db")) +
  labs(title = "Lineage Flow with Sample Identification",
       subtitle = "Left: Tumour ID | Right: Metastatic Region ID",
       y = "Number of Metastases", 
       fill = "Metastatic Phenotype") +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.grid = element_blank())

# Recommended: Save with high height to prevent ID overlap
ggsave("sankey_with_id_labels1.pdf", height = 50, width = 12, dpi = 300)
 

# Add a unique row ID to ensure every single metastasis is treated individually
plot_df <- plot_df %>% mutate(unique_flow_id = row_number())

ggplot(plot_df, aes(axis1 = primary_label, axis2 = met_label, y = 1, alluvium = unique_flow_id)) +
  # 1. Flows (Colored by Met Phenotype)
  geom_alluvium(aes(fill = met_label), width = 1/12, alpha = 0.4) +
  
  # 2. Categorical Boxes
  geom_stratum(width = 1/4, fill = "grey95", color = "grey70") +
  
  # 3. Box Labels
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, fontface = "bold") +
  
  # 4. LEFT TEXT (Tumour ID): Only print when at the first axis
  geom_text(stat = "alluvium", 
            aes(label = ifelse(after_stat(x) == 1, as.character(tumour_id), NA)), 
            x = 1, hjust = 1.2, size = 1.5, na.rm = TRUE) +
  
  # 5. RIGHT TEXT (Region ID): Only print when at the second axis
  geom_text(stat = "alluvium", 
            aes(label = ifelse(after_stat(x) == 2, as.character(regionID), NA)), 
            x = 2, hjust = -0.2, size = 1.5, na.rm = TRUE) +
  
  # Expand x-axis to make room for the labels outside the plot
  scale_x_discrete(limits = c("Primary Status", "Metastatic Status"), 
                   expand = c(.5, .5)) +
  
  scale_fill_manual(values = c("Have Keratin" = "#e74c3c", "No Keratin" = "#3498db")) +
  labs(title = "Lineage Flow with Sample Identification",
       subtitle = "Each flow line represents one distinct Metastatic Region",
       y = "Count of Metastases", 
       fill = "Metastatic Phenotype") +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.text.y = element_blank())

# At height 30-40, this will now be perfectly clear
ggsave("sankey_final_clean.pdf", height = 35, width = 12)




library(dplyr)
library(ggplot2)
library(ggalluvial)
library(tidyr)

# --- SECTION 1: Data Preparation (Keep original logic) ---
primary_summary <- primary_seeding %>%
  filter(seeding_region_label == "Seeding") %>%
  group_by(tumour_id) %>%
  summarise(
    n_kerat = sum(keratinization == "have_keratin", na.rm = TRUE),
    n_non   = sum(keratinization == "no_keratin", na.rm = TRUE),
    primary_label = case_when(
      n_kerat > 0 & n_non == 0 ~ "Primary: Keratin",
      n_non > 0 & n_kerat == 0 ~ "Primary: No Keratin",
      n_kerat > 0 & n_non > 0  ~ "Primary: Mixed",
      TRUE                     ~ "No data"
    )
  )

plot_df <- met_plot %>%
  left_join(primary_summary, by = "tumour_id") %>%
  filter(!is.na(met_keratinization)) %>%
  mutate(
    met_label = ifelse(met_keratinization == "met_have_keratin", "Have Keratin", "No Keratin"),
    primary_label = replace_na(as.character(primary_label), "Not in Primary File")
  )

# --- SECTION 2: Convert to Lodes (Long) Format ---
# This "unrolls" the data so we can label individual axis points accurately.
# It creates columns: 'x' (the axis), 'stratum' (the category), and 'alluvium' (the unique region ID)
plot_lodes <- to_lodes_form(plot_df, 
                            key = "Axis", 
                            value = "Status", 
                            id = "FlowID", 
                            axes = c("primary_label", "met_label"))

# --- SECTION 3: The Final Sankey Plot ---
ggplot(plot_lodes, aes(x = Axis, stratum = Status, alluvium = FlowID, y = 1)) +
  # 1. Main Alluvial Flows (Colored by Metastatic Phenotype)
  geom_alluvium(aes(fill = met_label), width = 1/12, alpha = 0.4) +
  
  # 2. Categorical Boxes (Strata)
  geom_stratum(width = 1/4, fill = "grey95", color = "grey70") +
  
  # 3. Categorical Labels (Inside the boxes)
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, fontface = "bold") +
  
  # 4. LEFT LABELS (Tumour ID)
  # We filter data to ONLY the first axis to prevent double-printing.
  geom_text(stat = "alluvium", 
            data = filter(plot_lodes, Axis == "primary_label"),
            aes(label = tumour_id), 
            hjust = 1.2, size = 1.2, color = "black") +
  
  # 5. RIGHT LABELS (Region ID)
  # We filter data to ONLY the second axis.
  geom_text(stat = "alluvium", 
            data = filter(plot_lodes, Axis == "met_label"),
            aes(label = regionID), 
            hjust = -0.2, size = 1.2, color = "black") +
  
  # Layout and Styling
  scale_x_discrete(limits = c("primary_label", "met_label"), 
                   labels = c("Primary Seeding Status", "Metastatic Status"),
                   expand = c(.5, .5)) +
  scale_fill_manual(values = c("Have Keratin" = "#e74c3c", "No Keratin" = "#3498db")) +
  labs(title = "Lineage Flow: Primary Seeding to Individual Metastatic Regions",
       subtitle = "Every distinct line represents one unique metastatic sample",
       y = "Number of Regions", fill = "Metastatic Status") +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.text.y = element_blank())

# Recommended: High height ensures the 59 region labels have vertical breathing room.
ggsave("sankey_final_distinct_labels.pdf", height = 30, width = 12)




library(dplyr)
library(ggplot2)
library(cowplot)
library(ggfittext)

# 1. Prepare the data: Aggregate counts to match your preferred 'stat="identity"' style
df_counts <- plot_df %>%
  group_by(primary_label, met_label) %>%
  summarise(count = n(), .groups = "drop") 

df_counts1 <- df_counts %>% filter(!primary_label == "Not in Primary File" & !primary_label == "No data")

# 2. Define your preferred color list
color_list <- c("Have Keratin" = "#A6CEE3", "No Keratin" = "#FB9A99")

# 3. Create the plot using your preferred formatting
g1 <- ggplot(df_counts, aes(x = primary_label, y = count, fill = met_label, label = count)) +
  # Use identity and fill, with the black outline from your preferred style
  geom_bar(stat = "identity", position = "fill", color = "black") +
  
  # Add the labels using geom_bar_text
  geom_bar_text(position = "fill", reflow = TRUE) + 
  
  # Apply the manual colors and cowplot theme
  scale_fill_manual(values = color_list, na.value = '#bababa', drop = FALSE) +
  theme_cowplot() +
  
  # Labels and Titles
  ylab("proportion") +
  labs(title = "Proportion of Metastatic Status by Primary Seeding Status") +
  
  # Detailed Theme adjustments to match your preferred style exactly
  theme(
    axis.title.x = element_blank(),  
    axis.title.y = element_text(size = 24),
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1), 
    axis.text.y = element_text(size = 16),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 18)
  )

# Display the plot
print(g1)
ggsave("seeding_vs_met_stacked.pdf")
# USED


# --- SECTION 7: PLOT 4 - Sankey / Alluvial Diagram ---
# This plot now visualizes the flow even for samples with missing primary data
ggplot(plot_df, aes(axis1 = primary_label, axis2 = met_label, y = 1)) +
  geom_alluvium(aes(fill = met_label), width = 1/12, alpha = 0.4) +
  geom_stratum(width = 1/4, fill = "grey95", color = "grey") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("Primary Seeding Status", "Metastatic Status"), expand = c(.1, .1)) +
  scale_fill_manual(values = c("Have Keratin" = "#e74c3c", "No Keratin" = "#3498db")) +
  labs(title = "Lineage Flow: Primary Seeding to Metastasis",
       y = "Number of Tumours", fill = "Metastatic Phenotype") +
  theme_minimal() +
  theme(legend.position = "bottom")



# --- SECTION 7: PLOT 4 - Sankey / Alluvial Diagram with Sample Labels ---
ggplot(plot_df, aes(axis1 = primary_label, axis2 = met_label, y = 1)) +
  # The flows (alluvia)
  geom_alluvium(aes(fill = met_label), width = 1/12, alpha = 0.4) +
  
  # The strata (the boxes/categories)
  geom_stratum(width = 1/4, fill = "grey95", color = "grey") +
  
  # Category labels (e.g., "Primary: Keratin")
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.5, fontface = "bold") +
  
  # NEW: Individual Sample Labels (Tumour ID + Region)
  # We place them in the middle of the flows. check_overlap helps if it's too crowded.
  # geom_text(stat = "alluvium", aes(label = regionID), 
  #           size = 1.8, color = "black", alpha = 0.8) + 
  
  geom_text(stat = "alluvium", 
            Axis = "primary_label",
            aes(label = tumour_id), 
            hjust = 1.2, size = 1.2, color = "black") +
  
  # 5. RIGHT LABELS (Region ID)
  # We filter data to ONLY the second axis.
  geom_text(stat = "alluvium", 
            Axis = "met_label",
            aes(label = regionID), 
            hjust = -0.2, size = 1.2, color = "black") +
  
  
  scale_x_discrete(limits = c("Primary Seeding Status", "Metastatic Status"), expand = c(.15, .15)) +
  scale_fill_manual(values = c("Have Keratin" = "#e74c3c", "No Keratin" = "#3498db")) +
  labs(title = "Lineage Flow: Primary Seeding to Metastasis (by Sample)",
       subtitle = "Labeled by [Tumour ID] ([Metastatic Region])",
       y = "Number of Tumours", 
       fill = "Metastatic Phenotype") +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
#THIS IS ALMOST OK TO USE
# Pro-tip: Save this plot with extra height to see all labels clearly:
 ggsave("sankey_with_labels.pdf", height = 15, width = 10)
 
 
 
 
 
 ggplot(plot_df, aes(axis1 = primary_label, axis2 = met_label, y = 1)) +
   # 1. The flows (alluvia)
   geom_alluvium(aes(fill = met_label), width = 1/12, alpha = 0.4) +
   
   # 2. The strata (the boxes/categories)
   geom_stratum(width = 1/4, fill = "grey95", color = "grey") +
   
   # 3. Category labels (e.g., "Primary: Keratin")
   geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.5, fontface = "bold") +
   
   # 4. LEFT LABELS (Tumour ID)
   geom_text(stat = "alluvium", 
             aes(label = tumour_id), 
             x = 0.9, hjust = 1, size = 1.2, color = "black") +
   
   # 5. NEW: Organ Metastasis Tile
   # Placed right next to Axis 2 (x = 2.08)
   geom_tile(stat = "alluvium", 
             aes(x = 2.08, fill = organ_mets_complete), 
             width = 0.04, height = 0.95) +
   
   # 6. RIGHT LABELS (Region ID)
   # Pushed slightly further right (x = 2.13) so they sit next to the new tile
   geom_text(stat = "alluvium", 
             aes(label = regionID), 
             x = 2.13, hjust = 0, size = 1.2, color = "black") +
   
   # Increased the expansion slightly so the new right-sided elements don't get cut off
   scale_x_discrete(limits = c("Primary Seeding Status", "Metastatic Status"), 
                    expand = c(.25, .25)) +
   
   # COMBINED COLOR SCALES: Put both sets of colors in the same list
   scale_fill_manual(values = c(
     # Original Met Labels
     "Have Keratin" = "#e74c3c", 
     "No Keratin"   = "#3498db",
     
     # New Organ Site Labels
     "lymphnode"    = "#c49a6c", 
     "brain"        = "#f15a29", 
     "pleura"       = "#fbb040",
     "peritoneal"   = "#9e1f63", 
     "lung"         = "#c2b59b",              
     "stump"        = "#00a79d", 
     "adrenal"      = "#f9ed32", 
     "neckLN"       = "#39b54a"
   )) +
   
   labs(title = "Lineage Flow: Primary Seeding to Metastasis (by Sample)",
        subtitle = "Labeled by: [Tumour ID]  ----->  [Organ Site Tile] [Metastatic Region]",
        y = "Number of Tumours", 
        fill = "Phenotype & Organ Site") +
   
   theme_minimal() +
   theme(legend.position = "bottom",
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank())
 
 # Pro-tip: Save this plot with extra height to see all labels and tiles clearly:
 ggsave("sankey_with_labels_and_tiles.pdf", height = 15, width = 10)
 
 
 
 
 # --- SECTION 7: PLOT 4 - Sankey / Alluvial Diagram with Sample Labels ---
 ggplot(plot_df %>% filter(!primary_label == "Not in Primary File" & !primary_label == "No data") , aes(axis1 = primary_label, axis2 = met_label, y = 1)) +
   # The flows (alluvia)
   geom_alluvium(aes(fill = met_label), width = 1/12, alpha = 0.4) +
   
   # The strata (the boxes/categories)
   geom_stratum(width = 1/4, fill = "grey95", color = "grey") +
   
   # Category labels (e.g., "Primary: Keratin")
   geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.5, fontface = "bold") +
   
   # NEW: Individual Sample Labels (Tumour ID + Region)
   # We place them in the middle of the flows. check_overlap helps if it's too crowded.
   # geom_text(stat = "alluvium", aes(label = regionID), 
   #           size = 1.8, color = "black", alpha = 0.8) + 
   
   geom_text(stat = "alluvium", 
             Axis = "primary_label",
             aes(label = tumour_id), 
             hjust = 1.2, size = 1.2, color = "black") +
   
   # 5. RIGHT LABELS (Region ID)
   # We filter data to ONLY the second axis.
   geom_text(stat = "alluvium", 
             Axis = "met_label",
             aes(label = organ_mets_complete, color = organ_mets_complete), 
             hjust = -0.2, size = 1.2) +
   
   scale_color_manual(values = c(
     "lymphnode"  = "#c49a6c", 
     "brain"      = "#f15a29", 
     "pleura"     = "#fbb040",
     "peritoneal" = "#9e1f63", 
     "lung"       = "#c2b59b",              
     "stump"      = "#00a79d", 
     "adrenal"    = "#f9ed32", 
     "neckLN"     = "#39b54a"
   )) +
   
   scale_x_discrete(limits = c("Primary Seeding Status", "Metastatic Status"), expand = c(.15, .15)) +
   scale_fill_manual(values = c("Have Keratin" = "#e74c3c", "No Keratin" = "#3498db")) +
   labs(title = "Lineage Flow: Primary Seeding to Metastasis (by Sample)",
        subtitle = "Labeled by [Tumour ID] ([Metastatic Region])",
        y = "Number of Tumours", 
        fill = "Metastatic Phenotype") +
   theme_minimal() +
   theme(legend.position = "bottom",
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank())
 #THIS IS ALMOST OK TO USE
 # Pro-tip: Save this plot with extra height to see all labels clearly:
 ggsave("sankey_with_labels_FINAL.pdf", height = 10, width = 7)
 
 
## plots for timing 
 
squamous_data <- read_csv("../df-alvin/500_entire_tumour_with_grading.csv")
entire_tumour_data <- readRDS("../data/Mets_paper_20221109_TRACERx421_all_tumour_df.rds") %>% dplyr::select(tumour_id_muttable_cruk , Histology_per_tumour_id_muttable)

load("../data/timingTable.rda") 
met_full <- met_keratin_plot %>% left_join(squamous_data[, c("cruk_tumourid", "keratinization", "tumour_budding", "single_cell_invasion")] , by = c("patient_id" = "cruk_tumourid")) %>% 
  left_join(timingTable, by = c("patient_id" = "tumour_id"))
table(met_full[, c("keratinization", "met_keratinization", "timing")])
table(met_full[, c("keratinization", "timing")])
 
timing_plot <- timingTable %>% left_join(entire_tumour_data, by = c("tumour_id" = "tumour_id_muttable_cruk"))
table(timing_plot[, c("Histology_per_tumour_id_muttable", "timing")]) 

#proportion of k in met case in our cohort
(met_full[, c("patient_id", "keratinization", "timing")]) %>% unique() %>% select(-patient_id) %>% table()
(met_full[, c("patient_id", "keratinization", "timing")]) %>% unique() %>% select(-patient_id, -timing) %>% table()
fisher.test(matrix(c(4,12,11,11), nrow = 2))

#proportion of budding in met case in our cohort
(met_full[, c("patient_id", "tumour_budding", "timing")]) %>% unique() %>% select(-patient_id) %>% table()
(met_full[, c("patient_id", "single_cell_invasion", "timing")]) %>% unique() %>% select(-patient_id) %>% table()
fisher.test(matrix(c(9,11,6,12), nrow = 2))
fisher.test(matrix(c(4,9,11,14), nrow = 2))

# plot 1
# 1. Prepare the dataset based on your unique patient pipeline
input_plot_df <- met_full %>% 
  dplyr::select(patient_id, keratinization, timing) %>% 
  unique()

# 2. Generate the table (swapping rows/columns to match template logic)
table1 <- table(input_plot_df$timing, input_plot_df$keratinization)

# 3. Restructure into a data frame matching the template layout
df_gathered1 <- data.frame(
  var1 = c(rownames(table1)[1], rownames(table1)[1], rownames(table1)[2], rownames(table1)[2]),
  var2 = c(colnames(table1)[1], colnames(table1)[2], colnames(table1)[1], colnames(table1)[2]),
  count = c(table1[1], table1[3], table1[2], table1[4])
)

# 4. Set factor levels for correct ordering
df_gathered1$var2 <- factor(df_gathered1$var2, levels = c("have_keratin", "no_keratin"))
df_gathered1$var1 <- factor(df_gathered1$var1, levels = c("early", "late"))

# 5. Run the Fisher-exact test
fisher_p <- fisher.test(matrix(c(table1[1], table1[2], table1[3], table1[4]), nrow = 2))

# 6. Define the color palette updated for your new variables
color_list <- c(
  "early"       = "#A6CEE3",
  "greencolour" = "#B2DF8A",
  "late"        = "#FB9A99"
)

# 7. Generate the stacked bar chart
g1 <- ggplot(data = df_gathered1, aes(x = var2, y = count, fill = var1, label = count)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  geom_bar_text(position = "fill", reflow = TRUE) + 
  scale_fill_manual(values = color_list, na.value = '#bababa', drop = F) +
  theme_cowplot() +
  ylab("proportion") +
  ggtitle(paste0(" Fisher-exact test p-value: ", format(round(as.numeric(fisher_p$p.value), 4), nsmall = 4))) +
  theme(
    axis.title.x = element_blank(),  
    axis.title.y = element_text(size = 24),
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1), 
    axis.text.y = element_text(size = 16),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 18)
  )

print(g1)
ggsave("../results/divergencetiming_keratinization.pdf")



# plot 1
# 1. Prepare the dataset based on your unique patient pipeline
input_plot_df <- met_full %>% 
  dplyr::select(patient_id, tumour_budding, timing) %>% 
  unique()

# 2. Generate the table (swapping rows/columns to match template logic)
table1 <- table(input_plot_df$timing, input_plot_df$tumour_budding)

# 3. Restructure into a data frame matching the template layout
df_gathered1 <- data.frame(
  var1 = c(rownames(table1)[1], rownames(table1)[1], rownames(table1)[2], rownames(table1)[2]),
  var2 = c(colnames(table1)[1], colnames(table1)[2], colnames(table1)[1], colnames(table1)[2]),
  count = c(table1[1], table1[3], table1[2], table1[4])
)

# 4. Set factor levels for correct ordering
df_gathered1$var2 <- factor(df_gathered1$var2, levels = c("have_budding", "no_budding"))
df_gathered1$var1 <- factor(df_gathered1$var1, levels = c("early", "late"))

# 5. Run the Fisher-exact test
fisher_p <- fisher.test(matrix(c(table1[1], table1[2], table1[3], table1[4]), nrow = 2))

# 6. Define the color palette updated for your new variables
color_list <- c(
  "no_budding"       = "#A6CEE3",
  "greencolour" = "#B2DF8A",
  "have_budding"        = "#FB9A99"
)

# 7. Generate the stacked bar chart
g1 <- ggplot(data = df_gathered1, aes(x = var1, y = count, fill = var2, label = count)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  geom_bar_text(position = "fill", reflow = TRUE) + 
  scale_fill_manual(values = color_list, na.value = '#bababa', drop = F) +
  theme_cowplot() +
  ylab("proportion") +
  ggtitle(paste0(" Fisher-exact test p-value: ", format(round(as.numeric(fisher_p$p.value), 4), nsmall = 4))) +
  theme(
    axis.title.x = element_blank(),  
    axis.title.y = element_text(size = 24),
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1), 
    axis.text.y = element_text(size = 16),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 18)
  )

print(g1)
ggsave("../results/divergencetiming_budding.pdf", height = 8, width = 10)


# plot 2

# 1. Filter and group the dataset based on the new histology logic
input_plot_df <- timing_plot %>% 
  dplyr::select(Histology_per_tumour_id_muttable, timing) %>% 
  dplyr::filter(Histology_per_tumour_id_muttable != "Collision LUAD and LUSC") %>%
  dplyr::mutate(
    histology_group = ifelse(Histology_per_tumour_id_muttable == "Squamous cell carcinoma", 
                             "LUSC histology", 
                             "non-LUSC histology")
  )

# 2. Generate the table 
table1 <- table(input_plot_df$timing, input_plot_df$histology_group)

# 3. Restructure into a data frame matching the template layout
df_gathered1 <- data.frame(
  var1 = c(rownames(table1)[1], rownames(table1)[1], rownames(table1)[2], rownames(table1)[2]),
  var2 = c(colnames(table1)[1], colnames(table1)[2], colnames(table1)[1], colnames(table1)[2]),
  count = c(table1[1], table1[3], table1[2], table1[4])
)

# 4. Set factor levels for correct ordering
df_gathered1$var2 <- factor(df_gathered1$var2, levels = c("LUSC histology", "non-LUSC histology"))
df_gathered1$var1 <- factor(df_gathered1$var1, levels = c("early", "late"))

# 5. Run the Fisher-exact test
fisher_p <- fisher.test(matrix(c(table1[1], table1[2], table1[3], table1[4]), nrow = 2))

# 6. Define the color palette updated for your fill variables
color_list <- c(
  "early"       = "#A6CEE3",
  "greencolour" = "#B2DF8A",
  "late"        = "#FB9A99"
)

# 7. Generate the stacked bar chart
g2 <- ggplot(data = df_gathered1, aes(x = var2, y = count, fill = var1, label = count)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  geom_bar_text(position = "fill", reflow = TRUE) + 
  scale_fill_manual(values = color_list, na.value = '#bababa', drop = F) +
  theme_cowplot() +
  ylab("proportion") +
  ggtitle(paste0(" Fisher-exact test p-value: ", format(round(as.numeric(fisher_p$p.value), 4), nsmall = 4))) +
  theme(
    axis.title.x = element_blank(),  
    axis.title.y = element_text(size = 24),
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1), 
    axis.text.y = element_text(size = 16),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 18)
  )

print(g2)
ggsave("../results/divergencetiming_histologies.pdf")



squamous_data_entire <- read_csv("../df-alvin/500_entire_tumour_with_grading.csv") %>% select(cruk_tumourid, keratinization)
met_table_full <- met_keratin_plot %>% left_join(squamous_data_entire, by = c("tumour_id" = "cruk_tumourid"))

FINAL_MET_TABLE <- met_table_full %>%
  group_by(tumour_id, met_keratinization) %>%
  summarise(count = n(), .groups = "drop") %>%
  # Pivot the unique met_keratinization values into individual columns
  pivot_wider(
    names_from = met_keratinization, 
    values_from = count, 
    values_fill = 0
  ) %>%
  # Rename columns to match your exact naming scheme
  rename(
    met_have_keratin_count = met_have_keratin,
    met_no_keratin_count = met_no_keratin
  ) %>%
  # Ensure strict column order layout
  select(tumour_id, met_have_keratin_count, met_no_keratin_count) %>% 
  left_join(squamous_data_entire, by = c("tumour_id" = "cruk_tumourid")) %>% 
  mutate(met_status = case_when(met_have_keratin_count > 0 ~ "at_least_one_keratinized",
                                met_no_keratin_count > 0 & met_have_keratin_count == 0 ~ "non_keratinized")) %>% 
  left_join(timingTable, by = c("tumour_id"))

#proportion of k in met case in TCGA cohort.
histology_TCGA <- read.csv("../df_TCGA/TCGA_histology_clean.csv") 
danaher_TCGA <- read.csv("../df_TCGA/tcga_lusc_danaher_tumor.csv") 
GO_TCGA <- read.csv("../df_TCGA/tcga_lusc_expression_gsva_GO.csv") 
mcp_TCGA <- read.csv("../df_TCGA/tcga_lusc_mcp_tumor.csv") 
prognosis_TCGA <- read.csv("../df_TCGA/TCGA_prognosis.csv", 
                           na.strings = c("#N/A", ""))

# > GO_TCGA$GO_patient_ID[duplicated(GO_TCGA$GO_patient_ID)]
# [1] "TCGA-22-1002" "TCGA-22-1005" "TCGA-22-1016" "TCGA-21-1076" "TCGA-22-1017" "TCGA-18-4083"
# [7] "TCGA-21-1070" "TCGA-21-1071" "TCGA-21-1082" "TCGA-21-1081"
# > mcp_TCGA$MCP_patient_ID[duplicated(mcp_TCGA$MCP_patient_ID)]
# [1] "TCGA-22-1002" "TCGA-22-1005" "TCGA-22-1016" "TCGA-21-1076" "TCGA-22-1017" "TCGA-18-4083"
# [7] "TCGA-21-1070" "TCGA-21-1071" "TCGA-21-1082" "TCGA-21-1081"

duplicate_ids <- GO_TCGA$GO_patient_ID[duplicated(GO_TCGA$GO_patient_ID)]

histology_TCGA$slide_patient[duplicated(histology_TCGA$slide_patient)]

#Average all numeric columns for duplicated patients
unique_GO_TCGA <- GO_TCGA %>%
  group_by(GO_patient_ID) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) %>%
  ungroup()

unique_mcp_TCGA <- mcp_TCGA %>%
  group_by(MCP_patient_ID) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) %>%
  ungroup()

unique_danaher_TCGA <- danaher_TCGA %>%
  group_by(patient_ID_danaher) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) %>%
  ungroup()

# GO_TCGA[GO_TCGA$GO_patient_ID %in% duplicate_ids, ] 
# mcp_TCGA[mcp_TCGA$MCP_patient_ID %in% duplicate_ids, ] 
# unique_mcp_TCGA[unique_mcp_TCGA$MCP_patient_ID %in% duplicate_ids, ]

#joining dfs
TCGA_df <- histology_TCGA %>%
  left_join(unique_danaher_TCGA,   by = c("slide_patient" = "patient_ID_danaher")) %>%
  left_join(prognosis_TCGA, by = c("slide_patient" = "prognosis_patient_barcode")) %>% 
  left_join(unique_GO_TCGA,        by = c("slide_patient" = "GO_patient_ID")) %>% 
  left_join(unique_mcp_TCGA,       by = c("slide_patient" = "MCP_patient_ID")) 

TCGA_df$slide_patient[duplicated(TCGA_df$slide_patient)]
# character(0)

TCGA_df$keratinization %>% table()  #459
TCGA_df$tumour_budding %>% table()  #455
#only cases with both budding and keratinization are used n = 455

TCGA_df %>% filter((tumour_budding != "")) %>% select(c("keratinization", "tumour_budding")) %>% table()
# 1. Prepare the dataset based on your new filtering logic
input_plot_df <- TCGA_df %>% 
  dplyr::filter(tumour_budding != "") %>% 
  dplyr::select(keratinization, tumour_budding)

# 2. Generate the table (oriented to match the template's layout)
table1 <- table(input_plot_df$tumour_budding, input_plot_df$keratinization)

# 3. Restructure into a data frame matching the template layout
df_gathered1 <- data.frame(
  var1 = c(rownames(table1)[1], rownames(table1)[1], rownames(table1)[2], rownames(table1)[2]),
  var2 = c(colnames(table1)[1], colnames(table1)[2], colnames(table1)[1], colnames(table1)[2]),
  count = c(table1[1], table1[3], table1[2], table1[4])
)

# 4. Set factor levels for correct ordering
df_gathered1$var2 <- factor(df_gathered1$var2, levels = c("have_keratin", "no_keratin"))
df_gathered1$var1 <- factor(df_gathered1$var1, levels = c("no_budding", "have_budding"))

# 5. Run the Fisher-exact test
fisher_p <- fisher.test(matrix(c(table1[1], table1[2], table1[3], table1[4]), nrow = 2))

# 6. Define the color palette matching the original budding template
color_list <- c(
  "no_budding"   = "#A6CEE3",
  "greencolour"  = "#B2DF8A",
  "have_budding" = "#FB9A99"
)

# 7. Generate the stacked bar chart
g1 <- ggplot(data = df_gathered1, aes(x = var2, y = count, fill = var1, label = count)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  geom_bar_text(position = "fill", reflow = TRUE) + 
  scale_fill_manual(values = color_list, na.value = '#bababa', drop = F) +
  theme_cowplot() +
  ylab("proportion") +
  ggtitle(paste0(" Fisher-exact test p-value: ", format(round(as.numeric(fisher_p$p.value), 4), nsmall = 4))) +
  theme(
    axis.title.x = element_blank(),  
    axis.title.y = element_text(size = 24),
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1), 
    axis.text.y = element_text(size = 16),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 18)
  )

print(g1)
ggsave("../results/TCGA_budding_keratinization.pdf")

