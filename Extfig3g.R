##Script for Tx421 LUSC paper
##by Alvin Cheung 20230228
##Acknowledgement to Tak Karasaki (takahiro.karasaki@crick.ac.uk) for previous script versions




#################################
### libraries and directories ###
#################################

library(tidyverse)
library(data.table)
library(readxl)
library(fst)
library(cowplot)
library(ggalluvial)
library(ggpubr)
library(nlme)
# setwd("")
workdir <- "../results/Supp2.LOH_count_across_chromosomes_for_ancestral_descendant_regions/" 
dir.create(workdir)

data_dir <- "../data/" 
dir.create(data_dir)

tmpdir <- paste0(workdir, "/tmp/")
dir.create(tmpdir)

#########################
### utility functions ###
#########################

# functions for growth pattern annotation
pattern_to_grade <- function(x){
  y = ifelse(x %in% c("have_keratin"),"low",ifelse(x %in% c("N.A."),"n.a.",ifelse(x %in% c("no_keratin"),"high","n.a.")))
  return(y)
}

pattern_to_grade.num <- function(x){
  y = ifelse(x %in% c("have_keratin"),1,ifelse(x %in% c("N.A."),NA_real_,ifelse(x %in% c("no_keratin"),3,NA_real_)))
  return(as.numeric(y))
}


# function to calculate euclidian distance of 2 vectors
calc_euc_dis <- function(x, y){sqrt(sum((x - y)^2)) }

# functions to add purity and add growth pattern for regional pairs 

add_purity <- function(df = df, # df : tumour_id (<- tumour_id_muttable), parent_region (e.g."region.1st" ), Region.1st (e.g. "R1"), Region.2nd
                       purity_df = purity_df){  #  purity_df : tumour_id (<- tumour_id_muttable), region (e.g. "R1), "purity"
  df_new <- df %>%
    mutate(region = Region.1st) %>%
    left_join(purity_df[,c("tumour_id", "region", "purity")], by = c("tumour_id","region")) %>%
    dplyr::rename(purity_1st = purity) %>%
    mutate(region = Region.2nd) %>%
    left_join(purity_df[,c("tumour_id", "region", "purity")], by = c("tumour_id","region")) %>%
    dplyr::rename(purity_2nd = purity)
  return(df_new)
}

add_growth_pattern <- function(df = df, # df : tumour_id (<- tumour_id_muttable), parent_region (e.g."region.1st" ), Region.1st (e.g. "R1"), Region.2nd
                               pattern_df = pattern_df){  #  num_gd_df : tumour_id (<- tumour_id_muttable), region (e.g. "R1), "num_gds"
  df_new <- df %>%
    mutate(region = Region.1st) %>%
    left_join(pattern_df[,c("tumour_id", "region", "regional_pattern.new")], by = c("tumour_id","region")) %>%
    dplyr::rename(pattern_1st = regional_pattern.new) %>%
    mutate(region = Region.2nd) %>%
    left_join(pattern_df[,c("tumour_id", "region", "regional_pattern.new")], by = c("tumour_id","region")) %>%
    dplyr::rename(pattern_2nd = regional_pattern.new) %>%
    mutate(grade_1st = pattern_to_grade(pattern_1st),
           grade_2nd = pattern_to_grade(pattern_2nd),
           grade.num_1st = pattern_to_grade.num(pattern_1st),
           grade.num_2nd = pattern_to_grade.num(pattern_2nd),)
  return(df_new)
}


#######################
### utility objects ###
#######################

# list of subtypes
subtype_list <- c("have_keratin","no_keratin")


################
## load data  ## 
################

mphase_cyto_df <- readRDS(paste0(data_dir, "20221108_mphase_cyto_df.RDS"))
mphase_arm_df <- readRDS(paste0(data_dir, "20221108_mphase_arm_df.RDS"))

sequenza_cyto_df <- readRDS(paste0(data_dir, "20221108_sequenza_cyto_df.RDS"))
sequenza_arm_df <- readRDS(paste0(data_dir, "20221108_sequenza_arm_df.RDS"))

trees <- readRDS(paste0(data_dir, "20221109_TRACERx421_phylogenetic_trees.rds"))

muttable_df <- read.fst(paste0(data_dir,"20221109_TRACERx421_mutation_table.fst"),
                        as.data.table = T)
muttable_region <- read.fst(paste0(data_dir,"20221123_TRACERx421_mutation_table_region.fst"),
                            as.data.table = T)

#LUAD_region_df_full <- readRDS(paste0(data_dir,"20221109_LUAD_input_plot_df_all_LUAD_region_df_full.RDS"))

#Alvin: regional SCC data, needs to have regional_pattern.new according to the pathologic feature to test
LUSC_region_df_full <- read.csv("../df-alvin/All_regions_500_regions-SCC-reviewed.csv")
LUSC_region_df_full <- LUSC_region_df_full %>%  
  mutate(regional_pattern.new = case_when(keratinization == "have_keratin" ~ "have_keratin", 
                                          keratinization == "no_keratin" ~ "no_keratin"))

#Alvin: Select histology here

LUSC_region_df <- LUSC_region_df_full  %>%
  filter(!is.na(Histology_per_tumour_id_muttable)) %>% 
  filter(Histology_per_tumour_id_muttable == "Squamous cell carcinoma") %>% 
  filter(included_500_regions == 1)

###############################
## organise regional pattern ##
###############################

regional_df <- LUSC_region_df %>%
  filter(regional_pattern.new %in% subtype_list) %>%
  arrange(sample_name_cruk)

## preparation to annotate growth pattern 
regional_df_join1 <- regional_df %>%
  dplyr::select(tumour_id_muttable_cruk, region, regional_pattern.new) %>%
  dplyr::rename(tumour_id = tumour_id_muttable_cruk,
                Region.1st = region)
regional_df_join2 <- regional_df %>%
  dplyr::select(tumour_id_muttable_cruk, region, regional_pattern.new) %>%
  dplyr::rename(tumour_id = tumour_id_muttable_cruk,
                Region.2nd = region)

# list of LUSC tumour_id (tumours with at least one regional pattern annotation). IMA included.
tumour_id_list <- unique(regional_df$tumour_id_muttable_cruk)
LUSC_tumour_id <- tumour_id_list 



##############################
##############################
## Step 1: data preparation ##
##############################
##############################

##############################################
## shared & private LOH : mphase (cytoband) ##
##############################################

in_dt <- mphase_cyto_df %>% 
  mutate(sample = sample_name_cruk) %>%
  distinct(sample, cytoband_name, .keep_all = T) %>%
  mutate(region = str_split_fixed(sample, pattern = "-", n =2)[,2])

tumour_list <- as.character(na.omit(unique(in_dt$tumour_id_muttable_cruk) )) # including non-LUSCs, CRUK0555_LN1 is missing but no problem (single region tumour)
num_sample <- length(unique(in_dt$sample)) 

## now calculate euc dist and shared & private LOH for all regional pairs 
# consider parallel LOH as private events 

all_vec <- vector() # all results will be saved in this vector (will be a matrix in the end)
joined_df_list <- list()  #Alvin: initialize an array of df
paired_results_list <- list()#Alvin: initialize an array of df

for (tum in tumour_list){
  print(tum)
  
  dt_tum <- in_dt %>%
    filter(tumour_id_muttable_cruk  == tum) 
  
  sample_list <- sort(unique(dt_tum$region))
  sample_num <- length(sample_list)
  
  
  # calculate distance between all regional pairs 
  
  if(sample_num > 1){  # sample_num : same as number of regions in each tumour
    for (i in 1:(sample_num - 1)){
      for (j in (1+i):sample_num){
        #i=1
        #j=2
        dt_i <- dt_tum %>%
          filter(region == sample_list[i]) %>%
          dplyr::select(cytoband_name,  loh_status, dominant_allele)
        dt_j <- dt_tum %>%
          filter(region == sample_list[j]) %>%
          dplyr::select(cytoband_name, loh_status, dominant_allele)
        
        pairwise_identifier <- paste0(tum, sample_list[i], sample_list[j])
        print(pairwise_identifier)     #Alvin: create array of df with index as pairwise_identifer
        
        # this time focus on LOH only 
        joined_df <- dt_i %>%
          inner_join(dt_j, by = "cytoband_name") %>%
          filter(!is.na(loh_status.x), !is.na(loh_status.y),!is.na(dominant_allele.x), !is.na(dominant_allele.y) ) 
        
        # use categorical loh status
        euc_dist_loh <- calc_euc_dis(joined_df$loh_status.x, joined_df$loh_status.y) 
        
        # shared / private loh 
        shared_loh <- nrow(joined_df[joined_df$loh_status.x == "1" & joined_df$loh_status.y == "1" & joined_df$dominant_allele.x == joined_df$dominant_allele.y,  ]) # parallel LOH (joined_df$dominant_allele.x != joined_df$dominant_allele.y) are not counted as shared LOH
        private_loh_x <- nrow(joined_df[(joined_df$loh_status.x == "1" & joined_df$loh_status.y == "0") | (joined_df$loh_status.x == "1" & joined_df$loh_status.y == "1" & joined_df$dominant_allele.x != joined_df$dominant_allele.y), ])  # any parallel LOH will be counted as private LOH in both regions
        private_loh_y <- nrow(joined_df[(joined_df$loh_status.x == "0" & joined_df$loh_status.y == "1") | (joined_df$loh_status.x == "1" & joined_df$loh_status.y == "1" & joined_df$dominant_allele.x != joined_df$dominant_allele.y), ])  # -> this means that any regional pairs with parallel LOH are treated as "siblings" in the downstream analysis
        
        paired_results <- c(tum, sample_list[i],sample_list[j],
                            euc_dist_loh, 
                            shared_loh,
                            private_loh_x, private_loh_y)
        
        joined_df_list[[pairwise_identifier]] <- joined_df   #Alvin: create array of df
        paired_results_list[[pairwise_identifier]] <- paired_results  #Alvin: create array of paired results
        all_vec <- rbind(all_vec,   paired_results )
      }
    }
  }
}

colnames(all_vec) <- c("tumour_id","Region.1st","Region.2nd",
                       "euc_dist_cyto_loh",
                       "shared_loh",
                       "private_loh_1st", "private_loh_2nd")    

all_res_df <- as.data.frame(all_vec, stringsAsFactors = F)  %>%
  mutate(compaired_regions = paste0(Region.1st, "-",Region.2nd))


##############Figure 2c
####assess total number of LOH events for every region
######################################################
LOH_count_for_every_region <- in_dt %>% group_by(Region_ID_cruk, loh_status) %>% summarize(total_LOH = n()) %>% ungroup
all_lusc_regions <- LUSC_region_df$Region_ID_cruk

#LUAD region list
all_luad_regions <- read_csv("../df-alvin/LUAD-all-regions-id.csv") %>% select(region_id) %>% unlist() 

LOH_count_for_every_region %>% mutate(histology = case_when(Region_ID_cruk %in% all_luad_regions ~ "LUAD", 
                                                            Region_ID_cruk %in% all_lusc_regions ~ "LUSC")) -> 
                        LOH_count_for_every_region1


##############Figure 2c
#plot figure 2c. LOH counts for LUAD Vs LUSC

LOH_count_for_every_region1 %>% filter(!is.na(histology)) %>% 
  filter(loh_status == 1) ->
  LOH_count_histology_plot

color_list <-  c("LUAD" = "#8DC63F",
                 "greencolour" = "#053061",
                 "LUSC"= "#FB9A99")

compare_means(total_LOH ~ histology, data = LOH_count_histology_plot, method = "wilcox.test")
ggboxplot(LOH_count_histology_plot, width=0.5,lwd=1, x = "histology", y = "total_LOH", fill = "histology", ) + 
  geom_jitter(width = 0.15, alpha = 0.4) +
  stat_compare_means(comparisons = list( c("LUAD", "LUSC")), label.y = c(450)) +
  theme_cowplot() +
  scale_fill_manual(values = color_list, na.value = '#bababa', drop = F) +
  scale_x_discrete(labels = c("LUAD", "LUSC")) +
  theme(  axis.title.y = element_text(size=26),
          axis.text.y = element_text(size=26),
          axis.title.x = element_text(size=26),
          axis.text.x = element_text(size=26),
          legend.position = "none",
          legend.title = element_blank()) 
# Extfig 3g - p value need to be modified by region correction below
#n number
LOH_count_histology_plot[LOH_count_histology_plot$histology == "LUAD", ] %>% nrow()
LOH_count_histology_plot[LOH_count_histology_plot$histology == "LUSC", ] %>% nrow()
table(LOH_count_histology_plot$histology)

# multi-region correction
LOH_count_histology_plot <- LOH_count_histology_plot %>% mutate(tumour_id = substr(Region_ID_cruk, 1, 11))
library(nlme)
effect_model <- lme(total_LOH ~ histology, random = ~ 1 |tumour_id, data = LOH_count_histology_plot, method = "ML")
null_model <- lme(total_LOH ~ 1, random = ~ 1 | tumour_id, data = LOH_count_histology_plot, method = "ML")
anova_res <- anova(null_model, effect_model)
print(anova_res["effect_model",]$`p-value`)


