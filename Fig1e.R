##Script for Tx421 LUSC paper
##Fig 1b

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
all_identifiers <- vector()

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
        
        pairwise_identifier <- paste0(tum, "-", sample_list[i], sample_list[j])
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
        all_identifiers <- c(all_identifiers, pairwise_identifier)
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

##############Figure 1g
####assess total number of truncal LOH events for every tumour
######################################################

#bind all df of one single tumour and then only keep the rows with all values being 1 and no 0

cruk_tumour_all <- paste0(sub("-.*", "", all_identifiers))


summary_truncal_LOH <- data.frame(cruk_tumour = c("testing"), 
                                  truncal_LOH = c(999), 
                                  number_with_NA = c(999))
for(patient in cruk_tumour_all){
print(patient)
cruk_tumour_grep <- paste0(patient, "-")
df_identifiers <- all_identifiers[grep(cruk_tumour_grep, all_identifiers)]
df_to_bind <- sapply(df_identifiers, function(x){return(joined_df_list[[x]][, c(2,4)])}) 
tmp2 <- do.call(cbind , df_to_bind)

# cbind(joined_df_list[["CRUK0005R1R2"]][, c(2,4)], joined_df_list[["CRUK0005R1R3"]][, c(2,4)]) -> tmp2
tmp3 <- tmp2[apply(tmp2!=0, 1, all),]  #row with no 0
# tmp4 <- tmp2[apply(tmp2==0, 1, any),]  #row with any 0
tmp4 <- tmp2[apply(is.na(tmp2), 1, any),]  # row with any NA
patient_LOH_info <- c(patient, as.numeric(nrow(tmp3)), as.numeric(nrow(tmp4)))  # nrow(tmp3) = number of truncal LOH for that tumour.; nrow(tmp4) = number of rows with any NA for that tumour
summary_truncal_LOH <- rbind(summary_truncal_LOH, patient_LOH_info)
}

summary_truncal_LOH1 <- summary_truncal_LOH[!duplicated(summary_truncal_LOH), ] %>% filter(cruk_tumour != "testing")
summary_truncal_LOH1$truncal_LOH <- as.numeric(summary_truncal_LOH1$truncal_LOH)
summary_truncal_LOH1$number_with_NA <- as.numeric(summary_truncal_LOH1$number_with_NA)



##############plot
######################################################
#boxplot for all regions of poorly-diff LUSC vs all regions of well-diff LUSC
grading_by_entire_tumour_data <- read_csv("../df-alvin/500_entire_tumour_with_grading.csv")
LOH_truncal_count <- left_join(summary_truncal_LOH1, grading_by_entire_tumour_data[, c("cruk_tumourid", "keratinization")], by = c("cruk_tumour" = "cruk_tumourid")) %>% 
        filter(!is.na(keratinization))
#only 127 cases of 136 cases have available mphase_cytoband data

color_list <-  c("have_keratin" = "#A6CEE3",
                 "greencolour" = "#B2DF8A",
                 "no_keratin"= "#FB9A99")

compare_means(truncal_LOH ~ keratinization, data = LOH_truncal_count, method = "wilcox.test")
ggboxplot(LOH_truncal_count, width=0.5,lwd=1, x = "keratinization", y = "truncal_LOH", fill = "keratinization", add = "jitter") + 
  stat_compare_means(comparisons = list( c("have_keratin", "no_keratin")), label.y = c(450)) +
  theme_cowplot() +
  scale_fill_manual(values = color_list, na.value = '#bababa', drop = F) +
  theme(  axis.title.y = element_text(size=26),
          axis.text.y = element_text(size=26),
          axis.title.x = element_text(size=26),
          axis.text.x = element_text(size=26),
          legend.title = element_blank()) 
#Fig 1e

#n number
LOH_truncal_count[LOH_truncal_count$keratinization == "have_keratin", ] %>% nrow()
LOH_truncal_count[LOH_truncal_count$keratinization == "no_keratin", ] %>% nrow()
table(LOH_truncal_count$keratinization)

