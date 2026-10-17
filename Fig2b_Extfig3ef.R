#Script for Tx421 LUSC paper
#Fig 2b
#Ext fig 3e, 3f
#Ancestral descendant analysis


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
workdir <- "../results/Fig.2f_keratinization_ancestraldescendant_assignment/" 
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
  y = ifelse(x %in% c("well_dx"),"low",ifelse(x %in% c("N.A."),"n.a.",ifelse(x %in% c("poor_dx"),"high","n.a.")))
  return(y)
}

pattern_to_grade.num <- function(x){
  y = ifelse(x %in% c("well_dx"),1,ifelse(x %in% c("N.A."),NA_real_,ifelse(x %in% c("poor_dx"),3,NA_real_)))
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
subtype_list <- c("well_dx","poor_dx")


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

#Alvin: Modify input file here
LUSC_region_df_full <- read_csv("../df-alvin/All_regions_500_regions-SCC_reviewed_With_mean_keratinMitosis.csv")

#Alvin: Select histology here

LUSC_region_df <- LUSC_region_df_full  %>%
  filter(!is.na(Histology_per_tumour_id_muttable)) %>% 
  filter(Histology_per_tumour_id_muttable == "Squamous cell carcinoma") %>% 
  filter(included_500_regions == 1)

#Alvin: Add regional pattern here

LUSC_region_df <- LUSC_region_df %>% 
  mutate(regional_pattern.new = case_when(keratinization == "have_keratin" ~ "well_dx", 
                                          keratinization == "no_keratin" ~ "poor_dx"))

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

# out_f <- paste0(tmpdir, "/LOH_all.csv")
# write.csv(all_res_df, out_f, row.names = F)

# read 
#all_res_df <- fread( paste0(data_dir, "/LOH_all.csv"))

## annotate growth pattern 

LOH_annotated <- all_res_df %>%
  mutate(shared_loh = as.numeric(shared_loh),
         private_loh_1st = as.numeric(private_loh_1st),
         private_loh_2nd = as.numeric(private_loh_2nd),
         euc_dist_cyto_loh = as.numeric(euc_dist_cyto_loh)) %>%
  left_join(regional_df_join1, by = c("tumour_id", "Region.1st") ) %>%
  left_join(regional_df_join2, by = c("tumour_id", "Region.2nd") ) %>%
  dplyr::rename(pattern_1st = regional_pattern.new.x,
                pattern_2nd = regional_pattern.new.y) %>%
  mutate(grade_1st = pattern_to_grade(pattern_1st),
         grade_2nd = pattern_to_grade(pattern_2nd),
         grade_num_1st = pattern_to_grade.num(pattern_1st),
         grade_num_2nd = pattern_to_grade.num(pattern_2nd))

# save (load) LOH table (cytoband)

out_f <- paste0(tmpdir, "/LOH_annotated.rds")
# saveRDS(LOH_annotated, out_f)


#########################################
## shared & private LOH : mphase (arm) ##
#########################################

in_dt <- mphase_arm_df %>%
  mutate(sample = sample_name_cruk) %>%
  distinct(sample, chr_arm, .keep_all = T) %>%
  mutate(region = str_split_fixed(sample, pattern = "-", n =2)[,2])

tumour_list <- as.character(na.omit(unique(in_dt$tumour_id_muttable_cruk) )) # including non-LUSCs, CRUK0555_LN1 is missing but no problem (single region tumour)
num_sample <- length(unique(in_dt$sample)) 

## now calculate euc dist and shared & private LOH for all regional pairs 

all_vec <- vector() # all results will be saved in this vector (will be a matrix in the end)
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
          dplyr::select(chr_arm,  loh_status)
        dt_j <- dt_tum %>%
          filter(region == sample_list[j]) %>%
          dplyr::select(chr_arm,  loh_status)
        
        # this time focus on LOH only 
        joined_df <- dt_i %>%
          inner_join(dt_j, by = "chr_arm") %>%
          filter(!is.na(loh_status.x), !is.na(loh_status.y) ) 
        
        # use categorical loh status
        euc_dist_loh <- calc_euc_dis(joined_df$loh_status.x, joined_df$loh_status.y) 
        
        # shared / private loh 
        shared_loh <- nrow(joined_df[joined_df$loh_status.x == "1" & joined_df$loh_status.y == "1", ])
        private_loh_x <- nrow(joined_df[joined_df$loh_status.x == "1" & joined_df$loh_status.y == "0", ])
        private_loh_y <- nrow(joined_df[joined_df$loh_status.x == "0" & joined_df$loh_status.y == "1", ])
        
        paired_results <- c(tum, sample_list[i],sample_list[j],
                            euc_dist_loh, 
                            shared_loh,
                            private_loh_x, private_loh_y)
        
        all_vec <- rbind(all_vec,   paired_results )
      }
    }
  }
}

colnames(all_vec) <- c("tumour_id","Region.1st","Region.2nd",
                       "euc_dist_arm_loh", 
                       "shared_loh",
                       "private_loh_1st", "private_loh_2nd")    

all_res_df <- as.data.frame(all_vec, stringsAsFactors = F)  %>%
  mutate(compaired_regions = paste0(Region.1st, "-",Region.2nd))

## annotate growth pattern 

LOH_arm_annotated <- all_res_df %>%
  mutate(shared_loh = as.numeric(shared_loh),
         private_loh_1st = as.numeric(private_loh_1st),
         private_loh_2nd = as.numeric(private_loh_2nd),
         euc_dist_arm_loh = as.numeric(euc_dist_arm_loh)) %>%
  left_join(regional_df_join1, by = c("tumour_id", "Region.1st") ) %>%
  left_join(regional_df_join2, by = c("tumour_id", "Region.2nd") ) %>%
  dplyr::rename(pattern_1st = regional_pattern.new.x,
                pattern_2nd = regional_pattern.new.y) %>%
  mutate(grade_1st = pattern_to_grade(pattern_1st),
         grade_2nd = pattern_to_grade(pattern_2nd),
         grade_num_1st = pattern_to_grade.num(pattern_1st),
         grade_num_2nd = pattern_to_grade.num(pattern_2nd))

# save (load) LOH table (arm)

out_f <- paste0(tmpdir, "/LOH_arm_annotated.rds")
# saveRDS(LOH_arm_annotated, out_f)


################################################
## shared & private LOH : sequenza (cytoband) ##
################################################

in_dt <- sequenza_cyto_df %>% 
  mutate(sample = sample_name_cruk) %>%
  distinct(sample, cytoband_name, .keep_all = T) %>%
  mutate(region = str_split_fixed(sample, pattern = "-", n =2)[,2])

tumour_list <- as.character(na.omit(unique(in_dt$tumour_id_muttable_cruk) )) # including non-LUSCs, CRUK0555_LN1 is missing but no problem (single region tumour)
num_sample <- length(unique(in_dt$sample)) 

## now calculate euc dist and shared & private LOH for all regional pairs 

all_vec <- vector() # all results will be saved in this vector (will be a matrix in the end)
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
        dt_i <- dt_tum %>%
          filter(region == sample_list[i]) %>%
          dplyr::select(cytoband_name,  loh_status)
        dt_j <- dt_tum %>%
          filter(region == sample_list[j]) %>%
          dplyr::select(cytoband_name,  loh_status)
        
        # this time focus on LOH only 
        joined_df <- dt_i %>%
          inner_join(dt_j, by = "cytoband_name") %>%
          filter(!is.na(loh_status.x), !is.na(loh_status.y) ) 
        
        # use categorical loh status
        euc_dist_loh <- calc_euc_dis(joined_df$loh_status.x, joined_df$loh_status.y) 
        
        # shared / private loh 
        shared_loh <- nrow(joined_df[joined_df$loh_status.x == "1" & joined_df$loh_status.y == "1", ])
        private_loh_x <- nrow(joined_df[joined_df$loh_status.x == "1" & joined_df$loh_status.y == "0", ])
        private_loh_y <- nrow(joined_df[joined_df$loh_status.x == "0" & joined_df$loh_status.y == "1", ])
        
        paired_results <- c(tum, sample_list[i],sample_list[j],
                            euc_dist_loh,
                            shared_loh,
                            private_loh_x, private_loh_y)
        
        all_vec <- rbind(all_vec,   paired_results )
      }
    }
  }
}

colnames(all_vec) <- c("tumour_id","Region.1st","Region.2nd",
                       "euc_dist_seq_cyto_loh",
                       "shared_loh",
                       "private_loh_1st", "private_loh_2nd")    

all_res_df <- as.data.frame(all_vec, stringsAsFactors = F)  %>%
  mutate(compaired_regions = paste0(Region.1st, "-",Region.2nd))

## annotate growth pattern 

LOH_annotated_seq <- all_res_df %>%
  mutate(shared_loh = as.numeric(shared_loh),
         private_loh_1st = as.numeric(private_loh_1st),
         private_loh_2nd = as.numeric(private_loh_2nd),
         euc_dist_seq_cyto_loh = as.numeric(euc_dist_seq_cyto_loh)) %>%
  left_join(regional_df_join1, by = c("tumour_id", "Region.1st") ) %>%
  left_join(regional_df_join2, by = c("tumour_id", "Region.2nd") ) %>%
  dplyr::rename(pattern_1st = regional_pattern.new.x,
                pattern_2nd = regional_pattern.new.y) %>%
  mutate(grade_1st = pattern_to_grade(pattern_1st),
         grade_2nd = pattern_to_grade(pattern_2nd),
         grade_num_1st = pattern_to_grade.num(pattern_1st),
         grade_num_2nd = pattern_to_grade.num(pattern_2nd))

# save (load) LOH table (cytoband)

out_f <- paste0(tmpdir, "/LOH_annotated_seq.rds")
# saveRDS(LOH_annotated_seq, out_f)


###########################################
## shared & private LOH : sequenza (arm) ##
###########################################

in_dt <- sequenza_arm_df %>%
  mutate(sample = sample_name_cruk) %>%
  distinct(sample, chr_arm, .keep_all = T) %>%
  mutate(region = str_split_fixed(sample, pattern = "-", n =2)[,2])

all_vec <- vector() # all results will be saved in this vector (will be a matrix in the end)
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
          dplyr::select(chr_arm,  loh_status)
        dt_j <- dt_tum %>%
          filter(region == sample_list[j]) %>%
          dplyr::select(chr_arm,  loh_status)
        
        # this time focus on LOH only 
        joined_df <- dt_i %>%
          inner_join(dt_j, by = "chr_arm") %>%
          filter(!is.na(loh_status.x), !is.na(loh_status.y) ) 
        
        # use categorical loh status
        euc_dist_loh <- calc_euc_dis(joined_df$loh_status.x, joined_df$loh_status.y) 
        
        
        # shared / private loh 
        shared_loh <- nrow(joined_df[joined_df$loh_status.x == "1" & joined_df$loh_status.y == "1", ])
        private_loh_x <- nrow(joined_df[joined_df$loh_status.x == "1" & joined_df$loh_status.y == "0", ])
        private_loh_y <- nrow(joined_df[joined_df$loh_status.x == "0" & joined_df$loh_status.y == "1", ])
        
        paired_results <- c(tum, sample_list[i],sample_list[j],
                            euc_dist_loh, 
                            shared_loh,
                            private_loh_x, private_loh_y)
        
        all_vec <- rbind(all_vec,   paired_results )
      }
    }
  }
}

colnames(all_vec) <- c("tumour_id","Region.1st","Region.2nd",
                       "euc_dist_seq_arm_loh", 
                       "shared_loh",
                       "private_loh_1st", "private_loh_2nd")    

all_res_df <- as.data.frame(all_vec, stringsAsFactors = F)  %>%
  mutate(compaired_regions = paste0(Region.1st, "-",Region.2nd))

## annotate growth pattern 

LOH_arm_annotated_seq <- all_res_df %>%
  mutate(shared_loh = as.numeric(shared_loh),
         private_loh_1st = as.numeric(private_loh_1st),
         private_loh_2nd = as.numeric(private_loh_2nd),
         euc_dist_seq_arm_loh = as.numeric(euc_dist_seq_arm_loh)) %>%
  left_join(regional_df_join1, by = c("tumour_id", "Region.1st") ) %>%
  left_join(regional_df_join2, by = c("tumour_id", "Region.2nd") ) %>%
  dplyr::rename(pattern_1st = regional_pattern.new.x,
                pattern_2nd = regional_pattern.new.y) %>%
  mutate(grade_1st = pattern_to_grade(pattern_1st),
         grade_2nd = pattern_to_grade(pattern_2nd),
         grade_num_1st = pattern_to_grade.num(pattern_1st),
         grade_num_2nd = pattern_to_grade.num(pattern_2nd))

## save LOH table (sequenza arm)

out_f <- paste0(tmpdir, "/LOH_arm_annotated_seq.rds")
# saveRDS(LOH_arm_annotated_seq, out_f)


########################################
## extract dominant mutation profile  ## 
########################################

## shared & private dominant mutations per regional pairs (mean cluster CCF>=95) 
mphase_tumours <- names(trees)
mphase_tumours_tmp <- setdiff(mphase_tumours, c("CRUK0721_Tumour1", "CRUK0555_Tumour1"))  # exclude single region tree

clusters <- list()
all_search_result  <- vector()

for(k in 1:length(mphase_tumours_tmp)){
  tumour = mphase_tumours_tmp[k]
  print(paste(k, tumour))
  
  # The following tree_cluster will not include mutations (clusters) removed during tree building.
  tree_cluster <- unique(sort(trees[[tumour]][["graph_pyclone"]][["Corrected_tree"]])) # same as above
  
  count_table <-  data.table(ClusterName = names(trees[[tumour]][["graph_pyclone"]][["edgelength"]]),
                             numSNVs = as.matrix(trees[[tumour]][["graph_pyclone"]][["edgelength"]])[,1])
  
  count_table_df <- count_table %>%
    mutate(numSNVs = as.numeric(numSNVs)) %>%
    filter(ClusterName %in% tree_cluster)
  
  ccf_table <- trees[[tumour]][["nested_pyclone"]][["ccf_cluster_table"]] 
  ccf_table <- ccf_table[rownames(ccf_table) %in% tree_cluster,]
  
  # account for tumour with only 1 pyclone cluster
  if(is.vector(ccf_table)){
    ccf_table <- as.matrix(rbind(ccf_table))
  }
  
  # take only clusters with mean ccf >=95% into account
  clusters[[tumour]] <-  apply(ccf_table, 2, function(j){
    return(rownames(ccf_table)[j >= 95])
  })
  
  # account for tumours with few clusters to be assessed
  if(is.matrix(clusters[[tumour]])){
    clusters[[tumour]] <- as.data.frame(clusters[[tumour]])
  }
  
  if(!is.list(clusters[[tumour]])){
    clusters[[tumour]] <- as.data.frame(rbind(clusters[[tumour]]))
  }
  
  # investigate parent-daughter related regions
  
  for(m in 1:(ncol(ccf_table)-1)){
    for(n in (m+1):ncol(ccf_table)){
      r1 <- colnames(ccf_table)[m]
      r2 <- colnames(ccf_table)[n]
      
      # shared clusters
      shared_clusters <- intersect(unlist(clusters[[tumour]][[r1]]), unlist(clusters[[tumour]][[r2]]))
      
      # private_clusters
      private_clusters_r1 <- setdiff(unlist(clusters[[tumour]][[r1]]), unlist(clusters[[tumour]][[r2]]))
      private_clusters_r2 <- setdiff(unlist(clusters[[tumour]][[r2]]), unlist(clusters[[tumour]][[r1]]))
      
      # counts of shared and private mutations
      
      shared_counts <- sum(count_table_df$numSNVs[count_table_df$ClusterName %in% shared_clusters])
      
      if(length(private_clusters_r1) >0){
        private_counts_r1 <- sum(count_table_df$numSNVs[count_table_df$ClusterName %in% private_clusters_r1])
      } else {
        private_counts_r1 <- 0
      }
      
      if(length(private_clusters_r2) >0){
        private_counts_r2 <- sum(count_table_df$numSNVs[count_table_df$ClusterName %in% private_clusters_r2])
      } else {
        private_counts_r2 <- 0
      }
      
      search_result <- c(r1,r2,shared_counts, private_counts_r1 ,private_counts_r2)
      all_search_result <- rbind(all_search_result, search_result )
    }
  }
}

all_search_result_df <- as.data.frame(all_search_result, stringsAsFactors = F)
colnames(all_search_result_df) <- c("Region.1st", "Region.2nd","shared_counts", "private_counts_r1" ,"private_counts_r2")
all_search_result_df$shared_counts <- as.numeric(all_search_result_df$shared_counts)
all_search_result_df$private_counts_r1 <- as.numeric(all_search_result_df$private_counts_r1)
all_search_result_df$private_counts_r2 <- as.numeric(all_search_result_df$private_counts_r2)

shared.private.numSNV_clusterCCF95 <- all_search_result_df

# save
# saveRDS(shared.private.numSNV_clusterCCF95, paste0(tmpdir, "/shared.private.numSNV_clusterCCF95.rds")) 


## dominant mutation burden per each region (CCF>=95) - for Monte Carlo permutation test

all_vec <- vector()
for(k in 1:length(mphase_tumours_tmp)){
  tumour = mphase_tumours_tmp[k]
  print(paste(k, tumour))
  
  # The following tree_cluster will include mutations (clusters) removed in tree building. 
  # To avoid potential bias introduced by tree building process, for the current analysis, do include such mutations.
  tree_cluster <- rownames(trees[[tumour]][["clonality_out"]][["clonality_table_corrected"]]) 
  
  ccf_table <- trees[[tumour]][["nested_pyclone"]][["ccf_cluster_table"]] 
  ccf_table <- ccf_table[rownames(ccf_table) %in% tree_cluster,]
  
  # account for tumour with only 1 pyclone cluster
  if(is.vector(ccf_table)){
    ccf_table <- as.matrix(rbind(ccf_table))
  }
  
  for(m in 1:ncol(ccf_table) ){ # for each region
    r1 <- colnames(ccf_table)[m] # get region name
    
    mut_r1 <- data.frame(PhyloCCF = trees[[tumour]]$ccf_table_pyclone_clean[[paste0(r1, "_PhyloCCF")]],
                         Cluster = trees[[tumour]]$ccf_table_pyclone_clean[["PycloneCluster"]],
                         stringsAsFactors = F) %>%
      filter(PhyloCCF >= 0.95 & Cluster %in% tree_cluster)
    
    tmb_r1 <- nrow(mut_r1)
    
    R1 <- gsub("\\.","-",r1)
    sample_vec <- c(R1,  tmb_r1 )
    all_vec <- rbind(all_vec, sample_vec )
  }
}

tmb_ccf95_df <- data.frame(all_vec, stringsAsFactors = F)
colnames(tmb_ccf95_df ) <- c("sample","tmb_95")

# saveRDS(tmb_ccf95_df, paste0(tmpdir, "/tmb_ccf95.RDS"))


#######################
## mutation distance ##
#######################

in_dt <- muttable_region %>%  # this will include single region tumours
  mutate(region = str_split_fixed(RegionID, pattern = "\\.", n = 2)[,2]) %>%
  filter(!grepl("LN", RegionID)) %>%
  filter(Is.present.region == TRUE) %>%
  mutate(mut_id = mutation_id)

# calculate distance between all regional pairs

tumour_list <- unique(in_dt$tumour_id)

all_vec <- vector() # all results will be saved in this vector (will be a matrix in the end)
for (tum in tumour_list){
  print(tum)
  
  dt_tum <- in_dt %>%
    filter(tumour_id == tum) %>%
    distinct(mut_id, region) %>%
    mutate(marker = 1)
  
  sample_list <- sort(unique(dt_tum$region))
  sample_num <- length(sample_list)
  
  # calculate Euclidean distance between all regional pairs
  
  if(sample_num > 1){  # sample_num : same as number of regions in each tumour
    for (i in 1:(sample_num - 1)){
      for (j in (1+i):sample_num){
        #i=1
        #j=2
        dt_i <- dt_tum %>%
          filter(region == sample_list[i]) %>%
          dplyr::select(mut_id, marker) %>%
          dplyr::rename(marker.x = marker)
        dt_j <- dt_tum %>%
          filter(region == sample_list[j]) %>%
          dplyr::select(mut_id, marker) %>%
          dplyr::rename(marker.y = marker)
        
        ## Euclidean distance
        
        joined_df <- dt_i %>%
          full_join(dt_j, by = "mut_id") %>%
          # distinct(mut_id, .keep_all = T) %>%
          mutate(marker.x = replace_na(marker.x, 0), # mutations in dt_j but not in dt_i will be marked as 0 (private mutation in sample_list[j])
                 marker.y = replace_na(marker.y, 0) ) # mutations in dt_i but not in dt_j will be marked as 0 (private mutation in sample_list[i])
        
        euc_dist <- calc_euc_dis(joined_df$marker.x, joined_df$marker.y) 
        
        dist_results <- c(tum, sample_list[i],sample_list[j],
                          euc_dist)
        all_vec <- rbind(all_vec,   dist_results )
      }
    }
  }
}

colnames(all_vec) <- c("tumour_id","Region.1st","Region.2nd",
                       "euc_dist_mut")

all_mut_df <- as.data.frame(all_vec, stringsAsFactors = F)  %>%
  mutate(compaired_regions = paste0(Region.1st, "-",Region.2nd)) 

mut_annotated <- all_mut_df %>%
  mutate(euc_dist_mut = as.numeric(euc_dist_mut)) %>%
  left_join(regional_df_join1, by = c("tumour_id", "Region.1st") ) %>%
  left_join(regional_df_join2, by = c("tumour_id", "Region.2nd") ) %>%
  dplyr::rename(pattern_1st = regional_pattern.new.x,
                pattern_2nd = regional_pattern.new.y) %>%
  mutate(grade_1st = pattern_to_grade(pattern_1st),
         grade_2nd = pattern_to_grade(pattern_2nd),
         grade_num_1st = pattern_to_grade.num(pattern_1st),
         grade_num_2nd = pattern_to_grade.num(pattern_2nd))

# save
# saveRDS(mut_annotated, paste0(tmpdir, "/mut_annotated.rds"))


###############################################################################
##  create data frame with all regional pairs with growth pattern annotation ##
###############################################################################

all_labels_vec <- vector()
for (tum in tumour_id_list){
  print(tum)
  
  dt_tum <- regional_df %>%
    filter(tumour_id_muttable_cruk == tum) %>%
    filter(regional_pattern.new %in% c("well_dx","poor_dx"))
  
  sample_list <- sort(unique(dt_tum$region))
  sample_num <- length(sample_list)
  
  if(sample_num > 1){  # sample_num : same as number of regions in each tumour
    for (i in 1:(sample_num - 1)){
      for (j in (1+i):sample_num){
        #i=1
        #j=2
        label_regional_patten <- ifelse(dt_tum$regional_pattern.new[dt_tum$region == sample_list[i]] == dt_tum$regional_pattern.new[dt_tum$region == sample_list[j]],
                                        "same","diff")
        
        labels_vec <- c(tum, paste(tum,sample_list[i],sample_list[j], sep = ":"),
                        sample_list[i],sample_list[j],
                        as.character(dt_tum$regional_pattern.new[dt_tum$region == sample_list[i]]), 
                        as.character(dt_tum$regional_pattern.new[dt_tum$region == sample_list[j]]) ,
                        label_regional_patten)
        all_labels_vec <- rbind(all_labels_vec,  labels_vec )
      }
    }
  }
}
colnames(all_labels_vec) <- c("tumour_id","comb_id",
                              "Region.1st","Region.2nd",
                              "pattern_1", "pattern_2",
                              "label_regional_patten")

all_patterns_df <- as.data.frame(all_labels_vec, stringsAsFactors = F) %>%
  mutate(compaired_regions = paste0(Region.1st, "-",Region.2nd))

# saveRDS(all_patterns_df,paste0(tmpdir, "all_patterns_df.rds") )

################################################################
## create data frame of minimum purity per each regional pair ##
################################################################

all_purity_vec <- vector()
for (tum in tumour_id_list){
  print(tum)
  
  dt_tum <- regional_df %>%   # regional_df only include regions with LUSC subtype 
    filter(tumour_id_muttable_cruk == tum) %>%
    filter(!is.na(purity)) %>%
    arrange(region)
  
  sample_list <- sort(unique(dt_tum$region))
  sample_num <- length(sample_list)
  
  if(sample_num > 1){  # sample_num : same as number of regions in each tumour
    for (i in 1:(sample_num - 1)){
      for (j in (1+i):sample_num){
        #i=1
        #j=2
        min_purity = min(dt_tum[which(dt_tum$region == sample_list[i]),][["purity"]], dt_tum[which(dt_tum$region == sample_list[j]),][["purity"]])
        
        purity_vec <- c(tum, paste(tum,sample_list[i],sample_list[j], sep = ":"),
                        sample_list[i],sample_list[j],
                        min_purity)
        all_purity_vec <- rbind(all_purity_vec,  purity_vec)
      }
    }
  }
}

colnames(all_purity_vec) <- c("tumour_id","comb_id",
                              "Region.1st","Region.2nd",
                              "min_purity")

all_purity_df <- as.data.frame(all_purity_vec, stringsAsFactors = F) %>%
  mutate(compaired_regions = paste0(Region.1st, "-",Region.2nd)) %>%
  mutate( min_purity = as.numeric(min_purity))

# save
# saveRDS(all_purity_df, paste0(tmpdir, "/all_purity_df.rds"))


####################################################
## combine all genomic distance (also add purity) ##
####################################################

# load 
# LOH_annotated <- readRDS(paste0(tmpdir, "/LOH_annotated.rds"))
# LOH_annotated_seq <- readRDS(paste0(tmpdir, "/LOH_annotated_seq.rds"))
# mut_annotated <- readRDS(paste0(tmpdir, "/mut_annotated.rds"))
# all_patterns_df <- readRDS(paste0(tmpdir, "/all_patterns_df.rds"))
# all_purity_df <- readRDS(paste0(tmpdir, "/all_purity_df.rds"))

all_dist_df <- all_patterns_df %>% 
  left_join(dplyr::select(mut_annotated, tumour_id,compaired_regions,euc_dist_mut ), by = c("tumour_id", "compaired_regions")) %>%
  left_join(dplyr::select(LOH_annotated,tumour_id,compaired_regions, euc_dist_cyto_loh), by = c("tumour_id", "compaired_regions")) %>%
  left_join(dplyr::select(LOH_annotated_seq,tumour_id,compaired_regions, euc_dist_seq_cyto_loh), by = c("tumour_id", "compaired_regions")) %>%
  left_join(dplyr::distinct(all_purity_df,tumour_id,compaired_regions, min_purity), by = c("tumour_id", "compaired_regions")) %>%
  mutate(subtype_pair = paste(pattern_1, pattern_2, sep = "-"))

# save 
# saveRDS(all_dist_df, paste0(tmpdir, "/all_dist_df.RDS") )


#########################################################
## Preliminary call for Ancestor-Descendant-like pairs ##
#########################################################

# load  
# LOH_annotated <- readRDS(paste0(tmpdir, "/LOH_annotated.rds"))
# LOH_arm_annotated <- readRDS(paste0(tmpdir, "/LOH_arm_annotated.rds"))
# LOH_annotated_seq <- readRDS(paste0(tmpdir, "/LOH_annotated_seq.rds"))
# LOH_arm_annotated_seq <- readRDS(paste0(tmpdir, "/LOH_arm_annotated_seq.rds"))
# shared.private.numSNV_clusterCCF95 <- readRDS(paste0(tmpdir, "/shared.private.numSNV_clusterCCF95.rds"))


purity_df <- LUSC_region_df %>% 
  dplyr::select(tumour_id_muttable_cruk, region, purity) %>%
  dplyr::rename(tumour_id = tumour_id_muttable_cruk )

## call A-D-like pairs using mphase LOH cytoband & arm, sequenza LOH cytoband & arm, and dominant mutation data 
# cutoffs - ancestor-like: < 1-5 %, descendant-like: > 7.5 - 15%

ratio_filter_p_list <- c(0.01, 0.02, 0.05)
ratio_filter_d_list <- c(0.075, 0.1, 0.15)

for(n in 1:length(ratio_filter_p_list)){
  for(m in 1:length(ratio_filter_d_list)){
    
    ratio_filter_p <- ratio_filter_p_list[n]
    ratio_filter_d <- ratio_filter_d_list[m]
    
    if(ratio_filter_p <= ratio_filter_d) {
      
      ### preliminary calling of ancestor (parent)-like and descendant (daughter)-like regions
      # only use private branch proportion (purity filter and arm filter will be applied later on)
      
      ## LOH cytoband (mphase)
      
      loh_cyto_branc.prop <- LOH_annotated %>%
        mutate(private_prop_1st = private_loh_1st/shared_loh,
               private_prop_2nd = private_loh_2nd/shared_loh) %>%
        mutate(parent_region = case_when( private_prop_1st > ratio_filter_d & private_prop_2nd < ratio_filter_p  ~ "region.2nd",
                                          private_prop_2nd > ratio_filter_d & private_prop_1st < ratio_filter_p  ~ "region.1st",
                                          TRUE ~ NA_character_)) %>%
        mutate(parent_R = case_when(parent_region == "region.1st" ~ Region.1st,
                                    parent_region == "region.2nd" ~ Region.2nd),
               daughter_R = case_when(parent_region == "region.1st" ~ Region.2nd,
                                      parent_region == "region.2nd" ~ Region.1st)) %>%
        mutate(parent_length = case_when(parent_region == "region.1st" ~ private_loh_1st,
                                         parent_region == "region.2nd" ~ private_loh_2nd),
               daughter_length = case_when(parent_region == "region.1st" ~ private_loh_2nd,
                                           parent_region == "region.2nd" ~ private_loh_1st)) %>%
        mutate(p.d_comb = paste0(tumour_id,":",compaired_regions,":",parent_region),
               p.d_pair = paste0(tumour_id,":",compaired_regions)) 
      
      ## LOH arm (mphase)
      
      loh_arm_branc.prop <- LOH_arm_annotated %>%
        mutate(private_prop_1st = private_loh_1st/shared_loh ,
               private_prop_2nd = private_loh_2nd/shared_loh) %>%
        mutate(parent_region = case_when( private_prop_1st > ratio_filter_d & private_prop_2nd < ratio_filter_p  ~ "region.2nd",
                                          private_prop_2nd > ratio_filter_d & private_prop_1st < ratio_filter_p  ~ "region.1st",
                                          TRUE ~ NA_character_)) %>%
        mutate(parent_R = case_when(parent_region == "region.1st" ~ Region.1st,
                                    parent_region == "region.2nd" ~ Region.2nd),
               daughter_R = case_when(parent_region == "region.1st" ~ Region.2nd,
                                      parent_region == "region.2nd" ~ Region.1st)) %>%
        mutate(parent_length = case_when(parent_region == "region.1st" ~ private_loh_1st,
                                         parent_region == "region.2nd" ~ private_loh_2nd),
               daughter_length = case_when(parent_region == "region.1st" ~ private_loh_2nd,
                                           parent_region == "region.2nd" ~ private_loh_1st)) %>%
        mutate(p.d_comb = paste0(tumour_id,":",compaired_regions,":",parent_region),
               p.d_pair = paste0(tumour_id,":",compaired_regions))  
      
      # LOH cyto (sequenza) 
      
      loh_cyto_seq_branc.prop <- LOH_annotated_seq %>%
        mutate(private_prop_1st = private_loh_1st/shared_loh,
               private_prop_2nd = private_loh_2nd/shared_loh ) %>%
        mutate(parent_region = case_when( private_prop_1st > ratio_filter_d & private_prop_2nd < ratio_filter_p  ~ "region.2nd",
                                          private_prop_2nd > ratio_filter_d & private_prop_1st < ratio_filter_p  ~ "region.1st",
                                          TRUE ~ NA_character_)) %>%
        mutate(parent_R = case_when(parent_region == "region.1st" ~ Region.1st,
                                    parent_region == "region.2nd" ~ Region.2nd),
               daughter_R = case_when(parent_region == "region.1st" ~ Region.2nd,
                                      parent_region == "region.2nd" ~ Region.1st)) %>%
        mutate(parent_length = case_when(parent_region == "region.1st" ~ private_loh_1st,
                                         parent_region == "region.2nd" ~ private_loh_2nd),
               daughter_length = case_when(parent_region == "region.1st" ~ private_loh_2nd,
                                           parent_region == "region.2nd" ~ private_loh_1st)) %>%
        mutate(p.d_comb = paste0(tumour_id,":",compaired_regions,":",parent_region),
               p.d_pair = paste0(tumour_id,":",compaired_regions)) 
      
      ## LOH arm (sequenza)
      
      loh_arm_seq_branc.prop <- LOH_arm_annotated_seq %>%
        mutate(private_prop_1st = private_loh_1st/shared_loh,
               private_prop_2nd = private_loh_2nd/shared_loh) %>%
        mutate(parent_region = case_when( private_prop_1st > ratio_filter_d & private_prop_2nd < ratio_filter_p  ~ "region.2nd",
                                          private_prop_2nd > ratio_filter_d & private_prop_1st < ratio_filter_p  ~ "region.1st",
                                          TRUE ~ NA_character_)) %>%
        mutate(parent_R = case_when(parent_region == "region.1st" ~ Region.1st,
                                    parent_region == "region.2nd" ~ Region.2nd),
               daughter_R = case_when(parent_region == "region.1st" ~ Region.2nd,
                                      parent_region == "region.2nd" ~ Region.1st)) %>%
        mutate(parent_length = case_when(parent_region == "region.1st" ~ private_loh_1st,
                                         parent_region == "region.2nd" ~ private_loh_2nd),
               daughter_length = case_when(parent_region == "region.1st" ~ private_loh_2nd,
                                           parent_region == "region.2nd" ~ private_loh_1st)) %>%
        mutate(p.d_comb = paste0(tumour_id,":",compaired_regions,":",parent_region),
               p.d_pair = paste0(tumour_id,":",compaired_regions))  
      
      ## mutation (mean cluster CCF>= 95% in tree building)
      
      mut_ccf95_branc.prop <- shared.private.numSNV_clusterCCF95 %>%
        filter(!Region.1st%in% c("CRUK0519_SU_T2.R1", "CRUK0519_SU_T2.R2"), !Region.2nd%in% c("CRUK0519_SU_T2.R1", "CRUK0519_SU_T2.R2")) %>%  # exclude intra-pulm mets (CRUK0495 T3 should also be removed, but this is LUSC, so we can ignore here)
        mutate(region.1st = gsub("\\.","-",Region.1st ),
               region.2nd = gsub("\\.","-",Region.2nd )) %>%
        mutate(private_prop_1st = private_counts_r1/shared_counts ,
               private_prop_2nd = private_counts_r2/shared_counts ) %>%
        mutate(parent_region = case_when( private_prop_1st > ratio_filter_d & private_prop_2nd < ratio_filter_p  ~ "region.2nd",
                                          private_prop_2nd > ratio_filter_d & private_prop_1st < ratio_filter_p  ~ "region.1st",
                                          TRUE ~ NA_character_)) %>%
        filter(!grepl("LN",Region.1st) & !grepl("LN",Region.2nd)) %>%
        mutate(Region.1st = str_split_fixed(Region.1st, pattern = "\\.", n =2)[,2],
               Region.2nd = str_split_fixed(Region.2nd, pattern = "\\.", n =2)[,2]) %>%
        mutate(parent = case_when(parent_region == "region.1st" ~ Region.1st,
                                  parent_region == "region.2nd" ~ Region.2nd),
               daughter = case_when(parent_region == "region.1st" ~ Region.2nd,
                                    parent_region == "region.2nd" ~ Region.1st)) %>%
        mutate(parent_length = case_when(parent_region == "region.1st" ~ private_counts_r1,
                                         parent_region == "region.2nd" ~ private_counts_r2),
               daughter_length = case_when(parent_region == "region.1st" ~ private_counts_r2,
                                           parent_region == "region.2nd" ~ private_counts_r1)) %>%
        mutate(compaired_regions = paste0(Region.1st,"-",Region.2nd)) %>%
        left_join(regional_df[,c("sample_name_cruk", "tumour_id_muttable_cruk")], by = c("region.1st" = "sample_name_cruk")) %>%
        dplyr::rename(tumour_id  = tumour_id_muttable_cruk )  %>%
        mutate(p.d_comb = paste0(tumour_id,":",compaired_regions,":",parent_region),
               p.d_pair = paste0(tumour_id,":",compaired_regions))  
      
      # combine all profiles 
      
      loh_cyto_join <- loh_cyto_branc.prop %>%
        dplyr::select(tumour_id,compaired_regions, p.d_comb, parent_R , daughter_R,  parent_length, daughter_length, shared_loh,private_loh_1st,private_loh_2nd,  private_prop_1st,private_prop_2nd) %>%
        dplyr::rename(parent_length_LOH_cyto = parent_length, 
                      daughter_length_LOH_cyto = daughter_length, 
                      shared_LOH_cyto = shared_loh,
                      p.d_comb_LOH_cyto = p.d_comb,
                      parent_R_LOH_cyto =  parent_R, 
                      daughter_R_LOH_cyto = daughter_R,
                      private_loh_1st_LOH_cyto = private_loh_1st,
                      private_loh_2nd_LOH_cyto = private_loh_2nd,
                      private_prop_1st_LOH_cyto = private_prop_1st,
                      private_prop_2nd_LOH_cyto = private_prop_2nd)
      
      loh_arm_join <- loh_arm_branc.prop %>%
        dplyr::select(tumour_id,compaired_regions, p.d_comb, parent_R , daughter_R,  parent_length, daughter_length, shared_loh,private_loh_1st,private_loh_2nd,  private_prop_1st,private_prop_2nd) %>%
        dplyr::rename(parent_length_LOH_arm = parent_length, 
                      daughter_length_LOH_arm = daughter_length, 
                      shared_LOH_arm = shared_loh,
                      p.d_comb_LOH_arm = p.d_comb,
                      parent_R_LOH_arm =  parent_R, 
                      daughter_R_LOH_arm = daughter_R ,
                      private_loh_1st_LOH_arm = private_loh_1st,
                      private_loh_2nd_LOH_arm = private_loh_2nd,
                      private_prop_1st_LOH_arm = private_prop_1st,
                      private_prop_2nd_LOH_arm = private_prop_2nd)
      
      loh_cyto_seq_join <- loh_cyto_seq_branc.prop %>%
        dplyr::select(tumour_id,compaired_regions, p.d_comb, parent_R , daughter_R,  parent_length, daughter_length, shared_loh,private_loh_1st,private_loh_2nd,  private_prop_1st,private_prop_2nd) %>%
        dplyr::rename(parent_length_LOH_cyto_seq = parent_length, 
                      daughter_length_LOH_cyto_seq = daughter_length, 
                      shared_LOH_cyto_seq = shared_loh,
                      p.d_comb_LOH_cyto_seq = p.d_comb,
                      parent_R_LOH_cyto_seq =  parent_R, 
                      daughter_R_LOH_cyto_seq = daughter_R,
                      private_loh_1st_LOH_cyto_seq = private_loh_1st,
                      private_loh_2nd_LOH_cyto_seq = private_loh_2nd,
                      private_prop_1st_LOH_cyto_seq = private_prop_1st,
                      private_prop_2nd_LOH_cyto_seq = private_prop_2nd)
      
      loh_arm_seq_join <- loh_arm_seq_branc.prop %>%
        dplyr::select(tumour_id,compaired_regions, p.d_comb, parent_R , daughter_R,  parent_length, daughter_length, shared_loh,private_loh_1st,private_loh_2nd,  private_prop_1st,private_prop_2nd) %>%
        dplyr::rename(parent_length_LOH_arm_seq = parent_length, 
                      daughter_length_LOH_arm_seq = daughter_length, 
                      shared_LOH_arm_seq = shared_loh,
                      p.d_comb_LOH_arm_seq = p.d_comb,
                      parent_R_LOH_arm_seq =  parent_R, 
                      daughter_R_LOH_arm_seq = daughter_R ,
                      private_loh_1st_LOH_arm_seq = private_loh_1st,
                      private_loh_2nd_LOH_arm_seq = private_loh_2nd,
                      private_prop_1st_LOH_arm_seq = private_prop_1st,
                      private_prop_2nd_LOH_arm_seq = private_prop_2nd)
      
      
      mut_ccf95_join <- mut_ccf95_branc.prop %>%
        dplyr::select(tumour_id,compaired_regions, p.d_comb, parent , daughter,  parent_length, daughter_length, shared_counts,private_counts_r1,private_counts_r2, private_prop_1st,private_prop_2nd ) %>%
        dplyr::rename(parent_length_mut_ccf95 = parent_length, 
                      daughter_length_mut_ccf95 = daughter_length, 
                      shared_mut_ccf95 = shared_counts,
                      p.d_comb_mut_ccf95 = p.d_comb,
                      parent_R_mut_ccf95 =  parent, 
                      daughter_R_mut_ccf95 = daughter,
                      private_mut_1st_mut_ccf95 = private_counts_r1,
                      private_mut_2nd_mut_ccf95 = private_counts_r2, 
                      private_prop_1st_mut_ccf95 = private_prop_1st ,
                      private_prop_2nd_mut_ccf95 = private_prop_2nd)  %>%
        distinct(tumour_id,compaired_regions, .keep_all = T)
      
      # join all 
      all_joined_comb_branch.prop <-  LOH_annotated %>%
        dplyr::select(tumour_id,Region.1st,Region.2nd, compaired_regions) %>% 
        left_join(loh_cyto_join, by = c("tumour_id","compaired_regions")) %>%
        left_join(loh_arm_join, by = c("tumour_id","compaired_regions")) %>%
        left_join(loh_cyto_seq_join, by = c("tumour_id","compaired_regions")) %>%
        left_join(loh_arm_seq_join, by = c("tumour_id","compaired_regions")) %>%
        left_join(mut_ccf95_join, by = c("tumour_id","compaired_regions"))
      
      # add purity and only keep LUSC tumours
      all_joined_comb_branch.prop <- add_purity(all_joined_comb_branch.prop, purity_df)  %>%
        filter(tumour_id %in% LUSC_tumour_id) 
      
      # save 
      saveRDS(all_joined_comb_branch.prop, paste0(tmpdir, "/all_joined_comb_branch.prop_ratio_filter_p",ratio_filter_p,"_ratio_filter_d",ratio_filter_d,".tree.ratio.RDS"))
    }
  }
}





#############################
## create pattern_df_random ## 
##############################

pattern_df <- LUSC_region_df[,c("tumour_id_muttable_cruk","region","regional_pattern.new")] %>%
  dplyr::rename(tumour_id = tumour_id_muttable_cruk) 

# split pattern_df into annotated region and n.a. region (do not include n.a. region for randomisation)
pattern_df_anno <- pattern_df %>% filter(regional_pattern.new %in% subtype_list) 
pattern_df_na <- pattern_df %>% filter(!regional_pattern.new %in% subtype_list) 

# randomise growth pattern within each patient for 1000 times 
set.seed(123)
permute_num <- 1000

pattern_df_random <- list() 
for(i in 1:permute_num){
  print(i)
  shuffled_tmp <- sapply(unique(pattern_df_anno$tumour_id), function(x){
    pattern_vec_tum <- pattern_df_anno$regional_pattern.new[which(pattern_df_anno$tumour_id == x)]
    pattern_vec_tum_rand <- sample(pattern_vec_tum)
    return(pattern_vec_tum_rand)
  })
  shuffled_pattern <- unlist(shuffled_tmp)
  
  pattern_df_anno_random <- cbind(pattern_df_anno[,c("tumour_id", "region")],
                                  data.frame(regional_pattern.new = shuffled_pattern))
  
  pattern_df_random[[i]] <- rbind(pattern_df_anno_random, pattern_df_na)
}

#save this 
# saveRDS(pattern_df_random , file = paste0(tmpdir, "/pattern_df_random.RDS"))



#######################
#######################
## Step 2: analyse! ##
#######################
#######################

# number of permutations for Monte-Carlo simulations
permute_num <- 1000

####################################################################################
## genomic distance comparison between "same" vs "different" growth pattern pairs ##  
####################################################################################

# load 
all_dist_df <- readRDS(paste0(tmpdir, "all_dist_df.RDS"))


## boxplot: compare genomic distance between "same" vs "different" growth pattern pairs  -  use only mixed pattern tumours (at least one pair of both same & diff)

group_column_list <- c("label_regional_patten") 
var_list <- c("euc_dist_mut",  "euc_dist_cyto_loh" ,  "euc_dist_seq_cyto_loh")
jitter_width <- 0.25

# use same tumours as loh and mut 

all_dist_df_common <- all_dist_df %>%
  filter(!is.na(euc_dist_cyto_loh), !is.na(euc_dist_mut))

# check sequenza version too
sum(is.na(all_dist_df_common$euc_dist_seq_cyto_loh))
# [1] 0

# plot
out_f <- paste0(workdir,"/Fig.2b.S3l.genomic_distance_mixed.pattern_common.pdf")
pdf(out_f, width = 6, height = 6)

for (group in group_column_list){
  for (var in  var_list) {
    
    tmp_df <- all_dist_df_common
    tmp_df$group <- factor(tmp_df[[group]])
    tmp_df$var <- tmp_df[[var]]
    tmp_df <- tmp_df %>%
      filter(!is.na(group)) %>%
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
          
          # focus on heterogenous pattern tumours         
          count_group <- df %>%
            distinct(tumour_id, group) %>%
            group_by(tumour_id) %>%
            summarise(count =n())  
          
          het_tum <- count_group$tumour_id[which(count_group$count > 1)]
          
          data <- df %>%
            filter(tumour_id %in% c(het_tum))
          
          if(length(unique(data$tumour_id))>0) {
            
            # calculate pvalue using lme()
            if(length(unique(data$tumour_id))>2) {  #Alvin: changed to >=2 will cause Error in MEEM(object, conLin, control$niterEM) : Overfitted model
              
              data_df <- data %>%
                dplyr::select(comb_id , tumour_id, group, var) %>%
                mutate(group = factor(group)) 
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
            num_tum <- length(unique(data$tumour_id))
            
            g <-  ggplot(data, aes(x = group, y = var)) +
              geom_jitter( size = 4, height = 0.0,width = jitter_width , fill = "gray40",
                           alpha = 0.5, colour="black", pch=21) +
              geom_boxplot(size = 1, width=0.7, alpha = 0, outlier.shape = NA) +
              labs( y = var, x = group,
                    title = paste0("N= ", num_tum, " tumours (",sum(data_df$group == Group1)," & ",sum(data_df$group == Group2)," regional pairs)")) +
              theme_cowplot()+
              theme(
                axis.title.x = element_text(size = 22),
                axis.title.y = element_text(size = 22),
                axis.text.x = element_text(size = 22),
                axis.text.y = element_text(size = 18),
                legend.position = "right",
                legend.text = element_text(size = 18),
                legend.title = element_text(size = 22)) + 
              annotate("text",size = 7, x = 1.5, y = (min(data$var) + y_range*1.1),
                       label = paste("p=",pval,"(lme)"))
            
            print(g)
            
          } #if(length(unique(data$tumour_id))>0)
        } #for(j in (i+1):length(group.var_list))
      } #for(i in 1:(length(group.var_list)-1)) {
    } #if(length(group.var_list) > 1) {
  } #for (var in var_list) {
}#for (group in group_column_list){
dev.off()


## only use regions with purity > 0.4
group_column_list <- c("label_regional_patten") 
var_list <- c("euc_dist_mut",  "euc_dist_cyto_loh")
jitter_width <- 0.25

purity.cutoff = 0.4

all_dist_df_common_purity <- all_dist_df %>%
  filter(!is.na(euc_dist_cyto_loh), !is.na(euc_dist_mut)) %>%
  filter(min_purity >= purity.cutoff)

out_f <- paste0(workdir,"/Fig.S3mn.genomic_distance_mixed.pattern_common_min.purity.",purity.cutoff,".pdf")
pdf(out_f, width = 6, height = 6)

for (group in group_column_list){
  for (var in  var_list) {
    
    tmp_df <- all_dist_df_common_purity 
    
    tmp_df$group <- factor(tmp_df[[group]])
    tmp_df$var <- tmp_df[[var]]
    tmp_df <- tmp_df %>%
      filter(!is.na(group)) %>%
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
          
          # focus on heterogenous pattern tumous
          
          count_group <- df %>%
            distinct(tumour_id, group) %>%
            group_by(tumour_id) %>%
            summarise(count =n())  
          
          het_tum <- count_group$tumour_id[which(count_group$count > 1)]
          
          data <- df %>%
            filter(tumour_id %in% c(het_tum))
          
          if(length(unique(data$tumour_id))>0) {
            
            # calculate pvalue using lme()
            if(length(unique(data$tumour_id))>2) {  
              
              data_df <- data %>%
                dplyr::select(comb_id , tumour_id, group, var) %>%
                mutate(group = factor(group)) #, levels = rev(levels(regional_pattern))))
              names(data_df) <- c("short_sample", "tumour_id", "group", "var")
              
              effect_model <- lme(var ~ group , random = ~ 1 |tumour_id, data = data_df, method = "ML")
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
            num_tum <- length(unique(data$tumour_id))
            
            # plot
            g <-  ggplot(data, aes(x = group, y = var)) +
              geom_jitter( size = 4, height = 0.0,width = jitter_width , fill = "gray40",
                           alpha = 0.5, colour="black", pch=21) +
              geom_boxplot(size = 1, width=0.7, alpha = 0, outlier.shape = NA) +
              labs( y = var, x = group,
                    title = paste0("N= ", num_tum, " tumours (",sum(data_df$group == Group1)," & ",sum(data_df$group == Group2)," regional pairs)")) +
              theme_cowplot()+
              theme(
                axis.title.x = element_text(size = 22),
                axis.title.y = element_text(size = 22),
                axis.text.x = element_text(size = 22),
                axis.text.y = element_text(size = 18),
                legend.position = "right",
                legend.text = element_text(size = 18),
                legend.title = element_text(size = 22)) + 
              annotate("text",size = 7, x = 1.5, y = (min(data$var) + y_range*1.1),
                       label = paste("p=",pval,"(lme)"))
            
            print(g)
            
          } #if(length(unique(data$tumour_id))>0)
        } #for(j in (i+1):length(group.var_list))
      } #for(i in 1:(length(group.var_list)-1)) {
    } #if(length(group.var_list) > 1) {
  } #for (var in var_list) {
}#for (group in group_column_list){
dev.off()



#######################################################################
## compare dominant TMB (CCF>=95%) between A-like vs D-like (by LOH) ##
#######################################################################

# organise dominant tmb, purity and growth pattern for annotation

tmb_ccf95_df <- readRDS(paste0(tmpdir, "/tmb_ccf95.RDS"))
tmb_ccf95_df$Region_ID_cruk <- paste0(substr(tmb_ccf95_df$sample, 1,8),"_",substr(tmb_ccf95_df$sample, 13,17) )
tmb_ccf95_df$TMB_CCF95_log = log10(as.numeric(as.character(tmb_ccf95_df$tmb_95)) + 1 )

y_var <- "TMB_CCF95_log"

add_df <- LUSC_region_df %>%
  dplyr::select(tumour_id_muttable_cruk, region ,Region_ID_cruk) %>%
  left_join(tmb_ccf95_df, by = "Region_ID_cruk") %>%
  mutate(y_var = TMB_CCF95_log) %>%
  filter(!is.na(y_var)) %>%
  dplyr::select(tumour_id_muttable_cruk, y_var, region) %>%
  dplyr::rename(region_add = region, 
                tumour_id = tumour_id_muttable_cruk)

pattern_df <- LUSC_region_df[,c("tumour_id_muttable_cruk","region","regional_pattern.new")] %>%
  dplyr::rename(tumour_id = tumour_id_muttable_cruk) 

## plot paired dot plot, with permutation test 

ratio_filter_p_list <- c(0.02)
ratio_filter_d_list <- c(0.1)

out_pdf <- paste0(workdir,"/Fig.S4a.dominant.mut.0.02_0.1_TMBccf95_simple.pdf")
pdf(out_pdf, width = 5, height = 6)
for(n in 1:length(ratio_filter_p_list)){
  for(m in 1:length(ratio_filter_d_list)){
    
    ratio_filter_p <- ratio_filter_p_list[n]
    ratio_filter_d <- ratio_filter_d_list[m]
    
    if(ratio_filter_p <= ratio_filter_d) {
      
      all_joined_comb_branch.prop <- readRDS(paste0(tmpdir, "/all_joined_comb_branch.prop_ratio_filter_p",ratio_filter_p,"_ratio_filter_d",ratio_filter_d,".tree.ratio.RDS")) 
      
      # apply purity cutoff (>= 0.15)
      all_joined_comb_branch.prop <- all_joined_comb_branch.prop %>%
        filter(purity_1st >= 0.15, purity_2nd >= 0.15)
      
      ## call ancestor-descendant relation using LOH profile
      
      # step 1 extract preliminary parent-daughter pairs based on cytoband LOH -> add arm filters : shared_LOH_arm >= 1 and  daughter's private LOH arm >= 1
      p.d_comb_LOH_cyto_filtered <-  all_joined_comb_branch.prop$p.d_comb_LOH_cyto[!grepl("NA",all_joined_comb_branch.prop$p.d_comb_LOH_cyto) 
                                                                                   & !is.na(all_joined_comb_branch.prop$p.d_comb_LOH_cyto)
                                                                                   & all_joined_comb_branch.prop$shared_LOH_arm >= 1
                                                                                   & ifelse(grepl("region.1st",all_joined_comb_branch.prop$p.d_comb_LOH_cyto), all_joined_comb_branch.prop$private_loh_2nd_LOH_arm, 
                                                                                            ifelse(grepl("region.2nd",all_joined_comb_branch.prop$p.d_comb_LOH_cyto), all_joined_comb_branch.prop$private_loh_1st_LOH_arm, 
                                                                                                   NA_real_)) >=1
      ]
      
      # add growth pattern
      
      p.d_comb_LOH_cyto_df <- add_growth_pattern(all_joined_comb_branch.prop[all_joined_comb_branch.prop$p.d_comb_LOH_cyto %in%  p.d_comb_LOH_cyto_filtered , ],
                                                 pattern_df) %>%
        mutate(p.d_comb = p.d_comb_LOH_cyto) %>%
        mutate(parent_region =  case_when(grepl("region.1st",p.d_comb) ~ "region.1st",
                                          grepl("region.2nd",p.d_comb) ~ "region.2nd")) %>%
        mutate(parent_grade = case_when(parent_region == "region.1st" ~ grade_1st,
                                        parent_region == "region.2nd" ~ grade_2nd),
               daughter_grade = case_when(parent_region == "region.1st" ~ grade_2nd,
                                          parent_region == "region.2nd" ~ grade_1st)) %>%
        mutate(parent_grade.num = case_when(parent_region == "region.1st" ~ grade.num_1st,
                                            parent_region == "region.2nd" ~ grade.num_2nd),
               daughter_grade.num = case_when(parent_region == "region.1st" ~ grade.num_2nd,
                                              parent_region == "region.2nd" ~ grade.num_1st)) %>%
        mutate(parent_pattern = case_when(parent_region == "region.1st" ~ pattern_1st,
                                          parent_region == "region.2nd" ~ pattern_2nd),
               daughter_pattern = case_when(parent_region == "region.1st" ~ pattern_2nd,
                                            parent_region == "region.2nd" ~ pattern_1st)) %>%
        mutate(grade_shift = case_when(parent_grade.num < daughter_grade.num ~ "up",
                                       parent_grade.num > daughter_grade.num ~ "down",
                                       parent_grade.num == daughter_grade.num ~ "same")) %>%
        mutate(grade_shift = factor(grade_shift, levels =c("up","down","same" ))) %>%
        mutate(parent_grade = factor(parent_grade, levels = c("high" ,"low")),
               daughter_grade = factor(daughter_grade, levels = c("high" ,"low")))
      
      
      p.d_comb_list <- list(p.d_comb_LOH_cyto_df)
      
      title_list <- c("LOH cytoband")
      
      for (k in 1:length(p.d_comb_list)){
        
        ## create paired dot plot 
        
        p.d_comb_df <- p.d_comb_list[[k]]
        
        # add y_var 
        p.d_comb_longer <- p.d_comb_df %>%
          dplyr::mutate(region_add = parent_R_LOH_cyto) %>%
          inner_join(add_df, by = c("tumour_id", "region_add")) %>%
          dplyr::rename(parent_var = y_var) %>%
          dplyr::mutate(region_add = daughter_R_LOH_cyto) %>%
          inner_join(add_df, by = c("tumour_id", "region_add")) %>%
          dplyr::rename(daughter_var = y_var) %>%
          pivot_longer(cols = parent_var:daughter_var, values_to = "y_var", names_to = "parent_daughter") %>%
          mutate(parent_daughter = gsub("_var", "",parent_daughter )) %>%
          mutate(parent_daughter = factor(parent_daughter, levels = c("parent", "daughter"))) 
        
        
        data_df.pre <- p.d_comb_longer %>%
          filter(!is.na(y_var) ) 
        
        paired_comb <- data_df.pre %>% group_by(p.d_comb) %>% summarise(count = n()) %>% filter(count == 2)
        keep_comb <- paired_comb$p.d_comb
        
        data_df <- data_df.pre %>% filter(p.d_comb %in%  keep_comb )
        
        
        data_df_wider <- data_df  %>%
          dplyr::select(tumour_id, p.d_comb, y_var, parent_daughter) %>%
          pivot_wider(names_from = "parent_daughter", values_from = "y_var") %>%
          mutate(delta_var = case_when(daughter - parent > 0 ~ "up",
                                       daughter - parent == 0  ~ "same",
                                       daughter - parent < 0 ~ "down")) 
        
        # add up/same/down to data_df
        up_p.d_comb <- data_df_wider$p.d_comb[data_df_wider$delta_var == "up"]
        same_p.d_comb <- data_df_wider$p.d_comb[data_df_wider$delta_var == "same"]
        down_p.d_comb <- data_df_wider$p.d_comb[data_df_wider$delta_var == "down"]
        
        data_df <- data_df %>%
          mutate(shift_var = case_when(p.d_comb %in% up_p.d_comb ~ "A<D",
                                       p.d_comb %in% same_p.d_comb ~ "A=D",
                                       p.d_comb %in% down_p.d_comb ~ "A>D") ) %>%
          mutate(shift_var = factor(shift_var , levels = c("A<D","A=D", "A>D")))
        
        
        # some numbers adding to plot
        num_tum <- length(unique(data_df_wider$tumour_id))
        num_regional_pair <- nrow(data_df_wider)
        y_range <- max(data_df$y_var) - min(data_df$y_var)
        
        
        ## permutation test : is up shift pairs (dominant tmb higher in descendant-like regions) significantly enriched than by chance?  
        
        # observed frequency 
        num_up <- sum(data_df_wider$delta_var  == "up", na.rm = T)
        num_down <- sum(data_df_wider$delta_var == "down", na.rm = T)
        num_same <-  sum(data_df_wider$delta_var  == "same", na.rm = T)
        
        # run 1000 times -> check empirical pvalue 
        
        all_delta_var_test_down <- vector()
        all_delta_var_test_up <- vector()
        
        set.seed(123)
        
        for(i in 1:permute_num){
          print(i)
          
          # For each permutation, shuffle y_var within each tumour
          
          df_anno <- add_df 
          
          shuffled_tmp <- sapply(unique(df_anno$tumour_id), function(x){
            vec_tum <- df_anno$y_var[which(df_anno$tumour_id == x)]
            if(length(vec_tum) > 1){
              vec_tum_rand <- sample(vec_tum)
            } else {
              vec_tum_rand <- vec_tum}
            return(vec_tum_rand)
          })
          shuffled_df <- unlist(shuffled_tmp)
          
          add_df_rand  <- cbind(df_anno[,c("tumour_id", "region_add")],
                                data.frame(y_var = shuffled_df)) %>%
            filter(!is.na(y_var)) 
          
          # use shuffled y_var data (add_df_rand) and count number of "up" shift 
          
          p.d_comb_longer_random <- p.d_comb_df %>%
            dplyr::mutate(region_add = parent_R_LOH_cyto) %>%
            inner_join(add_df_rand, by = c("tumour_id", "region_add")) %>%
            dplyr::rename(parent_var = y_var) %>%
            dplyr::mutate(region_add = daughter_R_LOH_cyto) %>%
            inner_join(add_df_rand, by = c("tumour_id", "region_add")) %>%
            dplyr::rename(daughter_var = y_var) %>%
            pivot_longer(cols = parent_var:daughter_var, values_to = "y_var", names_to = "parent_daughter") %>%
            mutate(parent_daughter = gsub("_var", "",parent_daughter )) %>%
            mutate(parent_daughter = factor(parent_daughter, levels = c("parent", "daughter"))) 
          
          data_df.pre_random <- p.d_comb_longer_random %>%
            filter(!is.na(y_var) ) 
          
          paired_comb_random <- data_df.pre_random %>% group_by(p.d_comb) %>% summarise(count = n()) %>% filter(count == 2)
          keep_comb_random <- paired_comb_random$p.d_comb
          
          data_df_random <- data_df.pre_random %>% filter(p.d_comb %in%  keep_comb_random )
          
          
          data_df_wider_random <- data_df_random  %>%
            dplyr::select(tumour_id, p.d_comb, y_var, parent_daughter) %>%
            pivot_wider(names_from = "parent_daughter", values_from = "y_var") %>%
            mutate(delta_var = case_when(daughter - parent > 0 ~ "up",
                                         daughter - parent == 0  ~ "same",
                                         daughter - parent < 0 ~ "down")) 
          
          #table(p.d_comb_random$grade_shift)
          num_up_rand <- sum(data_df_wider_random$delta_var == "up", na.rm = T)
          num_down_rand <- sum(data_df_wider_random$delta_var == "down", na.rm = T)
          num_same_rand <-  sum(data_df_wider_random$delta_var == "same", na.rm = T)
          
          delta_var_count <- c(num_up,num_down,num_same) 
          
          # compare with observed count 
          delta_var_test_down <- ifelse(num_down_rand <= num_down, 1 , 0)
          delta_var_test_up <- ifelse(num_up_rand >= num_up, 1 , 0)
          
          all_delta_var_test_down <- c(all_delta_var_test_down, delta_var_test_down)
          all_delta_var_test_up <- c(all_delta_var_test_up, delta_var_test_up)
        }
        
        p_emp_down <- signif( (sum(all_delta_var_test_down) +1)/(permute_num+1), 3)
        p_emp_up <- signif( (sum(all_delta_var_test_up) +1)/(permute_num+1), 3)
        
        #my_col <- createPalette(length(het_tum),  c("#ff0000", "#00ff00", "#0000ff"))
        title_string <- paste0(title_list[k], "\n",
                               "ratio_filter_p=",ratio_filter_p," ratio_filter_d=",ratio_filter_d ,"\n" 
                               ,"n=",num_tum," tumours (",num_regional_pair," regional_pairs)\n",
                               "P.emp(up)=",p_emp_up," (up:",num_up ,", down:",num_down,")"
        )
        
        
        ## dot plot 
        g1 <-  ggplot(data_df %>% mutate(parent_daughter_jitter = jitter(as.numeric(parent_daughter)), y_var_jitter = jitter(as.numeric(y_var), amount = 0) ), 
                      aes(x = parent_daughter_jitter, y = y_var_jitter))+
          geom_point(color = "black", fill = "gray40", alpha = 0.8, size = 8, pch=21)+  #color = organ, shape = sample_timing
          geom_line(aes(group = p.d_comb, color = shift_var))+
          scale_color_manual(name = "Comparison", values = c("A<D" = "#F8766D", "A=D" = "#00BA38", "A>D" = "#619CFF"), drop = T) +
          # stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
          #              geom = "crossbar", width = 0.5, color = "red",
          #              aes(x = as.numeric(parent_daughter), group = parent_daughter)) +
          ggtitle(title_string) +
          labs(x = "Regions",  y = gsub("_"," ", y_var) ) +
          theme_cowplot()+
          theme(plot.title = element_text(size = 14),
                axis.title.x = element_text(size = 16),
                axis.title.y = element_text(size = 16),
                axis.text.x = element_text(size = 16),
                axis.text.y = element_text(size = 16),
                legend.position = "right",
                legend.text = element_text(size = 16),
                legend.title = element_text(size = 18)) + 
          scale_x_continuous(breaks = c(1,2), 
                             labels=c("Ancestor\n-like (A)", "Descendant\n-like (D)")) 
        print(g1)
      }
    }
  }
}
dev.off()


###################################################################
## ancestor-descendant-like pairs (filtered) -> grade transition ##  
###################################################################
## exclude conflicting call (ancestor-descendant-like by LOH, but descendant-ancestor-like by mutation)

# create dataframe of regional growht pattern
pattern_df <- LUSC_region_df[,c("tumour_id_muttable_cruk","region","regional_pattern.new")] %>%
  dplyr::rename(tumour_id = tumour_id_muttable_cruk) 

# load pattern_df_random : growth pattern shuffled within each tumour x1000 - this can be used for Monte-Carlo simulation 
pattern_df_random <- readRDS(paste0(tmpdir, "/pattern_df_random.RDS"))

## hetero grade tumours 
multi_grade <- LUSC_region_df %>%
  mutate(regional_grade = pattern_to_grade(regional_pattern.new)) %>%
  filter(regional_grade %in% c("low","high")) %>%
  group_by(tumour_id_muttable_cruk, regional_grade) %>%
  summarise(count = n()) %>%
  group_by(tumour_id_muttable_cruk) %>%
  summarise(count = n()) %>%
  filter(count > 1) 

### LOH - include homogenous grade tumours 
ratio_filter_p_list <- c(0.02)
ratio_filter_d_list <- c(0.1)

out_pdf <- paste0(workdir,"/Fig.S4b.LOH_cyto_with_mut_filter_all_permute.pdf")
pdf(out_pdf, width = 7, height = 6)

for(n in 1:length(ratio_filter_p_list)){
  for(m in 1:length(ratio_filter_d_list)){
    
    ratio_filter_p <- ratio_filter_p_list[n]
    ratio_filter_d <- ratio_filter_d_list[m]
    
    if(ratio_filter_p <= ratio_filter_d) {
      
      all_joined_comb_branch.prop <- readRDS(paste0(tmpdir, "/all_joined_comb_branch.prop_ratio_filter_p",ratio_filter_p,"_ratio_filter_d",ratio_filter_d,".tree.ratio.RDS")) 
      
      # apply purity cutoff (>= 0.15) + only hetero grade tumours
      all_joined_comb_branch.prop <- all_joined_comb_branch.prop %>%
        filter(purity_1st >= 0.15, purity_2nd >= 0.15) 
      
      # step 1 extract p-d pairs based on cytoband LOH - add arm filters : shared_LOH_arm >= 1 and  daughter's private LOH arm >= 1
      p.d_comb_LOH_cyto <-  all_joined_comb_branch.prop$p.d_comb_LOH_cyto[!grepl("NA",all_joined_comb_branch.prop$p.d_comb_LOH_cyto) 
                                                                          & !is.na(all_joined_comb_branch.prop$p.d_comb_LOH_cyto)
                                                                          & all_joined_comb_branch.prop$shared_LOH_arm >= 1
                                                                          & ifelse(grepl("region.1st",all_joined_comb_branch.prop$p.d_comb_LOH_cyto), all_joined_comb_branch.prop$private_loh_2nd_LOH_arm, 
                                                                                   ifelse(grepl("region.2nd",all_joined_comb_branch.prop$p.d_comb_LOH_cyto), all_joined_comb_branch.prop$private_loh_1st_LOH_arm, 
                                                                                          NA_real_)) >=1
                                                                          #(all_joined_comb_branch.prop$private_loh_1st_LOH_arm + all_joined_comb_branch.prop$private_loh_2nd_LOH_arm) >= 1
      ]
      
      # step 2 exclude p-d pairs called as opposite in dominant mutation method (mean cluster CCF95)
      p.d_comb_mut_ccf95 <-  all_joined_comb_branch.prop$p.d_comb_mut_ccf95[!grepl("NA",all_joined_comb_branch.prop$p.d_comb_mut_ccf95)
                                                                            & !is.na(all_joined_comb_branch.prop$p.d_comb_mut_ccf95)]
      
      p.d_comb_mut_ccf95_reverse <- ifelse(grepl("region.1st",p.d_comb_mut_ccf95), gsub("region.1st", "region.2nd", p.d_comb_mut_ccf95), gsub("region.2nd", "region.1st",p.d_comb_mut_ccf95))
      
      p.d_comb_LOH_cyto_filtered <- setdiff(p.d_comb_LOH_cyto, p.d_comb_mut_ccf95_reverse)  #p.d_comb_mut_clonal_reverse,
      
      
      # add growth pattern
      
      p.d_comb_LOH_cyto_df <- add_growth_pattern(all_joined_comb_branch.prop[all_joined_comb_branch.prop$p.d_comb_LOH_cyto %in%  p.d_comb_LOH_cyto_filtered , ],
                                                 pattern_df) %>%
        mutate(p.d_comb = p.d_comb_LOH_cyto) %>%
        mutate(parent_region =  case_when(grepl("region.1st",p.d_comb) ~ "region.1st",
                                          grepl("region.2nd",p.d_comb) ~ "region.2nd")) %>%
        mutate(parent_grade = case_when(parent_region == "region.1st" ~ grade_1st,
                                        parent_region == "region.2nd" ~ grade_2nd),
               daughter_grade = case_when(parent_region == "region.1st" ~ grade_2nd,
                                          parent_region == "region.2nd" ~ grade_1st)) %>%
        mutate(parent_grade.num = case_when(parent_region == "region.1st" ~ grade.num_1st,
                                            parent_region == "region.2nd" ~ grade.num_2nd),
               daughter_grade.num = case_when(parent_region == "region.1st" ~ grade.num_2nd,
                                              parent_region == "region.2nd" ~ grade.num_1st)) %>%
        mutate(parent_pattern = case_when(parent_region == "region.1st" ~ pattern_1st,
                                          parent_region == "region.2nd" ~ pattern_2nd),
               daughter_pattern = case_when(parent_region == "region.1st" ~ pattern_2nd,
                                            parent_region == "region.2nd" ~ pattern_1st)) %>%
        mutate(grade_shift = case_when(parent_grade.num < daughter_grade.num ~ "up",
                                       parent_grade.num > daughter_grade.num ~ "down",
                                       parent_grade.num == daughter_grade.num ~ "same")) %>%
        mutate(grade_shift = factor(grade_shift, levels =c("up","down","same" ))) %>%
        mutate(parent_grade = factor(parent_grade, levels = c("high", "low")),
               daughter_grade = factor(daughter_grade, levels = c("high", "low")))
      
      p.d_comb_list <- list(p.d_comb_LOH_cyto_df) 
      
      title_list <- c("LOH cytoband")
      
      for (k in 1:length(p.d_comb_list)){
        
        ## create dot plots and alluvial plots and run permutation test
        
        p.d_comb_df <- p.d_comb_list[[k]]
        
        p.d_comb_freq <- p.d_comb_df %>%
          filter(!is.na(parent_grade.num), !is.na(daughter_grade.num) ) %>%
          group_by(parent_grade, daughter_grade, grade_shift) %>%
          summarise(Freq = n())
        
        p.d_comb_longer <- p.d_comb_df %>%
          filter(!is.na(parent_grade.num), !is.na(daughter_grade.num) ) %>%
          pivot_longer(cols = parent_grade.num:daughter_grade.num, values_to = "grade_score", names_to = "parent_daughter") %>%
          mutate(parent_daughter = gsub("_grade.num", "",parent_daughter )) %>%
          mutate(parent_daughter = factor(parent_daughter, levels = c("parent", "daughter"))) %>%
          mutate(regional_pattern = case_when(parent_daughter == "parent" ~ parent_pattern, 
                                              parent_daughter == "daughter" ~ daughter_pattern))
        
        data_df <- p.d_comb_longer 
        
        num_tum <- length(unique(data_df$tumour_id))
        num_regional_pair <- nrow(data_df)/2
        y_range <- max(data_df$grade_score) - min(data_df$grade_score)
        
        ## permutation test : is upward shift significantly enriched than by chance?  
        
        # observed frequency 
        num_up <- sum(p.d_comb_df$grade_shift == "up", na.rm = T)
        num_down <- sum(p.d_comb_df$grade_shift == "down", na.rm = T)
        num_same <-  sum(p.d_comb_df$grade_shift == "same", na.rm = T)
        
        p.d_comb_tmp <- p.d_comb_df %>%
          dplyr::select(-grade_1st, -grade_2nd, -grade.num_1st, -grade.num_2nd, -pattern_1st, -pattern_2nd,
                        -parent_grade, -daughter_grade, -parent_grade.num,-daughter_grade.num ,-parent_pattern, -daughter_pattern, -grade_shift) 
        
        # run 1000 times -> check empirical pvalue 
        permute_num = 1000
        
        all_grade_shift_test_down <- vector()
        all_grade_shift_test_up <- vector()
        for(i in 1:permute_num){
          print(i)
          
          p.d_comb_random <- add_growth_pattern(p.d_comb_tmp, pattern_df_random[[i]]) %>%  # now new regional patterns are assigned
            mutate(parent_grade = case_when(parent_region == "region.1st" ~ grade_1st,
                                            parent_region == "region.2nd" ~ grade_2nd),
                   daughter_grade = case_when(parent_region == "region.1st" ~ grade_2nd,
                                              parent_region == "region.2nd" ~ grade_1st)) %>%
            mutate(parent_grade.num = case_when(parent_region == "region.1st" ~ grade.num_1st,
                                                parent_region == "region.2nd" ~ grade.num_2nd),
                   daughter_grade.num = case_when(parent_region == "region.1st" ~ grade.num_2nd,
                                                  parent_region == "region.2nd" ~ grade.num_1st)) %>%
            mutate(parent_pattern = case_when(parent_region == "region.1st" ~ pattern_1st,
                                              parent_region == "region.2nd" ~ pattern_2nd),
                   daughter_pattern = case_when(parent_region == "region.1st" ~ pattern_2nd,
                                                parent_region == "region.2nd" ~ pattern_1st)) %>%
            mutate(grade_shift = case_when(parent_grade.num < daughter_grade.num ~ "up",
                                           parent_grade.num > daughter_grade.num ~ "down",
                                           parent_grade.num == daughter_grade.num ~ "same")) %>%
            mutate(grade_shift = factor(grade_shift, levels =c("up","down","same" ))) %>%
            mutate(parent_grade = factor(parent_grade, levels = c("high", "low")),
                   daughter_grade = factor(daughter_grade, levels = c("high", "low")))
          
          #table(p.d_comb_random$grade_shift)
          num_up_rand <- sum(p.d_comb_random$grade_shift == "up", na.rm = T)
          num_down_rand <- sum(p.d_comb_random$grade_shift == "down", na.rm = T)
          num_same_rand <-  sum(p.d_comb_random$grade_shift == "same", na.rm = T)
          
          grade_shift_count <- c(num_up,num_down,num_same) 
          
          # compare with observed count (fixed!)  sum up the number of simulations for which the statistic was greater (less) than or equal to the observed statistic,
          grade_shift_test_down <- ifelse(num_down_rand <= num_down, 1 , 0)
          grade_shift_test_up <- ifelse(num_up_rand >= num_up, 1 , 0)
          
          all_grade_shift_test_down <- c(all_grade_shift_test_down, grade_shift_test_down)
          all_grade_shift_test_up <- c(all_grade_shift_test_up, grade_shift_test_up)
        }
        
        p_emp_down <- signif( (sum(all_grade_shift_test_down) +1)/(permute_num + 1) , 3)
        p_emp_up <- signif( (sum(all_grade_shift_test_up) + 1)/(permute_num + 1) , 3)
        
        title_string <- paste0(title_list[k], "\n",
                               "ratio_filter_p=",ratio_filter_p," ratio_filter_d=",ratio_filter_d ,"\n" 
                               ,"n=",num_tum," tumours (",num_regional_pair," regional_pairs)\n",
                               "P.emp(up)=",p_emp_up,"(up:",num_up," same:",num_same," down:",num_down,")")
        
        ## alluvial plot 
        g1 <- ggplot(data =  p.d_comb_freq,
                     aes(axis1 = parent_grade, axis2 = daughter_grade,
                         y = Freq)) +
          scale_x_discrete(limits = c("Ancestor\n-like", "Descendant\n-like")) + #, expand = c(.2, .05)) +
          labs(x ="Regions", fill = "Transition") +
          geom_alluvium(aes(fill = grade_shift)) +
          scale_fill_manual(values = c("up" = "#F8766D", "same" = "#00BA38", "down" = "#619CFF"), drop = T) +
          geom_stratum() +
          geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 6) +
          theme_cowplot() +
          ggtitle(title_string) +
          theme(plot.title = element_text(size = 14),
                axis.title.x = element_text(size = 16),
                axis.title.y = element_blank(),
                axis.text.x = element_text(size = 16),
                axis.text.y = element_blank(),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_blank(),
                axis.line.x = element_blank(),
                axis.line.y = element_blank(),
                legend.position = "right",
                legend.text = element_text(size = 16),
                legend.title = element_text(size = 18))
        
        print(g1)
        
        ## dot plot 
        g2 <-  ggplot(data_df %>% mutate(parent_daughter_jitter = jitter(as.numeric(parent_daughter)), grade_score_jitter = jitter(as.numeric(grade_score), amount = 0) ), 
                      aes(x = parent_daughter_jitter, y = grade_score_jitter))+
          geom_point(color = "black", aes(fill = regional_pattern), alpha = 0.8, size = 8, pch=21)+ 
          geom_line(aes(group = p.d_comb))+
          stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                       geom = "crossbar", width = 0.5, color = "red",
                       aes(x = as.numeric(parent_daughter), group = parent_daughter)) +
          ggtitle(title_string) +
          labs(x = "Regions",  y = "Pattern grade score", fill = "Regional pattern") +
          theme_cowplot()+
          theme(plot.title = element_text(size = 14),
                axis.title.x = element_text(size = 16),
                axis.title.y = element_text(size = 16),
                axis.text.x = element_text(size = 16),
                axis.text.y = element_text(size = 16),
                legend.position = "right",
                legend.text = element_text(size = 16),
                legend.title = element_text(size = 18)) + 
          scale_fill_manual(values = c(  "well_dx" = "#A6CEE3", "papillary" = "#1F78B4", "acinar" = "#B2DF8A",
                                         "cribriform" = "#33A02C","micropapillary" = "#FB9A99", "poor_dx" = "#E31A1C"), drop = T) +
          scale_x_continuous(breaks = c(1,2), 
                             labels=c("Ancestor\n-like", "Descendant\n-like")) 
        
        print(g2)
      }
    }
  }
}
dev.off()


### LOH - include homogenous grade tumours , alluvial plot with all growth patterns
ratio_filter_p_list <- c(0.02)
ratio_filter_d_list <- c(0.1)

out_pdf <- paste0(workdir,"/Fig.S4b.LOH_cyto_with_mut_filter_each.pattern_alluvial.pdf")

pdf(out_pdf, width = 7, height = 6)

for(n in 1:length(ratio_filter_p_list)){
  for(m in 1:length(ratio_filter_d_list)){
    
    ratio_filter_p <- ratio_filter_p_list[n]
    ratio_filter_d <- ratio_filter_d_list[m]
    
    if(ratio_filter_p <= ratio_filter_d) {
      
      all_joined_comb_branch.prop <- readRDS(paste0(tmpdir, "/all_joined_comb_branch.prop_ratio_filter_p",ratio_filter_p,"_ratio_filter_d",ratio_filter_d,".tree.ratio.RDS")) 
      
      # apply purity cutoff (>= 0.15)
      all_joined_comb_branch.prop <- all_joined_comb_branch.prop %>%
        filter(purity_1st >= 0.15, purity_2nd >= 0.15)
      
      # step 1 extract p-d pairs based on cytoband LOH - add arm filters : shared_LOH_arm >= 1 and  daughter's private LOH arm >= 1
      p.d_comb_LOH_cyto <-  all_joined_comb_branch.prop$p.d_comb_LOH_cyto[!grepl("NA",all_joined_comb_branch.prop$p.d_comb_LOH_cyto) 
                                                                          & !is.na(all_joined_comb_branch.prop$p.d_comb_LOH_cyto)
                                                                          & all_joined_comb_branch.prop$shared_LOH_arm >= 1
                                                                          & ifelse(grepl("region.1st",all_joined_comb_branch.prop$p.d_comb_LOH_cyto), all_joined_comb_branch.prop$private_loh_2nd_LOH_arm, 
                                                                                   ifelse(grepl("region.2nd",all_joined_comb_branch.prop$p.d_comb_LOH_cyto), all_joined_comb_branch.prop$private_loh_1st_LOH_arm, 
                                                                                          NA_real_)) >=1
      ]
      
      # step 2 exclude p-d pairs called as opposite in mutation method (mean cluster CCF95)
      p.d_comb_mut_ccf95 <-  all_joined_comb_branch.prop$p.d_comb_mut_ccf95[!grepl("NA",all_joined_comb_branch.prop$p.d_comb_mut_ccf95) 
                                                                            & !is.na(all_joined_comb_branch.prop$p.d_comb_mut_ccf95)]
      
      p.d_comb_mut_ccf95_reverse <- ifelse(grepl("region.1st",p.d_comb_mut_ccf95), gsub("region.1st", "region.2nd", p.d_comb_mut_ccf95), gsub("region.2nd", "region.1st",p.d_comb_mut_ccf95))
      
      p.d_comb_LOH_cyto_filtered <- setdiff(p.d_comb_LOH_cyto,  p.d_comb_mut_ccf95_reverse ) 
      
      # add growth pattern
      
      p.d_comb_LOH_cyto_df <- add_growth_pattern(all_joined_comb_branch.prop[all_joined_comb_branch.prop$p.d_comb_LOH_cyto %in%  p.d_comb_LOH_cyto_filtered , ],
                                                 pattern_df) %>%
        mutate(p.d_comb = p.d_comb_LOH_cyto) %>%
        mutate(parent_region =  case_when(grepl("region.1st",p.d_comb) ~ "region.1st",
                                          grepl("region.2nd",p.d_comb) ~ "region.2nd")) %>%
        mutate(parent_grade = case_when(parent_region == "region.1st" ~ grade_1st,
                                        parent_region == "region.2nd" ~ grade_2nd),
               daughter_grade = case_when(parent_region == "region.1st" ~ grade_2nd,
                                          parent_region == "region.2nd" ~ grade_1st)) %>%
        mutate(parent_grade.num = case_when(parent_region == "region.1st" ~ grade.num_1st,
                                            parent_region == "region.2nd" ~ grade.num_2nd),
               daughter_grade.num = case_when(parent_region == "region.1st" ~ grade.num_2nd,
                                              parent_region == "region.2nd" ~ grade.num_1st)) %>%
        mutate(parent_pattern = case_when(parent_region == "region.1st" ~ pattern_1st,
                                          parent_region == "region.2nd" ~ pattern_2nd),
               daughter_pattern = case_when(parent_region == "region.1st" ~ pattern_2nd,
                                            parent_region == "region.2nd" ~ pattern_1st)) %>%
        mutate(grade_shift = case_when(parent_grade.num < daughter_grade.num ~ "up",
                                       parent_grade.num > daughter_grade.num ~ "down",
                                       parent_grade.num == daughter_grade.num ~ "same")) %>%
        mutate(grade_shift = factor(grade_shift, levels =c("up","down","same" ))) %>%
        mutate(parent_grade = factor(parent_grade, levels = c("high", "low")),
               daughter_grade = factor(daughter_grade, levels = c("high", "low")))
      
      p.d_comb_list <- list(p.d_comb_LOH_cyto_df) #,
      
      title_list <- c("LOH cytoband") #,
      
      for (k in 1:length(p.d_comb_list)){
        
        ## create dot plots and alluvial plots and run lme (although statistically not ideal...)
        
        p.d_comb_df <- p.d_comb_list[[k]]
        
        p.d_comb_freq <- p.d_comb_df %>%
          filter(!is.na(parent_grade.num), !is.na(daughter_grade.num) ) %>%
          group_by(parent_pattern, daughter_pattern, grade_shift) %>%
          summarise(Freq = n()) %>%
          mutate(parent_pattern = factor(parent_pattern, levels = rev(subtype_list)),
                 daughter_pattern = factor(daughter_pattern, levels = rev(subtype_list)))
        
        p.d_comb_longer <- p.d_comb_df %>%
          filter(!is.na(parent_grade.num), !is.na(daughter_grade.num) ) %>%
          pivot_longer(cols = parent_grade.num:daughter_grade.num, values_to = "grade_score", names_to = "parent_daughter") %>%
          mutate(parent_daughter = gsub("_grade.num", "",parent_daughter )) %>%
          mutate(parent_daughter = factor(parent_daughter, levels = c("parent", "daughter"))) %>%
          mutate(regional_pattern = case_when(parent_daughter == "parent" ~ parent_pattern, 
                                              parent_daughter == "daughter" ~ daughter_pattern))
        
        # save 
        # write.csv(p.d_comb_longer, paste0(out_pdf,"_",ratio_filter_p,"_",ratio_filter_d, ".full.csv"), row.names = F )
        
        # prep for plot 
        data_df <- p.d_comb_longer 
        
        num_tum <- length(unique(data_df$tumour_id))
        num_regional_pair <- nrow(data_df)/2
        
        title_string <- paste0(title_list[k], "\n",
                               "ratio_filter_p=",ratio_filter_p," ratio_filter_d=",ratio_filter_d ,"\n" 
                               ,"n=",num_tum," tumours (",num_regional_pair," regional_pairs)\n" ) 
        
        ## alluvial plot 
        g <- ggplot(data =  p.d_comb_freq,
                    aes(axis1 = parent_pattern, axis2 = daughter_pattern,
                        y = Freq)) +
          scale_x_discrete(limits = c("Ancestor\n-like", "Descendant\n-like")) + #, expand = c(.2, .05)) +
          labs(x ="Regions", fill = "Transition") +
          geom_alluvium(aes(fill = grade_shift)) +
          scale_fill_manual(values = c("up" = "#F8766D", "same" = "#00BA38", "down" = "#619CFF"), drop = T) +
          geom_stratum() +
          geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 6) +
          theme_cowplot() +
          ggtitle(title_string) +
          theme(plot.title = element_text(size = 14),
                axis.title.x = element_text(size = 16),
                axis.title.y = element_blank(),
                axis.text.x = element_text(size = 16),
                axis.text.y = element_blank(),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_blank(),
                axis.line.x = element_blank(),
                axis.line.y = element_blank(),
                legend.position = "right",
                legend.text = element_text(size = 16),
                legend.title = element_text(size = 18))
        
        print(g)
      }
    }
  }
}

dev.off()


### LOH  - use various cutoffs - only plot hetero grade tumours
ratio_filter_p_list <- c(0.01, 0.02, 0.05)
ratio_filter_d_list <- c(0.075, 0.1, 0.15)

#Fig 2b
out_pdf <- paste0(workdir,"/Fig.2d.S4c.LOH_cyto_with_mut_filter_hetero_permute_full.cutoffs.pdf")
pdf(out_pdf, width = 7, height = 6)

for(n in 1:length(ratio_filter_p_list)){
  for(m in 1:length(ratio_filter_d_list)){
    
    ratio_filter_p <- ratio_filter_p_list[n]
    ratio_filter_d <- ratio_filter_d_list[m]
    
    if(ratio_filter_p <= ratio_filter_d) {
      
      all_joined_comb_branch.prop <- readRDS(paste0(tmpdir, "/all_joined_comb_branch.prop_ratio_filter_p",ratio_filter_p,"_ratio_filter_d",ratio_filter_d,".tree.ratio.RDS")) 
      
      # apply purity cutoff (>= 0.15) + only hetero grade tumours
      all_joined_comb_branch.prop <- all_joined_comb_branch.prop %>%
        filter(purity_1st >= 0.15, purity_2nd >= 0.15) %>%
        filter(tumour_id %in% multi_grade$tumour_id_muttable_cruk) 
      
      # step 1 extract p-d pairs based on cytoband LOH - add arm filters : shared_LOH_arm >= 1 and  daughter's private LOH arm >= 1
      p.d_comb_LOH_cyto <-  all_joined_comb_branch.prop$p.d_comb_LOH_cyto[!grepl("NA",all_joined_comb_branch.prop$p.d_comb_LOH_cyto) 
                                                                          & !is.na(all_joined_comb_branch.prop$p.d_comb_LOH_cyto)
                                                                          & all_joined_comb_branch.prop$shared_LOH_arm >= 1
                                                                          & ifelse(grepl("region.1st",all_joined_comb_branch.prop$p.d_comb_LOH_cyto), all_joined_comb_branch.prop$private_loh_2nd_LOH_arm, 
                                                                                   ifelse(grepl("region.2nd",all_joined_comb_branch.prop$p.d_comb_LOH_cyto), all_joined_comb_branch.prop$private_loh_1st_LOH_arm, 
                                                                                          NA_real_)) >=1
                                                                          #(all_joined_comb_branch.prop$private_loh_1st_LOH_arm + all_joined_comb_branch.prop$private_loh_2nd_LOH_arm) >= 1
      ]
      
      # step 2 exclude p-d pairs called as opposite in dominant mutation method (mean cluster CCF95)
      p.d_comb_mut_ccf95 <-  all_joined_comb_branch.prop$p.d_comb_mut_ccf95[!grepl("NA",all_joined_comb_branch.prop$p.d_comb_mut_ccf95)
                                                                            & !is.na(all_joined_comb_branch.prop$p.d_comb_mut_ccf95)]
      
      p.d_comb_mut_ccf95_reverse <- ifelse(grepl("region.1st",p.d_comb_mut_ccf95), gsub("region.1st", "region.2nd", p.d_comb_mut_ccf95), gsub("region.2nd", "region.1st",p.d_comb_mut_ccf95))
      
      p.d_comb_LOH_cyto_filtered <- setdiff(p.d_comb_LOH_cyto, p.d_comb_mut_ccf95_reverse)  #p.d_comb_mut_clonal_reverse,
      
      
      # add growth pattern
      
      p.d_comb_LOH_cyto_df <- add_growth_pattern(all_joined_comb_branch.prop[all_joined_comb_branch.prop$p.d_comb_LOH_cyto %in%  p.d_comb_LOH_cyto_filtered , ],
                                                 pattern_df) %>%
        mutate(p.d_comb = p.d_comb_LOH_cyto) %>%
        mutate(parent_region =  case_when(grepl("region.1st",p.d_comb) ~ "region.1st",
                                          grepl("region.2nd",p.d_comb) ~ "region.2nd")) %>%
        mutate(parent_grade = case_when(parent_region == "region.1st" ~ grade_1st,
                                        parent_region == "region.2nd" ~ grade_2nd),
               daughter_grade = case_when(parent_region == "region.1st" ~ grade_2nd,
                                          parent_region == "region.2nd" ~ grade_1st)) %>%
        mutate(parent_grade.num = case_when(parent_region == "region.1st" ~ grade.num_1st,
                                            parent_region == "region.2nd" ~ grade.num_2nd),
               daughter_grade.num = case_when(parent_region == "region.1st" ~ grade.num_2nd,
                                              parent_region == "region.2nd" ~ grade.num_1st)) %>%
        mutate(parent_pattern = case_when(parent_region == "region.1st" ~ pattern_1st,
                                          parent_region == "region.2nd" ~ pattern_2nd),
               daughter_pattern = case_when(parent_region == "region.1st" ~ pattern_2nd,
                                            parent_region == "region.2nd" ~ pattern_1st)) %>%
        mutate(grade_shift = case_when(parent_grade.num < daughter_grade.num ~ "up",
                                       parent_grade.num > daughter_grade.num ~ "down",
                                       parent_grade.num == daughter_grade.num ~ "same")) %>%
        mutate(grade_shift = factor(grade_shift, levels =c("up","down","same" ))) %>%
        mutate(parent_grade = factor(parent_grade, levels = c("high","low")),
               daughter_grade = factor(daughter_grade, levels = c("high","low")))
      
      p.d_comb_list <- list(p.d_comb_LOH_cyto_df) 
      
      title_list <- c("LOH cytoband")
      
      for (k in 1:length(p.d_comb_list)){
        
        ## create dot plots and alluvial plots and run permutation test
        
        p.d_comb_df <- p.d_comb_list[[k]]
        
        p.d_comb_freq <- p.d_comb_df %>%
          filter(!is.na(parent_grade.num), !is.na(daughter_grade.num) ) %>%
          group_by(parent_grade, daughter_grade, grade_shift) %>%
          summarise(Freq = n())
        
        p.d_comb_longer <- p.d_comb_df %>%
          filter(!is.na(parent_grade.num), !is.na(daughter_grade.num) ) %>%
          pivot_longer(cols = parent_grade.num:daughter_grade.num, values_to = "grade_score", names_to = "parent_daughter") %>%
          mutate(parent_daughter = gsub("_grade.num", "",parent_daughter )) %>%
          mutate(parent_daughter = factor(parent_daughter, levels = c("parent", "daughter"))) %>%
          mutate(regional_pattern = case_when(parent_daughter == "parent" ~ parent_pattern, 
                                              parent_daughter == "daughter" ~ daughter_pattern))
        
        data_df <- p.d_comb_longer 
        
        num_tum <- length(unique(data_df$tumour_id))
        num_regional_pair <- nrow(data_df)/2
        y_range <- max(data_df$grade_score) - min(data_df$grade_score)
        
        ## permutation test : is upward shift significantly enriched than by chance?  
        
        # observed frequency 
        num_up <- sum(p.d_comb_df$grade_shift == "up", na.rm = T)
        num_down <- sum(p.d_comb_df$grade_shift == "down", na.rm = T)
        num_same <-  sum(p.d_comb_df$grade_shift == "same", na.rm = T)
        
        p.d_comb_tmp <- p.d_comb_df %>%
          dplyr::select(-grade_1st, -grade_2nd, -grade.num_1st, -grade.num_2nd, -pattern_1st, -pattern_2nd,
                        -parent_grade, -daughter_grade, -parent_grade.num,-daughter_grade.num ,-parent_pattern, -daughter_pattern, -grade_shift) 
        
        # run 1000 times -> check empirical pvalue 
        permute_num = 1000
        
        all_grade_shift_test_down <- vector()
        all_grade_shift_test_up <- vector()
        for(i in 1:permute_num){
          print(i)
          
          p.d_comb_random <- add_growth_pattern(p.d_comb_tmp, pattern_df_random[[i]]) %>%  # now new regional patterns are assigned
            mutate(parent_grade = case_when(parent_region == "region.1st" ~ grade_1st,
                                            parent_region == "region.2nd" ~ grade_2nd),
                   daughter_grade = case_when(parent_region == "region.1st" ~ grade_2nd,
                                              parent_region == "region.2nd" ~ grade_1st)) %>%
            mutate(parent_grade.num = case_when(parent_region == "region.1st" ~ grade.num_1st,
                                                parent_region == "region.2nd" ~ grade.num_2nd),
                   daughter_grade.num = case_when(parent_region == "region.1st" ~ grade.num_2nd,
                                                  parent_region == "region.2nd" ~ grade.num_1st)) %>%
            mutate(parent_pattern = case_when(parent_region == "region.1st" ~ pattern_1st,
                                              parent_region == "region.2nd" ~ pattern_2nd),
                   daughter_pattern = case_when(parent_region == "region.1st" ~ pattern_2nd,
                                                parent_region == "region.2nd" ~ pattern_1st)) %>%
            mutate(grade_shift = case_when(parent_grade.num < daughter_grade.num ~ "up",
                                           parent_grade.num > daughter_grade.num ~ "down",
                                           parent_grade.num == daughter_grade.num ~ "same")) %>%
            mutate(grade_shift = factor(grade_shift, levels =c("up","down","same" ))) %>%
            mutate(parent_grade = factor(parent_grade, levels = c("high", "low")),
                   daughter_grade = factor(daughter_grade, levels = c("high", "low")))
          
          #table(p.d_comb_random$grade_shift)
          num_up_rand <- sum(p.d_comb_random$grade_shift == "up", na.rm = T)
          num_down_rand <- sum(p.d_comb_random$grade_shift == "down", na.rm = T)
          num_same_rand <-  sum(p.d_comb_random$grade_shift == "same", na.rm = T)
          
          grade_shift_count <- c(num_up,num_down,num_same) 
          
          # compare with observed count (fixed!)  sum up the number of simulations for which the statistic was greater (less) than or equal to the observed statistic,
          grade_shift_test_down <- ifelse(num_down_rand <= num_down, 1 , 0)
          grade_shift_test_up <- ifelse(num_up_rand >= num_up, 1 , 0)
          
          all_grade_shift_test_down <- c(all_grade_shift_test_down, grade_shift_test_down)
          all_grade_shift_test_up <- c(all_grade_shift_test_up, grade_shift_test_up)
        }
        
        p_emp_down <- signif( (sum(all_grade_shift_test_down) +1)/(permute_num + 1) , 3)
        p_emp_up <- signif( (sum(all_grade_shift_test_up) + 1)/(permute_num + 1) , 3)
        
        title_string <- paste0(title_list[k], "\n",
                               "ratio_filter_p=",ratio_filter_p," ratio_filter_d=",ratio_filter_d ,"\n" 
                               ,"n=",num_tum," tumours (",num_regional_pair," regional_pairs)\n",
                               "P.emp(up)=",p_emp_up,"(up:",num_up," same:",num_same," down:",num_down,")")
        
        
        ## Fig 2b - alluvial plot 
        g1 <- ggplot(data =  p.d_comb_freq,
                     aes(axis1 = parent_grade, axis2 = daughter_grade,
                         y = Freq)) +
          scale_x_discrete(limits = c("Ancestor\n-like", "Descendant\n-like")) + #, expand = c(.2, .05)) +
          labs(x ="Regions", fill = "Transition") +
          geom_alluvium(aes(fill = grade_shift)) +
          scale_fill_manual(values = c("up" = "#F8766D", "same" = "#00BA38", "down" = "#619CFF"), drop = T) +
          geom_stratum() +
          geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 6) +
          theme_cowplot() +
          ggtitle(title_string) +
          theme(plot.title = element_text(size = 14),
                axis.title.x = element_text(size = 16),
                axis.title.y = element_blank(),
                axis.text.x = element_text(size = 16),
                axis.text.y = element_blank(),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_blank(),
                axis.line.x = element_blank(),
                axis.line.y = element_blank(),
                legend.position = "right",
                legend.text = element_text(size = 16),
                legend.title = element_text(size = 18))
        
        print(g1)
        
        
        ## dot plot 
        g2 <-  ggplot(data_df %>% mutate(parent_daughter_jitter = jitter(as.numeric(parent_daughter)), grade_score_jitter = jitter(as.numeric(grade_score), amount = 0) ), 
                      aes(x = parent_daughter_jitter, y = grade_score_jitter))+
          geom_point(color = "black", aes(fill = regional_pattern), alpha = 0.8, size = 8, pch=21)+ 
          geom_line(aes(group = p.d_comb))+
          stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                       geom = "crossbar", width = 0.5, color = "red",
                       aes(x = as.numeric(parent_daughter), group = parent_daughter)) +
          ggtitle(title_string) +
          labs(x = "Regions",  y = "Pattern grade score", fill = "Regional pattern") +
          theme_cowplot()+
          theme(plot.title = element_text(size = 14),
                axis.title.x = element_text(size = 16),
                axis.title.y = element_text(size = 16),
                axis.text.x = element_text(size = 16),
                axis.text.y = element_text(size = 16),
                legend.position = "right",
                legend.text = element_text(size = 16),
                legend.title = element_text(size = 18)) + 
          scale_fill_manual(values = c(  "well_dx" = "#A6CEE3", "papillary" = "#1F78B4", "acinar" = "#B2DF8A",
                                         "cribriform" = "#33A02C","micropapillary" = "#FB9A99", "poor_dx" = "#E31A1C"), drop = T) +
          scale_x_continuous(breaks = c(1,2), 
                             labels=c("Ancestor\n-like", "Descendant\n-like")) 
        
        print(g2)
        
        
      }
    }
  }
}
dev.off()

# Fig 2b - page 9 of the output pdf.



####################################
#==========plot private and public LOH for both LUSC and LUAD

lohscc <- LOH_annotated
loharmscc <- LOH_arm_annotated

LOH_annotated_LUAD <- read_rds("../data/LOH_annotated_LUAD.rds")
LOH_arm_annotated_LUAD <- read_rds("../data/LOH_arm_annotated_LUAD.rds")
lohadc <- LOH_annotated_LUAD
loharmadc <- LOH_arm_annotated_LUAD


#plot ADC LOH arranged
lohadc %>% filter(!(grade_1st == "n.a.") & !(grade_2nd == "n.a.")) %>% 
  mutate(pairwise_id = paste0(tumour_id,"_",compaired_regions)) %>% 
  select(pairwise_id, shared_loh, private_loh_1st, private_loh_2nd) %>% 
  arrange((private_loh_2nd + private_loh_1st)/(private_loh_1st + private_loh_2nd + shared_loh), private_loh_2nd/(private_loh_1st + private_loh_2nd + shared_loh)) -> lohadcuse

lohadcuse$indexnumber = seq.int(nrow(lohadcuse))

lohadc_plot <- melt(lohadcuse, id.vars = c("pairwise_id", "indexnumber")) %>% arrange(indexnumber) 

ggplot(lohadc_plot, aes(fill=variable, y=value, x=indexnumber)) + 
  geom_bar(position="fill", stat="identity") +
  ggtitle("LUAD, LOH_annotated, 814 pairwise regions") +
  scale_fill_manual(values = c(  "shared_loh" = "#A6CEE3", "private_loh_1st" = "#1F78B4", 
                                 "private_loh_2nd" = "#E31A1C"), drop = T) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 13, family = "ArialMT"),
        axis.title.y = element_text(size = 13, family = "ArialMT"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 13, family = "ArialMT"),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank())


#plot ADC arm LOH
loharmadcuse <- loharmadc %>% filter(!(grade_1st == "n.a.") & !(grade_2nd == "n.a.")) %>% 
  mutate(pairwise_id = paste0(tumour_id,"_",compaired_regions)) %>% 
  select(pairwise_id, shared_loh, private_loh_1st, private_loh_2nd) %>% 
  arrange((private_loh_2nd + private_loh_1st)/(private_loh_1st + private_loh_2nd + shared_loh), private_loh_2nd/(private_loh_1st + private_loh_2nd + shared_loh))

loharmadcuse$indexnumber = seq.int(nrow(loharmadcuse))

loharmadc_plot <- melt(loharmadcuse, id.vars = c("pairwise_id", "indexnumber")) %>% arrange(indexnumber) 


ggplot(loharmadc_plot, aes(fill=variable, y=value, x=indexnumber)) + 
  geom_bar(position="fill", stat="identity") +
  ggtitle("LUAD, LOH_chromosome_arm, 814 pairwise regions") +
  scale_fill_manual(values = c(  "shared_loh" = "#A6CEE3", "private_loh_1st" = "#1F78B4", 
                                 "private_loh_2nd" = "#E31A1C"), drop = T) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 13, family = "ArialMT"),
        axis.title.y = element_text(size = 13, family = "ArialMT"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 13, family = "ArialMT"),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank())
#Ext fig 3f - LOH_arm_plot


lohscc %>% filter(!(grade_1st == "n.a.") & !(grade_2nd == "n.a.")) %>% 
  mutate(pairwise_id = paste0(tumour_id,"_",compaired_regions)) %>% 
  select(pairwise_id, shared_loh, private_loh_1st, private_loh_2nd) %>% 
  arrange((private_loh_2nd + private_loh_1st)/(private_loh_1st + private_loh_2nd + shared_loh), private_loh_2nd/(private_loh_1st + private_loh_2nd + shared_loh)) -> lohsccuse

lohsccuse$indexnumber = seq.int(nrow(lohsccuse))

lohscc_plot <- melt(lohsccuse, id.vars = c("pairwise_id", "indexnumber")) %>% arrange(indexnumber) 

ggplot(lohscc_plot, aes(fill=variable, y=value, x=indexnumber)) + 
  geom_bar(position="fill", stat="identity") +
  ggtitle("LUSC, LOH_annotated, 607 pairwise regions") +
  scale_fill_manual(values = c(  "shared_loh" = "#A6CEE3", "private_loh_1st" = "#1F78B4", 
                                 "private_loh_2nd" = "#E31A1C"), drop = T) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 13, family = "ArialMT"),
        axis.title.y = element_text(size = 13, family = "ArialMT"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 13, family = "ArialMT"),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank())


#plot LUSC arm LOH
loharmsccuse <- loharmscc %>% filter(!(grade_1st == "n.a.") & !(grade_2nd == "n.a.")) %>% 
  mutate(pairwise_id = paste0(tumour_id,"_",compaired_regions)) %>% 
  select(pairwise_id, shared_loh, private_loh_1st, private_loh_2nd) %>% 
  arrange((private_loh_2nd + private_loh_1st)/(private_loh_1st + private_loh_2nd + shared_loh), private_loh_2nd/(private_loh_1st + private_loh_2nd + shared_loh))

loharmsccuse$indexnumber = seq.int(nrow(loharmsccuse))

loharmscc_plot <- melt(loharmsccuse, id.vars = c("pairwise_id", "indexnumber")) %>% arrange(indexnumber) 


ggplot(loharmscc_plot, aes(fill=variable, y=value, x=indexnumber)) + 
  geom_bar(position="fill", stat="identity") +
  ggtitle("LUSC, LOH_chromosome_arm, 607 pairwise regions") +
  scale_fill_manual(values = c(  "shared_loh" = "#A6CEE3", "private_loh_1st" = "#1F78B4", 
                                 "private_loh_2nd" = "#E31A1C"), drop = T) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 13, family = "ArialMT"),
        axis.title.y = element_text(size = 13, family = "ArialMT"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 13, family = "ArialMT"),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank())
#Ext fig 3f - LUSC LOH arm plot


##Keratinization plot for ancestor-descendant regions
#Ext fig 3e

purity_data <- all_purity_df
arm_filter_data <- readRDS("../results/Fig.2f_keratinization_ancestraldescendant_assignment/tmp/all_joined_comb_branch.prop_ratio_filter_p0.02_ratio_filter_d0.1.tree.ratio.RDS")
pattern_data <- all_patterns_df

purity_data$pairwise_id <- paste0(purity_data$tumour_id, "_", purity_data$compaired_regions)
arm_filter_data$pairwise_id <- paste0(arm_filter_data$tumour_id, "_", arm_filter_data$compaired_regions)
arm_filter_data %>% select(pairwise_id, shared_LOH_arm) -> 
  arm_filter_data_use

#create table for private LOH proportion - scc
lohscc %>% filter(!(grade_1st == "n.a.") & !(grade_2nd == "n.a.")) %>% 
  mutate(pairwise_id = paste0(tumour_id,"_",compaired_regions)) %>% 
  select(pairwise_id, shared_loh, private_loh_1st, private_loh_2nd) -> lohsccuse

lohsccuse1 <- lohsccuse %>% 
  mutate(private1proportion = private_loh_1st/(private_loh_1st + shared_loh), 
         private2proportion = private_loh_2nd/(private_loh_2nd + shared_loh))

lohsccuse1 %>% mutate(region1st = substr(pairwise_id, 10, 11), 
                      region2nd = substr(pairwise_id, 13, 14)) ->
  lohsccuse2

scc_pairwise_table0 <- left_join(lohsccuse2, arm_filter_data_use, by = c("pairwise_id" = "pairwise_id"))
scc_pairwise_table <- left_join(scc_pairwise_table0, purity_data, by = c("pairwise_id" = "pairwise_id"))

# filtering by min purity > 0.15 and shared_LOH_arm >= 1 and then analyse all regions with A<0.02 and D>0.1
#Criteria
# P1 = min_purity >= 0.15
# S1 = shared_LOH_arm >= 1
# AD1 = private1proportion < 0.02 & private2proportion > 0.1
# AD2 = private1proportion >0.1 & private2proportion < 0.02

scc_pairwise_table %>% mutate(valid = case_when(min_purity >= 0.15 & shared_LOH_arm >= 1 & private1proportion < 0.02 & private2proportion > 0.1 ~ "valid", 
                                                min_purity >= 0.15 & shared_LOH_arm >= 1 & private1proportion > 0.1 & private2proportion < 0.02 ~ "valid", 
                                                min_purity < 0.15 | shared_LOH_arm < 1 ~ "excluded", 
                                                TRUE ~ "cannot_assign")) ->
  scc_pairwise_filtered

scc_pairwise_filtered %>% mutate(AD_region1st = case_when(valid == "excluded" ~ "excluded",
                                                          valid == "cannot_assign" ~ "cannot_assign",
                                                          valid == "valid" & private1proportion > 0.1 ~ "D_like", 
                                                          valid == "valid" & private1proportion < 0.02 ~ "A_like"
)) %>% 
  mutate(AD_region2nd = case_when(valid == "excluded" ~ "excluded",
                                  valid == "cannot_assign" ~ "cannot_assign",
                                  valid == "valid" & private2proportion > 0.1 ~ "D_like", 
                                  valid == "valid" & private2proportion < 0.02 ~ "A_like"
  )) ->
  scc_ad_relationship

scc_ad_relationship %>% select(pairwise_id, comb_id, AD_region1st, AD_region2nd, region1st, region2nd) %>% 
  left_join(., pattern_data[, c("comb_id", "pattern_1", "pattern_2")], by = c("comb_id" = "comb_id")) ->
  scc_ad_pattern

scc_ad_pattern %>% filter(AD_region1st == "A_like" | AD_region1st == "D_like") ->
  scc_ad_selected

#observe A-D relationship between poor/well dx pairs
scc_ad_selected %>% filter(pattern_1 != pattern_2) -> scc_ad_wellpoorpaired
ancestral_descendant_table <- scc_ad_wellpoorpaired

GSVA_data <- read_csv("../df-alvin/gsva_GOBP_keratin_df.csv")

LUSC_region_df <- LUSC_region_df_full  %>%
  filter(!is.na(Histology_per_tumour_id_muttable)) %>% 
  filter(Histology_per_tumour_id_muttable == "Squamous cell carcinoma") %>% 
  filter(included_500_regions == 1)

scc_regionalstat <- LUSC_region_df[, c("tumour_id_muttable_cruk", "sample_name_cruk", "regional_pattern.new", "wGII", "all_TMB", "clonal_TMB", "subclonal_TMB", "FLOH", "purity", "ploidy")]

#match regions with data
ancestral_descendant_table0 <- ancestral_descendant_table %>% 
  mutate(cruk_id = substr(pairwise_id, 1, 8)) %>% 
  mutate(A_like_region = case_when(AD_region1st == "A_like" ~ region1st, 
                                   AD_region2nd == "A_like" ~ region2nd)) %>% 
  mutate(D_like_region = case_when(AD_region1st == "D_like" ~ region1st, 
                                   AD_region2nd == "D_like" ~ region2nd)) %>% 
  mutate(A_like_pattern = case_when(AD_region1st == "A_like" ~ pattern_1, 
                                    AD_region2nd == "A_like" ~ pattern_2)) %>% 
  mutate(D_like_pattern = case_when(AD_region1st == "D_like" ~ pattern_1, 
                                    AD_region2nd == "D_like" ~ pattern_2)) 


#check that the ancestral-descendant regions have no T2 but all are T1 before mutating region name
cruk_id <- ancestral_descendant_table0$cruk_id
lapply(cruk_id, function(x) {scc_regionalstat$sample_name_cruk[grep(x, scc_regionalstat$sample_name_cruk)]})

#mutate region name
ancestral_descendant_table <- ancestral_descendant_table0 %>% 
  mutate(A_like_region = paste0(cruk_id, "_T1-", A_like_region)) %>% 
  mutate(D_like_region = paste0(cruk_id, "_T1-", D_like_region))

table_for_compare<- {
  param = "GOBP_KERATINIZATION"
  left_join(ancestral_descendant_table , select(GSVA_data, sample_name, param), by = c("A_like_region" = "sample_name")) %>% 
    rename(A_like_value = param) %>% 
    left_join(., select(GSVA_data, sample_name, param), by = c("D_like_region" = "sample_name")) %>%
    rename(D_like_value = param) %>% 
    filter(!is.na(A_like_value) & !is.na(D_like_value))
}


#get the parameters
full_table <- function(param){left_join(ancestral_descendant_table , select(GSVA_data, sample_name, param), by = c("A_like_region" = "sample_name")) %>% 
    rename(A_like_value = param) %>% 
    left_join(., select(GSVA_data, sample_name, param), by = c("D_like_region" = "sample_name")) %>%
    rename(D_like_value = param) 
}

get_region_parameter <- function(param){left_join(ancestral_descendant_table , select(GSVA_data, sample_name, param), by = c("A_like_region" = "sample_name")) %>% 
    rename(A_like_value = param) %>% 
    left_join(., select(GSVA_data, sample_name, param), by = c("D_like_region" = "sample_name")) %>%
    rename(D_like_value = param) %>% 
    filter(!is.na(A_like_value) & !is.na(D_like_value) & A_like_pattern == "well_dx")
} # this only captures well dx to poor dx transition

get_region_parameter_reverse <- function(param){left_join(ancestral_descendant_table , select(GSVA_data, sample_name, param), by = c("A_like_region" = "sample_name")) %>% 
    rename(A_like_value = param) %>% 
    left_join(., select(GSVA_data, sample_name, param), by = c("D_like_region" = "sample_name")) %>%
    rename(D_like_value = param) %>% 
    filter(!is.na(A_like_value) & !is.na(D_like_value) & A_like_pattern == "poor_dx")
} # this only captures poor dx to well dx transition

get_region_parameter_all <- function(param){left_join(ancestral_descendant_table , select(GSVA_data, sample_name, param), by = c("A_like_region" = "sample_name")) %>% 
    rename(A_like_value = param) %>% 
    left_join(., select(GSVA_data, sample_name, param), by = c("D_like_region" = "sample_name")) %>%
    rename(D_like_value = param) %>% 
    filter(!is.na(A_like_value) & !is.na(D_like_value))
} # this captures all transition


tab1 <- get_region_parameter("GOBP_KERATINIZATION")
GOBP_KERATINIZATION <- get_region_parameter("GOBP_KERATINOCYTE_DIFFERENTIATION")
GOBP_KERATINIZATION_reverse <- get_region_parameter_reverse("GOBP_KERATINOCYTE_DIFFERENTIATION")

library(ggpubr)
library(dplyr)

# Add a dataset label to each data frame
df1 <- GOBP_KERATINIZATION %>% mutate(dataset = "Kerat")
df2 <- GOBP_KERATINIZATION_reverse %>% mutate(dataset = "Kerat_reverse")
ylabel = "GOBP_Keratinization"
plottitle = "pairwise regions"

# Combine the two data frames
combined_df <- bind_rows(df1, df2)

# Plot with ggpaired, mapping color or facet by dataset
ggpaired(combined_df, cond1 = "A_like_value", cond2 = "D_like_value", id = "pairwise_id",
         palette = c("darkgreen", "red"), line.color = "gray", line.size = 1,
         color = "condition", xlab = "Group", point.size = 3, ylab = ylabel) +
  facet_wrap(~dataset) +  # or use color = dataset if you want overlay instead of facets
  ggtitle(plottitle) +
  theme(legend.direction = 'vertical', legend.position = 'right', legend.key = element_blank(),
        legend.text = element_text(size = 13, family = "ArialMT"),
        axis.text.x = element_text(size = 13, family = "ArialMT"),
        axis.title.x = element_text(size = 13, family = "ArialMT"),
        axis.text.y = element_text(size = 13, family = "ArialMT"),
        axis.title.y = element_text(size = 13, family = "ArialMT"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 13, family = "ArialMT"),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank()) +
  stat_compare_means(method = "t.test", paired = TRUE, label.y = 0.5)
#Ext fig 3e


