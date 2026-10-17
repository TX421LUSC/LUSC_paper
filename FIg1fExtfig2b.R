## Script for Tx421 LUSC paper
#Fig 1f
#Ext fig 2b


library(tidyverse)
library(data.table)
library(readxl)
library(fst)
library(ggfittext)
library(Rediscover) 
library(RColorBrewer)
library(qvalue) 
library(cowplot)
library(ggpubr)
library(ggsci)
library(MASS)
library(car)
library(oddsratio)
library(scales) 


subtype_list_short <- c("have_keratin", "no_keratin")
data_dir <- "../data/"
mphase_sample_path <- paste0(data_dir, "/20221109_TRACERx421_scna_table.fst")
tumour_arm_path <- paste0(data_dir, "/20221108_tumour.arm.summary_cruk.RDS")
all_tumour_path <- paste0(data_dir, "/20221109_TRACERx421_all_tumour_df.rds")
mut_cn_driver_gd_path <- paste0(data_dir, "/20221109_TRACERx421_gd_driver_signature_events.tsv")
all_histology_data <- paste0(data_dir, "/full_histology_df.20220818.RDS")
cutoff_val <- 0.05

## LUSC subtype specific events ##
cutoff_num <- 10  # keep events with at least 10 tumours in either group

### load SCNA table ###
mphase_sample_df <- read.fst(mphase_sample_path, as.data.table=TRUE)
multi_region_tumours <- mphase_sample_df %>%
  filter(!grepl('LN', sample)) %>%
  distinct(tumour_id, sample) %>%
  group_by(tumour_id) %>%
  summarise(count = n()) %>%
  filter(count >= 2)

# load arm events 
tumour_arm <-  readRDS(tumour_arm_path)

### load histology data ### 
clin <-  readRDS(all_tumour_path)

####select cases with tumour histology	
clin <- clin[Histology_per_tumour_id_muttable == "Squamous cell carcinoma",]	
clin[, tumour_id := tumour_id_muttable_cruk]	
LUSC_tumour <- unique(clin$tumour_id)   #tumour_id of the selected histology	
clin_grading <- read.csv("../df-alvin/500_entire_tumour_with_grading.csv") #Alvin:%>% filter(differentiation != "N.A.")
clin %>% left_join(., clin_grading, by = c("tumour_id_muttable_cruk" = "cruk_tumourid")) ->
  clin_withgrade
clin0 <- clin

#Alvin: select pathologic feature here
clin <- clin_withgrade %>%
  mutate(subtype_grade = case_when(keratinization %in% c("have_keratin") ~ "kerat_group",
                                   keratinization %in% c("no_keratin") ~ "nonkerat_group"))
#Alvin: select pathologic feature here

LUSC_tumour <- unique(clin$tumour_id[clin$keratinization %in% c("have_keratin", "no_keratin")]) 
LUSC_tumour_cn <- intersect(LUSC_tumour , multi_region_tumours$tumour_id) # 224
LUSC_tumour <- LUSC_tumour_cn 
length(LUSC_tumour)
LUSC_df <- readRDS(all_histology_data) %>%
  mutate(First_WGD = First_GD) %>% 
  left_join(., clin_withgrade[, c("tumour_id_muttable_cruk", "keratinization")], by = "tumour_id_muttable_cruk") %>% 
  dplyr::select(tumour_id_mphase_cruk, keratinization, First_WGD, wGII, clonal_TMB) %>%
  mutate(genome.doubling = factor(First_WGD, levels = c("Clonal","Subclonal", "No GD"))) %>%
  dplyr::select(-First_WGD)

# load prebuilt object of driver events (generated for companion Tx421 paper "Natural history of ...")
events.ori <- fread(mut_cn_driver_gd_path, data.table = T)

###  Combine all events (drivers and passengers)  ###
arm_events <- tumour_arm #<- arm_summary 
arm_events[, gene := chr_arm]
arm_events[, is_gain := (gain.clonality.threshold == "clonal")]
arm_events[, is_loss.LOH := (loss.clonality.threshold == "clonal" |  loh.clonality.threshold == "clonal")]
arm_events[ (is_gain), `:=`(consequence = 'gain', is_clonal = TRUE)]
arm_events[ (is_loss.LOH), `:=`(consequence = 'loss.LOH', is_clonal = TRUE)]
arm_events <- arm_events[ !is_clonal %in% c('', FALSE) & !tumour_id == '',
                          .( tumour_id, gene, is_clonal,consequence,
                             is_driver = TRUE, type = 'arm') ]
mut_events <- events.ori[type == "mutation",
                         .(tumour_id, gene, 
                           is_clonal, 
                           consequence, 
                           is_driver , 
                           type)] %>%
  filter(type == "mutation")
CN_events <- events.ori[consequence %in% c("Amp", "Homdel"),
                        .(tumour_id, gene,
                          is_clonal,
                          consequence,
                          is_driver ,
                          type)]
events_all <- rbind(mut_events, CN_events)
events_all <- rbind(events_all, arm_events) 
events_all <- events_all %>%
  filter(tumour_id %in% LUSC_tumour) %>%
  as.data.table()
length(unique(events_all$tumour_id))

###  Count the frequency of truncal events ###
events_clonal <- events_all[is_clonal == TRUE, ]
driver.scna_df_count_LUSC_amp <- events_clonal %>%  #cn_drivers_tumour
  filter(consequence == "Amp"  & is_driver == TRUE) %>%  
  filter(tumour_id %in% LUSC_tumour) %>%
  distinct(tumour_id, gene) %>%
  group_by(gene) %>%
  summarise(count=n()) %>%
  arrange(desc(count)) %>%
  print()
freq_drivers.scna_amp_LUSC <- driver.scna_df_count_LUSC_amp$gene[driver.scna_df_count_LUSC_amp$count > length(LUSC_tumour)*cutoff_val]
driver_df_count_LUSC <- events_clonal %>% 
  filter(tumour_id %in% LUSC_tumour) %>%
  filter(is_driver == TRUE & type ==  "mutation") %>%
  distinct(tumour_id, gene) %>%
  group_by(gene) %>%
  summarise(count=n()) %>%
  arrange(desc(count)) %>%
  print() 
freq_drivers_LUSC <- driver_df_count_LUSC$gene[driver_df_count_LUSC$count > length(LUSC_tumour)*cutoff_val]
arm_df_all_histology_LUSC_gain <- events_clonal %>%
  filter(tumour_id %in% LUSC_tumour) %>%
  filter(type == "arm" & consequence == "gain") %>%
  group_by(gene) %>% 
  summarise(count = n()) %>%
  arrange(desc(count)) %>%
  print() 
arm_df_all_histology_LUSC_loss.LOH <- events_clonal %>%
  filter(tumour_id %in% LUSC_tumour) %>%
  filter(type == "arm" & consequence == "loss.LOH") %>%
  group_by(gene) %>% 
  summarise(count = n()) %>%
  arrange(desc(count)) %>%
  print() 
freq_gain_LUSC <- as.character(arm_df_all_histology_LUSC_gain$gene[arm_df_all_histology_LUSC_gain$count > length(LUSC_tumour)*cutoff_val])
freq_loss.LOH_LUSC <- as.character(arm_df_all_histology_LUSC_loss.LOH$gene[arm_df_all_histology_LUSC_loss.LOH$count > length(LUSC_tumour)*cutoff_val])

#########Frequency of events in nonkerat_group/kerat group ######	
subtype_grade_list <- c("kerat_group", "nonkerat_group")
freq_drivers_subtype <- list()
freq_gain_subtype <- list()
freq_loss.LOH_subtype <- list()
freq_drivers.scna_amp_subtype <- list()
for(subtype_grade in c("kerat_group", "nonkerat_group")){
  subtype_tumour <- clin$tumour_id[which(clin$subtype_grade == subtype_grade)]
  driver.scna_df_count_subtype_amp <- events_clonal %>%  #cn_drivers_tumour
    filter(consequence == "Amp"  & is_driver == TRUE) %>%  
    filter(tumour_id %in% LUSC_tumour) %>%
    filter(tumour_id %in% subtype_tumour) %>%
    distinct(tumour_id, gene) %>%
    group_by(gene) %>%
    summarise(count=n()) %>%
    arrange(desc(count)) %>%
    print()
  freq_drivers.scna_amp_subtype[[subtype_grade]] <- driver.scna_df_count_subtype_amp$gene[driver.scna_df_count_subtype_amp$count >= cutoff_num]
  driver_df_count_subtype <- events_clonal %>% 
    filter(tumour_id %in% LUSC_tumour) %>%
    filter(tumour_id %in% subtype_tumour) %>%
    filter(is_driver == TRUE & type ==  "mutation") %>%
    distinct(tumour_id, gene) %>%
    group_by(gene) %>%
    summarise(count=n()) %>%
    arrange(desc(count)) %>%
    print() 
  freq_drivers_subtype[[subtype_grade]] <- driver_df_count_subtype$gene[driver_df_count_subtype$count >= cutoff_num]
  arm_df_all_histology_subtype_gain <- events_clonal %>%
    filter(tumour_id %in% LUSC_tumour) %>%
    filter(tumour_id %in% subtype_tumour) %>%
    filter(type == "arm" & consequence == "gain") %>%
    group_by(gene) %>% 
    summarise(count = n()) %>%
    arrange(desc(count)) %>%
    print() 
  arm_df_all_histology_subtype_loss.LOH <- events_clonal %>%
    filter(tumour_id %in% LUSC_tumour) %>%
    filter(tumour_id %in% subtype_tumour) %>%
    filter(type == "arm" & consequence == "loss.LOH") %>%
    group_by(gene) %>% 
    summarise(count = n()) %>%
    arrange(desc(count)) %>%
    print() 
  freq_gain_subtype[[subtype_grade]] <- as.character(arm_df_all_histology_subtype_gain$gene[arm_df_all_histology_subtype_gain$count >= cutoff_num])
  freq_loss.LOH_subtype[[subtype_grade]] <- as.character(arm_df_all_histology_subtype_loss.LOH$gene[arm_df_all_histology_subtype_loss.LOH$count >= cutoff_num])
}
freq_drivers_subtype <- unique(unlist(freq_drivers_subtype))
freq_gain_subtype <-unique(unlist(freq_gain_subtype))
freq_loss.LOH_subtype <- unique(unlist(freq_loss.LOH_subtype))
freq_drivers.scna_amp_subtype <- unique(unlist(freq_drivers.scna_amp_subtype))

## select the events to be tested in Regression analysis and Rediscover ##
tumour_list <- LUSC_tumour # 224
top_drivers_use <- intersect(freq_drivers_LUSC, freq_drivers_subtype)
top_amp_use <- intersect(freq_drivers.scna_amp_LUSC, freq_drivers.scna_amp_subtype)
top_arm.gain_use <- intersect(freq_gain_LUSC, freq_gain_subtype)
top_arm.loss.LOH_use <- intersect(freq_loss.LOH_LUSC, freq_loss.LOH_subtype)
lapply(list(top_drivers_use, top_amp_use,top_arm.gain_use,top_arm.loss.LOH_use), length) 
top_drivers <- freq_drivers_LUSC
top_amp <- freq_drivers.scna_amp_LUSC
top_arm.gain <- freq_gain_LUSC
top_arm.loss.LOH <- freq_loss.LOH_LUSC
lapply(list(top_drivers, top_amp,top_arm.gain,top_arm.loss.LOH), length) 

### organise input data for logistic regression analysis ###
clonal_driver_df_top <- events_clonal %>%
  filter(is_driver == TRUE) %>%
  filter(type == "mutation") %>%
  filter(is_clonal == TRUE) %>%
  filter(gene %in% top_drivers) %>%
  filter(tumour_id %in% LUSC_tumour) 
clonal_amp_df_top <- events_clonal %>%  #cn_drivers_tumour
  filter(is_driver == TRUE) %>%
  filter(consequence == "Amp") %>%
  filter(is_clonal  == TRUE) %>%  
  filter(tumour_id %in% LUSC_tumour) %>%
  distinct(tumour_id, gene) %>%
  dplyr::rename(driver_gene = gene) %>%
  filter(driver_gene %in% top_amp) %>%
  mutate(scna_status = "Clonal_amp")
clonal_gain_df_top <- events_clonal %>%
  filter(consequence == "gain") %>%
  filter(gene %in% top_arm.gain) 
clonal_loss.LOH_df_top <- events_clonal %>%
  filter(consequence == "loss.LOH") %>%
  filter(gene %in% top_arm.loss.LOH) 
intersect(top_arm.gain, top_arm.loss.LOH)
alt_list <- unique(c(top_drivers, paste0(top_amp,"_cn"), paste0(top_arm.gain,"_gain"), paste0(top_arm.loss.LOH,"_loss.LOH"))) #
alt_list_use <- unique(c(top_drivers_use, paste0(top_amp_use,"_cn"), paste0(top_arm.gain_use,"_gain"), paste0(top_arm.loss.LOH_use,"_loss.LOH"))) #
clonal_amp_list <- clonal_amp_df_top %>%
  dplyr::select(tumour_id, driver_gene) %>%
  mutate(tumour_alt = paste0(tumour_id, "_",driver_gene,"_cn"))
clonal_amp_df <- data.frame(tumour_id = rep(tumour_list,length(alt_list)),  
                            alteration = sort(rep(alt_list, length(tumour_list)))) %>%
  mutate(tumour_alt = paste0(tumour_id, "_",alteration)) %>%
  mutate(status = case_when(tumour_alt %in% clonal_amp_list$tumour_alt ~ 1,
                            TRUE ~ 0))  %>%
  dplyr::select(-tumour_alt) %>%
  pivot_wider(names_from = tumour_id, values_from = status) 
AMP_clonal <- as.matrix(clonal_amp_df[,-1])
rownames(AMP_clonal) <- clonal_amp_df$alteration
clonal_mut_list <- clonal_driver_df_top %>%
  dplyr::select(tumour_id, gene) %>%
  mutate(tumour_alt = paste0(tumour_id, "_",gene))
clonal_mut_df <- data.frame(tumour_id = rep(tumour_list,length(alt_list)),  
                            alteration = sort(rep(alt_list, length(tumour_list)))) %>%
  mutate(tumour_alt = paste0(tumour_id, "_",alteration)) %>%
  mutate(status = case_when(tumour_alt %in% clonal_mut_list$tumour_alt ~ 1,
                            TRUE ~ 0))  %>%
  dplyr::select(-tumour_alt) %>%
  pivot_wider(names_from = tumour_id, values_from = status) #%>%
Mut_clonal <- as.matrix(clonal_mut_df[,-1])
rownames(Mut_clonal) <- clonal_mut_df$alteration
clonal_gain_list <- clonal_gain_df_top %>%
  dplyr::select(tumour_id, gene) %>%
  mutate(tumour_alt = paste0(tumour_id, "_",gene, "_gain"))
clonal_gain_df <- data.frame(tumour_id = rep(tumour_list,length(alt_list)),
                             alteration = sort(rep(alt_list, length(tumour_list)))) %>%
  mutate(tumour_alt = paste0(tumour_id, "_",alteration)) %>%
  mutate(status = case_when(tumour_alt %in% clonal_gain_list$tumour_alt ~ 1,
                            TRUE ~ 0))  %>%
  dplyr::select(-tumour_alt) %>%
  pivot_wider(names_from = tumour_id, values_from = status) #%>%
arm.gain_clonal <- as.matrix(clonal_gain_df[,-1])
rownames(arm.gain_clonal) <- clonal_gain_df$alteration
clonal_loss.LOH_list <- clonal_loss.LOH_df_top %>%
  dplyr::select(tumour_id, gene) %>%
  mutate(tumour_alt = paste0(tumour_id, "_",gene, "_loss.LOH"))
clonal_loss.LOH_df <- data.frame(tumour_id = rep(tumour_list,length(alt_list)),
                                 alteration = sort(rep(alt_list, length(tumour_list)))) %>%
  mutate(tumour_alt = paste0(tumour_id, "_",alteration)) %>%
  mutate(status = case_when(tumour_alt %in% clonal_loss.LOH_list$tumour_alt ~ 1,
                            TRUE ~ 0))  %>%
  dplyr::select(-tumour_alt) %>%
  pivot_wider(names_from = tumour_id, values_from = status) #%>%
arm.loss.LOH_clonal <- as.matrix(clonal_loss.LOH_df[,-1])
rownames(arm.loss.LOH_clonal) <- clonal_loss.LOH_df$alteration
anno_df <- data.frame(tumour_id = tumour_list) %>%
  left_join(LUSC_df, 
            by = c("tumour_id" = "tumour_id_mphase_cruk"))  #same as "tumour_id_muttable_cruk"
rownames(anno_df) <- anno_df$tumour_id
all(colnames(Mut_clonal) == anno_df$tumour_id) # [1] TRUE
all(rownames(arm.gain_clonal) == rownames(Mut_clonal)) # [1] TRUE
all(rownames(arm.loss.LOH_clonal) == rownames(arm.gain_clonal)) # [1] TRUE
all(rownames(AMP_clonal) == rownames(arm.loss.LOH_clonal)) # [1] TRUE
add_df <- anno_df %>%
  dplyr::select(tumour_id, genome.doubling, wGII, clonal_TMB) %>%
  mutate(GD_clonal = case_when(genome.doubling == "Clonal" ~ 1,
                               TRUE ~ 0)) %>%
  dplyr::select(-genome.doubling)
add_clonal <- t(matrix(add_df$GD_clonal, ncol = 1))
colnames(add_clonal) <- add_df$tumour_id
rownames(add_clonal) <- c("GD_clonal")
GD_analysed_tumours <- anno_df$tumour_id[!is.na(anno_df$genome.doubling)]
rownames(Mut_clonal)
arm.gain_rows <-  which(grepl("_gain",rownames(Mut_clonal) )) #
arm.loss.LOH_rows <-  which(grepl("_loss.LOH",rownames(Mut_clonal) )) #
scna_rows <- which(grepl("_cn", rownames(Mut_clonal)))
mut_rows <- setdiff(c(1:nrow(Mut_clonal)), c(arm.gain_rows,arm.loss.LOH_rows, scna_rows))
alt_clonal_gain <- arm.gain_clonal[arm.gain_rows ,]
rownames(alt_clonal_gain) <- paste0("gain_",gsub("_gain","",rownames(alt_clonal_gain)))
alt_clonal_loss.LOH <- arm.loss.LOH_clonal[arm.loss.LOH_rows,]
rownames(alt_clonal_loss.LOH) <- paste0("loss.LOH_",gsub("_loss.LOH","",rownames(alt_clonal_loss.LOH)))
alt_clonal_AMP <- AMP_clonal[scna_rows,]
rownames(alt_clonal_AMP) <- paste0(rownames(alt_clonal_AMP), "_AMP")
alt_clonal_Mut <- Mut_clonal[mut_rows, ]
rownames(alt_clonal_Mut) <- paste0(rownames(alt_clonal_Mut),"_Mut")
# step 4 combine and select tumours to be analysed -> make dataframe for regression analysis
alt_clonal <- rbind(alt_clonal_gain,alt_clonal_loss.LOH , alt_clonal_AMP,alt_clonal_Mut, add_clonal)   #alt_clonal_DEL,
alt_clonal <- alt_clonal[rowSums(alt_clonal) >0, ] #colnames(alt_clonal) %in% tumours_keep]
rownames(alt_clonal) <-  gsub("-","\\.",rownames(alt_clonal))
#Alvin: select pathologic feature here
in_df_clonal <- rownames_to_column(as.data.frame(t(alt_clonal)))
in_df_clonal <- in_df_clonal %>%
  dplyr::rename(tumour_id = rowname) %>%
  left_join(anno_df, by = "tumour_id") %>%
  mutate(nonkerat_group_vs_kerat_group = case_when(keratinization %in% c("no_keratin") ~ "nonkerat_group",
                                    keratinization %in% c("have_keratin") ~ "kerat_group")) %>%
  mutate(nonkerat_group_vs_kerat_group = factor( nonkerat_group_vs_kerat_group, levels = c("kerat_group","nonkerat_group")))
## frequency table (per subtype) ##
count_event_clonal <- sapply(subtype_list_short , function(x){
  tumour_list_subtype <- in_df_clonal$tumour_id[which(in_df_clonal$keratinization == x)]  ##
  alt_clonal_subtype <- alt_clonal[,colnames(alt_clonal) %in% tumour_list_subtype]
  rowSums(alt_clonal_subtype)
})
count_subtype <- sapply(subtype_list_short , function(x){   # tumour counts - now using same tumours for mutations and cn
  sum(in_df_clonal$keratinization == x, na.rm = T) })
count_event_clonal_mut <- count_event_clonal[grepl("_Mut",rownames(count_event_clonal) ),]
count_event_clonal_cn <- count_event_clonal[!grepl("_Mut",rownames(count_event_clonal) ),]
percent_event_clonal_mut <-  apply(count_event_clonal_mut, 1, function(x){round(x/count_subtype,2)})
percent_event_clonal_cn <-  apply(count_event_clonal_cn, 1, function(x){round(x/count_subtype,2)})
percent_event_clonal <- cbind(percent_event_clonal_mut,percent_event_clonal_cn)
## create table ##  # only clonal
freq_mut_LUSC_df <- driver_df_count_LUSC %>% # this is only including clonal
  mutate(subtype = "All") %>%
  mutate(category = "driver_mut") %>%
  filter(gene %in% freq_drivers_LUSC) %>%
  mutate(frequency = round(count/length(LUSC_tumour),2)) %>%
  mutate(event_all = paste0(gene, "_Mut")) %>%
  dplyr::select(subtype , event_all, frequency, category)
freq_amp_LUSC_df <- driver.scna_df_count_LUSC_amp %>%
  mutate(subtype = "All") %>%
  mutate(category = "driver_amp") %>%
  filter(gene %in% freq_drivers.scna_amp_LUSC) %>%
  mutate(frequency = round(count/length(LUSC_tumour),2)) %>%
  mutate(event_all = paste0(gene, "_cn_AMP")) %>%
  dplyr::select(subtype , event_all, frequency, category)
freq_arm.gain_LUSC_df <- arm_df_all_histology_LUSC_gain %>%
  mutate(subtype = "All") %>%
  mutate(category = "arm_gain") %>%
  filter(gene %in% freq_gain_LUSC) %>%
  mutate(frequency = round(count/length(LUSC_tumour),2)) %>%
  mutate(event_all = paste0("gain_",gene)) %>%
  dplyr::select(subtype , event_all, frequency, category)
freq_arm.loss.LOH_LUSC_df <- arm_df_all_histology_LUSC_loss.LOH %>%
  mutate(subtype = "All") %>%
  mutate(category = "arm_loss") %>%
  filter(gene %in% freq_loss.LOH_LUSC) %>%
  mutate(frequency = round(count/length(LUSC_tumour),2)) %>%
  mutate(event_all = paste0("loss.LOH_",gene)) %>%
  dplyr::select(subtype , event_all, frequency, category)
freq_GD_df <- LUSC_df %>%
  filter(tumour_id_mphase_cruk %in% LUSC_tumour) %>%
  filter(genome.doubling == "Clonal") %>%
  group_by(genome.doubling) %>%
  summarise(count = n()) %>%
  mutate(subtype = "All") %>%
  mutate(category = "other") %>%
  mutate(frequency = round(count/length(LUSC_tumour),2)) %>%
  mutate(event_all = "GD_clonal") %>%
  dplyr::select(subtype , event_all, frequency, category)
percent_event_all_long <- as.data.frame(percent_event_clonal) %>%
  rownames_to_column() %>%
  dplyr::rename(subtype = rowname) %>%
  pivot_longer(cols = -subtype, names_to = "event_all", values_to = "frequency") %>%
  #  mutate(subtype = factor(gsub("Mucinous","IMA", subtype), levels = c("well_dx", "poor_dx"))) %>%  ####
mutate(category = case_when(grepl("Mut",event_all) ~ "driver_mut",
                            grepl("AMP",event_all) ~ "driver_amp",
                            grepl("gain",event_all) ~ "arm_gain",
                            grepl("loss.LOH",event_all) ~ "arm_loss",
                            TRUE ~ "other"
)) %>%
  mutate(category = factor(category, levels =c("driver_mut", "driver_amp","arm_gain","arm_loss","other"))) %>%
  arrange(desc(category)) %>%
  mutate(event_all = factor(event_all, levels = unique(event_all))) 
percent_event_all_long_full <- rbind(percent_event_all_long, freq_mut_LUSC_df , freq_amp_LUSC_df, freq_arm.gain_LUSC_df , freq_arm.loss.LOH_LUSC_df, freq_GD_df ) %>% #freq_del_LUSC_df,
  mutate(subtype = gsub("All", "Total", subtype)) %>%
  mutate(subtype = factor(subtype, levels = c("well_dx", "poor_dx", "Total"))) %>%
  mutate(event_all = gsub("-",".", event_all))
# save
table(in_df_clonal$keratinization)
## Logistic regression analysis ##
variable_list <- c(paste0(top_drivers_use,"_Mut") , 
                   paste0(top_amp_use,"_cn_AMP"),
                   paste0("gain_",top_arm.gain_use), 
                   paste0("loss.LOH_", top_arm.loss.LOH_use) ) 
variable_list  <- gsub("-","\\.", variable_list)
## nonkerat_group vs kerat_group ##
dep_var <-"nonkerat_group_vs_kerat_group"   
expl_var.ori <- variable_list  
in_df_nona <- as.data.frame(lapply(in_df_clonal, factor))
in_df_nona$dep_var <- in_df_nona[[dep_var]]
in_df_nona <- in_df_nona %>% 
  filter(tumour_id %in% LUSC_tumour) %>%
  filter(keratinization %in% c("have_keratin", "no_keratin")) %>%  
  mutate(keratinization = factor(keratinization, levels = c("have_keratin", "no_keratin")))
table(in_df_nona$nonkerat_group_vs_kerat_group)
remove_var <- vector()
for(var in variable_list) {
  df <- in_df_nona
  df$var <- df[[var]]
  df$dep_var <- df[[dep_var]]
  if(any(table(df$dep_var, df$var) == 0) ) {
    remove_var <- c(remove_var , var)
  }
}
remove_var
expl_var <-  setdiff(expl_var.ori, remove_var) 
# run glm() and stepAIC
f <- as.formula(paste(dep_var,
                      paste(expl_var , collapse = " + "), 
                      sep = " ~ "))
res_glm <- glm(f , data = in_df_nona, family = "binomial") %>% stepAIC()
expl_var_final_H.vs.LM <- names(res_glm$xlevels)
summary(res_glm)
odds_ratio.ori <- or_glm(data = in_df_nona, model = res_glm) %>%
  mutate(predictor =  gsub("q1","q",gsub("p1","p",gsub("AMP1","AMP",gsub("Mut1","Mut",predictor)))))
odds_ratio <- odds_ratio.ori
colnames(odds_ratio) <- gsub("CI_","ci_",  colnames(odds_ratio))
coef_vec <- coef(summary(res_glm))[,1][-1]
coef_pval_vec <- coef(summary(res_glm))[,4][-1]
anova_table2 <- car::Anova(res_glm)  # type2 anova (no need to consider the order of variables)
anova_table2
anova_pval2 <- anova_table2$`Pr(>Chisq)`  
names(anova_pval2) <- rownames(anova_table2)
anova_Chi2 <- anova_table2$`LR Chisq`  
names(anova_Chi2) <- rownames(anova_table2)
res_df_LUSC <- odds_ratio %>%
  mutate(coef_pval = coef_pval_vec,
         Chisq = round(anova_Chi2, 1),
         pval2 = anova_pval2) %>%
  mutate(predictor = gsub("_clonal1","", predictor)) %>%
  mutate(col_coef = case_when( (grepl("ain", predictor) | grepl("AMP", predictor) | grepl("GD", predictor) ) & pval2 < 0.05  ~ "gain",
                               (grepl("oss", predictor) ) & pval2 < 0.05 ~ "loss",
                               grepl("Mut", predictor) & pval2 < 0.05 ~ "mut"),
         sig = case_when(coef_pval < 0.001 ~ "***",
                         coef_pval < 0.01 ~ "**",
                         coef_pval < 0.05 ~ "*"),
         sig2 = case_when(pval2 < 0.001 ~ "***",
                          pval2 < 0.01 ~ "**",
                          pval2 < 0.05 ~ "*")
  ) %>%
  mutate(predictor = gsub("_"," ",predictor))
res_df_LUSC_sig <-  res_df_LUSC %>%
  filter(pval2 < 0.05)

# Fig 1f
g1 <- ggplot(data = res_df_LUSC, aes(x = reorder(predictor, -oddsratio), y = oddsratio, ymin=`ci_low (2.5)`, ymax=`ci_high (97.5)`)) +
  geom_pointrange(aes(colour= col_coef) ) + 
  geom_hline(yintercept=1, lty=2) +
  theme_cowplot() +
  labs( y = "Odds ratio (95%CI)" , color = "Alteration\n(P<0.05)") + 
  scale_color_manual(values = c("gain" = "red", "loss" = "blue", "mut" = "darkorange"), na.value = "#bababa",  drop = F ) +
  geom_text(aes(label=sig2), size = 10) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  ggtitle("  kerat_group      nonkerat_group") + 
  coord_flip() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    axis.title.y=element_blank(),
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size=18),
    plot.title = element_text(size=18)
  ) 
print(g1)

### model adjusted for wGII ###
dep_var <-"nonkerat_group_vs_kerat_group"   
expl_var <-  expl_var_final_H.vs.LM
title_label <-  "  kerat_group      nonkerat_group" 
# run glm() and stepAIC
f <- as.formula(paste(dep_var,
                      paste(c(expl_var, "wGII") , collapse = " + "),   
                      sep = " ~ "))
in_df_nona1 <- in_df_nona %>%
  dplyr::select(-wGII, -clonal_TMB) %>%
  left_join(add_df, by = "tumour_id")
res_glm <- glm(f , data = in_df_nona1, family = "binomial") 
summary(res_glm)
odds_ratio.ori <- or_glm(data = in_df_nona1, model = res_glm,  incr = list(wGII = 0.1)) %>%
  mutate(predictor =  gsub("q1","q",gsub("p1","p",gsub("AMP1","AMP",gsub("Mut1","Mut",predictor)))))
odds_ratio <- odds_ratio.ori
colnames(odds_ratio) <- gsub("CI_","ci_",  colnames(odds_ratio))
coef_vec <- coef(summary(res_glm))[,1][-1]
coef_pval_vec <- coef(summary(res_glm))[,4][-1]
anova_table2 <- car::Anova(res_glm)  # type2 anova (no need to consider the order of variables)
anova_table2
anova_pval2 <- anova_table2$`Pr(>Chisq)`  
names(anova_pval2) <- rownames(anova_table2)
anova_Chi2 <- anova_table2$`LR Chisq`  
names(anova_Chi2) <- rownames(anova_table2)
res_df_LUSC <- odds_ratio %>%
  mutate(coef_pval = coef_pval_vec,
         Chisq = round(anova_Chi2, 1),
         pval2 = anova_pval2) %>%
  mutate(predictor = gsub("_clonal1","", predictor)) %>%
  mutate(col_coef = case_when( (grepl("ain", predictor) | grepl("AMP", predictor) | grepl("GD", predictor) ) & pval2 < 0.05  ~ "gain",
                               (grepl("oss", predictor) ) & pval2 < 0.05 ~ "loss",
                               grepl("Mut", predictor) & pval2 < 0.05 ~ "mut"),
         sig = case_when(coef_pval < 0.001 ~ "***",
                         coef_pval < 0.01 ~ "**",
                         coef_pval < 0.05 ~ "*"),
         sig2 = case_when(pval2 < 0.001 ~ "***",
                          pval2 < 0.01 ~ "**",
                          pval2 < 0.05 ~ "*")
  ) %>%
  mutate(predictor = gsub("_"," ",predictor))
res_df_LUSC_sig <-  res_df_LUSC %>%
  filter(pval2 < 0.05)

# Ext fig 2b
g1 <- ggplot(data = res_df_LUSC, aes(x = reorder(predictor, -oddsratio), y = oddsratio, ymin=`ci_low (2.5)`, ymax=`ci_high (97.5)`)) +
  geom_pointrange(aes(colour= col_coef) ) + 
  geom_hline(yintercept=1, lty=2) +
  theme_cowplot() +
  labs( y = "Odds ratio (95%CI)" , color = "Alteration\n(P<0.05)") + 
  scale_color_manual(values = c("gain" = "red", "loss" = "blue", "mut" = "darkorange","other" = "darkgreen"), na.value = "#bababa",  drop = F ) +
  geom_text(aes(label=sig2), size = 10) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  ggtitle(title_label) + 
  coord_flip() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    axis.title.y=element_blank(),
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 18), 
    axis.text.y = element_text(size=18),
    plot.title = element_text(size=18)
  ) 
print(g1)

### model adjusted for clonal TMB ###
dep_var <- "nonkerat_group_vs_kerat_group"   
expl_var <- expl_var_final_H.vs.LM
title_label <- "  kerat_group      nonkerat_group" 
f <- as.formula(paste(dep_var,
                      paste(c(expl_var, "clonal_TMB"), collapse = " + "),
                      sep = " ~ "))
in_df_nona1 <- in_df_nona %>%
  dplyr::select(-wGII, -clonal_TMB) %>%
  left_join(add_df, by = "tumour_id")
res_glm <- glm(f, data = in_df_nona1, family = "binomial")
summary(res_glm)
odds_ratio.ori <- or_glm(data = in_df_nona1, model = res_glm, incr = list(clonal_TMB = 100)) %>%
  mutate(predictor = gsub("q1", "q", gsub("p1", "p", gsub("AMP1", "AMP", gsub("Mut1", "Mut", predictor)))))
odds_ratio <- odds_ratio.ori
colnames(odds_ratio) <- gsub("CI_", "ci_", colnames(odds_ratio))
coef_vec <- coef(summary(res_glm))[, 1][-1]
coef_pval_vec <- coef(summary(res_glm))[, 4][-1]
anova_table2 <- car::Anova(res_glm)
anova_pval2 <- anova_table2$`Pr(>Chisq)`
names(anova_pval2) <- rownames(anova_table2)
anova_Chi2 <- anova_table2$`LR Chisq`
names(anova_Chi2) <- rownames(anova_table2)
res_df_LUSC <- odds_ratio %>%
  mutate(coef_pval = coef_pval_vec,
         Chisq = round(anova_Chi2, 1),
         pval2 = anova_pval2) %>%
  mutate(predictor = gsub("_clonal1", "", predictor)) %>%
  mutate(col_coef = case_when(
    (grepl("ain", predictor) | grepl("AMP", predictor) | grepl("GD", predictor)) & pval2 < 0.05 ~ "gain",
    grepl("oss", predictor) & pval2 < 0.05 ~ "loss",
    grepl("Mut", predictor) & pval2 < 0.05 ~ "mut"
  ),
  sig = case_when(
    coef_pval < 0.001 ~ "***",
    coef_pval < 0.01 ~ "**",
    coef_pval < 0.05 ~ "*"
  ),
  sig2 = case_when(
    pval2 < 0.001 ~ "***",
    pval2 < 0.01 ~ "**",
    pval2 < 0.05 ~ "*"
  )) %>%
  mutate(predictor = gsub("_", " ", predictor))
res_df_LUSC_sig <- res_df_LUSC %>%
  filter(pval2 < 0.05)
# Ext fig 2b
g1 <- ggplot(data = res_df_LUSC, aes(x = reorder(predictor, -oddsratio), y = oddsratio, ymin = `ci_low (2.5)`, ymax = `ci_high (97.5)`)) +
  geom_pointrange(aes(colour = col_coef)) + 
  geom_hline(yintercept = 1, lty = 2) +
  theme_cowplot() +
  labs(y = "Odds ratio (95%CI)", color = "Alteration\n(P<0.05)") + 
  scale_color_manual(values = c("gain" = "red", "loss" = "blue", "mut" = "darkorange", "other" = "darkgreen"), na.value = "#bababa", drop = FALSE) +
  geom_text(aes(label = sig2), size = 10) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  ggtitle(title_label) + 
  coord_flip() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 18), 
    axis.text.y = element_text(size = 18),
    plot.title = element_text(size = 18)
  )
print(g1)

