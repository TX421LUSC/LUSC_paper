library(tidyverse)
library(data.table)
library(readxl)
library(fst)
library(ggfittext)
library(Rediscover) 
library(qvalue) 
library(RColorBrewer)
library(cowplot)
library(ggpubr)
library(ggsci)
library(ggforce)
library(MASS)
library(car)
library(oddsratio)
library(scales) 


# chr arm list
arm_list <-  c("1p","1q", "2p","2q","3p","3q","4p","4q","5p","5q","6p","6q","7p","7q","8p","8q","9p","9q","10p","10q",
               "11p","11q", "12p","12q","13q","14q","15q","16p","16q","17p","17q","18p","18q","19p","19q","20p","20q",
               "21q","22q")


###############################
### set cutoff for filtering ###
################################

## set cutoff values for events to take into account
## frequency of truncal alteration (percentage of cases)
cutoff_val <- 0.1 
## LUSC subtype specific events ##
cutoff_num <- 10  # keep events with at least 10 tumours in either group


##########################
### path to input data ###
##########################
data_dir <- "../data/"

muttable_path <- paste0(data_dir, "/20221109_TRACERx421_mutation_table.fst")
mphase_sample_path <- paste0(data_dir, "/20221109_TRACERx421_scna_table.fst")
mphase_driver_path <- paste0(data_dir, "/mphase_cyto_gene_drivers_cruk.RDS")
mut_driver_path   <- paste0(data_dir,"/20200401_lung_drivers.csv")
cn_driver_genes_path  <- paste0(data_dir, "/cn_drivers.csv")
tumour_arm_path   <- paste0(data_dir, "/tumour.arm.summary.mphase_cruk.RDS")
all_tumour_path <- paste0(data_dir, "/20221109_TRACERx421_all_tumour_df.rds")
mut_cn_driver_gd_path <- paste0(data_dir, "/20221109_TRACERx421_gd_driver_signature_events.tsv")

#######################
### load input data ###
#######################

### load muttable ###
muttable_df <- read.fst(muttable_path, as.data.table=TRUE)

### load SCNA table ###
mphase_sample_df <- read.fst(mphase_sample_path, as.data.table=TRUE)

mphase_tumours <- mphase_sample_df  %>%
  distinct(tumour_id)

scnatable_df <- readRDS(mphase_driver_path)
scnatable_df <- as.data.table(scnatable_df)

# load arm events 
tumour_arm <-  readRDS(tumour_arm_path)

# tumours with multi-region scna data
multi_region_tumours <- mphase_sample_df %>%
  filter(!grepl('LN', sample)) %>%
  distinct(tumour_id, sample) %>%
  group_by(tumour_id) %>%
  summarise(count = n()) %>%
  filter(count >= 2)

### load drivers ###
cn_drivers_full <- data.table::fread( cn_driver_genes_path )
cn_drivers <- cn_drivers_full[, gene]

mut_drivers <- fread(mut_driver_path)
mut_drivers <- mut_drivers[(is_expressed_tx421),]

### load histology data ### 

clin <-  readRDS(all_tumour_path)

####for multiregion cases only
#clin <- clin[tumour_id_muttable_cruk %in% multi_region_tumours$tumour_id,]

######################################
####select cases with tumour histology
######################################

clin1 <- clin[Histology_per_tumour_id_muttable == "Squamous cell carcinoma",]

clin1[, tumour_id := tumour_id_muttable_cruk]

LUSC_tumour <- unique(clin1$tumour_id)

clin <- clin1
######################################
####separate cases into different cases, for LUSC only
######################################



# load prebuilt object of driver events 
events.ori <- fread(mut_cn_driver_gd_path, data.table = T)

### add chr arm events

arm_events <- tumour_arm #<- arm_summary 

arm_events[, gene := chr_arm]
arm_events[, is_gain := (gain.clonality.threshold == "clonal")]
arm_events[, is_loss.LOH := (loss.clonality.threshold == "clonal" |  loh.clonality.threshold == "clonal")]
arm_events[ (is_gain), `:=`(consequence = 'gain', is_clonal = TRUE)]
arm_events[ (is_loss.LOH), `:=`(consequence = 'loss.LOH', is_clonal = TRUE)]
arm_events <- arm_events[ !is_clonal %in% c('', FALSE) & !tumour_id == '',
                          .( tumour_id, gene, is_clonal,consequence,
                             is_driver = TRUE, type = 'arm') ]

# use prebuilt object 

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


###Alvin: separate into 2 groups of SCC
scc_entire_data <- read.csv("../df-alvin/500_entire_tumour_with_grading.csv")
keratinized_tumours <- scc_entire_data$cruk_tumourid[scc_entire_data$keratinization == "have_keratin"]
nonkeratinized_tumours <- scc_entire_data$cruk_tumourid[scc_entire_data$keratinization == "no_keratin"]


##############################################
###  Count the frequency of truncal events ###
##############################################

events_clonal <- events_all[is_clonal == TRUE, ]

driver_df_count_LUSC <- events_clonal %>% 
  filter(tumour_id %in% keratinized_tumours) %>%
  filter(is_driver == TRUE & type ==  "mutation") %>%
  distinct(tumour_id, gene) %>%
  group_by(gene) %>%
  summarise(count=n()) %>%
  mutate(proportion_altered = count/length(keratinized_tumours)) %>% 
  arrange(desc(count)) %>%
  print() 

#Ext fig 2e
ggplot(head(driver_df_count_LUSC, 12)) +
  geom_col(aes(proportion_altered, reorder(gene, proportion_altered)), fill = "#A6CEE3", width = 0.6) +
  theme(axis.text.x = element_text(size = 13, family = "ArialMT"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 20, family = "ArialMT"),
        axis.title.y = element_text(size = 13, family = "ArialMT"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 13, family = "ArialMT"),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank())



##############################################
###  Count the frequency of truncal events ###
##############################################

events_clonal <- events_all[is_clonal == TRUE, ]

driver_df_count_LUSC <- events_clonal %>% 
  filter(tumour_id %in% nonkeratinized_tumours) %>%
  filter(is_driver == TRUE & type ==  "mutation") %>%
  distinct(tumour_id, gene) %>%
  group_by(gene) %>%
  summarise(count=n()) %>%
  mutate(proportion_altered = count/length(nonkeratinized_tumours)) %>% 
  arrange(desc(count)) %>%
  print() 

#Ext fig 2e
ggplot(head(driver_df_count_LUSC, 12)) +
  geom_col(aes(proportion_altered, reorder(gene, proportion_altered)), fill = "#FB9A99", width = 0.6) +
  theme(axis.text.x = element_text(size = 13, family = "ArialMT"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 20, family = "ArialMT"),
        axis.title.y = element_text(size = 13, family = "ArialMT"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 13, family = "ArialMT"),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank())


#KRAS and EGFR mutated cases -- to be removed in publication
muttable_df %>% filter(Hugo_Symbol == "EGFR" | Hugo_Symbol == "KRAS") %>% dplyr::select(Hugo_Symbol, AAChange, NucleotideChange, func, exonic.func, patient_id) -> tmp1
events_clonal %>% filter(tumour_id %in% scc_entire_data$cruk_tumourid & (gene == "EGFR" | gene == "KRAS")) -> tmp2
tmp1 %>% left_join(tmp2, by = c("patient_id" = "tumour_id", "Hugo_Symbol" = "gene")) %>% filter(!is.na(is_clonal)) %>%  view()


