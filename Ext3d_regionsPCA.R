#Script for Tx421 LUSC paper
#Ext fig 3d
# PCA for pathologic features of tumour regions graded by entire tumour in LUSC

library(dplyr)
library(ggplot2)
library(reshape2)
library(ggbiplot)
library(tidyverse)
library(ggpubr)
library(cowplot)
#Script to obtain PCs from expression data (Surgery tumour samples only) and correlate with clinical/genomic covariates

# LOAD REGIONAL MEDIAN VST DATA
# FOR EVERY REGIONS
vst_su <- read.csv("../data/vst_transf_counts_surgery.csv",
                   stringsAsFactors = FALSE)

#Clean
rownames(vst_su) <- vst_su$X
vst_su$X <- NULL
colnames(vst_su) <- gsub(pattern = "\\.",
                         replacement = "-",
                         x = colnames(vst_su))

colnames(vst_su) <- gsub(pattern = "--.+",
                         replacement = "",
                         x = colnames(vst_su))

#Get all metadata---------------------------------------------
all_metadata <- read.csv("../data/all_metadata_RNAseq.csv", stringsAsFactors = FALSE)

all_metadata$region <- gsub(pattern = "LN00",
                            replacement = "LN0",
                            x = all_metadata$region)

#Get tumour ID and histology per region
tumour_id_region <- unique(all_metadata[, c("region", "tumour_id", "Histology")])

#only primary tumour regions
tumour_id_region <- tumour_id_region[grepl("SU_T", tumour_id_region$region), ]


# LOAD ENTIRE TUMOUR MEDIAN VST DATA
# load median vst per tumour 
# clean the vsd file column names
idlist <- read.csv("../df-alvin/id_conversion.csv") 
vsdfile <- readRDS("../data/vsd_median_alltumour_df.RDS") 

#this file vsd_median_alltumour_df.RDS was obtained using the exact same method which 
#resulted in the 20221108_vsd_median_df_luad.RDS df, already deposited in Zenodo. In fact,
#20221108_vsd_median_df_luad.RDS is a subset of this file, and this file contains all tumours LUAD+LUSC. 


vsdfileold <- vsdfile
vsdfile <- vsdfile %>% dplyr::select(-CRUK0495_SU_T2)
colnames(vsdfile) -> tumourlists
tumourlists[2:349] -> tumourlists1
sapply(tumourlists1, function(x) {unique(idlist$tumour_id_muttable_cruk[idlist$cruk_su_tumourid == x])[1]}) %>% unname -> tumourlist2
colnames(vsdfile)[2:349] <- tumourlist2  #WARNING:cases with _tumour1 and _tumour2 CRUK id may catch problems in this conversion
vsd_median <- vsdfile  
rm(idlist, vsdfile, tumourlists, vsdfileold)

squamous_entire_data <- read_csv("../df-alvin/500_entire_tumour_with_grading_with_wGIITMB_info.csv")
# join tumour_id_region df with pathologic features
tumour_id_region0 <- tumour_id_region
tumour_id_region$region_trim <- paste0(substr(tumour_id_region$region, 1, 3), str_sub(tumour_id_region$region, 5, -1))
tumour_id_region$ltx_tumour_id <- str_sub(gsub("_SU", "", tumour_id_region$region_trim), 1, -4)
tumour_id_region <- left_join(tumour_id_region, dplyr::select( squamous_entire_data, ltx_region_id, tumour_id_by.label, differentiation, keratinization, tumour_budding, single_cell_invasion, atypia), by = c("ltx_tumour_id" = "tumour_id_by.label"))

############

# Define parameters for PCA
param <- "keratinization"
param_poor_group <- "no_keratin"
param_well_group <- "have_keratin"

# Assign parameter to dataframe
tumour_id_region$param <- tumour_id_region[[param]]

# Get regions for each group
poor_regions <- tumour_id_region$region[tumour_id_region$param == param_poor_group & !is.na(tumour_id_region$param)]
well_regions <- tumour_id_region$region[tumour_id_region$param == param_well_group & !is.na(tumour_id_region$param)]

# Combine regions for analysis
all_regions <- c(well_regions, poor_regions)
region_list <- all_regions
vst_mat <- vst_su

# Keep only regions of interest and prepare data for PCA
to_pca_vst <- vst_mat[, colnames(vst_mat) %in% region_list]
to_pca_vst <- as.data.frame(t(to_pca_vst))
to_pca_vst$casenumber <- rownames(to_pca_vst)

# Annotate each case with its group
to_pca_vst1 <- to_pca_vst %>% mutate(
  grouping = case_when(
    casenumber %in% well_regions ~ param_well_group,
    casenumber %in% poor_regions ~ param_poor_group,
    TRUE ~ "allothers"
  )
)

# Perform PCA
pcs_region <- prcomp(dplyr::select(to_pca_vst1, -c(grouping, casenumber)))

# Print the distribution of groups
print(table(to_pca_vst1$grouping))

# Define color schemes
color_scheme <- c("#A6CEE3", "#FB9A99")
# Plot PCA with original color scheme
pcaplot1 <- ggbiplot(pcs_region, obs.scale = 1, var.scale = 1, choices = c(1, 2),
                     groups = to_pca_vst1$grouping, ellipse = FALSE, circle = TRUE,
                     var.axes = FALSE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top') +
  scale_color_manual(name = param, values = color_scheme) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 13, colour = "black"),
    axis.text.x = element_text(size = 13, colour = "black"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13),
    axis.line = element_line(colour = "black"),
    panel.border = element_blank(),
    panel.background = element_blank()
  )

# Check match between PCA row order and data
print(identical(rownames(pcs_region$x), to_pca_vst1$casenumber))  # Should be TRUE

# Show top genes with highest loadings on PC1
print(head(sort(abs(pcs_region$rotation[, "PC1"]), decreasing = TRUE)))

# Assign the plot to variable
p2 <- pcaplot1
print(p2)



##check the purity of the regions
input_plot_df_all <-  read_csv("../df-alvin/All_regions_500_regions-SCC_reviewed_With_max_budding_data.csv")

# histology data, on a regional  level
LUSC_input_plot_df_all <- input_plot_df_all %>%
  dplyr::select(ltx_region_id, tumour_budding, keratinization, purity) %>%
  mutate(tumour_id_muttable_cruk = ltx_region_id) %>% 
  filter(!is.na(keratinization))
purity_df <- tumour_id_region %>% left_join(LUSC_input_plot_df_all[, c("ltx_region_id", "purity")], by = c("region_trim" = "ltx_region_id"))


poor_regions <- tumour_id_region$region[tumour_id_region$keratinization == "no_keratin" & !is.na(tumour_id_region$keratinization)]
well_regions <- tumour_id_region$region[tumour_id_region$keratinization == "have_keratin" & !is.na(tumour_id_region$keratinization)]
all_regions <- c(well_regions, poor_regions)
region_list <- all_regions
vst_mat <- vst_su

#Keep only regions of interest and calculate PCs
to_pca_vst <- vst_mat[, colnames(vst_mat) %in% region_list]
to_pca_vst <- as.data.frame(t(to_pca_vst))
to_pca_vst$casenumber <- rownames(to_pca_vst)
to_pca_vst1 <- to_pca_vst %>% mutate(grouping = case_when(casenumber %in% well_regions ~ "have_keratin", 
                                                          casenumber %in% poor_regions ~ "no_keratin",
                                                          TRUE ~ "allothers"))
pcs_region <- prcomp(dplyr::select(to_pca_vst1, -c(grouping, casenumber)))
print(table(to_pca_vst1$grouping))
rownames(to_pca_vst1[to_pca_vst1$grouping == "have_keratin", ]) -> have_keratin_cases
rownames(to_pca_vst1[to_pca_vst1$grouping == "no_keratin", ]) -> no_keratin_cases

purity_df1 <- purity_df %>% mutate(includedInPCA = case_when(region %in% have_keratin_cases ~ "PCA_have_keratin", 
                                                             region %in% no_keratin_cases ~ "PCA_no_keratin"
)) %>% 
  mutate(have_purity_data = ifelse(!is.na(purity), "have_purity_data", "no_purity_data")) %>% 
  filter(!is.na(includedInPCA))

color_list <-  c("PCA_have_keratin" = "#A6CEE3",
                 "homo_have_keratin" = "#B2DF8A",
                 "PCA_no_keratin"= "#FB9A99")

g1 <- ggboxplot(purity_df1, width=0.5,lwd=1, x = "includedInPCA", y = "purity", fill = "includedInPCA", add = "jitter") + 
  stat_compare_means(comparisons = list( c("PCA_have_keratin", "PCA_no_keratin")), label.y = c(0.8)) +
  theme_cowplot() +
  scale_fill_manual(values = color_list, na.value = '#bababa', drop = F) +
  #  scale_x_discrete(labels = c("Well-differentiated", "Poorly-differentiated")) +
  theme(  axis.title.y = element_text(size=26),
          axis.text.y = element_text(size=26),
          axis.title.x = element_text(size=26),
          axis.text.x = element_text(size=26),
          legend.position = "none",
          legend.title = element_blank()) 
print(paste0("n number of have_keratin:", nrow(purity_df1[purity_df1$includedInPCA == "PCA_have_keratin",]) ))
print(paste0("n number of no_keratin:", nrow(purity_df1[purity_df1$includedInPCA == "PCA_no_keratin",]) ))
print(table(purity_df1[, c("includedInPCA", "have_purity_data")]))
print(g1)

head(sort(abs(pcs_region$rotation[, "PC1"]), decreasing = TRUE), n = 100) %>% view()
head(sort(abs(pcs_region$rotation[, "PC2"]), decreasing = TRUE), n = 100) %>% view()




