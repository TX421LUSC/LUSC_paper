##Script for Tx421 LUSC paper
##Fig 3c
##Ext fig 4d

# Description:
#  GSEA of regional level - budding 

#################################
### libraries and directories ###
#################################

library(tidyverse)
library(data.table)
library(cowplot)
library(ggfittext)
library(nlme)
library(GSEABase)
library(fgsea)
library(fst)

# setwd("") 
workdir <- "../results/new_results/"
dir.create(workdir)

tmpdir <- paste0(workdir,"/tmp/")
dir.create(tmpdir)

## load LUSC regional data 
input_plot_df_all <-  read_csv("../df-alvin/All_regions_500_regions-SCC_reviewed_With_max_budding_data.csv")

# histology data, on a regional  level
LUSC_input_plot_df_all <- input_plot_df_all %>%
  dplyr::select(sample_name_cruk, tumour_budding) %>%
  mutate(High.vs.LowMid =  case_when(tumour_budding %in% c("no_budding") ~ "LowMid",
                                     tumour_budding %in% c("have_budding") ~ "High")) %>%
  mutate(tumour_id_muttable_cruk = sample_name_cruk) %>% 
  mutate(High.vs.LowMid = factor(High.vs.LowMid , levels = c("LowMid","High"))) %>% 
  filter(!is.na(tumour_budding))

## load ploidy adjusted copy number 
Tx421_weigted_gene_mean.adj.cpn_df.ori <- readRDS("../data/20221108_Tx421_weigted_gene_mean.adj.cpn_df_cruk.RDS")

Tx421_weigted_gene_mean.adj.cpn_df.tmp <- as.data.frame(t(Tx421_weigted_gene_mean.adj.cpn_df.ori[,-1]))
colnames(Tx421_weigted_gene_mean.adj.cpn_df.tmp) <- Tx421_weigted_gene_mean.adj.cpn_df.ori$gene_name
Tx421_weigted_gene_mean.adj.cpn_df  <- cbind(data.frame(sample_name_cruk = colnames(Tx421_weigted_gene_mean.adj.cpn_df.ori[-1])),
                                             Tx421_weigted_gene_mean.adj.cpn_df.tmp)

# load median vst per tumour 
# clean the vsd file column names
idlist <- read.csv("../df-alvin/id_conversion.csv") 
vsdfile <- readRDS("../data/vsd_median_alltumour_df.RDS") 
vsdregionfile <- read_fst("../data/2022-10-18_vst_counts.fst")
vsdregionfile$gene <- vsdregionfile$gene_id 
vsdregionfile <-  vsdregionfile %>% dplyr::select(-gene_id) %>% 
                                    dplyr::select(gene, !(gene))  #move gene column to first column

vsdregionfileold <- vsdregionfile
colnames(vsdregionfile) -> tumourlists
tumourlists[2:1052] -> tumourlists1

#error caused by there is two columns with CRUK0495 in vsd_median. 
#both CRUK0495_SU_T1 and CRUK0495_SU_T2 converted to CRUK0495
#however, on checking, the tumour should be CRUK0495_T1 is squamous carcinoma, T2 is ?normal according to clinical data table.
#remvoe CRUK0495_T2 from raw vsdfile

# load Hallmark genesets 
pathway_file <- "../data/h.all.v7.2.symbols.gmt"
pathways_hm <- gmtPathways(pathway_file)

# load Hallmark geneset category
hallmark_process <- fread("../data/hallmark_category_2015CellSystems_Table1.csv")


##################################################
### Fig 1g: gsea - High vs LowMid grade tumour ###
##################################################
## run fgseaLabel() : more stringent than default fgsea()

vsd_median <- vsdregionfile
mat <- vsd_median[,colnames(vsd_median) %in% LUSC_input_plot_df_all$tumour_id_muttable_cruk[!is.na(LUSC_input_plot_df_all$High.vs.LowMid)]]
rownames(mat) <- vsd_median$gene


label_status <- sapply(colnames(mat), function(x){
  status <- LUSC_input_plot_df_all$High.vs.LowMid[LUSC_input_plot_df_all$tumour_id_muttable_cruk == x]
  return(status)
})
labels <- as.numeric(label_status)


# run fgsea 
fgseaRes_hm_High.vs.LowMid <- fgseaLabel(pathways_hm, mat, labels, nperm = 1000, minSize = 15, maxSize = 500)
#fgseaRes_hm_High.vs.LowMid$p.adj <- p.adjust(fgseaRes_hm_High.vs.LowMid$pval, method = "BH")

## plot 

fgseaResTidy <- fgseaRes_hm_High.vs.LowMid %>%
  as_tibble() %>%
  arrange(desc(NES))

in_df <- fgseaResTidy %>%
  mutate(pathway = gsub("HALLMARK_","",pathway)) %>%
  left_join(hallmark_process, by = c("pathway" = "Hallmark Name") ) %>%
  dplyr::rename(Category = `Process Category`) %>%
  mutate(Category = factor(Category , levels = rev(c("proliferation","DNA damage","signaling" ,
                                                 "metabolic" ,"pathway")))) %>%
  mutate(ID = gsub("_"," ",pathway)) %>%
  mutate(sig = ifelse(padj < 0.25, "padj < 0.05","padj > 0.05")) %>%
  mutate(idsig = ifelse(padj < 0.25, gsub("_"," ",ID),"")) %>%
  mutate(direction = ifelse(NES>0, "pos","neg")) %>%
  filter(padj < 0.25)

# color (fill) by category

g <- ggplot(in_df,  aes(x = reorder(ID, NES), y = NES, fill = Category)) +
  geom_bar(stat="identity", #fill = "gray40",  
           colour="black", pch=21) + 
  labs(x = "Geneset", y = "NES")  +
  theme_cowplot() +
  coord_flip() +
  scale_fill_brewer(palette = "YlGnBu",
                    guide = guide_legend(reverse = TRUE)) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    axis.title.y=element_blank(),
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size=14),
    plot.title = element_text(size=22) )

print(g)


