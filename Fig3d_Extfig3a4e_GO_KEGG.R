##Script for Tx421 LUSC paper
#Fig 3d
#Ext fig 3a, 4e

library(tidyverse)
library(data.table)
library(cowplot)
library(dplyr)
library(GSEABase)
library(fgsea)
library(fst)
library(DESeq2)
library(tidyverse)
library(tinyarray)
library(ggrepel)
library(FactoMineR)
library(pheatmap)
library(affy)
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(enrichplot)
library(ggupset)
library(ggrepel)
library(GOplot)


## convert the vst_count file to median 
rnaseq_count <- read_fst("../data/2022-10-17_rsem_counts_mat.fst")

tumourcolumns <- colnames(rnaseq_count)[grepl("SU_T", colnames(rnaseq_count))]
crukid <- unique(substr(tumourcolumns, 1, 14))

# collapse to median per tumour
vst_tumour <- rnaseq_count[, c(tumourcolumns, "gene_id")]
gene_list <- vst_tumour$gene_id
vst_tumour <- dplyr::select(vst_tumour, -gene_id)

vsd_median_alltumour <- sapply(crukid, function(x){
  samples <- colnames(vst_tumour)[grep(x, colnames(vst_tumour))]
  if (length(samples) == 1) {
    return(vst_tumour[, samples])
  } else {
    return(matrixStats::rowMedians(as.matrix(vst_tumour[, samples])))
  }
})
rownames(vsd_median_alltumour) <- gene_list

vsd_median_alltumour_df <- as.data.frame(vsd_median_alltumour) %>%
  rownames_to_column("gene")

# this df vsd_median_alltumour_df is different from the file vsd_median_alltumour_df.RDS -
# this contains raw read counts, while vsd_median_alltumour_df.RDS was normalized.

# histology data, on an entire tumour level
input_plot_df_all <- read_csv("../df-alvin/500_entire_tumour_with_grading_with_wGIITMB_info.csv")
LUSC_input_plot_df_all <- input_plot_df_all %>%
  dplyr::select(cruk_tumourid, keratinization) %>%
  mutate(High.vs.LowMid = case_when(
    keratinization %in% c("no_keratin") ~ "LowMid",
    keratinization %in% c("have_keratin") ~ "High"
  )) %>%
  mutate(tumour_id_muttable_cruk = cruk_tumourid) %>% 
  mutate(High.vs.LowMid = factor(High.vs.LowMid , levels = c("LowMid","High"))) %>% 
  filter(!is.na(keratinization))

# tumour ID conversion 
idlist <- read.csv("../df-alvin/id_conversion.csv") 
vsd_median <- vsd_median_alltumour_df %>% dplyr::select(-CRUK0495_SU_T2)

tumourlists <- colnames(vsd_median)
tumourlists1 <- tumourlists[2:length(tumourlists)]   # exclude "gene"
tumourlist2 <- sapply(tumourlists1, function(x) {
  unique(idlist$tumour_id_muttable_cruk[idlist$cruk_su_tumourid == x])[1]
}) %>% unname()

colnames(vsd_median)[2:length(tumourlists)] <- tumourlist2  

# round counts (base R, original script style)
rownames(vsd_median) <- vsd_median$gene
vsd_median <- vsd_median[-1]
vsd_median <- round(vsd_median)

# separate into two groups of cancer - high grade vs low grade
low.group <- LUSC_input_plot_df_all$cruk_tumourid[LUSC_input_plot_df_all$High.vs.LowMid == "LowMid"]
high.group <- LUSC_input_plot_df_all$cruk_tumourid[LUSC_input_plot_df_all$High.vs.LowMid == "High"]

exp4 <- cbind(
  vsd_median[, colnames(vsd_median) %in% high.group],
  vsd_median[, colnames(vsd_median) %in% low.group]
)

group_list = factor(c(rep("high", ncol(exp4[, colnames(exp4) %in% high.group])),
                      rep("low", ncol(exp4[, colnames(exp4) %in% low.group]))),
                    levels = c("high","low"))

################## DESeq2 for differential expressed genes ##############################
countData <- exp4
condition <- group_list   

dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), ~ condition)
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", levels(condition)))
RES <- as.data.frame(res)

# THE DEG matrix
RES$group <- case_when(
  RES$log2FoldChange > 1 & RES$pvalue < 0.05 ~ "Up",
  RES$log2FoldChange < -1 & RES$pvalue < 0.05 ~ "Down",
  abs(RES$log2FoldChange) <= 1 ~ "None",
  RES$pvalue >= 0.05 ~ "None"
)

RES.DEG  <- RES %>% filter(group %in% c("Up","Down"))
RES.up   <- RES %>% filter(group == "Up")
RES.down <- RES %>% filter(group == "Down")

# Choosing input for analysis: 
info <- RES.up

GO_database <- 'org.Hs.eg.db'
KEGG_database <- 'hsa'
gene <- bitr(rownames(info),fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)

GO <- enrichGO(gene$ENTREZID,
               OrgDb = GO_database,
               keyType = "ENTREZID",
               ont = "ALL",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               readable = TRUE)
GO@result[["FoldEnrichment"]] <- parse_ratio(GO@result[["GeneRatio"]]) / parse_ratio(GO@result[["BgRatio"]])

# Enriched GO terms (matrix)
go <- GO@result

# Figure
GO_top20 <- dotplot(GO, x="FoldEnrichment" ,showCategory = 20, orderBy = "FoldEnrichment")

KEGG2 <- enrichKEGG(gene$ENTREZID,
                    organism = KEGG_database,
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)
KEGG2@result[["FoldEnrichment"]] <- parse_ratio(KEGG2@result[["GeneRatio"]]) / parse_ratio(KEGG2@result[["BgRatio"]])
kegg <- KEGG2@result

barplot(GO, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
barplot(KEGG2, showCategory = 10, title = 'KEGG Pathway')
GO_top10 <- dotplot(GO, showCategory = 10, x="FoldEnrichment", orderBy = "FoldEnrichment") +
  facet_grid(ONTOLOGY~., scale="free") + 
  ggtitle("Keratinization") +
  scale_color_distiller(palette = "Spectral",direction = 1)
GO_top10



## ========================== TUMOUR BUDDING ANALYSIS ==========================
## convert the vst_count file to vsd median 
rnaseq_count <- read_fst("../data/2022-10-17_rsem_counts_mat.fst")
vst_df <- rnaseq_count
tumourcolumns <- colnames(vst_df)[grepl("SU_T", colnames(vst_df))]
crukid <- unique(substr(tumourcolumns, 1, 14))

vst_tumour <- vst_df[, c(tumourcolumns, "gene_id")]
gene_list <- vst_tumour$gene_id
vst_tumour <- vst_tumour %>% dplyr::select(-gene_id)

vsd_median_alltumour <- sapply(crukid, function(x){
  samples <- colnames(vst_tumour)[grep(x, colnames(vst_tumour))]
  if(sum(colnames(vst_tumour) %in% samples) == 1 ) {
    return(vst_tumour[, which(colnames(vst_tumour) %in% samples)])
  } else if (sum(colnames(vst_tumour) %in% samples) > 1) {
    return(rowMedians(as.matrix(vst_tumour[,which(colnames(vst_tumour) %in% samples)])))
  }
})
rownames(vsd_median_alltumour) <- gene_list

vsd_median_alltumour_df <- as.data.frame(vsd_median_alltumour) %>%
  rownames_to_column() %>%
  dplyr::rename(gene = rowname)

vsdfile <- vsd_median_alltumour_df 

# histology data
LUSC_input_plot_df_all <- input_plot_df_all %>%
  dplyr::select(cruk_tumourid, tumour_budding) %>%
  mutate(High.vs.LowMid =  case_when(
    tumour_budding %in% c("no_budding") ~ "LowMid",
    tumour_budding %in% c("have_budding") ~ "High")) %>%
  mutate(tumour_id_muttable_cruk = cruk_tumourid) %>% 
  mutate(High.vs.LowMid = factor(High.vs.LowMid , levels = c("LowMid","High"))) %>% 
  filter(!is.na(tumour_budding))

vsdfile <- vsdfile %>% dplyr::select(-CRUK0495_SU_T2)

tumourlists <- colnames(vsdfile)
tumourlists1 <- tumourlists[2:349]
tumourlist2 <- sapply(tumourlists1, function(x) {unique(idlist$tumour_id_muttable_cruk[idlist$cruk_su_tumourid == x])[1]}) %>% unname()
colnames(vsdfile)[2:349] <- tumourlist2

vsd_median <- vsdfile
exp3old <- vsd_median
exp3 <- exp3old
exp3[-1]  <- round(exp3old[-1])
rownames(exp3) <- exp3$gene
exp3 <- exp3[-1]

# separate groups
low.group <- LUSC_input_plot_df_all$cruk_tumourid[LUSC_input_plot_df_all$High.vs.LowMid == "LowMid"]
high.group <- LUSC_input_plot_df_all$cruk_tumourid[LUSC_input_plot_df_all$High.vs.LowMid == "High"]

exp4low <- exp3[,which(colnames(exp3)%in%low.group)]
exp4high <- exp3[,which(colnames(exp3)%in%high.group)]
exp4 <- cbind(exp4high,exp4low)

group_list = factor(c(rep("high",length(exp4high)),rep("low",length(exp4low))),
                    levels = c("high","low"))

################## DESeq2 for differential expressed genes (tumour budding) ##############################
countData <- exp4
condition <- group_list
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), ~ condition)
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition",levels(condition)))
RES<- as.data.frame(res)

# DEG classification
RES$group <- case_when(RES$log2FoldChange > 1 & RES$pvalue < 0.05 ~ "Up",
                       RES$log2FoldChange < -1 & RES$pvalue < 0.05 ~ "Down",
                       abs(RES$log2FoldChange) <= 1 ~ "None",
                       RES$pvalue >= 0.05 ~ "None")
RES.DEG <- RES[which(RES$group%in%c("Up","Down")),]
RES.up <- RES[which(RES$group%in%c("Up")),]
RES.down <- RES[which(RES$group%in%c("Down")),]

info <- RES.up

GO_database <- 'org.Hs.eg.db'
KEGG_database <- 'hsa'
gene <- bitr(rownames(info),fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)

GO<-enrichGO( gene$ENTREZID,
              OrgDb = GO_database,
              keyType = "ENTREZID",
              ont = "ALL",
              pvalueCutoff = 0.05,
              qvalueCutoff = 0.05,
              readable = T)
GO@result[["FoldEnrichment"]] <- parse_ratio(GO@result[["GeneRatio"]])/parse_ratio(GO@result[["BgRatio"]])
go <- GO@result

GO_top20 <- dotplot(GO, x="FoldEnrichment" ,showCategory = 20, orderBy = "FoldEnrichment")

KEGG2<-enrichKEGG(gene$ENTREZID,
                  organism = KEGG_database,
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)
KEGG2@result[["FoldEnrichment"]] <- parse_ratio(KEGG2@result[["GeneRatio"]])/parse_ratio(KEGG2@result[["BgRatio"]])
kegg <- KEGG2@result

barplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")
barplot(KEGG2,showCategory = 10,title = 'KEGG Pathway')

GO_top10 <- dotplot(GO,showCategory = 10, x="FoldEnrichment", orderBy = "FoldEnrichment")+
  facet_grid(ONTOLOGY~., scale="free") + 
  ggtitle("Tumour budding") +
  scale_color_distiller(palette = "Spectral",direction = 1)
GO_top10
#Ext fig 4e

KEGG_top20 <- dotplot(KEGG2, x="FoldEnrichment" ,showCategory = 20, orderBy = "FoldEnrichment")
KEGG_top10 <- dotplot(KEGG2,showCategory = 10, x="FoldEnrichment", orderBy = "FoldEnrichment")+
  ggtitle("Tumour budding") +
  scale_color_distiller(palette = "Spectral",direction = 1)
KEGG_top10
#Fig 3d
# The results have discrepancy with Fig 3d due to changes in KEGG version


## ========================== SINGLE CELL INVASION ==========================
## convert the vst_count file to vsd median 
rnaseq_count <- read_fst("../data/2022-10-17_rsem_counts_mat.fst")
vst_df <- rnaseq_count
tumourcolumns <- colnames(vst_df)[grepl("SU_T", colnames(vst_df))]
crukid <- unique(substr(tumourcolumns, 1, 14))

vst_tumour <- vst_df[, c(tumourcolumns, "gene_id")] 
gene_list <- vst_tumour$gene_id
vst_tumour <- vst_tumour %>% dplyr::select(-gene_id)

vsd_median_alltumour <- sapply(crukid, function(x){
  samples <- colnames(vst_tumour)[grep(x, colnames(vst_tumour))]
  if(sum(colnames(vst_tumour) %in% samples) == 1 ) {
    return(vst_tumour[, which(colnames(vst_tumour) %in% samples)])
  } else if (sum(colnames(vst_tumour) %in% samples) > 1) {
    return(rowMedians(as.matrix(vst_tumour[,which(colnames(vst_tumour) %in% samples)])))
  }
})
rownames(vsd_median_alltumour) <- gene_list

vsd_median_alltumour_df <- as.data.frame(vsd_median_alltumour) %>%
  rownames_to_column() %>%
  dplyr::rename(gene = rowname)

vsdfile <- vsd_median_alltumour_df 

# histology data
LUSC_input_plot_df_all <- input_plot_df_all %>%
  dplyr::select(cruk_tumourid, single_cell_invasion) %>%
  mutate(High.vs.LowMid =  case_when(single_cell_invasion %in% c("no_sci") ~ "LowMid",
                                     single_cell_invasion %in% c("have_sci") ~ "High")) %>%
  mutate(tumour_id_muttable_cruk = cruk_tumourid) %>% 
  mutate(High.vs.LowMid = factor(High.vs.LowMid , levels = c("LowMid","High"))) %>% 
  filter(!is.na(single_cell_invasion))

vsdfile <- vsdfile %>% dplyr::select(-CRUK0495_SU_T2)

tumourlists <- colnames(vsdfile)
tumourlists1 <- tumourlists[2:349]
tumourlist2 <- sapply(tumourlists1, function(x) {unique(idlist$tumour_id_muttable_cruk[idlist$cruk_su_tumourid == x])[1]}) %>% unname()
colnames(vsdfile)[2:349] <- tumourlist2

vsd_median <- vsdfile
exp3old <- vsd_median
exp3 <- exp3old
exp3[-1]  <- round(exp3old[-1])
rownames(exp3) <- exp3$gene
exp3 <- exp3[-1]

# separate groups
low.group <- LUSC_input_plot_df_all$cruk_tumourid[LUSC_input_plot_df_all$High.vs.LowMid == "LowMid"]
high.group <- LUSC_input_plot_df_all$cruk_tumourid[LUSC_input_plot_df_all$High.vs.LowMid == "High"]

exp4low <- exp3[,which(colnames(exp3)%in%low.group)]
exp4high <- exp3[,which(colnames(exp3)%in%high.group)]
exp4 <- cbind(exp4high,exp4low)

group_list = factor(c(rep("high",length(exp4high)),rep("low",length(exp4low))),
                    levels = c("high","low"))

################## DESeq2 for differential expressed genes (tumour budding) ##############################
countData <- exp4
condition <- group_list
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), ~ condition)
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition",levels(condition)))
RES<- as.data.frame(res)

# DEG classification
RES$group <- case_when(RES$log2FoldChange > 1 & RES$pvalue < 0.05 ~ "Up",
                       RES$log2FoldChange < -1 & RES$pvalue < 0.05 ~ "Down",
                       abs(RES$log2FoldChange) <= 1 ~ "None",
                       RES$pvalue >= 0.05 ~ "None")
RES.DEG <- RES[which(RES$group%in%c("Up","Down")),]
RES.up <- RES[which(RES$group%in%c("Up")),]
RES.down <- RES[which(RES$group%in%c("Down")),]

info <- RES.up

GO_database <- 'org.Hs.eg.db'
KEGG_database <- 'hsa'
gene <- bitr(rownames(info),fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)

GO<-enrichGO( gene$ENTREZID,
              OrgDb = GO_database,
              keyType = "ENTREZID",
              ont = "ALL",
              pvalueCutoff = 0.05,
              qvalueCutoff = 0.05,
              readable = T)
GO@result[["FoldEnrichment"]] <- parse_ratio(GO@result[["GeneRatio"]])/parse_ratio(GO@result[["BgRatio"]])
go <- GO@result

GO_top20 <- dotplot(GO, x="FoldEnrichment" ,showCategory = 20, orderBy = "FoldEnrichment")

KEGG2<-enrichKEGG(gene$ENTREZID,
                  organism = KEGG_database,
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)
KEGG2@result[["FoldEnrichment"]] <- parse_ratio(KEGG2@result[["GeneRatio"]])/parse_ratio(KEGG2@result[["BgRatio"]])
kegg <- KEGG2@result

barplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")
barplot(KEGG2,showCategory = 10,title = 'KEGG Pathway')

GO_top10 <- dotplot(GO,showCategory = 10, x="FoldEnrichment", orderBy = "FoldEnrichment")+
  facet_grid(ONTOLOGY~., scale="free") + 
  ggtitle("Single cell invasion") +
  scale_color_distiller(palette = "Spectral",direction = 1)
GO_top10
#Ext fig 4e

KEGG_top20 <- dotplot(KEGG2, x="FoldEnrichment" ,showCategory = 20, orderBy = "FoldEnrichment")
KEGG_top10 <- dotplot(KEGG2,showCategory = 10, x="FoldEnrichment", orderBy = "FoldEnrichment")+
  ggtitle("Single cell invasion") +
  scale_color_distiller(palette = "Spectral",direction = 1)
KEGG_top10
#Ext fig 4e
