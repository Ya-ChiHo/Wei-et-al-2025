library(Seurat)
options(Seurat.object.assay.version = 'v5')
library(Signac)
library(dplyr)
library(tidyverse)
library(stringr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(patchwork)
library(abind)
library(tidyverse)
library(VennDiagram)  
library(topGO)
library(fgsea)
library(msigdbr)
library(reshape2)
library(tibble)
library(circlize)

setwd("/home/yw786/palmer_scratch")

#these scripts apply fgsea and msigdbr packages to generate GSEA using RNA assay.
#shown are GSEA for CD4 TRM BACH2-high vs TRM BACH2-low cells
#similar steps were repeated for GSEA analyses on other cell subsets of interest

#additional files required:
#CD4.clean.rds generated from CD4 T cell data processing.R
#CD8.clean.rds generated from CD8 T cell data processing.R

#make GSEA sets####
m_df_H<- msigdbr(species = "Homo sapiens", category = "H")
m_df_H<- rbind(msigdbr(species = "Homo sapiens", category = "C2"), m_df_H)
m_df_H<- rbind(msigdbr(species = "Homo sapiens", category = "C7"), m_df_H)
m_df_H<- rbind(msigdbr(species = "Homo sapiens", category = "C5"), m_df_H)
fgsea_sets<- m_df_H %>% split(x = .$gene_symbol, f = .$gs_name)

#identify gsea in TRM Th1s vs TEM####
Idents(CD4.clean) <- "memory_celltype"; DefaultAssay(CD4.clean) <- "RNA.clean"

#build GSEA using genes upregulated in TRM-BACH2high
DEG.CD4.TRM_1_2.GSEA <- wilcoxauc(CD4.clean, group_by = "memory_celltype", seurat_assay = "RNA.clean")
DEG.CD4.TRM_1_2.GSEA <- subset(DEG.CD4.TRM_1_2.GSEA, DEG.CD4.TRM_1_2.GSEA$group == "TRM_1")

Rank.DEG.CD4.TRM_1_2.GSEA <- DEG.CD4.TRM_1_2.GSEA$logFC
names(Rank.DEG.CD4.TRM_1_2.GSEA) <- DEG.CD4.TRM_1_2.GSEA$feature

fgsea_DEG.CD4.TRM_1_2.GSEA <- fgsea(pathways = fgsea_sets, 
                                    stats = Rank.DEG.CD4.TRM_1_2.GSEA,
                                    eps   = 0.0, minSize=15, maxSize=500, scoreType = "pos")

#build GSEA using genes upregulated in TRM-BACH2low
DEG.CD4.TRM_2_1.GSEA <- wilcoxauc(CD4.clean, group_by = "memory_celltype", seurat_assay = "RNA.clean")
DEG.CD4.TRM_2_1.GSEA <- subset(DEG.CD4.TRM_2_1.GSEA, DEG.CD4.TRM_2_1.GSEA$group == "TRM_2")

Rank.DEG.CD4.TRM_2_1.GSEA <- DEG.CD4.TRM_2_1.GSEA$logFC
names(Rank.DEG.CD4.TRM_2_1.GSEA) <- DEG.CD4.TRM_2_1.GSEA$feature

fgsea_DEG.CD4.TRM_2_1.GSEA <- fgsea(pathways = fgsea_sets, 
                                    stats = Rank.DEG.CD4.TRM_2_1.GSEA,
                                    eps   = 0.0, minSize=15, maxSize=500, scoreType = "pos")

#identify leading edge genes in significantly enriched gene sets
#example:
fgsea_DEG.CD4.TEM_1_2.GSEA[which(fgsea_DEG.CD4.TEM_1_2.GSEA$pathway == 'GOBP_NEGATIVE_REGULATION_OF_ADAPTIVE_IMMUNE_RESPONSE'),]$leadingEdge
