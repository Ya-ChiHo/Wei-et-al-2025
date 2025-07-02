library(Seurat)
options(Seurat.object.assay.version = 'v5')
library(Signac)
library(ggplot2)
library(ggraph)
library(cowplot)
library(patchwork)
library(dplyr)
library(GenomicRanges)
library(future)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyverse)
library(GenomeInfoDb)
library(TFBSTools)
library(JASPAR2022)
library(harmony)
library(clustree)
library(presto)
library(reshape2)
library(presto)
library(tibble)

#see after CD4/CD8 T cell data processing scripts
#these scripts describe steps to merge processed CD4 T and CD8 T cell Seurat objects 

#additional files required:
#CD4.clean.rds generated from CD4 T cell data processing scripts
#CD8.clean.rds generated from CD8 T cell data processing scripts

#merge CD4 T cells from whole gut samples and CD3 isolated samples####
CD8 <- CreateSeuratObject(counts = CD8.clean[["RNA"]]$counts, assay =  "RNA")
CD8[["Antibody"]] <- CD8.clean[["Antibody"]]
CD8[["ATAC"]] <- CD8.clean[["ATAC"]]
CD8$nCount_ATAC <- CD8.clean$nCount_ATAC
CD8$nFeature_ATAC <- CD8.clean$nFeature_ATAC
CD8$percent.mt <- CD8.clean$percent.mt
CD8$HIV_DNA  <- CD8.clean$HIV_DNA 
CD8$HIV_RNA <- CD8.clean$HIV_RNA
CD8$HIV_RNA_DNA <- CD8.clean$HIV_RNA_DNA
CD8$nCount_Antibody <- CD8.clean$nCount_Antibody
CD8$nFeature_Antibody <- CD8.clean$nFeature_Antibody
CD8$major_celltype   <- CD8.clean$major_celltype  
CD8$Condition   <- CD8.clean$Condition  
CD8$Sample_ID   <- CD8.clean$Sample_ID 

CD4 <- CreateSeuratObject(counts = CD4.clean[["RNA"]]$counts, assay =  "RNA")
CD4[["Antibody"]] <- CD4.clean[["Antibody"]]
CD4[["ATAC"]] <- CD4.clean[["ATAC"]]
CD4$nCount_ATAC <- CD4.clean$nCount_ATAC
CD4$nFeature_ATAC <- CD4.clean$nFeature_ATAC
CD4$percent.mt <- CD4.clean$percent.mt
CD4$HIV_DNA  <- CD4.clean$HIV_DNA 
CD4$HIV_RNA <- CD4.clean$HIV_RNA
CD4$HIV_RNA_DNA <- CD4.clean$HIV_RNA_DNA
CD4$nCount_Antibody <- CD4.clean$nCount_Antibody
CD4$nFeature_Antibody <- CD4.clean$nFeature_Antibody
CD4$major_celltype   <- CD4.clean$major_celltype  
CD4$Condition   <- CD4.clean$Condition  
CD4$Sample_ID   <- CD4.clean$Sample_ID  

#remove duplicate from datasets and transfer cell type annotation labels
all.equal(colnames(CD4), colnames(CD8))
joint.bcs<- intersect(colnames(CD4), colnames(CD8))
joint.bcs <- as.data.frame(joint.bcs)
CD8$join_bc <- ifelse(colnames(CD8) %in% joint.bcs$joint.bcs, "TRUE", "FALSE")
Idents(CD8) <- "join_bc"
CD8 <- subset(CD8, idents = c("FALSE"))
CD8$join_bc <- NULL

CD4_CD8 <- merge(x = CD4, y = list(CD8))

rm(CD4); rm(CD8)

CD4.clean$Cells <- "CD4"

CD4_CD8$Cells <- ifelse(colnames(CD4_CD8) %in% colnames(CD4.clean), 
                        CD4.clean$Cells[match(colnames(CD4_CD8), colnames(CD4.clean))], "CD8")

CD4_CD8$celltype <- str_c(CD4_CD8$Cells, '_', CD4_CD8$major_celltype)

#perform MACS2 to recall peaks
frags.CD4 <- CD4.clean[["ATAC"]]@fragments
frags.CD8 <- CD8.clean[["ATAC"]]@fragments
fragments <- c(frags.CD4, frags.CD8)
CD4_CD8[["ATAC"]]@fragments <- fragments
MACS2peaks <- CallPeaks(CD4_CD8)
MACS2peaks <- keepStandardChromosomes(MACS2peaks, pruning.mode = "coarse")
MACS2_counts <- FeatureMatrix(fragments = Fragments(CD4_CD8), features = MACS2peaks, cells = colnames(CD4_CD8))
CD4_CD8[["MACS2ATAC"]] <- CreateChromatinAssay(counts = MACS2_counts, fragments = Fragments(CD4_CD8), annotation = annotation)

#normalize and scale data####
DefaultAssay(CD4_CD8) <- "RNA"
CD4_CD8 <- JoinLayers(CD4_CD8)
CD4_CD8 <- NormalizeData(CD4_CD8, normalization.method = 'LogNormalize')
CD4_CD8 <- FindVariableFeatures(CD4_CD8, selection.method = "vst", nfeatures = 2000)
CD4_CD8 <- ScaleData(CD4_CD8, verbose = TRUE)

DefaultAssay(CD4_CD8) <- "MACS2ATAC"
CD4_CD8 <- RunTFIDF(CD4_CD8)
CD4_CD8 <- FindTopFeatures(CD4_CD8, min.cutoff = 'q0')
CD4_CD8 <- RunSVD(CD4_CD8)

DefaultAssay(CD4_CD8) <- "Antibody"
CD4_CD8 <- JoinLayers(CD4_CD8)
CD4_CD8 <- NormalizeData(CD4_CD8, assay = "Antibody", normalization.method = "CLR", margin = 2)
CD4_CD8 <- ScaleData(CD4_CD8, features = rownames(CD4_CD8), verbose = FALSE)

#batch effect correction and WNN integration####
#perform batch effect correction on ATAC assay
CD4_CD8.atac_anchors <- FindIntegrationAnchors(object.list = SplitObject(CD4_CD8, split.by = "orig.ident"),
                                           anchor.features = rownames(subset(CD4_CD8, idents = c("B029A"))), 
                                           reduction = "rlsi", dims = 2:15)
CD4_CD8.ATAC_integration <- IntegrateEmbeddings(anchorset = CD4_CD8.atac_anchors, 
                                            reductions = CD4_CD8[["lsi"]],
                                            new.reduction.name = "integrated_lsi", 
                                            dims.to.integrate = 1:15, k.weight = 30) 

CD4_CD8@reductions$integrated_lsi <- CD4_CD8.ATAC_integration@reductions$integrated_lsi
CD4_CD8 <- RunUMAP(CD4_CD8, reduction = "integrated_lsi", dims = 2:15, 
               seed.use = 42,
               reduction.name = "umap.integrated.atac", 
               reduction.key = "atacIntegratedUMAP_")

#perform batch effect correction on RNA assay
DefaultAssay(CD4_CD8) <- "RNA"
CD4_CD8 <- NormalizeData(CD4_CD8, normalization.method = 'LogNormalize') %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>% RunPCA(npcs = 50) %>% RunHarmony("orig.ident", plot_convergence = TRUE) %>%
  RunUMAP(reduction = "harmony", reduction.name = "umap.integrated.rna", reduction.key = "rnaIntegratedUMAP_", 
          dims = 1:35, assay = 'RNA', seed.use = 42) %>%
  identity()

#perform weighted nearest neighbor (WNN) using RNA and ATAC 
DefaultAssay(CD4_CD8) <- "RNA"
CD4_CD8 <- FindMultiModalNeighbors(CD4_CD8, reduction.list = list("harmony", "integrated_lsi"), 
                               dims.list = list(1:35, 2:15), modality.weight.name = "WNN.weight")
CD4_CD8 <- RunUMAP(CD4_CD8, nn.name = "weighted.nn", reduction.name = "umap.wnn", reduction.key = "wnnUMAP_", seed.use = 42)

#build gene accessibility####
DefaultAssay(CD4_CD8) <- "ATAC"
gene.activities <- GeneActivity(CD4_CD8)
CD4_CD8[['activities']] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(CD4_CD8) <- "activities"
CD4_CD8 <- NormalizeData(object = CD4_CD8, assay = 'activities', 
                           normalization.method = 'LogNormalize', 
                           scale.factor = median(CD4_CD8$nCount_activities))
CD4_CD8 <- ScaleData(CD4_CD8)

#add chromVAR assay####
pfm <- getMatrixSet(x = JASPAR2022, opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))
main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
DefaultAssay(CD4_CD8) <- "MACS2ATAC"
CD4_CD8 <- AddMotifs(CD4_CD8, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)
motif.names <- CD4_CD8@assays$MACS2ATAC@motifs@motif.names

register(SerialParam(stop.on.error = TRUE, log = TRUE,
                     threshold = "INFO", logdir = NA_character_, progressbar = TRUE))
DefaultAssay(CD4_CD8) <- "MACS2ATAC"
CD4_CD8 <- RunChromVAR(CD4_CD8, genome = BSgenome.Hsapiens.UCSC.hg38)
