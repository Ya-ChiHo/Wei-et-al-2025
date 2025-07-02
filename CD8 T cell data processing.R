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

#see after preprocessing scripts on whole gut and CD3 isolated samples
#these scripts describe steps to merge and clean CD8 T cells from whole gut and CD3 isolated samples, re-processing RNA, protein, and ATAC assays, and adding TCR data
#see CD4 T cell data processing.R for processing CD4 T cells.

#additional files required:
#Combined.CD8.rds generated from pre-processing scripts
#CD8.total.rds generated from pre-processing scripts
#CD8_TCRB_shared_clone5copy.txt

#merge CD4 T cells from whole gut samples and CD3 isolated samples####
CD8.CD3 <- CreateSeuratObject(counts = Combined.CD8[["RNA"]]$counts, assay =  "RNA")
CD8.CD3[["ATAC"]] <- Combined.CD8[["ATAC"]]
CD8.CD3[["Antibody"]] <- Combined.CD8[["Antibody"]]
CD8.CD3$nCount_ATAC <- Combined.CD8$nCount_ATAC
CD8.CD3$nFeature_ATAC <- Combined.CD8$nFeature_ATAC
CD8.CD3$percent.mt <- Combined.CD8$percent.mt
CD8.CD3$nCount_Antibody <- Combined.CD8$nCount_Antibody
CD8.CD3$nFeature_Antibody <- Combined.CD8$nFeature_Antibody
CD8.CD3$atac_fragments  <- Combined.CD8$atac_fragments 
CD8.CD3$atac_TSS_fragments  <- Combined.CD8$atac_TSS_fragments 
CD8.CD3$nucleosome_signal   <- Combined.CD8$nucleosome_signal  
CD8.CD3$TSS.enrichment   <- Combined.CD8$TSS.enrichment  

CD8.total <- CreateSeuratObject(counts = Total.CD8[["RNA"]]$counts, assay =  "RNA")
CD8.total[["ATAC"]] <- Total.CD8[["ATAC"]]
CD8.total[["Antibody"]] <- Total.CD8[["Antibody"]]
CD8.total$nCount_ATAC <- Total.CD8$nCount_ATAC
CD8.total$nFeature_ATAC <- Total.CD8$nFeature_ATAC
CD8.total$percent.mt <- Total.CD8$percent.mt
CD8.total$HIV_DNA  <- Total.CD8$HIV_DNA 
CD8.total$HIV_RNA <- Total.CD8$HIV_RNA
CD8.total$HIV_RNA_DNA <- Total.CD8$HIV_RNA_DNA
CD8.total$nCount_Antibody <- Total.CD8$nCount_Antibody
CD8.total$nFeature_Antibody <- Total.CD8$nFeature_Antibody
CD8.total$atac_fragments  <- Total.CD8$atac_fragments 
CD8.total$atac_TSS_fragments  <- Total.CD8$atac_TSS_fragments 
CD8.total$nucleosome_signal   <- Total.CD8$nucleosome_signal  
CD8.total$TSS.enrichment   <- Total.CD8$TSS.enrichment  

CD8 <- merge(x = CD8.CD3, y = CD8.total)

#perform MACS2 to recall peaks
frags.CD8.CD3 <- CD8.CD3[["ATAC"]]@fragments
frags.CD8.total <- CD8.total[["ATAC"]]@fragments
fragments <- c(frags.CD8.CD3, frags.CD8.total)
CD8[["ATAC"]]@fragments <- fragments
MACS2peaks <- CallPeaks(CD8)
MACS2peaks <- keepStandardChromosomes(MACS2peaks, pruning.mode = "coarse")
MACS2_counts <- FeatureMatrix(fragments = Fragments(CD8), features = MACS2peaks, cells = colnames(CD8))
CD8[["MACS2ATAC"]] <- CreateChromatinAssay(counts = MACS2_counts, fragments = Fragments(CD8), annotation = annotation)

#normalize and scale data####
DefaultAssay(CD8) <- "RNA"
CD8 <- JoinLayers(CD8)
CD8 <- NormalizeData(CD8, normalization.method = 'LogNormalize')
CD8 <- FindVariableFeatures(CD8, selection.method = "vst", nfeatures = 2000)
CD8 <- ScaleData(CD8, verbose = TRUE)

DefaultAssay(CD8) <- "MACS2ATAC"
CD8 <- RunTFIDF(CD8)
CD8 <- FindTopFeatures(CD8, min.cutoff = 'q0')
CD8 <- RunSVD(CD8)

DefaultAssay(CD8) <- "Antibody"
CD8 <- JoinLayers(CD8)
CD8 <- NormalizeData(CD8, assay = "Antibody", normalization.method = "CLR", margin = 2)
CD8 <- ScaleData(CD8, features = rownames(CD8), verbose = FALSE)

#batch effect correction and WNN integration####
#perform batch effect correction on ATAC assay
CD8.atac_anchors <- FindIntegrationAnchors(object.list = SplitObject(CD8, split.by = "orig.ident"),
                                           anchor.features = rownames(subset(CD8, idents = c("B029A"))), 
                                           reduction = "rlsi", dims = 2:20)
CD8.ATAC_integration <- IntegrateEmbeddings(anchorset = CD8.atac_anchors, 
                                            reductions = CD8[["lsi"]],
                                            new.reduction.name = "integrated_lsi", 
                                            dims.to.integrate = 1:20, k.weight = 30) 

CD8@reductions$integrated_lsi <- CD8.ATAC_integration@reductions$integrated_lsi
CD8 <- RunUMAP(CD4, reduction = "integrated_lsi", dims = 2:20, 
               seed.use = 42,
               reduction.name = "umap.integrated.atac", 
               reduction.key = "atacIntegratedUMAP_")

#perform batch effect correction on RNA assay
DefaultAssay(CD8) <- "RNA"
CD8 <- NormalizeData(CD8, normalization.method = 'LogNormalize') %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>% RunPCA(npcs = 50) %>% RunHarmony("orig.ident", plot_convergence = TRUE) %>%
  RunUMAP(reduction = "harmony", reduction.name = "umap.integrated.rna", reduction.key = "rnaIntegratedUMAP_", 
          dims = 1:45, assay = 'RNA', seed.use = 42) %>%
  identity()

#perform weighted nearest neighbor (WNN) using RNA and ATAC 
DefaultAssay(CD8) <- "RNA"
CD4 <- FindMultiModalNeighbors(CD4, reduction.list = list("harmony", "integrated_lsi"), 
                               dims.list = list(1:40, 2:40), modality.weight.name = "WNN.weight")
CD4 <- RunUMAP(CD4, nn.name = "weighted.nn", reduction.name = "umap.wnn", reduction.key = "wnnUMAP_", seed.use = 42)

#generate and annotate cell clusters####
#Note that the UMAP algorithm is stochastic; therefore, results (including cluster structures and numbering) will vary between runs.
#As a result, a fully reproducible record of cluster annotation steps is not provided
CD8 <- FindClusters(CD8, graph.name = "wsnn", 
                    algorithm = 3, random.seed = 42, 
                    resolution = 1, verbose = TRUE)

#some clusters were further refined by subclustering, some clusters that were non-CD8 T cells (e.g., Combined.CD8 contained unresolved T cell clusters) were removed
#please refer to supplemental files which lists RNA/protein/TF markers used for each cell cluster annotation

#An example of annotation:
CD8$major_celltype <- CD8$seurat_clusters
Idents(CD8) <- "major_celltype"
CD8 <-RenameIdents(CD8, '12_0' = 'CD4 T', '21_2' = 'CD4 T', '3_1' = 'CD4 T', '20' = 'B')
CD8$major_celltype <- Idents(CD8)
CD8.clean <- subset(CD8, idents = c("B", "CD4 T"), invert = T)

#remove MT, RPS, IGH/L RNA expression noise from cleaned CD8 T cell data
DefaultAssay(CD8.clean) <- "RNA"
CD8.clean[["RNA.clean"]] <- CD8.clean[["RNA"]]
RNA.matrix <- CD8.clean[["RNA"]]$counts
RNA.matrix <- RNA.matrix[!grepl('MT', rownames(RNA.matrix)),]
RNA.matrix <- RNA.matrix[!grepl('IGH', rownames(RNA.matrix)),]
RNA.matrix <- RNA.matrix[!grepl('IGK', rownames(RNA.matrix)),]
RNA.matrix <- RNA.matrix[!grepl('IGL', rownames(RNA.matrix)),]
RNA.matrix <- RNA.matrix[!grepl('JCHAIN', rownames(RNA.matrix)),]
RNA.matrix <- RNA.matrix[!grepl('RPL', rownames(RNA.matrix)),]
RNA.matrix <- RNA.matrix[!grepl('RPS', rownames(RNA.matrix)),]
CD8.clean[["RNA.clean"]] <- CreateAssayObject(counts = RNA.matrix )

#Important: 
#the above MACS2 peak recall, data normalization and scaling (for MACS2ATAC, RNA.clean, Antobidy), ATAC and RNA integration, 
#WNN integration, and cluster identification processes were repeated for the CD8.clean data before cell type annotation.

#some clusters were further refined by subclustering
#An example of the final CD8 T cell anotation:
CD8.clean$major_celltype <- CD8.clean$seurat_clusters
Idents(CD8.clean) <- "major_celltype"
CD8.clean <-RenameIdents(CD8.clean, '13' = 'MAIT', '11' = 'TRM Tc17', '19' = 'Proliferating',
                         '10' = 'Naive', '15' = 'gdT', '14' = 'TRM_2', '7' = 'TRM_2', '17' = 'TRM_2', '4_0' = 'TRM_2',
                         '18' = 'TRM_2', '1' = 'TRM_2', '2' = 'TRM_1', '9' = 'NKT', '3' = 'TRM_3',
                         '16' = 'Memory', '8' = 'activated', '4' = 'TEM', '0' = 'TEM', '20' = 'TRM_2',
                         '12' = 'TRM_2', '5' = 'TEM', '0_3' = 'TEFF')

#build gene accessibility####
DefaultAssay(CD8.clean) <- "ATAC"
gene.activities <- GeneActivity(CD8.clean)
CD8.clean[['activities']] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(CD8.clean) <- "activities"
CD8.clean <- NormalizeData(object = CD8.clean, assay = 'activities', 
                           normalization.method = 'LogNormalize', 
                           scale.factor = median(CD8.clean$nCount_activities))
CD8.clean <- ScaleData(CD8.clean)

#add chromVAR assay####
pfm <- getMatrixSet(x = JASPAR2022, opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))
main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
DefaultAssay(CD8.clean) <- "MACS2ATAC"
CD8.clean <- AddMotifs(CD8.clean, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)
motif.names <- CD8.clean@assays$MACS2ATAC@motifs@motif.names

register(SerialParam(stop.on.error = TRUE, log = TRUE,
                     threshold = "INFO", logdir = NA_character_, progressbar = TRUE))
DefaultAssay(CD8.clean) <- "MACS2ATAC"
CD8.clean <- RunChromVAR(CD8.clean, genome = BSgenome.Hsapiens.UCSC.hg38)

#add demultiplexed sample ID to cells####
Total.CD8.sample_ID <- as.matrix(Total.CD8$sample_ID)
Combined.CD8.sample_ID <- as.matrix(Combined.CD8$sample_ID)
CD8.sample_ID <- rbind(Total.CD8.sample_ID, Combined.CD8.sample_ID, head = TRUE)
colnames(CD8.sample_ID) <- "ID"
CD8.sample_ID <- as.data.frame(CD8.sample_ID)
CD8.clean$Sample_ID <- ifelse(colnames(CD8.clean) %in% rownames(CD8.sample_ID), 
                              CD8.sample_ID$ID[match(colnames(CD8.clean), rownames(CD8.sample_ID))],
                              FALSE)
Idents(CD8.clean) <- "Sample_ID"
CD8.clean <- RenameIdents(CD8.clean, 'B008A' = 'B008', 'B008B' = 'B008', 'B012A' = 'B012', 'B012B' = 'B012',
                          'B015A' = 'B015', 'B015B' = 'B015', 'B017A' = 'B017', 'B017B' = 'B017',
                          'B023A' = 'B023', 'B023B' = 'B023', 'B027A' = 'B027', 'B027B' = 'B027',
                          'B029A' = 'B029', 'B029B' = 'B029', 'B035A' = 'B035', 'B035B' = 'B035',
                          'B037A' = 'B037', 'B037B' = 'B037', 'B040A' = 'B040', 'B040B' = 'B040',
                          'B351A' = 'HD351', 'B351B' = 'HD351', 'B351C' = 'HD351',
                          'B357A' = 'HD357', 'B357B' = 'HD357', 'B357C' = 'HD357',
                          'B360A' = 'HD360', 'B360B' = 'HD360', 'B360C' = 'HD360',
                          'B361A' = 'HD361', 'B361B' = 'HD361', 'B361C' = 'HD361',
                          'B363A' = 'HD363', 'B363B' = 'HD363', 'B363C' = 'HD363',
                          'BEAT008' = 'B008', 'BEAT012' = 'B012', 'BEAT015' = 'B015', 'BEAT017' = 'B017',
                          'BEAT023' = 'B023', 'BEAT027' = 'B027', 'BEAT029' = 'B029', 'BEAT035' = 'B035',
                          'BEAT037' = 'B037', 'BEAT040' = 'B040', 'BEAT361' = 'HD361')
CD8.clean$Sample_ID <- Idents(CD8.clean)

CD8.clean$Condition <- CD8.clean$Sample_ID
Idents(CD8.clean) <- "Condition"
CD8.clean <-RenameIdents(CD8.clean, 'B008A' = 'infected', 'B008B' = 'infected', 'B012A' = 'infected', 'B012B' = 'infected', 
                         'B015A' = 'infected', 'B015B' = 'infected', 'B017A' = 'infected', 'B017B' = 'infected', 
                         'B023A' = 'infected', 'B023B' = 'infected', 'B027A' = 'infected', 'B027B' = 'infected', 
                         'B029A' = 'infected', 'B029B' = 'infected', 'B035A' = 'infected', 'B035B' = 'infected', 
                         'B037A' = 'infected', 'B037B' = 'infected', 'B040A' = 'infected', 'B040B' = 'infected',
                         'B361A' = 'uninfected', 'B361B' = 'uninfected', 'B361C' = 'uninfected', 
                         'BEAT361' = 'uninfected', 'BEAT008' = 'infected', 'BEAT012' = 'infected', 'BEAT015' = 'infected', 
                         'BEAT017' = 'infected', 'BEAT023' = 'infected', 'BEAT027' = 'infected', 'BEAT029' = 'infected',
                         'BEAT035' = 'infected', 'BEAT037' = 'infected', 'BEAT040' = 'infected',
                         'B351A' = 'uninfected', 'B351B' = 'uninfected', 'B351C' = 'uninfected', 
                         'B357A' = 'uninfected', 'B357B' = 'uninfected', 'B357C' = 'uninfected', 
                         'B360A' = 'uninfected', 'B360B' = 'uninfected', 'B360C' = 'uninfected', 
                         'B363A' = 'uninfected', 'B363B' = 'uninfected', 'B363C' = 'uninfected') 
CD8.clean$Condition <- Idents(CD8.clean)
table(CD8.clean$Condition)

#add clones####
TCRB_shared_clone <- read.table(file = "~/CD8_TCRB_shared_clone5copy.txt", header = TRUE)
TCRB_shared_clone$is.TCR <- "TCRB"

CD8.clean$is.clone <- ifelse(colnames(CD8.clean) %in% TCRB_shared_clone$BC, 
                             TCRB_shared_clone$is.clone[match(colnames(CD8.clean), TCRB_shared_clone$BC)], FALSE)
CD8.clean$Clone_ID <- ifelse(colnames(CD8.clean) %in% TCRB_shared_clone$BC, 
                             TCRB_shared_clone$Clone_ID[match(colnames(CD8.clean), TCRB_shared_clone$BC)],FALSE)
CD8.clean$Clone_size <- ifelse(colnames(CD8.clean) %in% TCRB_shared_clone$BC, 
                               TCRB_shared_clone$Clone_size[match(colnames(CD8.clean), TCRB_shared_clone$BC)],FALSE)
CD8.clean$is.TCR <- ifelse(colnames(CD8.clean) %in% TCRB_shared_clone$BC, "TCRB", FALSE)



