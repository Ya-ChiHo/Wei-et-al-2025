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
#these scripts describe steps to merge and clean CD4 T cells from whole gut and CD3 isolated samples, re-processing RNA, protein, and ATAC assays, and adding TCR data
#see CD8 T cell data processing.R for processing CD8 T cells.

#additional files required:
#Combined.CD4.rds generated from pre-processing scripts
#CD4.total.rds generated from pre-processing scripts
#TCRB_shared_clone5copy.txt (for both CD4 and CD8 T cells)

#merge CD4 T cells from whole gut samples and CD3 isolated samples####
CD4.CD3 <- CreateSeuratObject(counts = Combined.CD4[["RNA"]]$counts, assay =  "RNA")
CD4.CD3[["ATAC"]] <- Combined.CD4[["ATAC"]]
CD4.CD3[["Antibody"]] <- Combined.CD4[["Antibody"]]
CD4.CD3$nCount_ATAC <- Combined.CD4$nCount_ATAC
CD4.CD3$nFeature_ATAC <- Combined.CD4$nFeature_ATAC
CD4.CD3$percent.mt <- Combined.CD4$percent.mt
CD4.CD3$HIV_DNA  <- Combined.CD4$HIV_DNA 
CD4.CD3$HIV_RNA <- Combined.CD4$HIV_RNA
CD4.CD3$HIV_RNA_DNA <- Combined.CD4$HIV_RNA_DNA
CD4.CD3$nCount_Antibody <- Combined.CD4$nCount_Antibody
CD4.CD3$nFeature_Antibody <- Combined.CD4$nFeature_Antibody
CD4.CD3$atac_fragments  <- Combined.CD4$atac_fragments 
CD4.CD3$atac_TSS_fragments  <- Combined.CD4$atac_TSS_fragments 
CD4.CD3$nucleosome_signal   <- Combined.CD4$nucleosome_signal  
CD4.CD3$TSS.enrichment   <- Combined.CD4$TSS.enrichment  

CD4.total <- CreateSeuratObject(counts = Total.CD4[["RNA"]]$counts, assay =  "RNA")
CD4.total[["ATAC"]] <- Total.CD4[["ATAC"]]
CD4.total[["Antibody"]] <- Total.CD4[["Antibody"]]
CD4.total$nCount_ATAC <- Total.CD4$nCount_ATAC
CD4.total$nFeature_ATAC <- Total.CD4$nFeature_ATAC
CD4.total$percent.mt <- Total.CD4$percent.mt
CD4.total$HIV_DNA  <- Total.CD4$HIV_DNA 
CD4.total$HIV_RNA <- Total.CD4$HIV_RNA
CD4.total$HIV_RNA_DNA <- Total.CD4$HIV_RNA_DNA
CD4.total$nCount_Antibody <- Total.CD4$nCount_Antibody
CD4.total$nFeature_Antibody <- Total.CD4$nFeature_Antibody
CD4.total$atac_fragments  <- Total.CD4$atac_fragments 
CD4.total$atac_TSS_fragments  <- Total.CD4$atac_TSS_fragments 
CD4.total$nucleosome_signal   <- Total.CD4$nucleosome_signal  
CD4.total$TSS.enrichment   <- Total.CD4$TSS.enrichment  

CD4 <- merge(x = CD4.CD3, y = CD4.total)

#perform MACS2 to recall peaks
frags.CD4.CD3 <- CD4.CD3[["ATAC"]]@fragments
frags.CD4.total <- CD4.total[["ATAC"]]@fragments
fragments <- c(frags.CD4.CD3, frags.CD4.total)
CD4[["ATAC"]]@fragments <- fragments
MACS2peaks <- CallPeaks(CD4)
MACS2peaks <- keepStandardChromosomes(MACS2peaks, pruning.mode = "coarse")
MACS2_counts <- FeatureMatrix(fragments = Fragments(CD4), features = MACS2peaks, cells = colnames(CD4))
CD4[["MACS2ATAC"]] <- CreateChromatinAssay(counts = MACS2_counts, fragments = Fragments(CD4), annotation = annotation)

#normalize and scale data####
DefaultAssay(CD4) <- "RNA"
CD4 <- JoinLayers(CD4)
CD4 <- NormalizeData(CD4, normalization.method = 'LogNormalize')
CD4 <- FindVariableFeatures(CD4, selection.method = "vst", nfeatures = 2000)
CD4 <- ScaleData(CD4, verbose = TRUE)

DefaultAssay(CD4) <- "MACS2ATAC"
CD4 <- RunTFIDF(CD4)
CD4 <- FindTopFeatures(CD4, min.cutoff = 'q0')
CD4 <- RunSVD(CD4)

DefaultAssay(CD4) <- "Antibody"
CD4 <- JoinLayers(CD4)
CD4 <- NormalizeData(CD4, assay = "Antibody", normalization.method = "CLR", margin = 2)
CD4 <- ScaleData(CD4, features = rownames(CD4), verbose = FALSE)

#batch effect correction and WNN integration####
#perform batch effect correction on ATAC assay
CD4.atac_anchors <- FindIntegrationAnchors(object.list = SplitObject(CD4, split.by = "orig.ident"),
                                           anchor.features = rownames(subset(CD4, idents = c("B029A"))), 
                                           reduction = "rlsi", dims = 2:20)
CD4.ATAC_integration <- IntegrateEmbeddings(anchorset = CD4.atac_anchors, 
                                            reductions = CD4[["lsi"]],
                                            new.reduction.name = "integrated_lsi", 
                                            dims.to.integrate = 1:20, k.weight = 30) 

CD4@reductions$integrated_lsi <- CD4.ATAC_integration@reductions$integrated_lsi
CD4 <- RunUMAP(CD4, reduction = "integrated_lsi", dims = 2:20, 
               seed.use = 42,
               reduction.name = "umap.integrated.atac", 
               reduction.key = "atacIntegratedUMAP_")

#perform batch effect correction on RNA assay
DefaultAssay(CD4) <- "RNA"
CD4 <- NormalizeData(CD4, normalization.method = 'LogNormalize') %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>% RunPCA(npcs = 50) %>% RunHarmony("orig.ident", plot_convergence = TRUE) %>%
  RunUMAP(reduction = "harmony", reduction.name = "umap.integrated.rna", reduction.key = "rnaIntegratedUMAP_", 
          dims = 1:40, assay = 'RNA', seed.use = 42) %>%
  identity()

#perform weighted nearest neighbor (WNN) using RNA and ATAC 
DefaultAssay(CD4) <- "RNA"
CD4 <- FindMultiModalNeighbors(CD4, reduction.list = list("harmony", "integrated_lsi"), 
                               dims.list = list(1:40, 2:40), modality.weight.name = "WNN.weight")
CD4 <- RunUMAP(CD4, nn.name = "weighted.nn", reduction.name = "umap.wnn", reduction.key = "wnnUMAP_", seed.use = 42)

#generate and annotate cell clusters####
#Note that the UMAP algorithm is stochastic; therefore, results (including cluster structures and numbering) will vary between runs.
#As a result, a fully reproducible record of cluster annotation steps is not provided
CD4 <- FindClusters(CD4, graph.name = "wsnn", 
                    algorithm = 3, random.seed = 42, 
                    resolution = 1, verbose = TRUE)

#some clusters were further refined by subclustering, some clusters that were non-CD4 T cells (e.g., Combined.CD4 contained unresolved T cell clusters) were removed
#please refer to supplemental files which lists RNA/protein/TF markers used for each cell cluster annotation

#An example of annotation:
CD4$major_celltype <- CD4$seurat_clusters
Idents(CD4) <- "major_celltype"
CD4 <-RenameIdents(CD4, '10' = 'CD8 T', '13' = 'CD8 T', '14' = 'CD8 T', '20' = 'CD8 T', '15' = 'CD8 T', 
                   '19' = 'B', '5' = 'CD8 T',
                   '0' = 'CD4 T', '1' = 'CD4 T', '2' = 'CD4 T', '3' = 'CD4 T', '4' = 'CD4 T', 
                   '6' = 'CD4 T', '7' = 'CD4 T', '8' = 'CD4 T', '9' = 'CD4 T', '11' = 'CD4 T', '12' = 'CD4 T',
                   '16' = 'CD4 T', '17' = 'CD4 T', '18' = 'CD4 T', '21' = 'CD4 T', '22' = 'CD4 T', 
                   '23' = 'CD4 T', '24' = 'CD4 T', '25' = 'CD4 T')
CD4$major_celltype <- Idents(CD4)
CD4.clean <- subset(CD4, idents = c("B", "CD8 T"), invert = T)

#remove MT, RPS, IGH/L RNA expression noise from cleaned CD4 T cell data
DefaultAssay(CD4.clean) <- "RNA"
CD4.clean[["RNA.clean"]] <- CD4.clean[["RNA"]]
RNA.matrix <- CD4.clean[["RNA"]]$counts
RNA.matrix <- RNA.matrix[!grepl('MT', rownames(RNA.matrix)),]
RNA.matrix <- RNA.matrix[!grepl('IGH', rownames(RNA.matrix)),]
RNA.matrix <- RNA.matrix[!grepl('IGK', rownames(RNA.matrix)),]
RNA.matrix <- RNA.matrix[!grepl('IGL', rownames(RNA.matrix)),]
RNA.matrix <- RNA.matrix[!grepl('JCHAIN', rownames(RNA.matrix)),]
RNA.matrix <- RNA.matrix[!grepl('RPL', rownames(RNA.matrix)),]
RNA.matrix <- RNA.matrix[!grepl('RPS', rownames(RNA.matrix)),]
CD4.clean[["RNA.clean"]] <- CreateAssayObject(counts = RNA.matrix )

#Important: 
#the above MACS2 peak recall, data normalization and scaling (for MACS2ATAC, RNA.clean, Antobidy), ATAC and RNA integration, 
#WNN integration, and cluster identification processes were repeated for the CD4.clean data before cell type annotation.

#some clusters were further refined by subclustering
#An example of the final CD4 T cell anotation:
CD4.clean$major_celltype <- CD4.clean$seurat_clusters
Idents(CD4.clean) <- "major_celltype"
CD4.clean <-RenameIdents(CD4.clean, '5' = 'Naive', '6' = 'Treg', '14' = 'Treg', '10' = 'iTreg',
                         '19' = 'Tfr', '17' = 'Proliferating', '18' = 'Th17',
                         '7' = 'Th1', 
                         '1' = 'TCM',
                         '12' = 'TRM', '0' = 'TRM-Th1-1', '16' = 'TRM-Th1-1',
                         '8' = 'TRM-Th1-1', 
                         '9' = 'TFH', '4' = 'TRM-Th1-2', 
                         '2' = 'TRM-Th17', '13' = 'TRM-Th17', '15' = 'TRM-Th17', 
                         '3' = 'TCM', '1' = 'TCM', '20' = 'TCM', '21' = 'TRM-Th1-2')
CD4.clean$major_celltype <- Idents(CD4.clean)

#build gene accessibility####
DefaultAssay(CD4.clean) <- "ATAC"
gene.activities <- GeneActivity(CD4.clean)
CD4.clean[['activities']] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(CD4.clean) <- "activities"
CD4.clean <- NormalizeData(object = CD4.clean, assay = 'activities', 
                           normalization.method = 'LogNormalize', 
                           scale.factor = median(CD4.clean$nCount_activities))
CD4.clean <- ScaleData(CD4.clean)

#add chromVAR assay####
pfm <- getMatrixSet(x = JASPAR2022, opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))
main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
DefaultAssay(CD4.clean) <- "MACS2ATAC"
CD4.clean <- AddMotifs(CD4.clean, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)
motif.names <- CD4.clean@assays$MACS2ATAC@motifs@motif.names

register(SerialParam(stop.on.error = TRUE, log = TRUE,
                     threshold = "INFO", logdir = NA_character_, progressbar = TRUE))
DefaultAssay(CD4.clean) <- "MACS2ATAC"
CD4.clean <- RunChromVAR(CD4.clean, genome = BSgenome.Hsapiens.UCSC.hg38)

#add demultiplexed sample ID to cells####
Total.CD4.sample_ID <- as.matrix(Total.CD4$sample_ID)
Combined.CD4.sample_ID <- as.matrix(Combined.CD4$sample_ID)
CD4.sample_ID <- rbind(Total.CD4.sample_ID, Combined.CD4.sample_ID, head = TRUE)
colnames(CD4.sample_ID) <- "ID"
CD4.sample_ID <- as.data.frame(CD4.sample_ID)
CD4.clean$Sample_ID <- ifelse(colnames(CD4.clean) %in% rownames(CD4.sample_ID), 
                              CD4.sample_ID$ID[match(colnames(CD4.clean), rownames(CD4.sample_ID))],
                              FALSE)
Idents(CD4.clean) <- "Sample_ID"
CD4.clean <- RenameIdents(CD4.clean, 'B008A' = 'B008', 'B008B' = 'B008', 'B012A' = 'B012', 'B012B' = 'B012',
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
CD4.clean$Sample_ID <- Idents(CD4.clean)

CD4.clean$Condition <- CD4.clean$Sample_ID
Idents(CD4.clean) <- "Condition"
CD4.clean <-RenameIdents(CD4.clean, 'B008A' = 'infected', 'B008B' = 'infected', 'B012A' = 'infected', 'B012B' = 'infected', 
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
CD4.clean$Condition <- Idents(CD4.clean)
table(CD4.clean$Condition)

#add clones####
TCRB_shared_clone <- read.table(file = "~/TCRB_shared_clone5copy.txt", header = TRUE)
TCRB_shared_clone$is.TCR <- "TCRB"

CD4.clean$is.clone <- ifelse(colnames(CD4.clean) %in% TCRB_shared_clone$BC, 
                                   TCRB_shared_clone$is.clone[match(colnames(CD4.clean), TCRB_shared_clone$BC)], FALSE)
CD4.clean$Clone_ID <- ifelse(colnames(CD4.clean) %in% TCRB_shared_clone$BC, 
                                   TCRB_shared_clone$Clone_ID[match(colnames(CD4.clean), TCRB_shared_clone$BC)],FALSE)
CD4.clean$Clone_size <- ifelse(colnames(CD4.clean) %in% TCRB_shared_clone$BC, 
                                     TCRB_shared_clone$Clone_size[match(colnames(CD4.clean), TCRB_shared_clone$BC)],FALSE)
CD4.clean$is.TCR <- ifelse(colnames(CD4.clean) %in% TCRB_shared_clone$BC, "TCRB", FALSE)



