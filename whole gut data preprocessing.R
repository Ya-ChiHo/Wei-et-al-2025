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
library(scDblFinder)

#see before preprocessing scripts on CD3 isolation samples
#these scripts describe steps to build and process RNA, protein, and ATAC assays for each whole gut sample (no CD3 isolation)
#these IDs contain participant samples that were multiplexed (by hashing antibody) and pooled: BT1, BT2, BT3, BT4
#these IDs contain each contain only 1 sample: B351A, B351B, B351C, B357A, B357B, B357C, B360A, B360B, B360C, B363A, B363B, B363C (e.g., B351A and B351B are the sample samples where cells were split for separate GEM emulsion and sequencing)
#hashtag info:
#BT1 contains samples B035 (hashtag 1), B037 (hashtag 2), B040 (hashtag 3)
#BT2 contains samples B008 (hashtag 4), B012 (hashtag 5), B015 (hashtag 6)
#BT3 contains samples B017 (hashtag 5), B023 (hashtag 6), B027 (hashtag 7), B029 (hashtag 8)
#BT4 contains samples B361 (hashtag 4), L158 (not needed, hashtag 5), L167 (not needed, hashtag 6), L202 (not needed, hashtag 7), L211 (not needed, hashtag 8)

#additional files needed:
#gut_dogma_total_demuxafy_results
#gut_HIV_cells_all

#build and merge RNA Seurat objects####
#repeat for each samples: BT1, BT2, BT3, BT4, B351A, B351B, B351C, B357A, B357B, B357C, B360A, B360B, B360C, B363A, B363B, B363C

BT1 <- Read10X("~/cellranger/BT1/outs/raw_feature_bc_matrix")
BT1_RNA <- CreateSeuratObject(counts = BT1$`Gene Expression`, assay =  "RNA")
BT1_RNA$orig.ident <- 'BT1'
BT1_RNA[["percent.mt"]] <- PercentageFeatureSet(BT1_RNA, pattern = "^MT-")

#merge all samples
Total_RNA <- merge(x = BT1_RNA, y = list(BT2_RNA, BT3_RNA, BT4_RNA, 
                                         B351A_RNA, B351B_RNA, B351C_RNA,
                                         B357A_RNA, B357B_RNA, B357C_RNA,
                                         B360A_RNA, B360B_RNA, B360C_RNA,
                                         B363A_RNA, B363B_RNA, B363C_RNA),
                   add.cell.ids = c("BT1", "BT2", "BT3", "BT4",
                                    "B351A", "B351B", "B351C",
                                    "B357A", "B357B", "B357C",
                                    "B360A", "B360B", "B360C",
                                    "B363A", "B363B", "B363C"))
Total_RNA <- JoinLayers(Total_RNA)

#remove poor quality cells by RNA QC metrics
Total_RNA <- subset(Total_RNA, nCount_RNA > 500 & nFeature_RNA > 240 & percent.mt < 25.00000)

#add protein assay to Seurat object####
#repeat for each samples: BT1, BT2, BT3, BT4, B351A, B351B, B351C, B357A, B357B, B357C, B360A, B360B, B360C, B363A, B363B, B363C

BT1_cite <- CreateSeuratObject(counts = BT1$`Antibody Capture`, assay =  "Antibody")

#merge protein seurat objects
Total_cite <- merge(x = BT1_cite, y = list(BT2_cite, BT3_cite, BT4_cite, 
                                           B351A_cite, B351B_cite, B351C_cite,
                                           B357A_cite, B357B_cite, B357C_cite,
                                           B360A_cite, B360B_cite, B360C_cite,
                                           B363A_cite, B363B_cite, B363C_cite),
                    add.cell.ids = c("BT1", "BT2", "BT3", "BT4",
                                     "B351A", "B351B", "B351C",
                                     "B357A", "B357B", "B357C",
                                     "B360A", "B360B", "B360C",
                                     "B363A", "B363B", "B363C"))

Total_cite <- JoinLayers(Total_cite)

#add protein assay to RNA Seurat objects
all.equal(colnames(Total_RNA), colnames(Total_cite))
joint.bcs<- intersect(colnames(Total_RNA), colnames(Total_cite))
joint.bcs <- as.data.frame(joint.bcs)
Total_cite$join_bc <- ifelse(colnames(Total_cite) %in% joint.bcs$joint.bcs, "TRUE", "FALSE")
Idents(Total_cite) <- "join_bc"
Total_cite <- subset(Total_cite, idents = c("TRUE"))

Total <- Total_RNA
Total[["Antibody"]] <- Total_cite[["Antibody"]] 
Total$nCount_Antibody  <- Total_cite$nCount_Antibody   
Total$nFeature_Antibody   <- Total_cite$nFeature_Antibody   

#make hashtag assay and demultiplex pooled, and remove doublets####
#repeat for each samples: BT1, BT2, BT3, BT4 only
#hashtag 5 = B351

#make a Seurat object containing the hashtag assay for each sample 
Idents(Total) <- "orig.ident"; DefaultAssay(Total) <- "Antibody"
BT1_hash <- subset(Total, idents = c("BT1"))
hashtag_BT1_counts <- BT1_hash@assays$Antibody$counts
hashtag_BT1_counts <- hashtag_BT1_counts[rownames(hashtag_BT1_counts) %in% c(
  "hashtag1-TotalA", "hashtag2-TotalA", "hashtag3-TotalA", "hashtag4-TotalA", "hashtag5-TotalA", "hashtag6-TotalA", "hashtag7-TotalA", "hashtag8-TotalA"),]

BT1_hash[["hashtag"]] <- CreateAssayObject(counts = hashtag_BT1_counts)
BT1_hash <- NormalizeData(BT1_hash, assay = "hashtag", normalization.method = 'CLR', margin = 2)

#demultiplex
BT1_hash <- MULTIseqDemux(BT1_hash, assay = "hashtag",  
                          autoThresh  = TRUE, maxiter = 10, qrange = seq(from = 0.1, to = 0.999, by = 0.001), verbose = TRUE)

#merge samples
Total_hash <- merge(x = BT1_hash, y = list(BT2_hash, BT3_hash, BT4_hash))
ncol(Total_hash)
table(Total_hash$MULTI_ID, Total_hash$orig.ident)

#also performed was freemuxlet, to demultiplex by SNPs and detect doublets, from https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/Demuxlet.html
#see document: gut_dogma_total_demuxafy_results
demuxafy_snp <- read.table("~/gut_dogma_total_demuxafy_results.txt", head = T)
rownames(demuxafy_snp) <- demuxafy_snp$BARCODE
Total$demuxafy_dropletType <- ifelse(colnames(Total) %in% demuxafy_snp$BARCODE, demuxafy_snp$DROPLET.TYPE[match(colnames(Total), demuxafy_snp$BARCODE)],FALSE)
Total$demuxafy_dropletSNP <- ifelse(colnames(Total) %in% demuxafy_snp$BARCODE, demuxafy_snp$SNG.BEST.GUESS[match(colnames(Total), demuxafy_snp$BARCODE)],FALSE)
Total$hashtag_snp_match <- str_c(Total$MULTI_ID, '_', Total$demuxafy_dropletSNP)
Total$sample_snp <- str_c(Total$orig.ident, '_', Total$demuxafy_dropletSNP)
table(Total$sample_snp, Total$MULTI_ID)
Total$sample_ID <- Total$sample_snp

Idents(Total) <- "sample_ID"
Total <-RenameIdents(Total, 'BT1_0' = 'BEAT035', 'BT1_2' = 'BEAT037', 'BT1_1' = 'BEAT040',
                     'BT2_0' = 'BEAT015', 'BT2_1' = 'BEAT012', 'BT2_2' = 'BEAT008',
                     'BT3_0' = 'BEAT029', 'BT3_1' = 'BEAT027', 'BT3_2' = 'BEAT017', 'BT3_3' = 'BEAT023',
                     'BT4_0' = 'BEAT361', 'BT4_1' = 'L167', 'BT4_2' = 'L202', 'BT4_3' = 'L158', 'BT4_4' = 'L211',
                     "B351A_0" = "B351A", "B351B_0" = "B351B", "B351C_0" = "B351C", 
                     "B357A_0" = "B357A", "B357B_0" = "B357B", "B357C_0" = "B357C", 
                     "B360A_0" = "B360A", "B360B_0" = "B360B", "B360C_0" = "B360C", 
                     "B363A_0" = "B363A", "B363B_0" = "B363B", "B363C_0" = "B363C")
Total$sample_ID <- Idents(Total)

#remove all samples with ID: L### (for a different study)
Idents(Total) <- "sample_ID"
Total <- subset(Total, idents = c("L167", "L202", "L158", "L211"), invert = T)

#in addition to freemuxlet and MULTIseqDemux, we will further remove doublets by scDblFInder
#repeat for each samples: BT1, BT2, BT3, BT4, B351A, B351B, B351C, B357A, B357B, B357C, B360A, B360B, B360C, B363A, B363B, B363C

bp <- MulticoreParam(3, RNGseed=1234)
Idents(Total) <- "orig.ident"
sce_BT1 <- subset(Total, idents = c("BT1"))
sce_BT1_RNA <- sce_BT1@assays$RNA@layers$counts
sce_BT1_RNA <- scDblFinder(sce_BT1_RNA, BPPARAM = bp)
sce_BT1$scDblFInder <- sce_BT1_RNA$scDblFinder.class
Total_sce <- merge(x = sce_BT1, y = list(sce_BT2, sce_BT3, sce_BT4,
                                         sce_B351A, sce_B351B, sce_B351C,
                                         sce_B357A, sce_B357B, sce_B357C,
                                         sce_B360A, sce_B360B, sce_B360C,
                                         sce_B363A, sce_B363B, sce_B363C))
Total$scDblFInder <- Total_sce$scDblFInder
Total$scdblfinder_demuxafy_singlets <- str_c(Total$scDblFInder, '_', Total$demuxafy_dropletType)

#we will now remove doublets detected by freemuxlet and MULTIseqDemux, and remove any cells with mismatch calling from the 2 methods
Idents(Total) <- "scdblfinder_demuxafy_singlets"
Total_singlets <- subset(Total, idents = c("singlet_SNG", "singlet_FALSE"))

#add ATAC Seurat to Seurat object####
#repeat for each samples: BT1, BT2, BT3, BT4, B351A, B351B, B351C, B357A, B357B, B357C, B360A, B360B, B360C, B363A, B363B, B363C

BT1.peaks <- read.table(file = "~a/cellranger_ARC/BEATrun1Total/outs/atac_peaks.bed",
                        col.names = c("chr", "start", "end"))
BT1.peaks.gr <- makeGRangesFromDataFrame(BT1.peaks)

#merge peaks
Total.peaks <- reduce(x = c(BT1.peaks.gr, BT2.peaks.gr, BT3.peaks.gr, BT4.peaks.gr,
                            B351A.peaks.gr, B351B.peaks.gr, B351C.peaks.gr,
                            B357A.peaks.gr, B357B.peaks.gr, B357C.peaks.gr,
                            B360A.peaks.gr, B360B.peaks.gr, B360C.peaks.gr,
                            B363A.peaks.gr, B363B.peaks.gr, B363C.peaks.gr))
peakwidths <- width(Total.peaks)
Total.peaks <- Total.peaks[peakwidths < 10000 & peakwidths > 20]
BT1.md <- read.table(
  file = "~/cellranger_ARC/BEATrun1Total/outs/per_barcode_metrics.csv",
  stringsAsFactors = FALSE, sep = ",", header = TRUE, row.names = 1)[-1, ]
BT1.frags <- CreateFragmentObject(
  path = "~/cellranger_ARC/BEATrun1Total/outs/atac_fragments.tsv.gz",
  cells = rownames(BT1.md))
BT1.counts <- FeatureMatrix( fragments = BT1.frags, features = Total.peaks, cells = rownames(BT1.md))
BT1.assay <- CreateChromatinAssay(BT1.counts, fragments = BT1.frags)
BT1.ATAC <- CreateSeuratObject(BT1.assay, assay = "ATAC", meta.data=BT1.md)
BT1.ATAC$orig.ident <- 'BT1'

#merge ATAC Seurat objects
Total.ATAC <- merge(x = BT1.ATAC, y = list(BT2.ATAC, BT3.ATAC, BT4.ATAC,
                                           B351A.ATAC, B351B.ATAC, B351C.ATAC,
                                           B357A.ATAC, B357B.ATAC, B357C.ATAC,
                                           B360A.ATAC, B360B.ATAC, B360C.ATAC,
                                           B363A.ATAC, B363B.ATAC, B363C.ATAC),
                    add.cell.ids = c("BT1", "BT2", "BT3", "BT4",
                                     "B351A", "B351B", "B351C",
                                     "B357A", "B357B", "B357C",
                                     "B360A", "B360B", "B360C",
                                     "B363A", "B363B", "B363C"))

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
Annotation(Total.ATAC) <- annotation

Total.ATAC <- NucleosomeSignal(Total.ATAC)
Total.ATAC <- TSSEnrichment(Total.ATAC)
Total.ATAC$high.tss <- ifelse(Total.ATAC$TSS.enrichment > 2, 'High', 'Low')

#add ATAC assay to RNA+protein Seurat object
all.equal(colnames(Total_singlets), colnames(Total.ATAC))
joint.bcs<- intersect(colnames(Total_singlets), colnames(Total.ATAC))
joint.bcs <- as.data.frame(joint.bcs)
Total.ATAC$join_bc <- ifelse(colnames(Total.ATAC) %in% joint.bcs$joint.bcs, "TRUE", "FALSE")
Idents(Total.ATAC) <- "join_bc"
Total.ATAC <- subset(Total.ATAC, idents = c("TRUE"))
Total_singlets[["ATAC"]] <- Total.ATAC[["ATAC"]] 
Total_singlets$atac_barcode <- Total.ATAC$atac_barcode   
Total_singlets$atac_raw_reads <- Total.ATAC$atac_raw_reads    
Total_singlets$atac_fragments <- Total.ATAC$atac_fragments    
Total_singlets$atac_TSS_fragments   <- Total.ATAC$atac_TSS_fragments    
Total_singlets$nCount_ATAC <- Total.ATAC$nCount_ATAC     
Total_singlets$nFeature_ATAC <- Total.ATAC$nFeature_ATAC  
Total_singlets$nucleosome_signal <- Total.ATAC$nucleosome_signal  
Total_singlets$TSS.enrichment <- Total.ATAC$TSS.enrichment  
Total_singlets$high.tss <- Total.ATAC$high.tss 

#remove poor quality cells by ATAC QC metrics
Total_singlets <- subset(Total_singlets, nucleosome_signal < 1 & TSS.enrichment > 2 & nCount_ATAC > 350 & nCount_ATAC < 20000 & nFeature_ATAC > 200)

#perform MACS2 to recall peaks
DefaultAssay(Total_singlets) <- "ATAC"
MACS2peaks <- CallPeaks(Total_singlets)
MACS2peaks <- keepStandardChromosomes(MACS2peaks, pruning.mode = "coarse")
MACS2_counts <- FeatureMatrix(fragments = Fragments(Total_singlets), features = MACS2peaks, cells = colnames(Total_singlets))

frag.list <- list(BT1.frags, BT2.frags, BT3.frags, BT4.frags,
                  B351A.frags, B351B.frags, B351C.frags,
                  B357A.frags, B357B.frags, B357C.frags,
                  B360A.frags, B360B.frags, B360C.frags,
                  B363A.frags, B363B.frags, B363C.frags)

Total_singlets[["MACS2ATAC"]] <- CreateChromatinAssay(counts = MACS2_counts, fragments = frag.list, annotation = annotation)

#add JASPAR2022 motifs to integrated ATAC object
pfm <- getMatrixSet(x = JASPAR2022, opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))
main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
DefaultAssay(Total_singlets) <- "MACS2ATAC"
Total_singlets <- AddMotifs(Total_singlets, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)

#batch effect correction and WNN integration####
#perform batch effect correction on ATAC assay
DefaultAssay(Total_singlets) <- "MACS2ATAC"
Total_singlets <- RunTFIDF(Total_singlets)
Total_singlets <- FindTopFeatures(Total_singlets, min.cutoff = 'q0')
Total_singlets <- RunSVD(Total_singlets)

Total.atac_anchors <- FindIntegrationAnchors(object.list = SplitObject(Total_singlets, split.by = "orig.ident"),
                                             anchor.features = rownames(subset(Total, idents = c("BT1"))), 
                                             reduction = "rlsi", dims = 2:20)
Total.ATAC.integrated <- IntegrateEmbeddings(anchorset = Total.atac_anchors, 
                                             reductions = Total_singlets[["lsi"]],
                                             new.reduction.name = "integrated_lsi", 
                                             dims.to.integrate = 1:20, k.weight = 50) 
DefaultAssay(Total.ATAC.integrated) <- "MACS2ATAC"
Total.ATAC.integrated <- RunUMAP(Total.ATAC.integrated, reduction = "integrated_lsi", dims = 2:20, 
                                 seed.use = 42,
                                 reduction.name = "umap.atac", reduction.key = "atacUMAP_")

#perform batch effect correction on RNA assay
DefaultAssay(Total.ATAC.integrated) <- "RNA"

Total.ATAC.integrated <- NormalizeData(Total.ATAC.integrated, normalization.method = 'LogNormalize') %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>% RunPCA(npcs = 35) %>% RunHarmony("orig.ident", plot_convergence = TRUE) %>%
  RunUMAP(reduction = "harmony", reduction.name = "umap.rna", reduction.key = "rnaUMAP_", dims = 1:35, assay = 'RNA', seed.use = 42) %>%
  identity()

#perform weighted nearest neighbor (WNN) using RNA and ATAC 
Total.ATAC.integrated <- FindMultiModalNeighbors(Total.ATAC.integrated, reduction.list = list("harmony", "integrated_lsi"), 
                                                 dims.list = list(1:35, 2:20), modality.weight.name = "WNN.weight")
Total.ATAC.integrated <- RunUMAP(Total.ATAC.integrated, nn.name = "weighted.nn", 
                                 reduction.name = "umap.wnn", reduction.key = "wnnUMAP_", seed.use = 42)

#normalize and scale data####
DefaultAssay(Total.ATAC.integrated) <- "RNA"
Total.ATAC.integrated <- NormalizeData(Total.ATAC.integrated, normalization.method = 'LogNormalize')
Total.ATAC.integrated <- ScaleData(Total.ATAC.integrated, features = rownames(Total.ATAC.integrated), verbose = TRUE)

DefaultAssay(Total.ATAC.integrated) <- "MACS2ATAC"
Total.ATAC.integrated <- RunTFIDF(Total.ATAC.integrated)
Total.ATAC.integrated <- FindTopFeatures(Total.ATAC.integrated, min.cutoff = 'q0')
Total.ATAC.integrated <- RunSVD(Total.ATAC.integrated)

DefaultAssay(Total.ATAC.integrated) <- "Antibody"
Total.ATAC.integrated <- JoinLayers(Total.ATAC.integrated)
Total.ATAC.integrated <- NormalizeData(Total.ATAC.integrated, assay = "Antibody", normalization.method = "CLR", margin = 2)
Total.ATAC.integrated <- ScaleData(Total.ATAC.integrated, features = rownames(Total.ATAC.integrated), verbose = FALSE)

#build gene accessibility####
DefaultAssay(Total.ATAC.integrated) <- "ATAC"
gene.activities <- GeneActivity(Total.ATAC.integrated)
Total.ATAC.integrated[['activities']] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(Total.ATAC.integrated) <- "activities"
Total.ATAC.integrated <- NormalizeData(object = Total.ATAC.integrated, assay = 'activities', 
                                       normalization.method = 'LogNormalize', 
                                       scale.factor = median(Total.ATAC.integrated$nCount_activities))
Total.ATAC.integrated <- ScaleData(Total.ATAC.integrated)

#add chromVAR assay####
register(SerialParam(stop.on.error = TRUE, log = TRUE,
                     threshold = "INFO", logdir = NA_character_, progressbar = TRUE))
DefaultAssay(Total.ATAC.integrated) <- "MACS2ATAC"
Total.ATAC.integrated <- RunChromVAR(Total.ATAC.integrated, genome = BSgenome.Hsapiens.UCSC.hg38)

#add HIV cells to object####
HIV_pos <- read.table("~/BEAT_dogma/gut_HIV_cells_all.txt", head = T)
Total.ATAC.integrated$HIV_DNA <- ifelse(colnames(Total.ATAC.integrated) %in% HIV_pos$BC, HIV_pos$HIV_DNApos[match(colnames(Total.ATAC.integrated), HIV_pos$BC)],FALSE)
Total.ATAC.integrated$HIV_RNA <- ifelse(colnames(Total.ATAC.integrated) %in% HIV_pos$BC, HIV_pos$HIV_RNApos[match(colnames(Total.ATAC.integrated), HIV_pos$BC)],FALSE)
Total.ATAC.integrated$HIV_RNA_DNA <- str_c(Total.ATAC.integrated$HIV_RNA, '_', Total.ATAC.integrated$HIV_DNA)

#add condition identification####
table(Total.ATAC.integrated$sample_ID)
Total.ATAC.integrated$Condition <- Total.ATAC.integrated$sample_ID
Idents(Total.ATAC.integrated) <- "Condition"
Total.ATAC.integrated <-RenameIdents(Total.ATAC.integrated, 'B351A' = 'uninfected', 'B351B' = 'uninfected', 'B351C' = 'uninfected',
                         'B357A' = 'uninfected', 'B357B' = 'uninfected', 'B357C' = 'uninfected',
                         'B360A' = 'uninfected', 'B360B' = 'uninfected', 'B360C' = 'uninfected',
                         'B363A' = 'uninfected', 'B363B' = 'uninfected', 'B363C' = 'uninfected',
                         'BEAT008' = 'infected', 'BEAT012' = 'infected', 'BEAT015' = 'infected',
                         'BEAT017' = 'infected', 'BEAT023' = 'infected', 'BEAT027' = 'infected',
                         'BEAT029' = 'infected', 'BEAT035' = 'infected', 'BEAT037' = 'infected',
                         'BEAT040' = 'infected', 'BEAT361' = 'uninfected')
Total.ATAC.integrated$Condition <- Idents(Total.ATAC.integrated)
table(Total.ATAC.integrated$Condition)
table(Total.ATAC.integrated$major_celltype, Total.ATAC.integrated$Condition)

#generate and annotate cell clusters####
#Note that the UMAP algorithm is stochastic; therefore, results (including cluster structures and numbering) will vary between runs.
#As a result, a fully reproducible record of cluster annotation steps is not provided

Total.ATAC.integrated <- FindClusters(Total.ATAC.integrated, graph.name = "wsnn", 
                                      algorithm = 3, random.seed = 42, 
                                      resolution = 1, verbose = TRUE)

#some clusters were further refined by subclustering, based on expression of key celltype defining markers, for example:
Idents(Total.ATAC.integrated) <- "seurat_clusters"
Total.ATAC.integrated <- FindSubCluster(Total.ATAC.integrated, cluster = "10", graph.name = "wsnn", 
                                        subcluster.name = "sub.seurat_clusters", resolution = 1, algorithm = 3)
Idents(Total.ATAC.integrated) <- "sub.seurat_clusters"
Total.ATAC.integrated <-RenameIdents(Total.ATAC.integrated, "10_5" = "24", "10_0" = "10", "10_1" = "10", "10_2" = "10",
                                     "10_3" = "10", "10_4" = "10", "10_6" = "10")
Total.ATAC.integrated$seurat_clusters <- Idents(Total.ATAC.integrated)

#please refer to supplemental files which lists RNA/protein/TF markers used for each cell cluster annotation
#an example of annotated cell types based on marker expressions:
Total.ATAC.integrated <-RenameIdents(Total.ATAC.integrated, '4' = 'CD4 T', '12' = 'CD4 T', '6' = 'CD4 T',
                         '3' = 'CD8 T', '11' = 'CD8 T', '10' = 'CD8 T', '24' = 'NK',
                         '0' = 'IgA Plasma', '13' = 'IgA Plasma', '1' = 'IgA Plasma', '5' = 'IgA Plasma', '7' = 'IgA Plasma',
                         '14' = 'IgG Plasma', '21' = 'MAST', '17' = 'Proliferating B',
                         '2' = 'memory B', '8' = 'Naive B',  '9' = 'Colonic subepithelial', 
                         '18' = 'Endothelial', '19' = 'Enteric Glia', '22' = 'Myofibroblast',
                         '23' = 'Myofibroblast', '20' = 'Epithelial', 
                         '16' = 'Myenteric ganglia', '15' = 'CD4 Monocytes')
Total.ATAC.integrated$major_celltype <- Idents(Total.ATAC.integrated)

#finally, we will subset out CD4 and CD8 T cells, to be merged with the CD3 selected cells from these same samples
Idents(Total_GUT) <- "major_celltype"
Total.CD8 <- subset(Total_GUT, idents = c("CD8 T"))
Total.CD4 <- subset(Total_GUT, idents = c("CD4 T"))


