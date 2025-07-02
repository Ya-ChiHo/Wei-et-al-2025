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

#see after preprocessing scripts on whole gut samples
#these scripts describe steps to build and process RNA, protein, and ATAC assays for CD3 isolated samples (same samples as whole gut data)
#these IDs contain participant samples that were multiplexed (by hashing antibody) and pooled: BHDA, BHDB, BHDC
#these IDs contain each contain only 1 sample: B008A, B008B, B012A, B012B, B015A, B015B, B017A, B017B, B023A, B023B, B027A, B027B, B029A, B029B, B035A, B035B, B037A, B037B, B040A, B040B  (e.g., B008A and B008B are the sample samples where cells were split for separate GEM emulsion and sequencing)
#hashtag info:
#BHDA contains samples B361A (hashtag 4), L158A (not needed, hashtag 5), L167A (not needed, hashtag 6), L202A (not needed, hashtag 7), L211A (not needed, hashtag 8)
#BHDB contains samples B361B (hashtag 4), L158B (not needed, hashtag 5), L167B (not needed, hashtag 6), L202B (not needed, hashtag 7), L211B (not needed, hashtag 8)
#BHDC contains samples B361C (hashtag 4), L158C (not needed, hashtag 5), L167C (not needed, hashtag 6), L202C (not needed, hashtag 7), L211C (not needed, hashtag 8)

#additional files needed:
#gut_dogma_CD3uninfected_demuxafy_results
#gut_HIV_cells_all

#build and merge RNA Seurat objects####
#repeat for each samples: B008A, B008B, B012A, B012B, B015A, B015B, B017A, B017B, B023A, B023B, B027A, B027B, B029A, B029B, B035A, B035B, B037A, B037B, B040A, B040B, BHDA, BHDB, BHDC

B008A <- Read10X("~/cellranger/B008A/outs/raw_feature_bc_matrix")
B008A_RNA <- CreateSeuratObject(counts = B008A$`Gene Expression`, assay =  "RNA")
B008A_RNA$orig.ident <- 'B008A'
B008A_RNA[["percent.mt"]] <- PercentageFeatureSet(B008A_RNA, pattern = "^MT-")

#merge all samples
Combined_RNA <- merge(x = B008A_RNA, y = list(B008B_RNA, B012A_RNA, B012B_RNA,
                                              B015A_RNA, B015B_RNA, B017A_RNA, B017B_RNA,
                                              B023A_RNA, B023B_RNA, B027A_RNA, B027B_RNA,
                                              B029A_RNA, B029B_RNA, B035A_RNA, B035B_RNA, 
                                              B037A_RNA, B037B_RNA, B040A_RNA, B040B_RNA, 
                                              BHDA_RNA, BHDB_RNA, BHDC_RNA),
                      add.cell.ids = c("B008A", "B008B", "B012A", "B012B",
                                       "B015A", "B015B", "B017A", "B017B",
                                       "B023A", "B023B", "B027A", "B027B",
                                       "B029A", "B029B", "B035A", "B035B", 
                                       "B037A", "B037B", "B040A", "B040B",
                                       "BHDA", "BHDB", "BHDC"))
Combined_RNA <- JoinLayers(Combined_RNA)

#remove poor quality cells by RNA QC metrics
Combined_RNA <- subset(Combined_RNA, nCount_RNA > 500 & nFeature_RNA > 240 & percent.mt < 25.000001)

#add protein assay to Seurat object####
#repeat for each samples: B008A, B008B, B012A, B012B, B015A, B015B, B017A, B017B, B023A, B023B, B027A, B027B, B029A, B029B, B035A, B035B, B037A, B037B, B040A, B040B, BHDA, BHDB, BHDC
B008A_cite <- CreateSeuratObject(counts = B008A$`Antibody Capture`, assay =  "Antibody")

#merge protein seurat objects
Combined_cite <- merge(x = B008A_cite, y = list(B008B_cite, B012A_cite, B012B_cite,
                                                B015A_cite, B015B_cite, B017A_cite, B017B_cite,
                                                B023A_cite, B023B_cite, B027A_cite, B027B_cite,
                                                B029A_cite, B029B_cite, B035A_cite, B035B_cite, 
                                                B037A_cite, B037B_cite, B040A_cite, B040B_cite,
                                                BHDA_cite, BHDB_cite, BHDC_cite),
                       add.cell.ids = c("B008A", "B008B", "B012A", "B012B", "B015A", "B015B",
                                        "B017A", "B017B", "B023A", "B023B", "B027A", "B027B",
                                        "B029A", "B029B","B035A", "B035B", "B037A", "B037B", "B040A", "B040B",
                                        "BHDA", "BHDB", "BHDC"))
Combined_cite <- JoinLayers(Combined_cite)

#add protein assay to RNA Seurat objects
all.equal(colnames(Combined), colnames(Combined_cite))
joint.bcs<- intersect(colnames(Combined), colnames(Combined_cite))
joint.bcs <- as.data.frame(joint.bcs)
Combined_cite$join_bc <- ifelse(colnames(Combined_cite) %in% joint.bcs$joint.bcs, "TRUE", "FALSE")
Idents(Combined_cite) <- "join_bc"
Combined_cite <- subset(Combined_cite, idents = c("TRUE"))

Combined[["Antibody"]] <- Combined_cite[["Antibody"]] 
Combined$nCount_Antibody  <- Combined_cite$nCount_Antibody   
Combined$nFeature_Antibody   <- Combined_cite$nFeature_Antibody   

#add ATAC Seurat to Seurat object####
B008A.peaks <- read.table(file = "~/cellranger_ARC/BEAT008A/outs/atac_peaks.bed",
                          col.names = c("chr", "start", "end"))
B008A.peaks.gr <- makeGRangesFromDataFrame(B008A.peaks)

#merge peaks
combined.peaks <- GenomicRanges::reduce(x = c(B008A.peaks.gr, B008B.peaks.gr, B012A.peaks.gr, B012B.peaks.gr,
                                              B015A.peaks.gr, B015B.peaks.gr, B017A.peaks.gr, B017B.peaks.gr,
                                              B023A.peaks.gr, B023B.peaks.gr, B027A.peaks.gr, B027B.peaks.gr,
                                              B029A.peaks.gr, B029B.peaks.gr, B035A.peaks.gr, B035B.peaks.gr, 
                                              B037A.peaks.gr, B037B.peaks.gr, B040A.peaks.gr, B040B.peaks.gr,
                                              BHDA.peaks.gr, BHDB.peaks.gr, BHDC.peaks.gr))
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths < 10000 & peakwidths > 20]
B008A.md <- read.table(
  file = "~/cellranger_ARC/BEAT008A/outs/per_barcode_metrics.csv",
  stringsAsFactors = FALSE, sep = ",", header = TRUE, row.names = 1)[-1, ]
B008A.frags <- CreateFragmentObject(
  path = "~/cellranger_ARC/BEAT008A/outs/atac_fragments.tsv.gz",
  cells = rownames(B008A.md))
B008A.counts <- FeatureMatrix( fragments = B008A.frags, features = combined.peaks, cells = rownames(B008A.md))
B008A.assay <- CreateChromatinAssay(B008A.counts, fragments = B008A.frags)
B008A.ATAC <- CreateSeuratObject(B008A.assay, assay = "ATAC", meta.data=B008A.md)
B008A.ATAC$orig.ident <- 'B008A'

#merge ATAC Seurat objects
Combined.ATAC <- merge(x = B008A.ATAC, y = list(B008B.ATAC, B012A.ATAC, B012B.ATAC,
                                                B015A.ATAC, B015B.ATAC, B017A.ATAC, B017B.ATAC,
                                                B023A.ATAC, B023B.ATAC, B027A.ATAC, B027B.ATAC,
                                                B029A.ATAC, B029B.ATAC, B035A.ATAC, B035B.ATAC, 
                                                B037A.ATAC, B037B.ATAC, B040A.ATAC, B040B.ATAC,
                                                BHDA.ATAC, BHDB.ATAC, BHDC.ATAC),
                       add.cell.ids = c("B008A", "B008B", "B012A", "B012B", "B015A", "B015B",
                                        "B017A", "B017B", "B023A", "B023B", "B027A", "B027B",
                                        "B029A", "B029B", "B035A", "B035B", "B037A", "B037B", "B040A", "B040B",
                                        "BHDA", "BHDB", "BHDC"))
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
Annotation(Combined.ATAC) <- annotation

pfm <- getMatrixSet(x = JASPAR2022, opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))
main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
DefaultAssay(Combined.ATAC) <- "MACS2ATAC"
Combined.ATAC <- AddMotifs(Combined.ATAC, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)
motif.names <- Combined.ATAC@assays$MACS2ATAC@motifs@motif.names

Combined.ATAC <- NucleosomeSignal(Combined.ATAC)
Combined.ATAC <- TSSEnrichment(Combined.ATAC)
Combined.ATAC$high.tss <- ifelse(Combined.ATAC$TSS.enrichment > 2, 'High', 'Low')

#add ATAC assay to RNA+protein Seurat object
all.equal(colnames(Combined), colnames(Combined.ATAC))
joint.bcs<- intersect(colnames(Combined), colnames(Combined.ATAC))
joint.bcs <- as.data.frame(joint.bcs)
Combined.ATAC$join_bc <- ifelse(colnames(Combined.ATAC) %in% joint.bcs$joint.bcs, "TRUE", "FALSE")
Idents(Combined.ATAC) <- "join_bc"
Combined.ATAC <- subset(Combined.ATAC, idents = c("TRUE"))
Combined[["ATAC"]] <- Combined.ATAC[["ATAC"]] 
Combined$atac_barcode <- Combined.ATAC$atac_barcode   
Combined$atac_raw_reads <- Combined.ATAC$atac_raw_reads    
Combined$atac_fragments <- Combined.ATAC$atac_fragments    
Combined$atac_TSS_fragments   <- Combined.ATAC$atac_TSS_fragments    
Combined$nCount_ATAC <- Combined.ATAC$nCount_ATAC     
Combined$nFeature_ATAC <- Combined.ATAC$nFeature_ATAC  
Combined$nucleosome_signal <- Combined.ATAC$nucleosome_signal  
Combined$TSS.enrichment <- Combined.ATAC$TSS.enrichment  
Combined$high.tss <- Combined.ATAC$high.tss 

#remove poor quality cells by ATAC QC metrics
Combined <- subset(Combined, nucleosome_signal < 1 & TSS.enrichment > 2 & nCount_ATAC > 350 & nCount_ATAC < 20000 & nFeature_ATAC > 200)

#add HIV cells to object####
HIV_pos <- read.table("~/gut_HIV_cells_all.txt", head = T)
Combined$HIV_DNA <- ifelse(colnames(Combined) %in% HIV_pos$BC, HIV_pos$HIV_DNApos[match(colnames(Combined), HIV_pos$BC)],FALSE)
Combined$HIV_RNA <- ifelse(colnames(Combined) %in% HIV_pos$BC, HIV_pos$HIV_RNApos[match(colnames(Combined), HIV_pos$BC)],FALSE)
Combined$HIV_RNA_DNA <- str_c(Combined$HIV_RNA, '_', Combined$HIV_DNA)

#remove doublets####
#repeat for each samples: B008A, B008B, B012A, B012B, B015A, B015B, B017A, B017B, B023A, B023B, B027A, B027B, B029A, B029B, B035A, B035B, B037A, B037B, B040A, B040B, BHDA, BHDB, BHDC

bp <- MulticoreParam(3, RNGseed=1234)
Idents(Combined) <- "orig.ident"; DefaultAssay(Combined) <- "RNA"
B008A_RNA1 <- subset(Combined, idents = c("B008A"))
sce_B008A_RNA <- B008A_RNA1@assays$RNA@layers$counts
sce_B008A_RNA <- scDblFinder(sce_B008A_RNA, BPPARAM = bp)
B008A_RNA1$scDblFInder <- sce_B008A_RNA$scDblFinder.class
Combined_scdblfinder <- merge(x = B008A_RNA1, y = list(B008B_RNA1, B012A_RNA1, B012B_RNA1,
                                                       B015A_RNA1, B015B_RNA1, 
                                                       B017A_RNA1, B017B_RNA1, B023A_RNA1, B023B_RNA1,
                                                       B027A_RNA1, B027B_RNA1, B029A_RNA1, B029B_RNA1, 
                                                       B035A_RNA1, B035B_RNA1, B037A_RNA1, B037B_RNA1, 
                                                       B040A_RNA1, B040B_RNA1, BHDA_RNA1, BHDB_RNA1, BHDC_RNA1))

Combined$scDblFInder <- Combined_scdblfinder$scDblFInder
Combined <- subset(Combined, idents = c("singlet"))

#make hashtag assay and demultiplex pooled samples####
#these steps are for the uninfected donor CD3 isolated samples (BHDA, BHDB, BHDC), which were each pooled, not for HIV+ CD3 isolated samples, which were not pooled
Idents(Combined1) <- "orig.ident"; DefaultAssay(Combined1) <- "Antibody"
Combined_hashtag_counts <- Combined1@assays$Antibody$counts
Combined_hashtag_counts <- Combined_hashtag_counts[rownames(Combined_hashtag_counts) %in% c(
  "hashtag1-TotalA", "hashtag2-TotalA", "hashtag3-TotalA", "hashtag4-TotalA", "hashtag5-TotalA", "hashtag6-TotalA", "hashtag7-TotalA", "hashtag8-TotalA"),]
Combined1[["hashtag"]] <- CreateAssayObject(counts = Combined_hashtag_counts)
Combined1 <- NormalizeData(Combined1, assay = "hashtag", normalization.method = 'CLR', margin = 2)

#demultiplex
Combined1 <- MULTIseqDemux(Combined1, assay = "hashtag",  
                           autoThresh  = TRUE, maxiter = 10, qrange = seq(from = 0.1, to = 0.999, by = 0.001), verbose = TRUE)

#also performed was freemuxlet, to demultiplex by SNPs and detect doublets, from https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/Demuxlet.html
#see document: gut_dogma_CD3uninfected_demuxafy_results
demuxafy_snp <- read.table("~/gut_dogma_CD3uninfected_demuxafy_results.txt", head = T)
Combined$demuxafy_dropletType <- ifelse(colnames(Combined) %in% demuxafy_snp$BARCODE, demuxafy_snp$DROPLET.TYPE[match(colnames(Combined), demuxafy_snp$BARCODE)],FALSE)
Combined$demuxafy_dropletSNP <- ifelse(colnames(Combined) %in% demuxafy_snp$BARCODE, demuxafy_snp$SNG.BEST.GUESS[match(colnames(Combined), demuxafy_snp$BARCODE)],FALSE)
Combined$hashtag_snp_match <- str_c(Combined$MULTI_ID, '_', Combined$demuxafy_dropletSNP)
Combined$sample_snp <- str_c(Combined$orig.ident, '_', Combined$demuxafy_dropletSNP)
Combined$sample_ID <- Combined$sample_snp

Idents(Combined) <- "sample_ID"
Combined <-RenameIdents(Combined, 'BHDA_0' = 'B361A', 'BHDA_1' = 'L202A', 'BHDA_2' = 'L211A', 'BHDA_3' = 'L167A', 'BHDA_4' = 'L158A',
                        'BHDB_0' = 'L167B', 'BHDB_1' = 'B361B', 'BHDB_2' = 'L202B', 'BHDB_3' = 'L211B', 'BHDB_4' = 'L158B', 
                        'BHDC_0' = 'B361C', 'BHDC_1' = 'L167C', 'BHDC_2' = 'L202C', 'BHDC_3' = 'L211C', 'BHDC_4' = 'L158C',
                        'BHDB_FALSE' = 'B361B', 'BHDC_FALSE' = 'B361C',
                        'B008A_FALSE' = 'B008A', 'B008B_FALSE' = 'B008B', 'B012A_FALSE' = 'B012A', 'B012B_FALSE' = 'B012B',
                        'B015A_FALSE' = 'B015A', 'B015B_FALSE' = 'B015B', 'B017A_FALSE' = 'B017A', 'B017B_FALSE' = 'B017B',
                        'B023A_FALSE' = 'B023A', 'B023B_FALSE' = 'B023B', 'B027A_FALSE' = 'B027A', 'B027B_FALSE' = 'B027B',
                        'B029A_FALSE' = 'B029A', 'B029B_FALSE' = 'B029B', 'B035A_FALSE' = 'B035A', 'B035B_FALSE' = 'B035B',
                        'B037A_FALSE' = 'B037A', 'B037B_FALSE' = 'B037B', 'B040A_FALSE' = 'B040A', 'B040B_FALSE' = 'B040B')
Combined$sample_ID <- Idents(Combined)

#remove all samples with ID: L### (for a different study)
Idents(Combined) <- "sample_ID"
Combined <- subset(Combined, idents = c("L167A", "L167B", "L167C", "L202A", "L202B", "L202C", 
                                        "L158A", "L158B", "L158C", "L211A", "L211B", "L211C"), invert = T)

#we will now remove doublets detected by freemuxlet and MULTIseqDemux, and remove any cells with mismatch calling from the 2 methods
Combined$scdblfinder_demuxafy_singlets <- str_c(Combined$scDblFInder, '_', Combined$demuxafy_dropletType)
Idents(Combined) <- "scdblfinder_demuxafy_singlets"
Combined <- subset(Combined, idents = c("singlet_SNG", "singlet_FALSE"))

#perform MACS2 to recall peaks####
DefaultAssay(Combined) <- "ATAC"
MACS2peaks <- CallPeaks(Combined)
MACS2peaks <- keepStandardChromosomes(MACS2peaks, pruning.mode = "coarse")
MACS2_counts <- FeatureMatrix(fragments = Fragments(Combined), features = MACS2peaks, cells = colnames(Combined))

frag.list <- list(B008A.frags, B008B.frags, B012A.frags, B012B.frags,
                  B015A.frags, B015B.frags, B017A.frags, B017B.frags,
                  B023A.frags, B023B.frags, B027A.frags, B027B.frags,
                  B029A.frags, B029B.frags, B035A.frags, B035B.frags, 
                  B037A.frags, B037A.frags, B040A.frags, B040B.frags, 
                  BHDA.frags, BHDB.frags, BHDC.frags)

Combined[["MACS2ATAC"]] <- CreateChromatinAssay(counts = MACS2_counts, fragments = frag.list, annotation = annotation)

#batch effect correction and WNN integration####
#perform batch effect correction on ATAC assay
DefaultAssay(Combined) <- "MACS2ATAC"
Combined <- RunTFIDF(Combined)
Combined <- FindTopFeatures(Combined, min.cutoff = 'q0')
Combined <- RunSVD(Combined)

Combined.atac_anchors <- FindIntegrationAnchors(object.list = SplitObject(Combined, split.by = "orig.ident"),
                                                anchor.features = rownames(subset(Combined, idents = c("B029A"))), 
                                                reduction = "rlsi", dims = 2:30)
Combined.ATAC.integrated <- IntegrateEmbeddings(anchorset = Combined.atac_anchors, 
                                                reductions = Combined[["lsi"]],
                                                new.reduction.name = "integrated_lsi", 
                                                dims.to.integrate = 1:30, k.weight = 50) 

Combined@reductions$integrated_lsi  <- Combined.ATAC.integrated@reductions$integrated_lsi
Combined <- RunUMAP(Combined, reduction = "integrated_lsi", dims = 2:25, 
                    seed.use = 42,
                    reduction.name = "umap.atac", reduction.key = "atacUMAP_", assay = 'MACS2ATAC')

#perform batch effect correction on RNA assay
DefaultAssay(Combined) <- "RNA"
Combined <- NormalizeData(Combined, normalization.method = 'LogNormalize') %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>% RunPCA(npcs = 50) %>% RunHarmony("orig.ident", plot_convergence = TRUE) %>%
  RunUMAP(reduction = "harmony", reduction.name = "umap.rna", reduction.key = "rnaUMAP_", 
          dims = 1:30, assay = 'RNA', seed.use = 42) %>%
  identity()

#perform weighted nearest neighbor (WNN) using RNA and ATAC 
Combined <- FindMultiModalNeighbors(Combined, reduction.list = list("harmony", "integrated_lsi"), 
                                    dims.list = list(1:20, 2:20), modality.weight.name = "WNN.weight")
Combined <- RunUMAP(Combined, nn.name = "weighted.nn", 
                    reduction.name = "umap.wnn", reduction.key = "wnnUMAP_", seed.use = 42)

#normalize and scale data####
DefaultAssay(Combined) <- "RNA"
Combined <- NormalizeData(Combined, normalization.method = 'LogNormalize')
Combined <- ScaleData(Combined, features = rownames(Combined), verbose = TRUE)

DefaultAssay(Combined) <- "MACS2ATAC"
Combined <- RunTFIDF(Combined)
Combined <- FindTopFeatures(Combined, min.cutoff = 'q0')
Combined <- RunSVD(Combined)

DefaultAssay(Combined) <- "Antibody"
VariableFeatures(Combined) <- antibody_features
Combined <- NormalizeData(Combined, assay = "Antibody", normalization.method = "CLR", margin = 2)
Combined <- ScaleData(Combined, features = antibody_features, verbose = FALSE)

#build gene accessibility####
DefaultAssay(Combined) <- "ATAC"
gene.activities <- GeneActivity(Combined)
Combined[['activities']] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(Combined) <- "activities"
Combined <- NormalizeData(object = Combined, assay = 'activities', 
                          normalization.method = 'LogNormalize', 
                          scale.factor = median(Combined$nCount_activities))
Combined <- ScaleData(Combined)

#add chromVAR assay####
register(SerialParam(stop.on.error = TRUE, log = TRUE,
                     threshold = "INFO", logdir = NA_character_, progressbar = TRUE))
DefaultAssay(Combined) <- "MACS2ATAC"
Combined <- RunChromVAR(Combined, genome = BSgenome.Hsapiens.UCSC.hg38)

#generate and annotate cell clusters####
#Note that although these samples were CD3 selected, isolation is not 100% and other celltypes remain, thus cells were annotated and non-T cell clusters were removed
#Note that the UMAP algorithm is stochastic; therefore, results (including cluster structures and numbering) will vary between runs.
#As a result, a fully reproducible record of cluster annotation steps is not provided

Combined <- FindClusters(Combined, graph.name = "wsnn", 
                         algorithm = 3, random.seed = 42, 
                         resolution = 1, verbose = TRUE)

#As the first step, we annotated cells to remove non-T cells from CD3 selected samples. 
#please refer to supplemental files which lists RNA/protein/TF markers used for each cell cluster annotation
#For example:
Combined$major_celltype <- Combined$seurat_clusters
Idents(Combined) <- "major_celltype"
Combined <-RenameIdents(Combined, '9' = 'CD4 T', '12' = 'CD4 T', '3' = 'CD4 T', '2' = 'CD4 T', '7' = 'CD4 T',
                        '10' = 'Naive CD4 T', '9' = 'Treg', '4' = 'Th17', '14' = 'T', 
                        '15' = 'CD8 T', '5' = 'CD8 T', '11' = 'CD8 T', '0' = 'CD8 T', '8' = 'CD8 T', '17' = 'CD8 T',
                        '24' = 'Proliferating T', '13' = 'B naive', '6' = 'B memory', '25' = 'B proliferating', '18' = 'B proliferating',
                        '1' = 'Plasma', '16' = 'Colonic subepithelial', '26' = 'Myenteric ganglia', '20' = 'CD14 Monocytes',
                        '19' = 'Endothelial', '22' = 'MAST', '21' = 'Enteric glia')
Combined$major_celltype <- Idents(Combined)

#finally, we will subset out T cells, to be merged with T cells subsetted out from whole gut cells from these same samples
#Note that some clusters (T, Proliferating) contain both CD4 and CD8 T cells. These will be further refined in the next steps.
Idents(Combined) <- "major_celltype"
Combined.CD4 <- subset(Combined, idents = c("Naive CD4 T", "CD4 T", "Th17", "T", "Proliferating T"))
Combined.CD8 <- subset(Combined, idents = c("CD8 T", "T", "Proliferating T"))


