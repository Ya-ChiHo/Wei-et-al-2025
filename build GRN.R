library(Seurat)
options(Seurat.object.assay.version = 'v5')
library(Signac)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(stringr)
library(cicero)
library(data.table)
library(Matrix)
library(GenomicRanges)
library(magrittr)
library(SummarizedExperiment)
library(optparse)
library(yaml)
library(Rcpp)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(chromVAR)
library(doParallel)
library(BuenColors)
library(cisTopic)
library(parallel)
library(FNN)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(networkD3)
library(Hmisc)
library(qlcMatrix)
library(reshape2)
library(scales)
library(Matrix)
library(topGO)
library(fgsea)
library(msigdbr)
library(tibble)

#see after T cell data processing scripts 
#these scripts apply cisTopic and FigR packages to generate the Gene Regulatory Networks using ATAC and RNA assays.
#shown are GRNs built for CD4 T cells in PLWH samples (CD4.infected)
#the same steps were repeated for CD4 T cells in HIV- samples (CD4.uninfected), and for CD8 T cells (CD8.infected, CD8.uninfected)

#additional files required:
#CD4.clean.rds
#CD8.clean.rds

#subset out PLWH vs HIV- donor cells####
Idents(CD4.clean) <- "Condition"
CD4.infected <- subset(CD4.clean, idents = c("infected"))
CD4.infected <- NormalizeData(CD4.infected, normalization.method = 'LogNormalize')
CD4.infected <- FindVariableFeatures(CD4.infected, selection.method = "vst", nfeatures = 2000)
CD4.infected <- ScaleData(CD4.infected, verbose = TRUE)

CD4.uninfected <- subset(CD4.clean, idents = c("uninfected"))
CD4.uninfected <- NormalizeData(CD4.uninfected, normalization.method = 'LogNormalize')
CD4.uninfected <- FindVariableFeatures(CD4.uninfected, selection.method = "vst", nfeatures = 2000)
CD4.uninfected <- ScaleData(CD4.uninfected, verbose = TRUE)

Idents(CD8.clean) <- "Condition"
CD8.infected <- subset(CD8.clean, idents = c("infected"))
CD8.infected <- NormalizeData(CD8.infected, normalization.method = 'LogNormalize')
CD8.infected <- FindVariableFeatures(CD8.infected, selection.method = "vst", nfeatures = 2000)
CD8.infected <- ScaleData(CD8.infected, verbose = TRUE)

CD8.uninfected <- subset(CD8.clean, idents = c("uninfected"))
CD8.uninfected <- NormalizeData(CD8.uninfected, normalization.method = 'LogNormalize')
CD8.uninfected <- FindVariableFeatures(CD8.uninfected, selection.method = "vst", nfeatures = 2000)
CD8.uninfected <- ScaleData(CD8.uninfected, verbose = TRUE)

#make cisTopic####
CD4.infected$barcode <- colnames(CD4.infected)
CD4.infected.MACS2ATAC.counts <- CD4.infected[["MACS2ATAC"]]@counts
rownames(CD4.infected.MACS2ATAC.counts) <- str_replace(rownames(CD4.infected.MACS2ATAC.counts), "-", ":")

barcode <- colnames(CD4.infected)
barcode <- as.data.frame(barcode)
barcode$is__cell_barcode <- "1"
barcode$total <- CD4.infected$nCount_MACS2ATAC
barcode$TSS_fragments <- CD4.infected$atac_TSS_fragments

cisTopicObject <- createcisTopicObject(CD4.infected.MACS2ATAC.counts, barcode, project.name='CD4', 
                                       keepCountsMatrix = FALSE,
                                       min.cells = 1, min.regions = 1, is.acc = 1)
cisTopicObject <- runCGSModels(cisTopicObject, topic=c(5, 10, 20, 30, 40), seed=987, 
                               nCores=20, burnin = 120, iterations = 150, addModels=FALSE)
logLikelihoodByIter(cisTopicObject, select=c(5, 10, 20, 30, 40))
cisTopicObject <- selectModel(cisTopicObject, select= 40)
cisTopicObject <- runUmap(cisTopicObject, target='cell')
topic.mat <- modelMatSelection(cisTopicObject, 'cell', 'Probability')
topic.mat <- t(topic.mat)
topic.mat <- as.matrix(topic.mat)

set.seed(123)
cellkNN <- get.knn(topic.mat,k = 30)$nn.index
dim(cellkNN)
rownames(cellkNN) <- barcode$barcode

#make ATAC matrix
CD4.infected.ATAC <- CD4.infected[["MACS2ATAC"]]@counts
CD4.infected.ATAC <- as.data.frame(as.matrix(CD4.infected.ATAC))
#create rows at the end of count data frame "chr", "start", "end"
CD4.infected.ATAC$seqnames <- sapply(strsplit(rownames(CD4.infected.ATAC),"-"), `[`, 1)
CD4.infected.ATAC$start <- sapply(strsplit(rownames(CD4.infected.ATAC),"-"), `[`, 2)
CD4.infected.ATAC$end <- sapply(strsplit(rownames(CD4.infected.ATAC),"-"), `[`, 3)
#remove all non-Chr defined regions
CD4.infected.ATAC <- subset(CD4.infected.ATAC, grepl('chr', rownames(CD4.infected.ATAC)))
#make summarized experiment object
CD4.infected.ATAC.se <- makeSummarizedExperimentFromDataFrame(CD4.infected.ATAC)
counts(CD4.infected.ATAC.se) <- assay(CD4.infected.ATAC.se)
assay(CD4.infected.ATAC.se) <- as(assay(CD4.infected.ATAC.se), 'sparseMatrix')

#make RNA matrix
CD4.infected.RNA <- CD4.infected[["RNA.clean"]]@data
# Remove genes with zero expression across all cells
CD4.infected.RNA <- CD4.infected.RNA[Matrix::rowSums(CD4.infected.RNA)!=0,]
CD4.infected.RNA <- as.matrix(CD4.infected.RNA)

#determine DORCs####
#make peaks to gene associations
cisCorr <- FigR::runGenePeakcorr(ATAC.se = CD4.infected.ATAC.se,
                                 RNAmat = CD4.infected.RNA,
                                 genome = "hg38", 
                                 nCores = 60,
                                 p.cut = NULL, # Set this to NULL and we can filter later
                                 n_bg = 100)
cisCorr.filt <- cisCorr %>% filter(pvalZ <= 0.05)
dorcGenes <- dorcJPlot(dorcTab = cisCorr.filt,
                       cutoff = 2, # No. sig peaks needed to be called a DORC
                       labelTop = 20,
                       returnGeneList = TRUE, # Set this to FALSE for just the plot
                       force=2)
dorcMat <- getDORCScores(ATAC.se = CD4.infected.ATAC.se, # Has to be same SE as used in previous step
                         dorcTab = cisCorr.filt,
                         geneList = dorcGenes,
                         nCores = 4)
dorcMat.s <- smoothScoresNN(NNmat = cellkNN[,1:30], mat = dorcMat,nCores = 4)
RNAmat.s <- smoothScoresNN(NNmat = cellkNN[,1:30], mat = CD4.infected.RNA, nCores = 10)

#filter normalized RNA matrix to only contain DORC genes
RNAmat.filtered.s <- RNAmat.s[rownames(RNAmat.s) %in% rownames(dorcMat.s), ]
RNAmat.filtered.s <- as.data.frame(RNAmat.filtered.s)
dorcMat.filtered.s <- as.data.frame(dorcMat.s)
dorcMat.filtered.s <- dorcMat.filtered.s[match(rownames(RNAmat.filtered.s), rownames(dorcMat.filtered.s)),]

match_RNA_access <- vector(mode = "list", length = nrow(RNAmat.filtered.s))
cor.match_RNA_access <- vector(mode = "list", length = nrow(RNAmat.filtered.s))
for(i in 1:nrow(RNAmat.filtered.s)){
  match_RNA_access[[i]] <- rbind(RNAmat.filtered.s[i,], dorcMat.filtered.s[i,])
  cor.match_RNA_access[[i]] <- corSparse(t(match_RNA_access[[i]]))
  cor.match_RNA_access[[i]] <- cor.match_RNA_access[[i]][1,2]
}

cor.match_RNA_access <- as.matrix(unlist(cor.match_RNA_access))
rownames(cor.match_RNA_access) <- rownames(RNAmat.filtered.s)
cor.match_RNA_access <- as.data.frame(cor.match_RNA_access)
ComplexHeatmap::Heatmap(cor.match_RNA_access, row_names_gp = grid::gpar(fontsize = 5))

#TF-gene association heatmap####
figR.d <- runFigRGRN(ATAC.se = CD4.infected.ATAC.se, # Must be the same input as used in runGenePeakcorr()
                     dorcTab = cisCorr.filt, # Filtered peak-gene associations
                     genome = "hg38",
                     dorcMat = dorcMat.s,
                     rnaMat = RNAmat.s, 
                     nCores = 4)

#make FigR matrix 
motifs <- as.data.frame(unlist(motif.names))
motifs$`unlist(motif.names)` <- base::toupper(motifs$`unlist(motif.names)`)
figR.d1 <- figR.d[figR.d$DORC %in% rownames(cor.match_RNA_access_top), ]
figR.d1 <- figR.d1[figR.d1$Motif %in% motifs$`unlist(motif.names)`, ]

#plot heatmap for all dorc gene - TF associations
plotfigRHeatmap(figR.d = figR.d1,
                score.cut = 1.5,
                column_names_gp = gpar(fontsize=6), # from ComplexHeatmap
                show_row_dend = FALSE # from ComplexHeatmap
)

rankDrivers(figR.d1,rankBy = "meanScore",interactive = FALSE)
rankDrivers(figR.d1,score.cut = 2,rankBy = "nTargets",interactive = FALSE)

figR.d1 %>% 
  ggplot(aes(Corr.log10P,Enrichment.log10P,color=Score)) + 
  ggrastr::geom_point_rast(size=0.01,shape=16) + 
  theme_classic() + 
  scale_color_gradientn(colours = jdb_palette("solar_extra"),
                        limits=c(-3,3),oob = scales::squish,breaks=scales::breaks_pretty(n=3))

#add module score for GRN genes that are activated by BACH2####
p.infected.all.BACH2 <- p.infected.all@matrix
p.infected.all.BACH2 <- as.data.frame(p.infected.all.BACH2)
p.infected.all.BACH2 <- p.infected.all.BACH2[p.infected.all.BACH2$BACH2 > 0, ] 

DefaultAssay(CD4.clean) <- "RNA.clean"
CD4.clean <- AddModuleScore(CD4.clean, features = list(rownames(p.infected.all.BACH2)), nbin = 10,
                            name = 'GRN.bach2.upreg.genes.module')
FeaturePlot(CD4.clean, feature = c("GRN.bach2.upreg.genes.module1"),  reduction = "umap.wnn",
            min.cutoff = 'q20', max.cutoff = 'q80')
