#Figure 6A
#require objects from CD4 T cell data processing.R 
DimPlot(CD4, reduction = "umap.integrated.atac", group.by = "HIV_RNA_DNA", pt.size = 0.8, cols= c("grey", "red", "blue", "purple"), order = c("TRUE_TRUE", "FALSE_TRUE", "TRUE_FALSE", "FALSE_FALSE"), raster=FALSE) 

#Figure 6D
#require objects from CD4 T cell data processing.R 
CD4.clean$HIV_pos <- CD4.clean$HIV_RNA_DNA
Idents(CD4.clean) <- "HIV_pos"
CD4.clean <-RenameIdents(CD4.clean, "HIV DNA and RNA positive" = "HIV positive", "HIV DNA positive only" = "HIV positive",
                         "HIV RNA positive only" = "HIV positive")
CD4.clean$HIV_pos <- Idents(CD4.clean)
Idents(CD4.infected) <- "HIV_pos"; DefaultAssay(CD4.infected) <- 'chromvar'

motif.names <- CD4.infected@assays$MACS2ATAC@motifs@motif.names

DEM.downsample.HIV_Pos_vs_Neg <- vector(mode = "list", length = 1000)
DEM.downsample.HIV_Pos_vs_Neg_bootstrap <- vector(mode = "list", length = 1000)
for(i in 1:length(DEM.downsample.HIV_Pos_vs_Neg)){
  DEM.downsample.HIV_Pos_vs_Neg[[i]] <- FindMarkers(object = CD4.infected, ident.1 = "HIV positive", ident.2 = "HIV negative", 
                                                    min.pct = 0, logfc.threshold = 0.25, 
                                                    max.cells.per.ident = 99, random.seed = i, only.pos = FALSE, 
                                                    mean.fxn = rowMeans, fc.name = "avg_diff",
                                                    test.use = 'LR')
  DEM.downsample.HIV_Pos_vs_Neg[[i]]$BH <- p.adjust(DEM.downsample.HIV_Pos_vs_Neg[[i]]$p_val, 
                                                                method = "BH", n = nrow(DEM.downsample.HIV_Pos_vs_Neg[[i]]))
  merged_ID <- unique(c(rownames(DEM.downsample.HIV_Pos_vs_Neg[[i]])))
  DEM.downsample.HIV_Pos_vs_Neg[[i]] <- DEM.downsample.HIV_Pos_vs_Neg[[i]][merged_ID,]
  DEM.downsample.HIV_Pos_vs_Neg[[i]] <- DEM.downsample.HIV_Pos_vs_Neg[[i]][order(rownames(DEM.downsample.HIV_Pos_vs_Neg[[i]])),]
  DEM.downsample.HIV_Pos_vs_Neg_bootstrap[[i]] <- data.frame(DEM.downsample.HIV_Pos_vs_Neg[[i]][,1])
  downsample.HIV_RNA_Pos_vs_Neg_bootstrap1 <- t(matrix(unlist(DEM.downsample.HIV_Pos_vs_Neg_bootstrap), nrow = length(merged_ID)))
  Cauchy <- ACAT(matrix(downsample.HIV_RNA_Pos_vs_Neg_bootstrap1, ncol = ncol(downsample.HIV_RNA_Pos_vs_Neg_bootstrap1)))
  matrix <- cbind(DEM.downsample.HIV_Pos_vs_Neg[[1]], Cauchy)
  DEM.HIV_Pos_vs_Neg_bootstrap <- matrix[, -c(1,5,6)]
}
DEM.HIV_Pos_vs_Neg_bootstrap$BH <- p.adjust(DEM.HIV_Pos_vs_Neg_bootstrap$Cauchy, method = "BH", n = nrow(DEM.HIV_Pos_vs_Neg_bootstrap))
DEM.HIV_Pos_vs_Neg_bootstrap$motifname <-rownames(DEM.HIV_Pos_vs_Neg_bootstrap)
DEM.HIV_Pos_vs_Neg_bootstrap$gene <- ConvertMotifID(motif.names, id = rownames(DEM.HIV_Pos_vs_Neg_bootstrap))

#Figure 6E
#require objects from CD4 T cell data processing.R 
Idents(CD4.infected) <- "HIV_pos"; DefaultAssay(CD4.infected) <- 'RNA.clean'

DEG.downsample.HIV_Pos_vs_Neg <- vector(mode = "list", length = 1000)
DEG.downsample.HIV_Pos_vs_Neg_bootstrap <- vector(mode = "list", length = 1000)
for(i in 1:length(DEG.downsample.HIV_Pos_vs_Neg)){
  DEG.downsample.HIV_Pos_vs_Neg[[i]] <- FindMarkers(object = CD4.infected, ident.1 = "HIV positive", ident.2 = "HIV negative", 
                                                        min.pct = 0.3, logfc.threshold = 0.33, max.cells.per.ident = 99, random.seed = i,
                                                        only.pos = FALSE, test.use = 'wilcox')
  DEG.downsample.HIV_Pos_vs_Neg[[i]]$BH <- p.adjust(DEG.downsample.HIV_Pos_vs_Neg[[i]]$p_val, 
                                                        method = "BH", n = nrow(DEG.downsample.HIV_Pos_vs_Neg[[i]]))
  merged_ID <- unique(c(rownames(DEG.downsample.HIV_Pos_vs_Neg[[i]])))
  DEG.downsample.HIV_Pos_vs_Neg[[i]] <- DEG.downsample.HIV_Pos_vs_Neg[[i]][merged_ID,]
  DEG.downsample.HIV_Pos_vs_Neg[[i]] <- DEG.downsample.HIV_Pos_vs_Neg[[i]][order(rownames(DEG.downsample.HIV_Pos_vs_Neg[[i]])),]
  DEG.downsample.HIV_Pos_vs_Neg_bootstrap[[i]] <- data.frame(DEG.downsample.HIV_Pos_vs_Neg[[i]][,1])
  downsample.HIV_DNA_Pos_vs_Neg_bootstrap1 <- t(matrix(unlist(DEG.downsample.HIV_Pos_vs_Neg_bootstrap), nrow = length(merged_ID)))
  Cauchy <- ACAT(matrix(downsample.HIV_DNA_Pos_vs_Neg_bootstrap1, ncol = ncol(downsample.HIV_DNA_Pos_vs_Neg_bootstrap1)))
  matrix <- cbind(DEG.downsample.HIV_Pos_vs_Neg[[1]], Cauchy)
  DEG.HIV_Pos_vs_Neg_bootstrap <- matrix[, -c(1,5,6)]
}
DEG.HIV_Pos_vs_Neg_bootstrap <- DEG.HIV_Pos_vs_Neg_bootstrap[!grepl('IGL', rownames(DEG.HIV_Pos_vs_Neg_bootstrap)),]
DEG.HIV_Pos_vs_Neg_bootstrap <- DEG.HIV_Pos_vs_Neg_bootstrap[!grepl('IGK', rownames(DEG.HIV_Pos_vs_Neg_bootstrap)),]
DEG.HIV_Pos_vs_Neg_bootstrap$BH <- p.adjust(DEG.HIV_Pos_vs_Neg_bootstrap$Cauchy, method = "BH", n = nrow(DEG.HIV_Pos_vs_Neg_bootstrap))
DEG.HIV_Pos_vs_Neg_bootstrap$gene <- rownames(DEG.HIV_Pos_vs_Neg_bootstrap)
rm(downsample.HIV_DNA_Pos_vs_Neg_bootstrap1); rm(Cauchy); rm(i); rm(merged_ID)
rm(DEG.downsample.HIV_Pos_vs_Neg); rm(DEG.downsample.HIV_Pos_vs_Neg_bootstrap); rm(matrix)

#Figure 6F
#require objects from CD4 T cell data processing.R 
#require a list of antibodies in Antibody_stats.CD4.infected with mean expression > mean expression of their specific isotype controls (Z > 2; two-sample Z-test), see Methods in associated manuscript.

DefaultAssay(TRM.infected) <- 'Antibody'; Idents(TRM.infected) <- "HIV_pos"

DEA.antibody.downsample.HIV_Pos_vs_Neg.TRM <- vector(mode = "list", length = 1000)
DEA.antibody.downsample.HIV_Pos_vs_Neg_bootstrap.TRM <- vector(mode = "list", length = 1000)
for(i in 1:length(DEA.antibody.downsample.HIV_Pos_vs_Neg.TRM)){
  DEA.antibody.downsample.HIV_Pos_vs_Neg.TRM[[i]] <- FindMarkers(object = TRM.infected, ident.1 = "HIV positive", ident.2 = "HIV negative", 
                                                             min.pct = 0.15, logfc.threshold = 0.3, only.pos = FALSE, 
                                                             max.cells.per.ident = 80, random.seed = i, features =  rownames(Antibody_stats.CD4.infected))
  DEA.antibody.downsample.HIV_Pos_vs_Neg.TRM[[i]]$BH <- p.adjust(DEA.antibody.downsample.HIV_Pos_vs_Neg.TRM[[i]]$p_val, 
                                                             method = "BH", n = nrow(DEA.antibody.downsample.HIV_Pos_vs_Neg.TRM[[i]]))
  merged_ID <- unique(c(rownames(DEA.antibody.downsample.HIV_Pos_vs_Neg.TRM[[i]])))
  DEA.antibody.downsample.HIV_Pos_vs_Neg.TRM[[i]] <- DEA.antibody.downsample.HIV_Pos_vs_Neg.TRM[[i]][merged_ID,]
  DEA.antibody.downsample.HIV_Pos_vs_Neg.TRM[[i]] <- DEA.antibody.downsample.HIV_Pos_vs_Neg.TRM[[i]][order(rownames(DEA.antibody.downsample.HIV_Pos_vs_Neg.TRM[[i]])),]
  DEA.antibody.downsample.HIV_Pos_vs_Neg_bootstrap.TRM[[i]] <- data.frame(DEA.antibody.downsample.HIV_Pos_vs_Neg.TRM[[i]][,1])
  downsample.HIV_RNA_Pos_vs_Neg_bootstrap1 <- t(matrix(unlist(DEA.antibody.downsample.HIV_Pos_vs_Neg_bootstrap.TRM), nrow = length(merged_ID)))
  Cauchy <- ACAT(matrix(downsample.HIV_RNA_Pos_vs_Neg_bootstrap1, ncol = ncol(downsample.HIV_RNA_Pos_vs_Neg_bootstrap1)))
  matrix <- cbind(DEA.antibody.downsample.HIV_Pos_vs_Neg.TRM[[1]], Cauchy)
  DEA.antibody.HIV_Pos_vs_Neg_TRM.bootstrap <- matrix[, -c(1,5,6)]
}
DEA.antibody.HIV_Pos_vs_Neg_TRM.bootstrap$BH <- p.adjust(DEA.antibody.HIV_Pos_vs_Neg_TRM.bootstrap$Cauchy, method = "BH", n = nrow(DEA.antibody.HIV_Pos_vs_Neg_TRM.bootstrap))
DEA.antibody.HIV_Pos_vs_Neg_TRM.bootstrap$gene <- rownames(DEA.antibody.HIV_Pos_vs_Neg_TRM.bootstrap)

#Figure 6G, H
#require objects from CD4 T cell data processing.R 
CD4.clean$HIV_pos_condition <- str_c(CD4.clean$HIV_pos, '_', CD4.clean$Condition)
Idents(CD4.clean) <- "HIV_pos_condition"
CoveragePlot(CD4.clean, region = "ITGAE",
             ranges.title = "MACS2",extend.upstream = 2000,extend.downstream = 2000,
             ranges = StringToGRanges("chr17-3765415-3765695"))
CoveragePlot(CD4.clean, region = "CCR6",
             ranges.title = "MACS2",extend.upstream = 2000,extend.downstream = 2000)