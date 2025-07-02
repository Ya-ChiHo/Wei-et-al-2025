library(Seurat)
options(Seurat.object.assay.version = 'v5')
library(Signac)
library(ggplot2)
library(ggraph)
library(ggrepel)
library(cowplot)
library(patchwork)
library(ComplexHeatmap)
library(dplyr)
library(RColorBrewer)
library(circlize)

#see after T cell data processing scripts 
#these scripts apply diffusion heatmap on differentially expressed genes/TFs in cells ordered by BACH2 TF accessibility in the TRM subsets of CD4 T and CD8 T

#additional files required:
#CD4.clean.rds
#CD8.clean.rds

#subset out TRM cells####
Idents(CD4.clean) <- "major_celltype"
TRM.Th1 <- subset(CD4.clean, ident = c("TRM-Th1-1", "TRM-Th1-2"))
Idents(TRM.Th1) <- "Condition"
TRM.Th1.infected <- subset(TRM.Th1, idents = c("infected"))
TRM.Th1.uninfected <- subset(TRM.Th1, idents = c("uninfected"))

Idents(CD8.clean) <- "major_celltype"
TRM.CD8 <- subset(CD8.clean, ident = c("CD8-TRM"))

Idents(TRM.CD8) <- "Condition"
TRM.CD8.infected <- subset(TRM, idents = c("infected"))
TRM.CD8.uninfected <- subset(TRM, idents = c("uninfected"))


#add BACH2 activity as cell score####
#repeat for cells in TRM.Th1.uninfected, TRM.CD8.infected, TRM.CD8.uninfected

p.infected <- as.data.frame(rownames(TRM.Th1.infected[["chromvar"]]$data))
#check which row is BACH2 (MA1101.2), in my case row 388 
q.infected <- as.data.frame(TRM.Th1.infected[["chromvar"]]$data[388,])
q.infected$BACH2_activity <- q.infected$`TRM.Th1.infected[["chromvar"]]$data[388, ]`
TRM.Th1.infected$BACH2_activity <- q.infected$BACH2_activity

#bin cells ordered by BACH2 into groups####
#repeat for cells in TRM.Th1.uninfected, TRM.CD8.infected, TRM.CD8.uninfected

#bin cells into 20 groups by BACH2 activity
TRM.Th1.infected$BACH2_bin20 <- cut(TRM.Th1.infected$BACH2_activity, 20)

Idents(TRM.Th1.infected) <- "BACH2_bin20"
table(TRM.Th1.infected$BACH2_bin20)

#make the groups somewhat equal in cell size
TRM.Th1.infected$BACH2_bin <- TRM.Th1.infected$BACH2_bin20
Idents(TRM.Th1.infected) <- "BACH2_bin"
TRM.Th1.infected <-RenameIdents(TRM.Th1.infected, '(-8.01,-7.2]' = '(-8.01,-3.23]', '(-7.2,-6.41]' = '(-8.01,-3.23]',
                                '(-6.41,-5.61]' = '(-8.01,-3.23]', '(-5.61,-4.82]' = '(-8.01,-3.23]', 
                                '(-4.82,-4.03]' = '(-8.01,-3.23]', '(-4.03,-3.23]' = '(-8.01,-3.23]',
                                '(-3.23,-2.44]' = '(-3.23,-1.64]', '(-2.44,-1.64]' = '(-3.23,-1.64]',
                                '(-1.64,-0.851]' = '(-1.64,0.737]', '(-0.851,-0.0569]' = '(-1.64,0.737]', '(-0.0569,0.737]' = '(-1.64,0.737]',
                                '(3.12,3.91]' = '(2.32,7.9]', '(2.32,3.12]' = '(2.32,7.9]',
                                '(3.91,4.71]' = '(2.32,7.9]', '(4.71,5.5]' = '(2.32,7.9]', '(5.5,6.29]' = '(2.32,7.9]',
                                '(6.29,7.09]' = '(2.32,7.9]', '(7.09,7.9]' = '(2.32,7.9]')
TRM.Th1.infected$BACH2_bin <- Idents(TRM.Th1.infected)

TRM.Th1.infected.downsample <- subset(x = TRM.Th1.infected, downsample = 800)

#identify differential expressed genes (RNA) for each bin####
#repeat for cells in TRM.Th1.uninfected, TRM.CD8.infected, TRM.CD8.uninfected

DefaultAssay(TRM.Th1.infected.downsample) <- "RNA.clean"; Idents(TRM.Th1.infected.downsample) <- "BACH2_bin"
TRM.Th1.infected.bin.RNA  <- FindAllMarkers(TRM.Th1.infected.downsample, 
                                            only.pos = TRUE, min.pct = 0.1,
                                            logfc.threshold = 0.5)
TRM.Th1.infected.bin.RNA$BH <- p.adjust(TRM.Th1.infected.bin.RNA $p_val, method = "BH", n = nrow(TRM.Th1.infected.bin.RNA ))
TRM.Th1.infected.bin.RNA <- TRM.Th1.infected.bin.RNA[TRM.Th1.infected.bin.RNA$BH < 0.05,]
TRM.Th1.infected.bin.RNA <- TRM.Th1.infected.bin.RNA[!duplicated(TRM.Th1.infected.bin.RNA$gene),]

#plot diffusion heatmap for DEGs (RNA)####
#repeat for cells in TRM.Th1.uninfected, TRM.CD8.infected, TRM.CD8.uninfected

pt.matrix <- as.matrix(TRM.Th1.infected.downsample[["RNA.clean"]]@data
                       [match(rownames(TRM.Th1.infected.bin.RNA), rownames(TRM.Th1.infected.downsample[["RNA.clean"]]@data)),
                         order(TRM.Th1.infected.downsample$BACH2_activity)])

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))

ht <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  #km = 6,
  row_title_rot                = 0,
  cluster_rows                 = FALSE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  use_raster = TRUE)


#heatmap plot
ht

#dynamic gene expression plot by cells ordered by BACH2 activity, e.g., IFNG
plot(pt.matrix[c("IFNG"),], xlab = "cells ordered by BACH2 accessibility", ylab = "IFNG gene expression", 
     main = "IFNG activity expression dynamic change", col = "red") 


#identify differential accessible TFs (ATAC) for each bin####
#repeat for cells in TRM.Th1.uninfected, TRM.CD8.infected, TRM.CD8.uninfected

DefaultAssay(TRM.Th1.infected.downsample) <- "chromvar"; Idents(TRM.Th1.infected.downsample) <- "BACH2_bin"
TRM.Th1.infected.bin.TFs <- FindAllMarkers(TRM.Th1.infected.downsample, only.pos = TRUE, min.pct = 0,
                                           logfc.threshold = 0.3, 
                                           mean.fxn = rowMeans, fc.name = "avg_diff")
TRM.Th1.infected.bin.TFs$BH <- p.adjust(TRM.Th1.infected.bin.TFs$p_val, 
                                        method = "BH", n = nrow(TRM.Th1.infected.bin.TFs))
TRM.Th1.infected.bin.TFs$motifs <- ConvertMotifID(motif.names, id = rownames(TRM.Th1.infected.bin.TFs))
TRM.Th1.infected.bin.TFs <- na.omit(TRM.Th1.infected.bin.TFs)
TRM.Th1.infected.bin.TFs <- TRM.Th1.infected.bin.TFs[!duplicated(TRM.Th1.infected.bin.TFs$motifs),]

#plot diffusion heatmap for differential accessible TFs (ATAC)####
#repeat for cells in TRM.Th1.uninfected, TRM.CD8.infected, TRM.CD8.uninfected

pt.matrix <- as.matrix(TRM.Th1.infected.downsample[["chromvar"]]@data
                       [match(rownames(TRM.Th1.infected.bin.TFs), rownames(TRM.Th1.infected.downsample[["chromvar"]]@data)),
                         order(TRM.Th1.infected.downsample$BACH2_activity)])
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- ConvertMotifID(motif.names, id = rownames(pt.matrix))

ht <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  #km = 6,
  row_title_rot                = 0,
  cluster_rows                 = FALSE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  use_raster = TRUE)

#heatmap plot
ht

#dynamic TF accessibility plot by cells ordered by BACH2 activity, e.g., IRF4
plot(pt.matrix[c("IRF4"),], xlab = "BACH2 accessibility rank", ylab = "IRF4 motif accessibility", 
     main = "IRF4 activity expression dynamic change", col = "red") 


