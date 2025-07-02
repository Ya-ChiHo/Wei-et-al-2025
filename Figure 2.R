#Figure 2A
#require objects from CD8 T cell data processing.R
DimPlot(CD8.clean, reduction = "umap.wnn", group.by = "major_celltype", label = TRUE, label.size = 4, raster=FALSE)

#Figure 2C
#require objects from CD8 T cell data processing.R and build GRN.R
rankDrivers(figR.d1, score.cut = 2,rankBy = "nTargets",interactive = FALSE)

#Figure 2D
#require objects from CD8 T cell data processing.R 
FeaturePlot(CD8.clean, feature = c("MA1101.2"),  reduction = "umap.wnn",
            min.cutoff = 'q20', max.cutoff = 'q80')

#Figure 2E
#require objects from CD8 T cell data processing.R and build GRN.R
FeaturePlot(CD8.clean, feature = c("GRN.upreg.bach2.genes.module1"),  reduction = "umap.wnn",
            min.cutoff = 'q20', max.cutoff = 'q80')


#Figure 2F
#require objects from CD8 T cell data processing.R 
Idents(CD8.clean) <- "major_celltype"
TRM <- subset(CD8.clean, idents = c("TRM_1", "TRM_2", "TEM"))
Idents(TRM) <- "Condition"
TRM.infected <- subset(TRM, idents = c("infected"))
TRM.uninfected <- subset(TRM, idents = c("uninfected"))

DefaultAssay(TRM.infected) <- "chromvar"; Idents(TRM.infected) <- "memory_celltype"
mylevels<- c("TRM_1", "TRM_2", "TEM")
Idents(TRM.infected) <- factor(Idents(TRM.infected), levels = mylevels)
DEM.TRM.infected <- FindAllMarkers(TRM.infected, only.pos = TRUE, min.pct = 0,
                                   logfc.threshold = 0.1, test.use = "wilcox_limma",
                                   mean.fxn = rowMeans, fc.name = "avg_diff")
DEM.TRM.infected$BH <- p.adjust(DEM.TRM.infected$p_val, method = "BH", n = nrow(DEM.TRM.infected))
DEM.TRM.infected$motifs <- ConvertMotifID(motif.names, id = rownames(DEM.TRM.infected))
DEM.TRM.infected <- na.omit(DEM.TRM.infected)
DEM.TRM.infected <- DEM.TRM.infected[!duplicated(DEM.TRM.infected$motifs),]
DEM.TRM.infected <- DEM.TRM.infected[DEM.TRM.infected$BH < 0.05,]
filtered_DEM.TRM.infected <- subset(DEM.TRM.infected, grepl("BACH|JUN|FOS|BATF|MAF|TBX|EOMES|FOXO|RUNX|ROR|NFKB|REL|IRF|STAT|ELK|KLF|NR4A|RUNX|NFAT", motifs))
filtered_DEM.TRM.infected <- filtered_DEM.TRM.infected[order(filtered_DEM.TRM.infected$cluster, filtered_DEM.TRM.infected$motifs), ]
filtered_DEM.TRM.infected <- filtered_DEM.TRM.infected[!filtered_DEM.TRM.infected$motifs %in%
                                                         c("HOXD12::ELK1", "HOXB2::ELK1", "ETV5::FOXO1"), ]
Idents(TRM.infected) <- "memory_celltype"; DefaultAssay(TRM.infected) <- "chromvar"
TRM.infected@assays$chromvar$counts <- TRM.infected@assays$chromvar$data
AVG.TRM.infected.TF <- AverageExpression(TRM.infected, return.seurat = T, group.by = 'memory_celltype')
AVG.TRM.infected.TF$orig.ident <- colnames(AVG.TRM.infected.TF)
DefaultAssay(AVG.TRM.infected.TF) <- "chromvar"; Idents(AVG.TRM.infected.TF) <- "orig.ident"
mylevels<- c("TRM-1", "TRM-2", "TEM")
Idents(AVG.TRM.infected.TF) <- factor(Idents(AVG.TRM.infected.TF), levels = mylevels)

DefaultAssay(TRM.uninfected) <- "chromvar"; Idents(TRM.uninfected) <- "memory_celltype"
mylevels<- c("TRM_1", "TRM_2", "TEM")
Idents(TRM.uninfected) <- factor(Idents(TRM.uninfected), levels = mylevels)
DEM.TRM.uninfected <- FindAllMarkers(TRM.uninfected, only.pos = TRUE, min.pct = 0,
                                     logfc.threshold = 0.1, test.use = "wilcox_limma",
                                     mean.fxn = rowMeans, fc.name = "avg_diff")
DEM.TRM.uninfected$BH <- p.adjust(DEM.TRM.uninfected$p_val, method = "BH", n = nrow(DEM.TRM.uninfected))
DEM.TRM.uninfected$motifs <- ConvertMotifID(motif.names, id = rownames(DEM.TRM.uninfected))
DEM.TRM.uninfected <- na.omit(DEM.TRM.uninfected)
DEM.TRM.uninfected <- DEM.TRM.uninfected[!duplicated(DEM.TRM.uninfected$motifs),]
DEM.TRM.uninfected <- DEM.TRM.uninfected[DEM.TRM.uninfected$BH < 0.05,]
filtered_DEM.TRM.uninfected <- subset(DEM.TRM.uninfected, grepl("BACH|JUN|FOS|BATF|MAF|TBX|EOMES|FOXO|RUNX|ROR|NFKB|REL|IRF|STAT|ELK|KLF|NR4A|RUNX|NFAT", motifs))
filtered_DEM.TRM.uninfected <- filtered_DEM.TRM.uninfected[order(filtered_DEM.TRM.uninfected$cluster, filtered_DEM.TRM.uninfected$motifs), ]
filtered_DEM.TRM.uninfected <- filtered_DEM.TRM.uninfected[!filtered_DEM.TRM.uninfected$motifs %in%
                                                             c("HOXD12::ELK1", "HOXB2::ELK1", "ETV5::FOXO1"), ]
Idents(TRM.uninfected) <- "memory_celltype"; DefaultAssay(TRM.uninfected) <- "chromvar"
TRM.uninfected@assays$chromvar$counts <- TRM.uninfected@assays$chromvar$data
AVG.TRM.uninfected.TF <- AverageExpression(TRM.uninfected, return.seurat = T, group.by = 'memory_celltype')
AVG.TRM.uninfected.TF$orig.ident <- colnames(AVG.TRM.uninfected.TF)
DefaultAssay(AVG.TRM.uninfected.TF) <- "chromvar"; Idents(AVG.TRM.uninfected.TF) <- "orig.ident"
mylevels<- c("TRM-1", "TRM-2", "TEM")
Idents(AVG.TRM.uninfected.TF) <- factor(Idents(AVG.TRM.uninfected.TF), levels = mylevels)

DoHeatmap(object = AVG.TRM.infected.TF, 
          features = filtered_DEM.TRM.infected_gene, slot = "scale.data", draw.lines = FALSE,
          disp.min = -1, disp.max = 1, assay = 'chromvar') +  scale_y_discrete(
            labels = rev(ConvertMotifID(motif.names, id = filtered_DEM.TRM.infected_gene))) + scale_fill_gradientn(colors = c("grey", "orange", "red"))
DoHeatmap(object = AVG.TRM.uninfected.TF, 
          features = filtered_DEM.TRM.uninfected_gene, slot = "scale.data", draw.lines = FALSE,
          disp.min = -1, disp.max = 1, assay = 'chromvar') +  scale_y_discrete(
            labels = rev(ConvertMotifID(motif.names, id = filtered_DEM.TRM.uninfected_gene))) + scale_fill_gradientn(colors = c("grey", "orange", "red"))


#Figure 2G
#require objects from CD8 T cell data processing.R 
DefaultAssay(TRM.infected) <- "RNA.clean"; Idents(TRM.infected) <- "memory_celltype"
mylevels<- c("TRM_1", "TRM_2", "TEM")
Idents(TRM.infected) <- factor(Idents(TRM.infected), levels = mylevels)
DEG.TRM.infected  <- FindAllMarkers(TRM.infected, only.pos = TRUE, min.pct = 0.1, test.use = "wilcox_limma",
                                    logfc.threshold = 0.25)
DEG.TRM.infected$BH <- p.adjust(DEG.TRM.infected $p_val, method = "BH", n = nrow(DEG.TRM.infected ))
DEG.TRM.infected <- na.omit(DEG.TRM.infected)
DEG.TRM.infected <- DEG.TRM.infected[!duplicated(DEG.TRM.infected$gene),]
DEG.TRM.infected <- DEG.TRM.infected[DEG.TRM.infected$BH < 0.05,]
filtered_DEG.TRM.infected <- subset(DEG.TRM.infected, 
                                    grepl("BACH|JUN|FOS|CD44|TCF7|CCR|CXCR|CXCL|CCL|CD69|ICOS|TNF|IFN|IL|CTLA|ENTPD|LAG3|TOX|ITGA|KLR|
                                          ISG|IFI|LCK|TIGIT|CCL|ZEB|RORA|TGF|FAS|LCP2|LAT|FYN|HIF1A|CASP|ATG|LEF1|NFAT|TGF|ID2|AHR|
                                          CD27|CD28|ZAP70|CD40LG|IKZF|BIRC|GZM|MAF", gene))
filtered_DEG.TRM.infected <- filtered_DEG.TRM.infected[order(filtered_DEG.TRM.infected$cluster, filtered_DEG.TRM.infected$gene), ]
filtered_DEG.TRM.infected <- filtered_DEG.TRM.infected[!filtered_DEG.TRM.infected$gene %in%
                                                         c("MALAT1", "LATS1", "LCLAT1", "RORA-AS1", "ENTPD1−AS1", "CD44-AS1", "ENTPD1-AS1", "FASTKD2", "IFNG−AS1"), ]
Idents(TRM.infected) <- "memory_celltype"; DefaultAssay(TRM.infected) <- "RNA.clean"
AVG.TRM.infected.RNA <- AverageExpression(TRM.infected, return.seurat = T, group.by = 'memory_celltype')
AVG.TRM.infected.RNA$orig.ident <- colnames(AVG.TRM.infected.RNA)
DefaultAssay(AVG.TRM.infected.RNA) <- "RNA.clean"; Idents(AVG.TRM.infected.RNA) <- "orig.ident"
mylevels<- c("TRM-1", "TRM-2", "TEM")
Idents(AVG.TRM.infected.RNA) <- factor(Idents(AVG.TRM.infected.RNA), levels = mylevels)

DefaultAssay(TRM.uninfected) <- "RNA.clean"; Idents(TRM.uninfected) <- "memory_celltype"
mylevels<- c("TRM_1", "TRM_2", "TEM")
Idents(TRM.uninfected) <- factor(Idents(TRM.uninfected), levels = mylevels)
DEG.TRM.uninfected  <- FindAllMarkers(TRM.uninfected, only.pos = TRUE, min.pct = 0.1, test.use = "wilcox_limma",
                                      logfc.threshold = 0.25)
DEG.TRM.uninfected$BH <- p.adjust(DEG.TRM.uninfected $p_val, method = "BH", n = nrow(DEG.TRM.uninfected ))
DEG.TRM.uninfected <- na.omit(DEG.TRM.uninfected)
DEG.TRM.uninfected <- DEG.TRM.uninfected[!duplicated(DEG.TRM.uninfected$gene),]
DEG.TRM.uninfected <- DEG.TRM.uninfected[DEG.TRM.uninfected$BH < 0.05,]
filtered_DEG.TRM.uninfected <- subset(DEG.TRM.uninfected, 
                                      grepl("BACH|JUN|FOS|CD44|TCF7|CCR|CXCR|CXCL|CCL|CD69|ICOS|TNF|IFN|IL|CTLA|ENTPD|LAG3|TOX|ITGA|KLR|
                                          ISG|IFI|LCK|TIGIT|CCL|ZEB|RORA|TGF|FAS|LCP2|LAT|FYN|HIF1A|CASP|ATG|LEF1|NFAT|TGF|ID2|AHR|
                                          CD27|CD28|ZAP70|CD40LG|IKZF|BIRC|GZM|MAF", gene))
filtered_DEG.TRM.uninfected <- filtered_DEG.TRM.uninfected[order(filtered_DEG.TRM.uninfected$cluster, filtered_DEG.TRM.uninfected$gene), ]
filtered_DEG.TRM.uninfected <- filtered_DEG.TRM.uninfected[!filtered_DEG.TRM.uninfected$gene %in%
                                                             c("MALAT1", "LATS1", "LCLAT1", "RORA-AS1", "ENTPD1−AS1", "CD44-AS1", "ENTPD1-AS1", "FASTKD2", "IFNG−AS1"), ]
Idents(TRM.uninfected) <- "memory_celltype"; DefaultAssay(TRM.uninfected) <- "RNA.clean"
AVG.TRM.uninfected.RNA <- AverageExpression(TRM.uninfected, return.seurat = T, group.by = 'memory_celltype')
AVG.TRM.uninfected.RNA$orig.ident <- colnames(AVG.TRM.uninfected.RNA)
DefaultAssay(AVG.TRM.uninfected.RNA) <- "RNA.clean"; Idents(AVG.TRM.uninfected.RNA) <- "orig.ident"
mylevels<- c("TRM-1", "TRM-2", "TEM")
Idents(AVG.TRM.uninfected.RNA) <- factor(Idents(AVG.TRM.uninfected.RNA), levels = mylevels)

DoHeatmap(object = AVG.TRM.infected.RNA, 
          features = filtered_DEG.TRM.infected_gene, slot = "scale.data", draw.lines = FALSE,
          disp.min = -1, disp.max = 1, assay = 'RNA.clean') + scale_fill_gradientn(colors = c("grey", "orange", "red"))
DoHeatmap(object = AVG.TRM.uninfected.RNA, 
          features = filtered_DEG.TRM.uninfected$gene, slot = "scale.data", draw.lines = FALSE,
          disp.min = -1, disp.max = 1, assay = 'RNA.clean') + scale_fill_gradientn(colors = c("grey", "orange", "red"))


#Figure 2H
#require objects from CD8 T cell data processing.R 
#require a list of antibodies in Antibody_stats.TEM.infected and Antibody_stats.TEM.uninfected with mean expression > mean expression of their specific isotype controls (Z > 2; two-sample Z-test), see Methods in associated manuscript.
DefaultAssay(TRM.infected) <- "Antibody"
Idents(TRM.infected) <- "memory_celltype"
mylevels<- c("TRM_1", "TRM_2", "TEM")
Idents(TRM.infected) <- factor(Idents(TRM.infected), levels = mylevels)
DEA.TRM.infected <- FindAllMarkers(TRM.infected, min.pct = 0.5, logfc.threshold = 0.1, test.use = "wilcox_limma",
                                   only.pos = TRUE, pseudocount.use = 0.000001, features = rownames(Antibody_stats.TRM.infected))
DEA.TRM.infected$BH <- p.adjust(DEA.TRM.infected$p_val, method = "BH", n = nrow(DEA.TRM.infected))
DEA.TRM.infected$gene <- rownames(DEA.TRM.infected)
DEA.TRM.infected <- DEA.TRM.infected[DEA.TRM.infected$BH < 0.05,]
DEA.TRM.infected <- DEA.TRM.infected[order(DEA.TRM.infected$cluster, DEA.TRM.infected$gene), ]

Idents(TRM.infected) <- "memory_celltype"; DefaultAssay(TRM.infected) <- "Antibody"
AVG.TRM.infected.Antibody <- AverageExpression(TRM.infected, return.seurat = T, group.by = 'memory_celltype')
AVG.TRM.infected.Antibody$orig.ident <- colnames(AVG.TRM.infected.Antibody)
DefaultAssay(AVG.TRM.infected.Antibody) <- "RNA.clean"; Idents(AVG.TRM.infected.Antibody) <- "orig.ident"
mylevels<- c("TRM-1", "TRM-2", "TEM")
Idents(AVG.TRM.infected.Antibody) <- factor(Idents(AVG.TRM.infected.Antibody), levels = mylevels)


DefaultAssay(TRM.uninfected) <- "Antibody"
Idents(TRM.uninfected) <- "memory_celltype"
mylevels<- c("TRM_1", "TRM_2", "TEM")
Idents(TRM.uninfected) <- factor(Idents(TRM.uninfected), levels = mylevels)
DEA.TRM.uninfected <- FindAllMarkers(TRM.uninfected, min.pct = 0.5, logfc.threshold = 0.1, test.use = "wilcox_limma",
                                     only.pos = TRUE, pseudocount.use = 0.000001, features = rownames(Antibody_stats.TRM.uninfected))
DEA.TRM.uninfected$BH <- p.adjust(DEA.TRM.uninfected$p_val, method = "BH", n = nrow(DEA.TRM.uninfected))
DEA.TRM.uninfected$gene <- rownames(DEA.TRM.uninfected)
DEA.TRM.uninfected <- DEA.TRM.uninfected[DEA.TRM.uninfected$BH < 0.05,]
DEA.TRM.uninfected <- DEA.TRM.uninfected[order(DEA.TRM.uninfected$cluster, DEA.TRM.uninfected$gene), ]

Idents(TRM.uninfected) <- "memory_celltype"; DefaultAssay(TRM.uninfected) <- "Antibody"
AVG.TRM.uninfected.Antibody <- AverageExpression(TRM.uninfected, return.seurat = T, group.by = 'memory_celltype')
AVG.TRM.uninfected.Antibody$orig.ident <- colnames(AVG.TRM.uninfected.Antibody)
DefaultAssay(AVG.TRM.uninfected.Antibody) <- "Antibody"; Idents(AVG.TRM.uninfected.Antibody) <- "orig.ident"
mylevels<- c("TRM-1", "TRM-2", "TEM")
Idents(AVG.TRM.uninfected.Antibody) <- factor(Idents(AVG.TRM.uninfected.Antibody), levels = mylevels)

DoHeatmap(object = AVG.TRM.infected.Antibody, 
          features = DEA.TRM.infected_gene, slot = "scale.data", draw.lines = FALSE,
          disp.min = -1, disp.max = 1, assay = 'Antibody') + scale_fill_gradientn(colors = c("grey", "orange", "red"))
DoHeatmap(object = AVG.TRM.uninfected.Antibody, 
          features = DEA.TRM.uninfected_gene, slot = "scale.data", draw.lines = FALSE,
          disp.min = -1, disp.max = 1, assay = 'Antibody') + scale_fill_gradientn(colors = c("grey", "orange", "red"))

#Figure 2I
#require objects from CD8 T cell data processing.R and build GSEA.R
Idents(CD8.clean) <- "memory_celltype"; DefaultAssay(CD8.clean) <- "RNA.clean"

DEG.CD8.TRM_1_3.GSEA <- wilcoxauc(CD8.clean, group_by = "memory_celltype", seurat_assay = "RNA.clean")
DEG.CD8.TRM_1_3.GSEA <- subset(DEG.CD8.TRM_1_3.GSEA, DEG.CD8.TRM_1_3.GSEA$group == "TRM_1")
Rank.DEG.CD8.TRM_1_3.GSEA <- DEG.CD8.TRM_1_3.GSEA$logFC
names(Rank.DEG.CD8.TRM_1_3.GSEA) <- DEG.CD8.TRM_1_3.GSEA$feature
fgsea_DEG.CD8.TRM_1_3.GSEA <- fgsea(pathways = fgsea_sets, 
                                    stats = Rank.DEG.CD8.TRM_1_3.GSEA,
                                    eps   = 0.0, minSize=15, maxSize=500, scoreType = "pos")

DEG.CD8.TRM_3_1.GSEA <- wilcoxauc(CD8.clean, group_by = "memory_celltype", seurat_assay = "RNA.clean")
DEG.CD8.TRM_3_1.GSEA <- subset(DEG.CD8.TRM_3_1.GSEA, DEG.CD8.TRM_3_1.GSEA$group == "TRM_3")
Rank.DEG.CD8.TRM_3_1.GSEA <- DEG.CD8.TRM_3_1.GSEA$logFC
names(Rank.DEG.CD8.TRM_3_1.GSEA) <- DEG.CD8.TRM_3_1.GSEA$feature
fgsea_DEG.CD8.TRM_3_1.GSEA <- fgsea(pathways = fgsea_sets, 
                                    stats = Rank.DEG.CD8.TRM_3_1.GSEA,
                                    eps   = 0.0, minSize=15, maxSize=500, scoreType = "pos")


DEG.CD8.TEM.GSEA <- wilcoxauc(CD8.clean, group_by = "memory_celltype", seurat_assay = "RNA.clean")
DEG.CD8.TEM.GSEA <- subset(DEG.CD8.TEM.GSEA, DEG.CD8.TEM.GSEA$group == "TEM")
Rank.DEG.CD8.TEM.GSEA <- DEG.CD8.TEM.GSEA$logFC
names(Rank.DEG.CD8.TEM.GSEA) <- DEG.CD8.TEM.GSEA$feature
fgsea_DEG.CD8.TEM.GSEA <- fgsea(pathways = fgsea_sets, 
                                    stats = Rank.DEG.CD8.TEM.GSEA,
                                    eps   = 0.0, minSize=15, maxSize=500, scoreType = "pos")

a1 <- plotEnrichment(fgsea_sets[["GOBP_NEGATIVE_REGULATION_OF_IMMUNE_EFFECTOR_PROCESS"]], 
                     Rank.DEG.CD8.TRM_1_3.GSEA) 
b1 <- plotEnrichment(fgsea_sets[["GOBP_NEGATIVE_REGULATION_OF_IMMUNE_EFFECTOR_PROCESS"]], 
                     Rank.DEG.CD8.TRM_3_1.GSEA) 
c1 <- plotEnrichment(fgsea_sets[["GOBP_NEGATIVE_REGULATION_OF_IMMUNE_EFFECTOR_PROCESS"]], 
                     Rank.DEG.CD8.TEM.GSEA)

ggplot(a1$data, aes(x, y)) + geom_line(color = "red", size = 2) + 
  geom_line(data = b1$data, color = "blue", size = 2) +
  geom_line(data = c1$data, color = "green", size = 2) + labs(title = "GOBP_NEGATIVE_REGULATION_OF_IMMUNE_EFFECTOR_PROCESS") + xlab("rank") + ylab("enrichment score") + geom_hline(
    yintercept=0) + geom_vline(xintercept=0) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                     panel.background = element_blank())

#Figure 2J
#require objects from CD8 T cell data processing.R and build GSEA.R
a1 <- plotEnrichment(fgsea_sets[["PID_AP1_PATHWAY"]], 
                     Rank.DEG.CD4.TEM_1_2.GSEA) 
b1 <- plotEnrichment(fgsea_sets[["PID_AP1_PATHWAY"]], 
                     Rank.DEG.CD4.TEM_2_1.GSEA) 
c1 <- plotEnrichment(fgsea_sets[["PID_AP1_PATHWAY"]], 
                     Rank.DEG.CD4.TEM.GSEA) 

ggplot(a1$data, aes(x, y)) + geom_line(color = "red", size = 2) + 
  geom_line(data = b1$data, color = "blue", size = 2) +
  geom_line(data = c1$data, color = "green", size = 2) + labs(title = "PID_AP1_PATHWAY") + xlab("rank") + ylab("enrichment score") + geom_hline(
    yintercept=0) + geom_vline(xintercept=0) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                     panel.background = element_blank())
