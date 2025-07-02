#Figure 1A
#require objects from whole gut data preprocessing.R
DimPlot(CD4.clean, reduction = "umap.wnn", group.by = "major_celltype", label = TRUE, label.size = 4, raster=FALSE)

#Figure 1C
#require objects from CD4 T cell data processing.R
DimPlot(CD4.clean, reduction = "umap.wnn", group.by = "major_celltype", label = TRUE, label.size = 8)

#Figure 1E
#require objects from CD4 T cell data processing.R and build GRN.R
rankDrivers(figR.d1,score.cut = 2,rankBy = "nTargets",interactive = FALSE)

#Figure 1F
#require objects from CD4 T cell data processing.R 
FeaturePlot(CD4.clean, feature = c("MA1101.2"),  reduction = "umap.wnn",
            min.cutoff = 'q20', max.cutoff = 'q80')

#Figure 1G
#require objects from CD4 T cell data processing.R and build GRN.R
FeaturePlot(CD4.clean, feature = c("GRN.bach2.upreg.genes.module1"),  reduction = "umap.wnn",
            min.cutoff = 'q20', max.cutoff = 'q80')

#Figure 1H
#require objects from CD4 T cell data processing.R 
Idents(CD4.clean) <- "memory_celltype"
TEM <- subset(CD4.clean, idents = c("TRM_1", "TRM_2", "TEM"))
Idents(TEM) <- "Condition"
TEM.infected <- subset(TEM, idents = c("infected"))
TEM.uninfected <- subset(TEM, idents = c("uninfected"))

DefaultAssay(TEM.infected) <- "chromvar"; Idents(TEM.infected) <- "memory_celltype"
mylevels<- c("TRM_1", "TRM_2", "TEM")
Idents(TEM.infected) <- factor(Idents(TEM.infected), levels = rev(mylevels))

DEM.TEM.infected <- FindAllMarkers(TEM.infected, only.pos = TRUE, min.pct = 0,
                                logfc.threshold = 0.1, test.use = "wilcox_limma",
                                mean.fxn = rowMeans, fc.name = "avg_diff")
DEM.TEM.infected$BH <- p.adjust(DEM.TEM.infected$p_val, method = "BH", n = nrow(DEM.TEM.infected))
DEM.TEM.infected$motifs <- ConvertMotifID(motif.names, id = rownames(DEM.TEM.infected))
DEM.TEM.infected <- na.omit(DEM.TEM.infected)
DEM.TEM.infected <- DEM.TEM.infected[!duplicated(DEM.TEM.infected$motifs),]
DEM.TEM.infected <- DEM.TEM.infected[DEM.TEM.infected$BH < 0.05,]

filtered_DEM.TEM.infected <- subset(DEM.TEM.infected, grepl("BACH|JUN|FOS|BATF|MAF|TBX|EOMES|FOXO|RUNX|ROR|NFKB|REL|IRF|STAT|ELK|KLF|NR4A", motifs))
custom_order <- c("BACH", "JUN", "FOS", "BATF", "MAF", "NFKB", "FOXO", "RUNX", "TBX", "EOMES", "ROR",
                  "KLF", "ELK", "NR4A", "IRF", "REL")
filtered_DEM.TEM.infected <- filtered_DEM.TEM.infected %>%
  mutate(motif_group = case_when(
    grepl("^BACH", motifs) ~ "BACH",
    grepl("^JUN", motifs) ~ "JUN",
    grepl("^FOS", motifs) ~ "FOS",
    grepl("^BATF", motifs) ~ "BATF",
    grepl("^MAF", motifs) ~ "MAF",
    grepl("^NFKB", motifs) ~ "NFKB",
    grepl("^FOXO", motifs) ~ "FOXO",
    grepl("^RUNX", motifs) ~ "RUNX",
    grepl("^TBX", motifs) ~ "TBX",
    grepl("^EOMES", motifs) ~ "EOMES",
    grepl("^ROR", motifs) ~ "ROR",
    grepl("^KLF", motifs) ~ "KLF",
    grepl("^ELK", motifs) ~ "ELK",
    grepl("^NR4A", motifs) ~ "NR4A",
    grepl("^IRF", motifs) ~ "IRF",
    grepl("^REL", motifs) ~ "REL")) %>%
  mutate(motif_group = factor(motif_group, levels = custom_order, ordered = TRUE)) %>%
  arrange(motif_group, motifs)
filtered_DEM.TEM.infected <- filtered_DEM.TEM.infected[order(-xtfrm(filtered_DEM.TEM.infected$cluster), 
                                                             filtered_DEM.TEM.infected$motifs), ]
Idents(TEM.infected) <- "memory_celltype"; DefaultAssay(TEM.infected) <- "chromvar"
TEM.infected@assays$chromvar$counts <- TEM.infected@assays$chromvar$data
AVG.TEM.infected.TF <- AverageExpression(TEM.infected, return.seurat = T, group.by = 'memory_celltype')
AVG.TEM.infected.TF$orig.ident <- colnames(AVG.TEM.infected.TF)
DefaultAssay(AVG.TEM.infected.TF) <- "chromvar"; Idents(AVG.TEM.infected.TF) <- "orig.ident"

DefaultAssay(TEM.uninfected) <- "chromvar"; Idents(TEM.uninfected) <- "memory_celltype"
mylevels<- c("TRM_1", "TRM_2", "TEM")
Idents(TEM.uninfected) <- factor(Idents(TEM.uninfected), levels = rev(mylevels))

DEM.TEM.uninfected <- FindAllMarkers(TEM.uninfected, only.pos = TRUE, min.pct = 0,
                                   logfc.threshold = 0.1, test.use = "wilcox_limma",
                                   mean.fxn = rowMeans, fc.name = "avg_diff")
DEM.TEM.uninfected$BH <- p.adjust(DEM.TEM.uninfected$p_val, method = "BH", n = nrow(DEM.TEM.uninfected))
DEM.TEM.uninfected$motifs <- ConvertMotifID(motif.names, id = rownames(DEM.TEM.uninfected))
DEM.TEM.uninfected <- na.omit(DEM.TEM.uninfected)
DEM.TEM.uninfected <- DEM.TEM.uninfected[!duplicated(DEM.TEM.uninfected$motifs),]
DEM.TEM.uninfected <- DEM.TEM.uninfected[DEM.TEM.uninfected$BH < 0.05,]

filtered_DEM.TEM.uninfected <- subset(DEM.TEM.uninfected, grepl("BACH|JUN|FOS|BATF|MAF|TBX|EOMES|FOXO|RUNX|ROR|NFKB|REL|IRF|STAT|ELK|KLF|NR4A", motifs))
custom_order <- c("BACH", "JUN", "FOS", "BATF", "MAF", "NFKB", "FOXO", "RUNX", "TBX", "EOMES", "ROR",
                  "KLF", "ELK", "NR4A", "IRF", "REL")

filtered_DEM.TEM.uninfected <- filtered_DEM.TEM.uninfected %>%
  mutate(motif_group = case_when(
    grepl("^BACH", motifs) ~ "BACH",
    grepl("^JUN", motifs) ~ "JUN",
    grepl("^FOS", motifs) ~ "FOS",
    grepl("^BATF", motifs) ~ "BATF",
    grepl("^MAF", motifs) ~ "MAF",
    grepl("^NFKB", motifs) ~ "NFKB",
    grepl("^FOXO", motifs) ~ "FOXO",
    grepl("^RUNX", motifs) ~ "RUNX",
    grepl("^TBX", motifs) ~ "TBX",
    grepl("^EOMES", motifs) ~ "EOMES",
    grepl("^ROR", motifs) ~ "ROR",
    grepl("^KLF", motifs) ~ "KLF",
    grepl("^ELK", motifs) ~ "ELK",
    grepl("^NR4A", motifs) ~ "NR4A",
    grepl("^IRF", motifs) ~ "IRF",
    grepl("^REL", motifs) ~ "REL")) %>%
  mutate(motif_group = factor(motif_group, levels = custom_order, ordered = TRUE)) %>%
  arrange(motif_group, motifs)

filtered_DEM.TEM.uninfected <- filtered_DEM.TEM.uninfected[order(-xtfrm(filtered_DEM.TEM.uninfected$cluster), 
                                                                 filtered_DEM.TEM.uninfected$motifs), ]
Idents(TEM.uninfected) <- "memory_celltype"; DefaultAssay(TEM.uninfected) <- "chromvar"
TEM.uninfected@assays$chromvar$counts <- TEM.uninfected@assays$chromvar$data
AVG.TEM.uninfected.TF <- AverageExpression(TEM.uninfected, return.seurat = T, group.by = 'memory_celltype')
AVG.TEM.uninfected.TF$orig.ident <- colnames(AVG.TEM.uninfected.TF)
DefaultAssay(AVG.TEM.uninfected.TF) <- "chromvar"; Idents(AVG.TEM.uninfected.TF) <- "orig.ident"

DoHeatmap(object = AVG.TEM.infected.TF, 
          features = filtered_DEM.TEM.infected_gene, slot = "scale.data", draw.lines = FALSE,
          disp.min = -1, disp.max = 1, assay = 'chromvar') +  scale_y_discrete(
            labels = rev(ConvertMotifID(motif.names, id = filtered_DEM.TEM.infected_gene))) + scale_fill_gradientn(colors = c("grey", "orange", "red"))
DoHeatmap(object = AVG.TEM.uninfected.TF, 
          features = filtered_DEM.TEM.uninfected_gene, slot = "scale.data", draw.lines = FALSE,
          disp.min = -1, disp.max = 1, assay = 'chromvar') +  scale_y_discrete(
            labels = rev(ConvertMotifID(motif.names, id = filtered_DEM.TEM.uninfected_gene))) + scale_fill_gradientn(colors = c("grey", "orange", "red"))

#Figure 1I
#require objects from CD4 T cell data processing.R 
DefaultAssay(TEM.infected) <- "RNA.clean"; Idents(TEM.infected) <- "memory_celltype"
mylevels<- c("TRM_1", "TRM_2", "TEM")
Idents(TEM.infected) <- factor(Idents(TEM.infected), levels = rev(mylevels))

DEG.TEM.infected  <- FindAllMarkers(TEM.infected, only.pos = TRUE, min.pct = 0.1, test.use = "wilcox_limma",
                                 logfc.threshold = 0.25)
DEG.TEM.infected$BH <- p.adjust(DEG.TEM.infected $p_val, method = "BH", n = nrow(DEG.TEM.infected ))
DEG.TEM.infected <- na.omit(DEG.TEM.infected)
DEG.TEM.infected <- DEG.TEM.infected[!duplicated(DEG.TEM.infected$gene),]
DEG.TEM.infected <- DEG.TEM.infected[DEG.TEM.infected$BH < 0.05,]
filtered_DEG.TEM.infected <- subset(DEG.TEM.infected, 
                                    grepl("BACH|JUN|FOS|CD44|TCF7|CCR|CXCR|CXCL|CCL|CD69|ICOS|TNF|IFN|IL|CTLA|ENTPD|LAG3|TOX|ITGA|KLRB1|
                                          ISG|IFI|LCK|TIGIT|CCL|ZEB|RORA|TGF|FAS|LCP2|LAT|FYN|HIF1A|CASP|ATG|LEF1|
                                          CD27|CD28|ZAP70|CD40LG|IKZF|BIRC", gene))
filtered_DEG.TEM.infected <- filtered_DEG.TEM.infected[order(filtered_DEG.TEM.infected$cluster, filtered_DEG.TEM.infected$gene), ]
filtered_DEG.TEM.infected <- filtered_DEG.TEM.infected[!filtered_DEG.TEM.infected$gene %in%
                                                         c("MALAT1", "LATS1", "LCLAT1", "RORA-AS1", "ENTPD1−AS1", "IFNG-AS1"), ]
filtered_DEG.TEM.infected <- filtered_DEG.TEM.infected[order(-xtfrm(filtered_DEG.TEM.infected$cluster), 
                                                             filtered_DEG.TEM.infected$gene), ]
Idents(TEM.infected) <- "memory_celltype"; DefaultAssay(TEM.infected) <- "RNA.clean"
AVG.TEM.infected.RNA <- AverageExpression(TEM.infected, return.seurat = T, group.by = 'memory_celltype')
AVG.TEM.infected.RNA$orig.ident <- colnames(AVG.TEM.infected.RNA)
DefaultAssay(AVG.TEM.infected.RNA) <- "RNA.clean"; Idents(AVG.TEM.infected.RNA) <- "orig.ident"

DefaultAssay(TEM.uninfected) <- "RNA.clean"; Idents(TEM.uninfected) <- "memory_celltype"
mylevels<- c("TRM_1", "TRM_2", "TEM")
Idents(TEM.uninfected) <- factor(Idents(TEM.uninfected), levels = rev(mylevels))
DEG.TEM.uninfected  <- FindAllMarkers(TEM.uninfected, only.pos = TRUE, min.pct = 0.1, test.use = "wilcox_limma",
                                    logfc.threshold = 0.25)
DEG.TEM.uninfected$BH <- p.adjust(DEG.TEM.uninfected $p_val, method = "BH", n = nrow(DEG.TEM.uninfected ))
DEG.TEM.uninfected <- na.omit(DEG.TEM.uninfected)
DEG.TEM.uninfected <- DEG.TEM.uninfected[!duplicated(DEG.TEM.uninfected$gene),]
DEG.TEM.uninfected <- DEG.TEM.uninfected[DEG.TEM.uninfected$BH < 0.05,]
filtered_DEG.TEM.uninfected <- subset(DEG.TEM.uninfected, 
                                    grepl("BACH|JUN|FOS|CD44|TCF7|CCR|CXCR|CXCL|CCL|CD69|ICOS|TNF|IFN|IL|CTLA|ENTPD|LAG3|TOX|ITGA|KLRB1|
                                          ISG|IFI|LCK|TIGIT|CCL|ZEB|RORA|TGF|FAS|LCP2|LAT|FYN|HIF1A|CASP|ATG|LEF1|
                                          CD|ZAP70|CD40LG|IKZF|BIRC", gene))
filtered_DEG.TEM.uninfected <- filtered_DEG.TEM.uninfected[order(filtered_DEG.TEM.uninfected$cluster, filtered_DEG.TEM.uninfected$gene), ]
filtered_DEG.TEM.uninfected <- filtered_DEG.TEM.uninfected[!filtered_DEG.TEM.uninfected$gene %in%
                                                         c("MALAT1", "LATS1", "LCLAT1", "RORA-AS1", "ENTPD1−AS1"), ]
filtered_DEG.TEM.uninfected <- filtered_DEG.TEM.uninfected[order(-xtfrm(filtered_DEG.TEM.uninfected$cluster), 
                                                                 filtered_DEG.TEM.uninfected$gene), ]
Idents(TEM.uninfected) <- "memory_celltype"; DefaultAssay(TEM.uninfected) <- "RNA.clean"
AVG.TEM.uninfected.RNA <- AverageExpression(TEM.uninfected, return.seurat = T, group.by = 'memory_celltype')
AVG.TEM.uninfected.RNA$orig.ident <- colnames(AVG.TEM.uninfected.RNA)
DefaultAssay(AVG.TEM.uninfected.RNA) <- "RNA.clean"; Idents(AVG.TEM.uninfected.RNA) <- "orig.ident"

DoHeatmap(object = AVG.TEM.infected.RNA, 
          features = filtered_DEG.TEM.infected_gene, slot = "scale.data", draw.lines = FALSE,
          disp.min = -1, disp.max = 1, assay = 'RNA.clean') + scale_fill_gradientn(colors = c("grey", "orange", "red"))
DoHeatmap(object = AVG.TEM.uninfected.RNA, 
          features = filtered_DEG.TEM.uninfected_gene, slot = "scale.data", draw.lines = FALSE,
          disp.min = -1, disp.max = 1, assay = 'RNA.clean') + scale_fill_gradientn(colors = c("grey", "orange", "red"))

#Figure 1J
#require objects from CD4 T cell data processing.R 
#require a list of antibodies in Antibody_stats.TEM.infected and Antibody_stats.TEM.uninfected with mean expression > mean expression of their specific isotype controls (Z > 2; two-sample Z-test), see Methods in associated manuscript.

DefaultAssay(TEM.infected) <- 'Antibody'; Idents(TEM.infected) <- "memory_celltype"
mylevels<- c("TRM_1", "TRM_2", "TEM")
Idents(TEM.infected) <- factor(Idents(TEM.infected), levels = rev(mylevels))

DEA.TEM.infected <- FindAllMarkers(TEM.infected, min.pct = 0.1, logfc.threshold = 0.1, test.use = "wilcox_limma",
                                             only.pos = TRUE, pseudocount.use = 0.000001, features = rownames(Antibody_stats.TEM.infected))
DEA.TEM.infected$BH <- p.adjust(DEA.TEM.infected$p_val, method = "BH", n = nrow(DEA.TEM.infected))
DEA.TEM.infected$gene <- rownames(DEA.TEM.infected)
DEA.TEM.infected <- DEA.TEM.infected[DEA.TEM.infected$BH < 0.05,]
DEA.TEM.infected <- DEA.TEM.infected[order(DEA.TEM.infected$cluster, DEA.TEM.infected$gene), ]
DEA.TEM.infected <- DEA.TEM.infected[order(-xtfrm(DEA.TEM.infected$cluster), 
                                           DEA.TEM.infected$gene), ]
Idents(TEM.infected) <- "memory_celltype"; DefaultAssay(TEM.infected) <- "Antibody"
AVG.TEM.infected.Antibody <- AverageExpression(TEM.infected, return.seurat = T, group.by = 'memory_celltype')
AVG.TEM.infected.Antibody$orig.ident <- colnames(AVG.TEM.infected.Antibody)
DefaultAssay(AVG.TEM.infected.Antibody) <- "RNA.clean"; Idents(AVG.TEM.infected.Antibody) <- "orig.ident"

DefaultAssay(TEM.uninfected) <- 'Antibody'; Idents(TEM.uninfected) <- "memory_celltype"
mylevels<- c("TRM_1", "TRM_2", "TEM")
Idents(TEM.uninfected) <- factor(Idents(TEM.uninfected), levels = rev(mylevels))
DEA.TEM.uninfected <- FindAllMarkers(TEM.uninfected, min.pct = 0.3, logfc.threshold = 0.3, test.use = "wilcox_limma",
                                   only.pos = TRUE, pseudocount.use = 0.000001, features = rownames(Antibody_stats.TEM.uninfected))
DEA.TEM.uninfected$BH <- p.adjust(DEA.TEM.uninfected$p_val, method = "BH", n = nrow(DEA.TEM.uninfected))
DEA.TEM.uninfected$gene <- rownames(DEA.TEM.uninfected)
DEA.TEM.uninfected <- DEA.TEM.uninfected[DEA.TEM.uninfected$BH < 0.05,]
DEA.TEM.uninfected <- DEA.TEM.uninfected[order(DEA.TEM.uninfected$cluster, DEA.TEM.uninfected$gene), ]
DEA.TEM.uninfected <- DEA.TEM.uninfected[order(-xtfrm(DEA.TEM.uninfected$cluster), 
                                               DEA.TEM.uninfected$gene), ]
Idents(TEM.uninfected) <- "memory_celltype"; DefaultAssay(TEM.uninfected) <- "Antibody"
AVG.TEM.uninfected.Antibody <- AverageExpression(TEM.uninfected, return.seurat = T, group.by = 'memory_celltype')
AVG.TEM.uninfected.Antibody$orig.ident <- colnames(AVG.TEM.uninfected.Antibody)
DefaultAssay(AVG.TEM.uninfected.Antibody) <- "Antibody"; Idents(AVG.TEM.uninfected.Antibody) <- "orig.ident"

DoHeatmap(object = AVG.TEM.infected.Antibody, 
          features = DEA.TEM.infected_gene, slot = "scale.data", draw.lines = FALSE,
          disp.min = -1, disp.max = 1, assay = 'Antibody') + scale_fill_gradientn(colors = c("grey", "orange", "red"))
DoHeatmap(object = AVG.TEM.uninfected.Antibody, 
          features = DEA.TEM.uninfected_gene, slot = "scale.data", draw.lines = FALSE,
          disp.min = -1, disp.max = 1, assay = 'Antibody') + scale_fill_gradientn(colors = c("grey", "orange", "red"))

#Figure 1K
#require objects from CD4 T cell data processing.R and build GSEA.R
Idents(CD4.clean) <- "memory_celltype"; DefaultAssay(CD4.clean) <- "RNA.clean"
DEG.CD4.TEM_1_2.GSEA <- wilcoxauc(CD4.clean, group_by = "memory_celltype", seurat_assay = "RNA.clean")
DEG.CD4.TEM_1_2.GSEA <- subset(DEG.CD4.TEM_1_2.GSEA, DEG.CD4.TEM_1_2.GSEA$group == "TRM_1")
Rank.DEG.CD4.TEM_1_2.GSEA <- DEG.CD4.TEM_1_2.GSEA$logFC
names(Rank.DEG.CD4.TEM_1_2.GSEA) <- DEG.CD4.TEM_1_2.GSEA$feature
fgsea_DEG.CD4.TEM_1_2.GSEA <- fgsea(pathways = fgsea_sets, 
                                      stats = Rank.DEG.CD4.TEM_1_2.GSEA,
                                      eps   = 0.0, minSize=15, maxSize=500, scoreType = "pos")
DEG.CD4.TEM_2_1.GSEA <- wilcoxauc(CD4.clean, group_by = "memory_celltype", seurat_assay = "RNA.clean")
DEG.CD4.TEM_2_1.GSEA <- subset(DEG.CD4.TEM_2_1.GSEA, DEG.CD4.TEM_2_1.GSEA$group == "TRM_2")
Rank.DEG.CD4.TEM_2_1.GSEA <- DEG.CD4.TEM_2_1.GSEA$logFC
names(Rank.DEG.CD4.TEM_2_1.GSEA) <- DEG.CD4.TEM_2_1.GSEA$feature
fgsea_DEG.CD4.TEM_2_1.GSEA <- fgsea(pathways = fgsea_sets, 
                                            stats = Rank.DEG.CD4.TEM_2_1.GSEA,
                                            eps   = 0.0, minSize=15, maxSize=500, scoreType = "pos")
DEG.CD4.TEM.GSEA <- wilcoxauc(CD4.clean, group_by = "memory_celltype", seurat_assay = "RNA.clean")
DEG.CD4.TEM.GSEA <- subset(DEG.CD4.TEM.GSEA, DEG.CD4.TEM.GSEA$group == "TEM")
Rank.DEG.CD4.TEM.GSEA <- DEG.CD4.TEM.GSEA$logFC
names(Rank.DEG.CD4.TEM.GSEA) <- DEG.CD4.TEM.GSEA$feature
fgsea_DEG.CD4.TEM.GSEA <- fgsea(pathways = fgsea_sets, 
                                    stats = Rank.DEG.CD4.TEM.GSEA,
                                    eps   = 0.0, minSize=15, maxSize=500, scoreType = "pos")

a1 <- plotEnrichment(fgsea_sets[["GOBP_NEGATIVE_REGULATION_OF_ADAPTIVE_IMMUNE_RESPONSE"]], 
                     Rank.DEG.CD4.TEM_1_2.GSEA) 
b1 <- plotEnrichment(fgsea_sets[["GOBP_NEGATIVE_REGULATION_OF_ADAPTIVE_IMMUNE_RESPONSE"]], 
                     Rank.DEG.CD4.TEM_2_1.GSEA) 
c1 <- plotEnrichment(fgsea_sets[["GOBP_NEGATIVE_REGULATION_OF_ADAPTIVE_IMMUNE_RESPONSE"]], 
                     Rank.DEG.CD4.TEM.GSEA) 

ggplot(a1$data, aes(x, y)) + geom_line(color = "red", size = 2) + 
  geom_line(data = b1$data, color = "blue", size = 2) +
  geom_line(data = c1$data, color = "green", size = 2) + labs(title = "GOBP_NEGATIVE_REGULATION_OF_ADAPTIVE_IMMUNE_RESPONSE") + xlab("rank") + ylab("enrichment score") + geom_hline(
    yintercept=0) + geom_vline(xintercept=0) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                     panel.background = element_blank())

#Figure 1L
#require objects from CD4 T cell data processing.R and build GSEA.R
a1 <- plotEnrichment(fgsea_sets[["PID_TCR_CALCIUM_PATHWAY"]], 
                     Rank.DEG.CD4.TEM_1_2.GSEA) 
b1 <- plotEnrichment(fgsea_sets[["PID_TCR_CALCIUM_PATHWAY"]], 
                     Rank.DEG.CD4.TEM_2_1.GSEA) 
c1 <- plotEnrichment(fgsea_sets[["PID_TCR_CALCIUM_PATHWAY"]], 
                     Rank.DEG.CD4.TEM.GSEA) 

ggplot(a1$data, aes(x, y)) + geom_line(color = "red", size = 2) + 
  geom_line(data = b1$data, color = "blue", size = 2) +
  geom_line(data = c1$data, color = "green", size = 2) + labs(title = "PID_TCR_CALCIUM_PATHWAY") + xlab("rank") + ylab("enrichment score") + geom_hline(
    yintercept=0) + geom_vline(xintercept=0) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                     panel.background = element_blank())

