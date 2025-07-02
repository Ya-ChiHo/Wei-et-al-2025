#Figure 5A
#require objects from CD8 T cell data processing.R and antigen specificity information (see supplemental files from associated manuscript)
DimPlot(CD8.clean.5copy.clone, reduction = "umap.wnn", group.by = "antigen.specificity", split.by = "Condition",
        pt.size = 0.8, raster=FALSE,
        cols= c("grey", "red", "blue"), order = c("CMV", "HIV", "FALSE"))

#Figure 5E
#require objects from CD8 T cell data processing.R and antigen specificity information (see supplemental files from associated manuscript)
Idents(CD8.clean.5copy) <- "is.TCR"
CD8.clean.5copy.TCR <- subset(CD8.clean.5copy, idents = c("TCRB"))
Idents(CD8.clean.5copy) <- "Condition"
CD8.clean.5copy.TCR.infected <- subset(CD8.clean.5copy, idents = c("infected"))

Idents(CD8.clean.5copy.TCR.infected) <- "antigen.specificity"; DefaultAssay(CD8.clean.5copy.TCR.infected) <- 'RNA.clean'
DEG.HIV_vs_CMV <- FindMarkers(CD8.clean.5copy.TCR.infected, ident.1 = "HIV", ident.2 = "CMV",
                                    min.pct = 0.33, logfc.threshold = 0.33,
                                    only.pos = FALSE, test.use = 'wilcox')
DEG.HIV_vs_CMV$BH <- p.adjust(DEG.HIV_vs_CMV$p_val, method = "BH", n = nrow(DEG.HIV_vs_CMV))
DEG.HIV_vs_CMV$gene <- rownames(DEG.HIV_vs_CMV)

#Figure 5G, H
#require objects from CD8 T cell data processing.R and antigen specificity information (see supplemental files from associated manuscript)
DefaultAssay(CD8.clean.5copy.TCR.infected) <- "MACS2ATAC"
Idents(CD8.clean.5copy.TCR.infected) <- "antigen.spec.celltype"

CD8.clean.5copy.TCR.infected.sample <- subset(CD8.clean.5copy.TCR.infected, 
                                            idents = c("HIV", "CMV", "Naive", "All_others))
Idents(CD8.clean.5copy.TCR.infected.sample) <- "antigen.spec.celltype"
CoveragePlot(CD8.clean.5copy.TCR.infected.sample, region = "LEF1",ranges.title = "MACS2", extend.upstream = 500,extend.downstream = 500)
CoveragePlot(CD8.clean.5copy.TCR.infected.sample, region = "ENTPD1",ranges.title = "MACS2", extend.upstream = 500,extend.downstream = 500)