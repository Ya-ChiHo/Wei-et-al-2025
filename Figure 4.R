#Figure 4B
#require objects from CD4 T cell data processing.R
DimPlot(CD4.clean.5copy.TCR, reduction = "umap.wnn", group.by = "is.clone", split.by = "Condition",
        pt.size = 0.5, cols= c("grey", "red"), order = c("TRUE", "FALSE"), raster=FALSE) 

#Figure 4G
#require objects from CD8 T cell data processing.R
DimPlot(CD8.clean.5copy.TCR, reduction = "umap.wnn", group.by = "is.clone", split.by = "Condition",
        pt.size = 0.5, cols= c("grey", "red"), order = c("TRUE", "FALSE"), raster=FALSE) 

#Figure 4L
#require objects from CD8 T cell data processing.R
Idents(CD8.clean.5copy.TCR) <- "Condition"
CD8.clean.5copy.TCR.infected <- subset(CD8.clean.5copy.TCR, idents = c("infected"))

DefaultAssay(CD8.clean.5copy.TCR.infected) <- "chromvar"; Idents(CD8.clean.5copy.TCR.infected) <- "is.clone"
DEM.infected.clone_vs_notclone <- FindMarkers(CD8.clean.5copy.TCR.infected, ident.1 = "TRUE", ident.2 = "FALSE",
                                              min.pct = 0, logfc.threshold = 0.1, 
                                              only.pos = FALSE, mean.fxn = rowMeans, fc.name = "avg_diff")
DEM.infected.clone_vs_notclone$BH <- p.adjust(DEM.infected.clone_vs_notclone$p_val, method = "BH", n = nrow(DEM.infected.clone_vs_notclone))
DEM.infected.clone_vs_notclone$gene <- ConvertMotifID(motif.names, id = rownames(DEM.infected.clone_vs_notclone))

#Figure 4M
#require objects from CD8 T cell data processing.R
DefaultAssay(CD8.clean.5copy.TCR.infected) <- "RNA.clean"; Idents(CD8.clean.5copy.TCR.infected) <- "is.clone"
DEG.infected.clone_vs_notclone <- FindMarkers(CD8.clean.5copy.TCR.infected, ident.1 = "TRUE", ident.2 = "FALSE",
                                              only.pos = FALSE, test.use = 'wilcox', min.pct = 0.25, logfc.threshold = 0.25)
DEG.infected.clone_vs_notclone$BH <- p.adjust(DEG.infected.clone_vs_notclone$p_val, method = "BH", n = nrow(DEG.infected.clone_vs_notclone))
DEG.infected.clone_vs_notclone$gene <- rownames(DEG.infected.clone_vs_notclone)

#Figure 4N
#require objects from CD8 T cell data processing.R
DefaultAssay(CD8.clean.5copy.TCR.infected) <- "Antibody"; Idents(CD8.clean.5copy.TCR.infected) <- "is.clone"
table(CD8.clean.5copy.TCR.infected$is.clone)
DEA.infected.clone_vs_notclone <- FindMarkers(CD8.clean.5copy.TCR.infected, ident.1 = "TRUE", ident.2 = "FALSE",
                                              only.pos = FALSE, test.use = 'wilcox', min.pct = 0.1, logfc.threshold = 0.1,
                                              features = rownames(Antibody_stats))
DEA.infected.clone_vs_notclone$BH <- p.adjust(DEA.infected.clone_vs_notclone$p_val, method = "BH", n = nrow(DEA.infected.clone_vs_notclone))
DEA.infected.clone_vs_notclone$gene <- rownames(DEA.infected.clone_vs_notclone)