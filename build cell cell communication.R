library(Seurat)
options(Seurat.object.assay.version = 'v5')
library(Signac)
library(ggplot2)
library(ggraph)
library(cowplot)
library(patchwork)
library(dplyr)
library(tidyverse)
library(GenomicRanges)
library(future)
library(EnsDb.Hsapiens.v86)
library(TFBSTools)
library(scSHC)
library(stringr)
library(scriabin)
library(pbapply)
library(scriabin)
library(ComplexHeatmap)
library(cowplot)
library(conflicted)
library(flashClust)
library(ade4)
library(glmGamPoi)
library(CelliD)

conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::intersect)
conflicts_prefer(dplyr::union)
conflicts_prefer(dplyr::setdiff)
scriabin::load_nichenet_database()

#these scripts apply scriabin package to identify ligand-receptor sender and receiver cells using RNA assay.
#a merged (CD4 T cell and CD8 T cell) Seurat Object is required. Please see "T cell data processing" for pre-processing steps (including re-normalization and scaling, ATAC and RNA batch effect correction, WNN integration)
#cell subset annotations for the CD4 T and CD8 T merged object were transferred from each CD4 T cell and CD8 T cell objects.

#modified functions from Scriabin####
#functions modified from Scriabin (https://github.com/BlishLab/scriabin) 
InteractionPrograms1 <- function (object, assay = "alra", slot = "data", species = "human", 
                                  database = "OmniPath", ligands = NULL, recepts = NULL, iterate.threshold = 300, 
                                  n.iterate = NULL, specific = F, ranked_genes = NULL, return.mat = T, 
                                  softPower = NULL, r2_cutoff = 0.6, min.size = 5, plot.mods = F, 
                                  tree.cut.quantile = 0.4, threads = NULL, cell_types = NULL, 
                                  min.cell = 3) 
{
  if (database == "custom") {
    if (is.null(ligands) | is.null(recepts)) {
      stop("To use custom database, please supply equidimensional character vectors of ligands and recepts")
    }
    message("Using custom database")
    ligands <- ligands
    recepts <- recepts
    lit.put <- data.frame(pair.name = paste(ligands, recepts, 
                                            sep = "_"), source_genesymbol = ligands, target_genesymbol = recepts)
  }
  else if ((!is.null(ligands) | !is.null(recepts)) & database != 
           "custom") {
    stop("To use custom ligand or receptor lists, set database = 'custom'")
  }
  else {
    lit.put <- scriabin::LoadLR(species = species, database = database)
    ligands <- as.character(lit.put[, "source_genesymbol"])
    recepts <- as.character(lit.put[, "target_genesymbol"])
  }
  ligands.use <- intersect(ligands, rownames(object@assays[[assay]]))
  recepts.use <- intersect(recepts, rownames(object@assays[[assay]]))
  genes.use = union(ligands.use, recepts.use)
  lit.put <- lit.put %>% dplyr::filter(source_genesymbol %in% 
                                         ligands.use) %>% dplyr::filter(target_genesymbol %in% 
                                                                          recepts.use)
  ligands <- as.character(lit.put[, "source_genesymbol"])
  recepts <- as.character(lit.put[, "target_genesymbol"])
  if (specific) {
    message("Only considering genes in per-cell gene signature")
    ranked_names <- lapply(ranked_genes, function(x) {
      names(x)
    })
    ranked_mat <- as.matrix(reshape2::dcast(reshape2::melt(t(bind_rows(ranked_names))), 
                                            formula = value ~ Var1) %>% column_to_rownames("value"))
    genes.use <- intersect(genes.use, rownames(ranked_mat))
    ranked_mat <- ranked_mat[genes.use, ]
    cell.exprs <- GetAssayData(object, assay = assay, slot = slot)[genes.use, 
    ]
    cell.exprs[is.na(ranked_mat)] <- 0
    cell.exprs <- as.data.frame(cell.exprs) %>% rownames_to_column(var = "gene")
  }
  else {
    cell.exprs <- GetAssayData(object, assay = assay, slot = slot)[genes.use, 
    ]
  }
  ligands.df <- data.frame(ligands)
  ligands.df$id <- 1:nrow(ligands.df)
  recepts.df <- data.frame(recepts)
  recepts.df$id <- 1:nrow(recepts.df)
  if (ncol(object) > iterate.threshold) {
    message("\nIteratively generating interaction matrix")
    if (is.null(n.iterate)) {
      n.rep <- round(sqrt(ncol(object)/iterate.threshold))
    }
    else {
      n.rep = n.iterate
    }
    message(paste("Will perform", n.rep, "iterations to approximate TOM"))
    if (is.null(cell_types)) {
      warning("We recommend setting a cell_types parameter so that all cell types are included in each sequence of TOM generation")
    }
    mat_list <- lapply(seq_along(1:n.rep), function(z) {
      if (!is.null(cell_types)) {
        sub_prop <- (object@meta.data %>% rownames_to_column("cell"))[, 
                                                                      c("cell", cell_types)]
        colnames(sub_prop) <- c("cell", "var")
        cells <- sub_prop %>% group_by(var) %>% dplyr::mutate(prop = round(iterate.threshold * 
                                                                             n()/nrow(.))) %>% dplyr::mutate(prop = ifelse(prop < 
                                                                                                                             min.cell, min.cell, prop)) %>% group_by(var) %>% 
          dplyr::sample_n(10) %>% pull(cell)
        cell.exprs.sub <- as.data.frame(cell.exprs[, 
                                                   cells]) %>% rownames_to_column(var = "gene")
      }
      else {
        cell.exprs.sub <- as.data.frame(cell.exprs[, 
                                                   sample(colnames(cell.exprs), iterate.threshold)]) %>% 
          rownames_to_column(var = "gene")
      }
      cell.exprs.rec <- merge(recepts.df, cell.exprs.sub, 
                              by.x = "recepts", by.y = "gene", all.x = T)
      cell.exprs.rec <- cell.exprs.rec[order(cell.exprs.rec$id), 
      ]
      cell.exprs.lig <- merge(ligands.df, cell.exprs.sub, 
                              by.x = "ligands", by.y = "gene", all.x = T)
      cell.exprs.lig <- cell.exprs.lig[order(cell.exprs.lig$id), 
      ]
      a <- as.matrix(cell.exprs.lig[, 3:ncol(cell.exprs.lig)])
      a[is.na(a)] <- 0
      b <- as.matrix(cell.exprs.rec[, 3:ncol(cell.exprs.rec)])
      b[is.na(b)] <- 0
      m <- sqrt(as.sparse((pbsapply(1:nrow(a), function(i) tcrossprod(a[i, 
      ], b[i, ])))))
      colnames(m) <- paste(cell.exprs.lig$ligands, cell.exprs.rec$recepts, 
                           sep = "=")
      cna <- rep(colnames(a), ncol(a))
      cnb <- rep(colnames(a), each = ncol(a))
      rownames(m) <- paste(cna, cnb, sep = "=")
      m <- m[, Matrix::colSums(m) > 0]
      m_cor <- 0.5 + (0.5 * corSparse(m))
      rownames(m_cor) <- colnames(m)
      colnames(m_cor) <- colnames(m)
      m_cor
    })
    lr_union <- unique(unlist(lapply(mat_list, rownames)))
    mat_list_2 <- pblapply(mat_list, function(m) {
      missing <- lr_union[lr_union %notin% rownames(m)]
      if (length(missing) == 0) {
        return(m[lr_union, lr_union])
      }
      else if (length(missing) == 1) {
        m <- rbind(m, missing = NA)
        rownames(m)[nrow(m)] <- missing
        m <- cbind(m, missing = NA)
        colnames(m)[ncol(m)] <- missing
        return(m[lr_union, lr_union])
      }
      else if (length(missing) > 1) {
        add_rows <- matrix(NA, nrow = length(missing), 
                           ncol = ncol(m))
        rownames(add_rows) <- missing
        m <- rbind(m, add_rows)
        add_cols <- matrix(NA, nrow = nrow(m), ncol = length(missing))
        colnames(add_cols) <- missing
        m <- cbind(m, add_cols)
        return(m[lr_union, lr_union])
      }
    })
    m_cor <- apply(simplify2array(mat_list_2), 1:2, function(x) {
      mean(x, na.rm = T)
    })
  }
  else {
    message(paste("\nGenerating Interaction Matrix..."))
    cell.exprs.sub <- as.data.frame(cell.exprs) %>% rownames_to_column(var = "gene")
    cell.exprs.rec <- merge(recepts.df, cell.exprs.sub, by.x = "recepts", 
                            by.y = "gene", all.x = T)
    cell.exprs.rec <- cell.exprs.rec[order(cell.exprs.rec$id), 
    ]
    cell.exprs.lig <- merge(ligands.df, cell.exprs.sub, by.x = "ligands", 
                            by.y = "gene", all.x = T)
    cell.exprs.lig <- cell.exprs.lig[order(cell.exprs.lig$id), 
    ]
    a <- as.matrix(cell.exprs.lig[, 3:ncol(cell.exprs.lig)])
    a[is.na(a)] <- 0
    b <- as.matrix(cell.exprs.rec[, 3:ncol(cell.exprs.rec)])
    b[is.na(b)] <- 0
    m <- sqrt(as.sparse((pbsapply(1:nrow(a), function(i) tcrossprod(a[i, 
    ], b[i, ])))))
    colnames(m) <- paste(cell.exprs.lig$ligands, cell.exprs.rec$recepts, 
                         sep = "=")
    cna <- rep(colnames(object), ncol(object))
    cnb <- rep(colnames(object), each = ncol(object))
    rownames(m) <- paste(cna, cnb, sep = "=")
    m <- m[, Matrix::colSums(m) > 0]
    m_cor <- 0.5 + (0.5 * corSparse(m))
    rownames(m_cor) <- colnames(m)
    colnames(m_cor) <- colnames(m)
  }
  if (!is.null(threads)) {
    enableWGCNAThreads(nThreads = threads)
  }
  if (is.null(softPower)) {
    message("Automatically selecting softPower . . .")
    sp_det <- pickSoftThreshold.fromSimilarity(m_cor, RsquaredCut = r2_cutoff)
    softPower = sp_det$powerEstimate
    if (is.na(softPower) | softPower > 3) {
      warning("No appropriate softPower found to reach minimum scale free topology fit. Proceeding without soft thresholding, interpret results with caution")
      softPower = 1
    }
  }
  message("Identifying modules")
  adj <- m_cor^softPower
  tom <- TOMsimilarity(adj, TOMType = "signed")
  colnames(tom) <- rownames(tom) <- rownames(m_cor)
  geneTree = flashClust(as.dist(1 - tom), method = "complete")
  dynamicMods = cutreeDynamic(dendro = geneTree, method = "tree", 
                              minClusterSize = min.size, cutHeight = quantile(geneTree$height, 
                                                                              0.5))
  dynamicColors = labels2colors(dynamicMods)
  if (plot.mods) {
    plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", 
                        dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, 
                        guideHang = 0.05, main = "Gene dendrogram and module colors")
  }
  module_colors = setdiff(unique(dynamicColors), "grey")
  modules = lapply(module_colors, function(x) {
    colnames(m_cor)[which(dynamicColors == x)]
  })
  names(modules) = module_colors
  module_melt <- reshape2::melt(modules)
  colors <- scriabin::mapvalues(colnames(m_cor), from = module_melt$value, 
                                to = module_melt$L1, warn_missing = F)
  colors[colors %notin% unique(module_melt$L1)] <- "grey"
  Alldegrees1 = intramodularConnectivity(adj, colors)
  names(modules) <- paste0("IP-", 1:length(modules))
  if (return.mat) {
    return(list(cor_mat = m_cor, tom = tom, modules = modules, 
                connectivity = Alldegrees1))
  }
  else {
    return(list(modules = modules, connectivity = im_results))
  }
}

FindAllInteractionPrograms1 <- function (seu, group.by = NULL, sim_threshold = 0.15, ...) 
{
  if (is.null(group.by)) {
    message("Grouping by current idents")
    seu$mod_grouping <- Idents(seu)
    group.by = "mod_grouping"
  }
  seu_split <- SplitObject(seu, split.by = group.by)
  q_mods <- lapply(seu_split, function(x) {
    InteractionPrograms1(object = x, return.mat = T, ...)
  })
  mod_list <- unlist(lapply(q_mods, function(x) {
    x[[3]]
  }), recursive = F)
  tom_list <- lapply(q_mods, function(x) {
    x[[2]]
  })
  m_cor_list <- lapply(q_mods, function(x) {
    x[[1]]
  })
  con_list <- lapply(q_mods, function(x) {
    x[[4]]
  })
  names(con_list) <- names(tom_list) <- names(m_cor_list) <- names(q_mods)
  message("Merging similar modules")
  lrsum <- t(do.call(rbind, lapply(1:length(names(mod_list)), 
                                   function(x) {
                                     unique(unlist(mod_list)) %in% mod_list[[x]]
                                   })))
  lrsum[lrsum == T] <- 1
  colnames(lrsum) <- names(mod_list)
  rownames(lrsum) <- unique(unlist(mod_list))
  d <- as.matrix(dist.binary(t(lrsum), method = 1))
  d <- 1 - d^2
  merge_num = 1
  while (max(d[d < 1]) > sim_threshold) {
    tomerge <- rownames(which(d == max(d[d < 1]), arr.ind = T))
    new_mod <- unique(unlist(mod_list[tomerge]))
    mod_list <- mod_list[names(mod_list) %notin% tomerge]
    mod_list[[length(mod_list) + 1]] <- new_mod
    names(mod_list)[length(mod_list)] <- paste0("merge_", 
                                                merge_num)
    merge_num = merge_num + 1
    lrsum <- t(do.call(rbind, lapply(1:length(names(mod_list)), 
                                     function(x) {
                                       unique(unlist(mod_list)) %in% mod_list[[x]]
                                     })))
    lrsum[lrsum == T] <- 1
    colnames(lrsum) <- names(mod_list)
    rownames(lrsum) <- unique(unlist(mod_list))
    d <- as.matrix(dist.binary(t(lrsum), method = 1))
    d <- 1 - d^2
  }
  names(mod_list) <- paste0("IP-", 1:length(mod_list))
  return(list(cor_mat = m_cor_list, tom = tom_list, modules = mod_list, 
              connectivity = con_list))
}

IPFeaturePlot1 <- function (seu, ip, cols = c("grey90", "blue", "orangered3"), 
                            order = T) 
{
  seu$lig_feature <- seu[["IPligands"]]@data[ip, ]
  seu$rec_feature <- seu[["IPreceptors"]]@data[ip, ]
  p <- FeaturePlot(seu, features = c("lig_feature", "rec_feature"), 
                   blend = T, combine = F, cols = cols, order = order, reduction = 'umap.wnn')
  p[[3]] + NoLegend() + labs(x = NULL, y = NULL, title = NULL) + 
    theme(aspect.ratio = 1, axis.ticks = element_blank(), 
          axis.text = element_blank())
}

#identify interaction programs IPs####
#perform ALRA-based denoising
Idents(CD4_CD8) <- "Condition"
DefaultAssay(CD4_CD8) <- "RNA.clean"
CD4_CD8 <- SCTransform(CD4_CD8, verbose = F) %>%
  RunPCA(verbose = F) %>%
  RunUMAP(dims = 1:30, verbose = F)
CD4_CD8 <- SeuratWrappers::RunALRA(CD4_CD8)

#find interaction programs
Idents(CD4_CD8) <- "major_celltype"; DefaultAssay(CD4_CD8) <- "alra"
Interaction_programs_compare <- FindAllInteractionPrograms1(CD4_CD8, 
                                                            group.by = "Condition", cell_types = "major_celltype",
                                                            assay = "alra")
Interaction_programs_compare_sig <- InteractionProgramSignificance(Interaction_programs_compare)
CD4_CD8 <- ScoreInteractionPrograms(CD4_CD8, mods = Interaction_programs_compare_sig)

#identify IPs that are enriched in HIV+ samples
infection_ip_ligands <- FindMarkers(CD4_CD8, group.by = "Condition", ident.1 = "infected", assay = "IPligands")
infection_ip_receptors <- FindMarkers(CD4_CD8, group.by = "Condition", ident.1 = "infected", assay = "IPreceptors")

#identify sender and receiver cells for significant IPs####
#repeat for each enriched IPs (from infection_ip_ligands and infection_ip_receptors)
infection_ip_ligands
infection_ip_receptors

#visualize sender and receiver cells for significant IP modules of ligand-receptor pairs, e.g., IP-2
IPFeaturePlot1(CD4_CD8, ip = "IP-2")
#visualize connectivity score of each ligand-receptor pairs within IP in HIV+ and HIV- samples
moi_IP2 <- reshape2::melt(Interaction_programs_compare_sig %>% dplyr::filter(name=="IP-2") %>%
                            select("lr_pair",contains("connectivity"))) %>% arrange(-value)
moi_IP2$lr_pair <- factor(moi_IP2$lr_pair, levels = unique(moi_IP2$lr_pair))
ggplot(moi_IP2, aes(x = lr_pair, y = value, color = variable)) + 
  geom_point() + theme_cowplot() + ggpubr::rotate_x_text() + labs(x = NULL, 
                      y = "Intramodular\nconnectivity") + ggtitle("IP-2")

#alluvial plot of NicheNet predicted ligand targets####
#repeat for each major_celltype
#repeat for any ligands_of_interest

DefaultAssay(CD4_CD8.infected) <- "RNA.clean"; DefaultAssay(CD4_CD8.uninfected) <- "RNA.clean"
variant_genes.infected <- IDVariantGenes(CD4_CD8.infected, group.by = "celltype1", assay = "RNA.clean")
gene_signature.infected <- GenerateCellSignature1(CD4_CD8.infected, variant_genes = variant_genes.infected)
active_ligands.infected <- RankActiveLigands(CD4_CD8.infected, signature_matrix = gene_signature.infected)

#alluvial plot for TRM-Th17:
PlotLigandTargetAlluvium(CD4_CD8.infected, signature_matrix = gene_signature.infected,
                         active_ligands = active_ligands.infected, 
                         receiver_cells = colnames(CD4_CD8.infected)[CD4_CD8.infected$major_celltype=="TRM-Th17"],
                         ligands_of_interest = c("TNF", "IFNG"))
