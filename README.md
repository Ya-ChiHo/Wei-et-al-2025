# Wei-et-al-2025

R Scripts related to the Wei et al 2025 manuscript.

Please see STAR Methods in Wei et al 2025 manuscript for software and package versions.

These scripts are minimally commented to be functional and reproducible, they are not tutorials. Please see the Methods section in the manuscript for additional details and non-R related analyses.

Scripts in 'CD3 isolation data preprocessing.R' and 'whole gut data preprocessing.R' describe steps taken to process Cellranger data to generate Seurat objects: Building and merging ATAC, RNA, and CITE Seurat objects, recall peaks by MACS2, doublet removal, low quality cell filter, labeling HIV-1-infected cells (identified through Bowtie2 and STAR), data normalization and batch effect removal, building gene accessibility data, adding chromVAR transcription factor motif scores, and celltype annotation.

Each sample was split into two aliquots: one underwent CD3-positive magnetic selection to enrich for T lymphocytes (CD3 isolation data), while the other remained untreated (whole gut data) prior to processing for DOGMA-seq.

CD4 T cell data processing.R: steps to merge CD4 T cells from whole gut and CD3 isolated samples after pre-processing, re-perform data normalization and batch effect removal, re-build gene accessibility, and re-adding chromVAR transcription factor motif scores, adding celltype annotations, and adding TCR information for subsetted cells.
CD8 T cell data processing.R: steps to merge CD8 T cells from whole gut and CD3 isolated samples after pre-processing, re-perform data normalization and batch effect removal, re-build gene accessibility, and re-adding chromVAR transcription factor motif scores, adding celltype annotations, and adding TCR information for subsetted cells.

build GSEA.R - minimal example scripts to build required R objects for Gene Set Enrichment Analysis (by msigdbr and fgsea): GSEA gene sets, intermediates for GSEA (e.g., gene ranks by differential expression, gene set normalized enrichment scores, leading edge genes, etc.).

build GRN.R - minimal example scripts to build required R objects for Gene Regulatory Network Analyses (by cisTopic and FigR): constructing cisTopics, determing DORCs, and making FigR matrix.

build diffusion map - minimal example scripts to bin cells by BACH2 TF accessibility (chromVAR), determining differentially expressed genes and transcription factors from BACH2-low to BACH2-high cells, and performing diffusion map.

build cell cell communication.R - minimal example scripts to build required R objects for ligand-receptor interaction analyses (by Scriabin): performing ALRA-based denoising, identifying significant IPs, and determining NichNet-predicted ligand targets.

All scripts in the 'Figures.R' require processed Seurat objects generated from either 'CD4 T cell data processing.R' or 'CD8 T cell data processing.R' and R objects generated from the 'build ---.R' scripts. Scripts in 'Figures' generate the plots as shown in manuscript main figures and, when relevant, the functions/data processing steps required to generate these plots. Scripts in 'Figures.R' may be modified to reproduce all supplemental figures. If there are particular aspects of the analyses you would like to see that are not here, or have any other questions, please email Yulong.Wei@yale.edu or Ya-Chi.Ho@yale.edu.

Raw reads and processed data (e.g., cellranger filtered feature bc matrices, hashtag keys, fragment.tsv etc.) are hosted on GEO:GSE299348 
For all other required processed data (e.g., HIV+ cell barcodes, TCR igblast results, cell barcode and clonotype) see this repository and the supplemental files in the associated manuscript.

Comment on Seurat v5 and Signac for Seurat v5 (2025-07-01): Note that some functions have been removed/replaced, while default parameters changed in others (e.g., FIndMarkers), and data structure different (e.g., layers), in Seurat v5 vs Seurat v4. Version difference may lead to some scripts breaking and may lead to slight graphical differences (e.g., volcanoplot because of change in FC in FindMarkers), but do not otherwise affect the main results (e.g., key differentially expressed genes).
