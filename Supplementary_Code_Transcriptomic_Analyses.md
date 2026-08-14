# Supplementary Code for Transcriptomic Analyses

(Object names and paths generalized for publication)

Organization: single-cell RNA-seq/regulon analyses, Visium HD spatial transcriptomics analyses, and other RNA-seq-related/transcriptome analyses. Headings were standardized by analysis step.


# Part I. Single-cell RNA-seq and regulon analyses


## 1. R setup, shared utility functions, and assumed Seurat objects


```r
## Required packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(Matrix)
library(SCENIC)
library(RcisTarget)
library(gplots)
library(VennDiagram)
library(ggridges)
library(ggrepel)
library(igraph)
 
############################################################
# 0. Utility functions shared across figures
############################################################
 
## 0-1. Gene filtering for removing pseudogenes and background noise 
gene_filter_pattern <- paste(
  "mt-", "MT-", "MMU", "SNO", "Gm", "Rn6", "RNA", "7SK",
  "SNORD", "SNORA", "SCARNA", "B3g", "Vmn", "Mir", "Rik",
  "Snora", "Snord", "LOC", "Rn4", "OTTMUSG", "Scarna",
  "RNU", "RMR", "Rpl", "Rps", "Rnu", "Rmr", "RPl", "RPS",
  "AW", "BC", "Malat1", "H2-", "B2m", "Actb",
  "AA", "AB", "AC", "AF", "AI", "AL", "AU", "AV", "AY",
  "Hist",
  sep = "|"
)
 
## 0-2. Spearman correlation heatmap function (Ward.D2 clustering)
make_spearman_heatmap <- function(mat,
                                  palette_low = "#0571b0",
                                  palette_mid = "#f7f7f7",
                                  palette_high = "#ca0020",
                                  breaks_q_low = 0.10,
                                  breaks_q_high = 0.91,
                                  step = 0.06,
                                  cex_row = 0.3,
                                  cex_col = 0.3,
                                  dendrogram = c("row", "both")) {
  random.matrix  <- matrix(runif(500, min = -1, max = 1), nrow = 50)
  quantile.range <- quantile(random.matrix, probs = seq(0, 1, 0.01))
  palette.breaks <- seq(
    quantile.range[paste0(100 * breaks_q_low, "%")],
    quantile.range[paste0(100 * breaks_q_high, "%")],
    step
  )
  color.palette  <- colorRampPalette(c(palette_low, palette_mid, palette_high))(
    length(palette.breaks) - 1
  )
 
  clustFunction <- function(x) {
    hclust(
      as.dist(1 - cor(t(as.matrix(x)), method = "spearman")),
      method = "ward.D2"
    )
  }
 
  gplots::heatmap.2(
    x          = mat,
    col        = color.palette,
    breaks     = palette.breaks,
    trace      = "none",
    symm       = TRUE,
    hclustfun  = clustFunction,
    dendrogram = match.arg(dendrogram),
    cexRow     = cex_row,
    cexCol     = cex_col,
    key        = TRUE
  )
}
 
## 0-3. Function to create binary annotation vectors (TF presence) for pheatmap
create_annotations <- function(tfs, data_matrix) {
  all_colnames <- colnames(data_matrix)
  annotation_vector <- as.numeric(all_colnames %in% tfs)
  names(annotation_vector) <- all_colnames
  return(annotation_vector)
}
 
## 0-4. Signature score calculation (Adam Haber / Itay Tirosh 2D scoring)
get_controls <- function(counts, gene.list, verbose = FALSE,
                         control.genes.per.gene = 10) {
  if (verbose) {
    cat(sprintf(
      "Finding %s background genes based on similarity to given gene set [%s genes]\n",
      control.genes.per.gene * length(gene.list), length(gene.list)
    ))
  }
  cat("Summarizing data\n")
 
  summary <- data.frame(
    gene       = rownames(counts),
    mean.expr  = Matrix::rowMeans(counts),
    fract.zero = Matrix::rowMeans(counts == 0),
    stringsAsFactors = FALSE
  )
 
  summary$mean.expr.s  <- scale(summary$mean.expr)
  summary$fract.zero.s <- scale(summary$fract.zero)
 
  actual.genes     <- summary[summary$gene %in% gene.list, ]
  background.genes <- summary[!summary$gene %in% gene.list, ]
 
  controls <- c()
 
  get_closest_genes <- function(i) {
    background.genes$dist <- sqrt(
      (background.genes$mean.expr.s - actual.genes$mean.expr.s[i])^2 +
        (background.genes$fract.zero.s - actual.genes$fract.zero.s[i])^2
    )
    ordered <- background.genes$gene[order(background.genes$dist)]
    ordered <- ordered[!ordered %in% controls]
    closest <- head(ordered, n = control.genes.per.gene)
    return(closest)
  }
 
  for (i in seq_len(nrow(actual.genes))) {
    closest  <- get_closest_genes(i)
    controls <- unique(c(controls, closest))
  }
 
  if (verbose) {
    cat(sprintf("Control gene selection complete. %s genes found.\n",
                length(controls)))
  }
  return(controls)
}
 
calculate_signature_score <- function(count_matrix, gene_list) {
  control_gene    <- get_controls(counts = count_matrix, gene.list = gene_list)
  signature_score <- colMeans(count_matrix[gene_list, , drop = FALSE], na.rm = TRUE) -
    colMeans(count_matrix[control_gene, , drop = FALSE], na.rm = TRUE)
  return(signature_score)
}
 
############################################################
# Assumed Seurat objects
############################################################
# seurat_all       : integrated EC object (4 organs)
# seurat_adipose   : adipose ECs
# seurat_heart     : heart ECs
# seurat_muscle    : skeletal muscle ECs
# seurat_liver     : liver ECs
#
# For SCENIC:
# seurat_capillary : capillary/sinusoidal EC subset
#
# For aging/diet regulon ratios:
# adipose_aging, adipose_diet, heart_aging, heart_diet,
# muscle_aging, muscle_diet, liver_aging, liver_diet
############################################################
```


## 2. Tissue-specific endothelial signatures and regulon inference


```r
############################################################
 
############################
 
# DEG-based correlation heatmap across tissues
############################
 
# Filter genes by name pattern
gene_names <- rownames(seurat_all)
gene_names <- gene_names[!grepl(gene_filter_pattern, gene_names)]
seurat_filtered <- seurat_all[gene_names, ]
 
Idents(seurat_filtered) <- "Tissue"
 
# Find top 100 markers per tissue cluster
deg_markers <- FindAllMarkers(
  seurat_filtered,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25
)
 
deg_top100 <- deg_markers %>%
  group_by(cluster) %>%
  dplyr::top_n(n = 100, wt = avg_log2FC)
 
# Spearman correlation matrix of SCT expression
cor_mat <- cor(
  method = "spearman",
  log2(t(as.matrix(
    seurat_filtered@assays[["SCT"]]@data[deg_top100$gene, ]
  )) + 1)
)
 
options(repr.plot.width = 12, repr.plot.height = 12)
make_spearman_heatmap(cor_mat, dendrogram = "row")
 
 
############################
 
# Regulon–target presence heatmap with organ annotation
############################
 
# SCENIC regulon-target table
regulon_targets <- read.csv(
  "/path/to/SCENIC/output/Step2_regulonTargetsInfo_target.csv",
  header = TRUE
)
# regulon_targets: columns 'TF', 'gene'
 
# TF order (from RSS plot)
TF_order <- c(
  "Etv1","Maf","Meis1","Gata4","Rarb","Nr2f1","Nr2f2","Irf8","Hdac1",
  "Creb3l2","Nr5a2","Esrrb","Nfic","Tfec","Zeb1","Mafg","Nfkb1","Mlxip",
  "Bach1","Lhx6","Dbp","Tef","Hoxb7","Xrcc4","Sox13","Srf","Zfp143","Ar",
  "Klf16","Nr1d1","Wt1","Klf12","Gata2","Klf13","Mycn","Maff","Creb5",
  "Mxi1","Prdm1","Crem","Tcf4","Klf4","Klf6","Foxp1","Jund","Klf2",
  "Hoxb8","Hoxa7","Hoxd8","Nr1h3","Cebpa","Hoxb6","Klf10","Pparg","Tbx3",
  "Stat3","Atf4","Hes1","Cebpd","Junb","Elf4","Klf11","Fosl2","Smad1",
  "Jun","Egr1"
)
 
# genes_reversed: ordered gene vector used as rownames (provided elsewhere)
gene_df <- setNames(
  as.data.frame(matrix(0L, nrow = length(genes_reversed), ncol = length(TF_order))),
  TF_order
)
rownames(gene_df) <- genes_reversed
 
# Fill matrix with 1 if (TF, gene) pair exists
for (i in seq_len(nrow(regulon_targets))) {
  tf   <- regulon_targets$TF[i]
  gene <- regulon_targets$gene[i]
  if (tf %in% TF_order && gene %in% genes_reversed) {
    gene_df[gene, tf] <- 1L
  }
}
 
heat_data <- gene_df
 
# TF sets for each organ (character vectors)
# Adipose_TFs, Heart_TFs, Muscle_TFs, Liver_TFs should be defined elsewhere
annotation_df <- data.frame(Group = NA_character_, row.names = colnames(heat_data))
annotation_df$Group[colnames(heat_data) %in% Adipose_TFs] <- "Adipose"
annotation_df$Group[colnames(heat_data) %in% Heart_TFs]   <- "Heart"
annotation_df$Group[colnames(heat_data) %in% Muscle_TFs]  <- "Muscle"
annotation_df$Group[colnames(heat_data) %in% Liver_TFs]   <- "Liver"
 
ann_colors <- list(
  Group = c(
    Adipose = "#F8766D",
    Heart   = "#7CAE00",
    Muscle  = "#00BFC4",
    Liver   = "#C77CFF"
  )
)
 
pheatmap::pheatmap(
  heat_data,
  cluster_rows      = FALSE,
  cluster_cols      = TRUE,
  col               = c("gray95", "orangered"),
  annotation_col    = annotation_df,
  annotation_colors = ann_colors,
  fontsize_row      = 4,
  fontsize_col      = 17,
  fontsize          = 20,
  cellwidth         = 14,
  angle_col         = 90,
  annotation_height = 2,
  clustering_distance_cols = "euclidean",
  clustering_method        = "ward.D"
)
 
 
############################
 
# SCENIC pipeline and RSS
############################
 
seurat_capillary <- readRDS(
  "/path/to/CapillarySinusoidal_subset_SeuratV5.rds"
)
 
exprMat <- seurat_capillary@assays$SCT@counts
exprMat <- as.matrix(exprMat)
 
cellInfo <- seurat_capillary@meta.data
saveRDS(cellInfo, file = "int/cellInfo_4organs.Rds")
 
colVars <- list(
  Tissue = c(
    Adipose = "#F8766D",
    Heart   = "#7CAE00",
    Muscle  = "#00BFC4",
    Liver   = "#C77CFF"
  )
)
colVars$Tissue <- colVars$Tissue[
  intersect(names(colVars$Tissue), cellInfo$Tissue)
]
saveRDS(colVars, file = "int/colVars_4organs.Rds")
 
org   <- "mgi"
dbDir <- "cisTarget_databases"
myDatasetTitle <- "FourCapillary_16groups_Adipose_aging"
 
data(defaultDbNames)
defaultDbNames$mgi[1] <- "mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"
defaultDbNames$mgi[2] <- "mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather"
 
dbs <- defaultDbNames[[org]]
data(list = "motifAnnotations_mgi_v9", package = "RcisTarget")
 
scenicOptions <- initializeScenic(
  org          = org,
  dbDir        = dbDir,
  dbs          = dbs,
  datasetTitle = myDatasetTitle,
  nCores       = 10
)
scenicOptions@settings$dbs <- dbs
 
genesKept <- geneFiltering(
  exprMat,
  scenicOptions    = scenicOptions,
  minCountsPerGene = 3 * 0.01 * ncol(exprMat),
  minSamples       = ncol(exprMat) * 0.01
)
 
exprMat_filtered <- exprMat[genesKept, ]
exprMat_filtered <- log2(exprMat_filtered + 1)
 
runGenie3(exprMat_filtered, scenicOptions, nTrees = 500)
 
# Reload (typical SCENIC workflow)
# scenicOptions <- readRDS("int/scenicOptions.Rds")
 
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores  <- 10
scenicOptions@settings$seed    <- 123
 
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
 
regulons <- loadInt(scenicOptions, "regulons")
 
regulonAUC_obj <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC_mtx <- getAUC(regulonAUC_obj)
 
AUCmat_extended <- regulonAUC_mtx[grep("_extended", rownames(regulonAUC_mtx)), ]
rownames(AUCmat_extended) <- gsub(" \\(\\d+g\\)", "", rownames(AUCmat_extended))
rownames(AUCmat_extended) <- gsub("_extended", "", rownames(AUCmat_extended))
 
regulonAUC <- AUCmat_extended
 
cellInfo_rss <- data.frame(seuratCluster = Idents(seurat_capillary))
rss <- calcRSS(
  AUC            = regulonAUC,
  cellAnnotation = cellInfo_rss[colnames(regulonAUC), , drop = FALSE]
)
rssPlot <- plotRSS(rss)
 
 
############################################################
```


## 3. TF activity visualization, aging DEGs, coexpression annotation, and signature scoring


```r
############################################################
 
## TF expression and regulon-activity scatter/bubble plots
############################
 
adipose_TF <- c(
  "Pparg","Atf4","Junb","Cebpd","Elf4","Fosl2","Tbx3",
  "Nr1h3","Hoxd8","Hoxb8","Hoxa7","Klf10","Klf11"
)
 
# Scatter: single gene example (Pparg)
DefaultAssay(seurat_all) <- "SCT"
exp_data <- FetchData(seurat_all, vars = c("Pparg", "Tissue"))
 
DefaultAssay(seurat_all) <- "regulonAUC_extended"
auc_data <- FetchData(seurat_all, vars = "Pparg")
 
scatter_data <- data.frame(
  Pparg_SCT                 = exp_data$Pparg,
  Pparg_regulonAUC_extended = auc_data$Pparg,
  Tissue                    = exp_data$Tissue
)
 
ggplot(scatter_data,
       aes(x = Pparg_SCT, y = Pparg_regulonAUC_extended, color = Tissue)) +
  geom_point(alpha = 0.5) +
  labs(x = "Pparg Expression (SCT)", y = "Pparg regulonAUC") +
  theme_minimal() +
  scale_color_manual(values = c(
    Adipose = "#F8766D",
    Heart   = "#7CAE00",
    Muscle  = "#00BFC4",
    Liver   = "#C77CFF"
  ))
 
# Bubble plot: avg expression vs avg regulonAUC across tissues
tissue_list <- list(
  Adipose = seurat_adipose,
  Heart   = seurat_heart,
  Muscle  = seurat_muscle,
  Liver   = seurat_liver
)
 
bubble_results <- list()
for (nm in names(tissue_list)) {
  obj_t <- tissue_list[[nm]]
 
  DefaultAssay(obj_t) <- "SCT"
  avg_expr <- sapply(adipose_TF, function(g) {
    mean(FetchData(obj_t, vars = g)[, g], na.rm = TRUE)
  })
 
  DefaultAssay(obj_t) <- "regulonAUC_extended"
  avg_reg <- sapply(adipose_TF, function(g) {
    mean(FetchData(obj_t, vars = g)[, g], na.rm = TRUE)
  })
 
  bubble_results[[nm]] <- data.frame(
    Gene              = adipose_TF,
    Average_Expression = avg_expr,
    Average_Regulon   = avg_reg,
    Tissue            = nm
  )
}
 
bubble_data <- do.call(rbind, bubble_results)
 
color_palette <- c(
  Adipose = "#F8766D",
  Heart   = "#7CAE00",
  Muscle  = "#00BFC4",
  Liver   = "#C77CFF"
)
 
ggplot(bubble_data,
       aes(x = Average_Expression, y = Average_Regulon,
           size = Average_Expression, label = Gene, fill = Tissue)) +
  geom_point(shape = 21, color = "black", alpha = 0.7) +
  scale_size(range = c(3, 10)) +
  scale_fill_manual(values = color_palette) +
  labs(
    x    = "Average Expression (SCT)",
    y    = "Average regulonAUC",
    size = "Expression"
  ) +
  geom_text(size = 3, vjust = -1) +
  theme_minimal() +
  theme(
    panel.background  = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.line         = element_line(color = "black")
  )
 
 
## Aging-associated DEGs and Venn diagram
############################
 
# Example: Adipose aging DEGs
deg_adipose <- FindMarkers(
  seurat_adipose,
  group.by        = "Age",
  ident.1         = "Aged",
  ident.2         = "Young",
  logfc.threshold = 0.25,
  recorrect_umi   = FALSE
)
 
deg_adipose_filtered <- deg_adipose[deg_adipose$avg_log2FC >= 0, ]
deg_adipose_sorted   <- deg_adipose_filtered[order(deg_adipose_filtered$p_val_adj), ]
top100_deg_adipose   <- head(deg_adipose_sorted, 100)
 
write.csv(
  top100_deg_adipose,
  file      = "top100_DEGs_Adipose.csv",
  row.names = TRUE
)
 
# top_genes_*: vectors of gene names from each organ
venn.plot <- venn.diagram(
  x = list(
    Adipose = top_genes_adipose,
    Heart   = top_genes_heart,
    Muscle  = top_genes_muscle,
    Liver   = top_genes_liver
  ),
  category.names = c("Adipose", "Heart", "Muscle", "Liver"),
  filename       = NULL,
  output         = TRUE,
  height         = 800,
  width          = 800,
  resolution     = 300,
  compression    = "lzw",
  lwd            = 2,
  col            = "black",
  fill           = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"),
  alpha          = 0.50,
  cex            = 2,
  fontface       = "bold",
  fontfamily     = "sans",
  cat.cex        = 1.5,
  cat.col        = "black",
  cat.pos        = 0,
  cat.dist       = 0.07,
  cat.fontface   = "bold",
  cat.default.pos  = "text",
  cat.default.dist = 0.035,
  cat.default.fontface = "plain"
)
grid::grid.draw(venn.plot)
 
 
## TF presence annotation on coexpression heatmap
############################
 
# average_matrix: coexpression matrix (e.g., CoexWeight)
ann_adipose <- create_annotations(top_TFs_adipose, average_matrix)
ann_heart   <- create_annotations(top_TFs_heart,   average_matrix)
ann_muscle  <- create_annotations(top_TFs_muscle,  average_matrix)
ann_liver   <- create_annotations(top_TFs_liver,   average_matrix)
 
annotations <- data.frame(
  Adipose = ann_adipose,
  Heart   = ann_heart,
  Muscle  = ann_muscle,
  Liver   = ann_liver
)
 
ann_colors_tf <- list(
  Adipose = c("gray92", "#F8766D"),
  Heart   = c("gray92", "#7CAE00"),
  Muscle  = c("gray92", "#00BFC4"),
  Liver   = c("gray92", "#C77CFF")
)
 
pheatmap::pheatmap(
  average_matrix,
  cluster_rows      = FALSE,
  cluster_cols      = TRUE,
  col               = my_palette,
  breaks            = seq(0, 0.02, length.out = 101),
  annotation_col    = annotations,
  annotation_colors = ann_colors_tf,
  fontsize_row      = 4,
  fontsize_col      = 18,
  fontsize          = 17,
  cellwidth         = 15,
  angle_col         = 90,
  annotation_height = 3,
  main              = "Aging_CoexWeight"
)
 
 
## Signature-score calculation and DotPlot
############################
 
# Example: apply signature scoring to integrated object seurat_all
normalized_matrix <- as.data.frame(seurat_all[["SCT"]]@data)
 
# List_1–4: gene sets for signature scoring
sig1 <- calculate_signature_score(normalized_matrix, List_1)
seurat_all <- AddMetaData(seurat_all, sig1, "Aging_1")
 
sig2 <- calculate_signature_score(normalized_matrix, List_2)
seurat_all <- AddMetaData(seurat_all, sig2, "Aging_2")
 
sig3 <- calculate_signature_score(normalized_matrix, List_3)
seurat_all <- AddMetaData(seurat_all, sig3, "Aging_3")
 
sig4 <- calculate_signature_score(normalized_matrix, List_4)
seurat_all <- AddMetaData(seurat_all, sig4, "Aging_4")
 
# DotPlot grouped by 'Type'
DotPlot(
  seurat_all,
  features = Lists,     # vector of signature columns or marker genes
  group.by = "Type",
  cols     = c("purple4", "orange")
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
 
 
############################################################
```


# Part II. Spatial transcriptomics analyses


## 1. R setup and sample definitions


```r
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(BPCells)
  library(Matrix)
  library(hdf5r)
  library(arrow)
  library(harmony)
  library(dplyr)
  library(ggplot2)
})
 
heart_outs <- list(
  Heart_Aged_1  = "/work/Visium/Visium_session1/Heart_Aged_1/outs",
  Heart_Young_1 = "/work/Visium/Visium_session1/Heart_Young_1/outs",
  Heart_Aged_2  = "/work/Visium/Visium_session2/Heart_Aged_2/outs",
  Heart_Young_2 = "/work/Visium/Visium_session2/Heart_Young_2/outs"
)
 
adipose_outs <- list(
  Adipose_Aged_1  = "/work/Visium/Visium_session1/Adipose_Aged_1/outs",
  Adipose_Young_1 = "/work/Visium/Visium_session1/Adipose_Young_1/outs",
  Adipose_Aged_2  = "/work/Visium/Visium_session2/Adipose_Aged_2/outs",
  Adipose_Young_2 = "/work/Visium/Visium_session2/Adipose_Young_2/outs"
)
```


## 2. Reference construction for spatial cell-type annotation


```r
# Tissue-specific references used for spatial cell-type annotation.
# Public mouse aging atlas datasets were obtained from CZ CELLxGENE Discover
# as h5ad files, including Tabula Muris Senis heart, adipose tissue,
# and gonadal fat pad datasets. Cell category columns were converted to strings
# before conversion/merging into Seurat objects.
 
# Heart reference:
#   Tabula Muris Senis heart dataset + present-study heart single-cell RNA-seq dataset
#   -> normalization/integration
#   -> RCTD reference object used by spacexr
reference_rctd_heart <- readRDS("reference_rctd_heart.rds")
 
# Adipose reference:
#   Tabula Muris Senis adipose tissue and gonadal fat pad datasets
#   + previously analyzed adipose single-nucleus RNA-seq dataset
#     deposited under BioProject accession PRJDB17809
#   -> normalization/integration
#   -> Seurat reference object used for anchor-based label transfer
adipose_reference <- readRDS("adipose_merged_obj_reference.rds")
 
# The heart reference was passed to create.RCTD()/run.RCTD().
# The adipose reference was passed to FindTransferAnchors()/TransferData().
```


## 3. Heart Visium HD 8-um bin processing in Seurat


```r
# Input: Space Ranger v4.0.1 binned output, square_008um/filtered_feature_bc_matrix.h5
h5 <- lapply(heart_outs, function(o) {
  file.path(o, "binned_outputs/square_008um/filtered_feature_bc_matrix.h5")
})
stopifnot(all(file.exists(unlist(h5))))
 
odir <- "/work/Visium/Heart_combined_samples/on_disk_mats"
dir.create(odir, recursive = TRUE, showWarnings = FALSE)
 
mat_list <- lapply(names(h5), function(nm) {
  m <- BPCells::open_matrix_10x_hdf5(path = h5[[nm]], feature_type = "Gene Expression")
  out <- file.path(odir, nm)
  BPCells::write_matrix_dir(mat = m, dir = out, overwrite = TRUE)
  BPCells::open_matrix_dir(out)
})
names(mat_list) <- names(h5)
 
obj_list <- lapply(names(mat_list), function(nm) {
  so <- CreateSeuratObject(counts = mat_list[[nm]], project = nm, assay = "RNA",
                           min.cells = 0, min.features = 0)
  so$sample <- nm
  RenameCells(so, add.cell.id = nm)
})
names(obj_list) <- names(mat_list)
 
heart <- Reduce(function(a, b) merge(a, b), obj_list)
heart <- subset(heart, subset = nCount_RNA >= 25)
 
DefaultAssay(heart) <- "RNA"
heart <- NormalizeData(heart, scale.factor = median(heart$nCount_RNA), verbose = FALSE)
heart <- FindVariableFeatures(heart, verbose = FALSE)
heart <- SketchData(heart, ncells = 50000, method = "LeverageScore", sketched.assay = "sketch")
 
DefaultAssay(heart) <- "sketch"
heart <- NormalizeData(heart, verbose = FALSE)
heart <- FindVariableFeatures(heart, verbose = FALSE)
heart <- ScaleData(heart, verbose = FALSE)
heart <- RunPCA(heart, npcs = 50, verbose = FALSE)
heart <- IntegrateLayers(heart, method = HarmonyIntegration, orig = "pca",
                         new.reduction = "harmony", dims = 1:30)
heart <- RunUMAP(heart, reduction = "harmony", dims = 1:30, reduction.name = "umap")
heart <- FindNeighbors(heart, reduction = "harmony", dims = 1:30)
heart <- FindClusters(heart, resolution = 0.5)
 
heart <- ProjectIntegration(heart, sketched.assay = "sketch", assay = "RNA", reduction = "harmony")
heart <- ProjectData(heart, assay = "RNA", full.reduction = "harmony.full",
                     sketched.assay = "sketch", sketched.reduction = "harmony.full",
                     umap.model = "umap", dims = 1:30,
                     refdata = list(cluster_full = "seurat_clusters"))
saveRDS(heart, "heart_mergeobj_preProject_full_sketchHarmony.rds")
```


## 4. Heart gene-symbol conversion and RCTD-based annotation


```r
convert_ensembl_to_symbol <- function(seurat_obj, assay = "RNA") {
  library(biomaRt)
  mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  ens_ids <- rownames(seurat_obj[[assay]])
  annot <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),
                 filters = "ensembl_gene_id", values = ens_ids, mart = mart)
  annot <- annot[annot$mgi_symbol != "", ]
  annot <- annot[!duplicated(annot$ensembl_gene_id), ]
  rownames(seurat_obj[[assay]]) <- make.unique(ifelse(
    ens_ids %in% annot$ensembl_gene_id,
    annot$mgi_symbol[match(ens_ids, annot$ensembl_gene_id)],
    ens_ids
  ))
  seurat_obj
}
 
heart <- readRDS("heart_mergeobj_preProject_full_sketchHarmony.rds")
heart <- convert_ensembl_to_symbol(heart, assay = "RNA")
saveRDS(heart, "heart_mergeobj_preProject_full_sketchHarmony_ENSMUSG16kept.rds")
 
# Cluster marker-based naming
DefaultAssay(heart) <- "RNA"
Idents(heart) <- heart$cluster_full
heart <- JoinLayers(heart)
markers <- FindAllMarkers(heart, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
 
heart$cluster_naming <- NA_character_
cl <- as.character(heart$cluster_full)
heart$cluster_naming[cl == "0"]  <- "Atrial cardiac muscle"
heart$cluster_naming[cl == "16"] <- "Platelet"
heart$cluster_naming[cl == "10"] <- "Red blood cell"
heart$cluster_naming[cl == "23"] <- "Airway epithelium"
heart$cluster_naming[cl == "17"] <- "Airway structural tissue"
heart$cluster_naming[cl == "26"] <- "Cardiac neuron"
heart$cluster_naming[cl == "8"]  <- "Endocardium"
heart$cluster_naming[cl == "12"] <- "Epicardium"
heart$cluster_naming[cl == "5"]  <- "Adipose tissue"
saveRDS(heart, "heart_mergeobj_preProject_full_sketchHarmony_ENSMUSG16kept_clusterNaming.rds")
 
# RCTD annotation per sample
library(spacexr)
library(org.Mm.eg.db)
reference_rctd <- readRDS("reference_rctd_heart.rds")
heart_raw <- readRDS("heart_mergeobj_preProject_full_sketchHarmony.rds")
heart_raw <- JoinLayers(heart_raw, assay = "RNA", layers = "counts", new = "counts_joined")
counts <- GetAssayData(heart_raw, assay = "RNA", layer = "counts_joined")
 
run_rctd_one_heart <- function(sample, pos_path, nCount_cutoff = 30) {
  cells <- colnames(heart_raw)[heart_raw$sample == sample & heart_raw$nCount_RNA >= nCount_cutoff]
  pos <- read_parquet(pos_path)
  coords_raw <- data.frame(x = pos$pxl_col_in_fullres, y = pos$pxl_row_in_fullres,
                           row.names = pos$barcode)
  cells_base <- sub("^.*_(s_008um_)", "\\1", cells)
  idx <- match(cells_base, rownames(coords_raw))
  stopifnot(all(!is.na(idx)))
  coords_sample <- coords_raw[idx, , drop = FALSE]
  rownames(coords_sample) <- cells
 
  counts_use <- counts[, cells, drop = FALSE]
  ens_keys <- unique(rownames(counts_use))
  map_df <- AnnotationDbi::select(org.Mm.eg.db, keys = ens_keys, keytype = "ENSEMBL",
                                  columns = "SYMBOL")
  map_df <- map_df[!is.na(map_df$SYMBOL) & map_df$SYMBOL != "", c("ENSEMBL", "SYMBOL")]
  map_df <- map_df[!duplicated(map_df$ENSEMBL), ]
  ens2sym <- setNames(map_df$SYMBOL, map_df$ENSEMBL)
  ens_in <- intersect(rownames(counts_use), names(ens2sym))
  M_ens <- counts_use[ens_in, , drop = FALSE]
  sym_ids <- ens2sym[ens_in]
  symbols <- unique(sym_ids)
  A <- sparseMatrix(i = match(sym_ids, symbols), j = seq_along(sym_ids), x = 1,
                    dims = c(length(symbols), length(sym_ids)))
  M_sym <- A %*% M_ens
  rownames(M_sym) <- symbols
 
  spatialRNA <- SpatialRNA(coords = coords_sample, counts = M_sym, nUMI = Matrix::colSums(M_sym))
  myRCTD <- create.RCTD(spatialRNA = spatialRNA, reference = reference_rctd,
                        max_cores = 8, keep_reference = TRUE)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")
  saveRDS(myRCTD, file = paste0("RCTD_fullref_cutoff30_", sample, ".rds"))
}
```


## 5. Heart merged cell labels, XY coordinates, and quality filtering


```r
# Merge RCTD predictions, treat rejected spots as Unknown, and attach to Seurat metadata.
paths <- c(
  Heart_Young_1 = "RCTD_fullref_cutoff30_Heart_Young_1.rds",
  Heart_Young_2 = "RCTD_fullref_cutoff30_Heart_Young_2.rds",
  Heart_Aged_1  = "RCTD_fullref_cutoff30_Heart_Aged_1.rds",
  Heart_Aged_2  = "RCTD_fullref_cutoff30_Heart_Aged_2.rds"
)
 
read_one_rctd <- function(sample, fp) {
  r <- readRDS(fp)
  cellid <- colnames(r@spatialRNA@counts)
  W <- r@results$weights
  pred_type <- colnames(W)[max.col(W, ties.method = "first")]
  spot_class <- as.character(r@results$results_df$spot_class)
  pred_type2 <- pred_type
  pred_type2[spot_class == "reject"] <- "Unknown"
  data.frame(cellid = cellid, sample = sample, spot_class = spot_class,
             rctd_type = pred_type, rctd_type_unknown = pred_type2,
             stringsAsFactors = FALSE)
}
 
rctd_merged <- do.call(rbind, Map(read_one_rctd, names(paths), paths))
stopifnot(!anyDuplicated(rctd_merged$cellid))
saveRDS(rctd_merged, "RCTD_merged_cutoff30_rejectUnknown.rds")
 
obj <- readRDS("heart_mergeobj_preProject_full_sketchHarmony_ENSMUSG16kept_clusterNaming.rds")
spotclass_vec <- setNames(rctd_merged$spot_class, rctd_merged$cellid)
type_vec <- setNames(rctd_merged$rctd_type_unknown, rctd_merged$cellid)
obj$rctd_spot_class <- unname(spotclass_vec[colnames(obj)])
obj$rctd_type <- unname(type_vec[colnames(obj)])
 
# Final labels combine cluster-based labels and RCTD labels, followed by manual curation.
# Valve-region fibroblasts were reassigned as Cardiac valve from manually selected valve cell IDs.
obj_use <- subset(obj, subset = !is.na(Merged_cellname))
obj_use$Merged_cellname <- factor(obj_use$Merged_cellname)
saveRDS(obj_use, "heart_mergeobj_preProject_full_sketchHarmony_ENSMUSG16kept_clusterNaming_RCTDmerged_MergedCellname_ASSIGNEDonly.rds")
 
# Per-sample Visium HD 8-um objects were filtered and XY coordinates were added from tissue_positions.parquet.
add_xy_heart <- function(obj, pos_path) {
  pos <- read_parquet(pos_path)
  coords_raw <- data.frame(x = pos$pxl_col_in_fullres, y = pos$pxl_row_in_fullres,
                           row.names = pos$barcode)
  cells <- colnames(obj)
  cells_base <- sub("^.*_(s_008um_)", "\\1", cells)
  idx <- match(cells_base, rownames(coords_raw))
  stopifnot(all(!is.na(idx)))
  obj$x <- coords_raw$x[idx]
  obj$y <- coords_raw$y[idx]
  obj
}
 
obj_hd_aged1_50 <- subset(obj_hd_aged1_lab, subset = nCount_Spatial.008um >= 50)
obj_hd_aged1_50 <- add_xy_heart(obj_hd_aged1_50,
  "Heart_Aged_1/outs/binned_outputs/square_008um/spatial/tissue_positions.parquet")
saveRDS(obj_hd_aged1_50, "Heart_Aged_1_visiumHD8um_ASSIGNEDonly_withXY.rds")
```


## 6. Adipose Visium HD cell-segmentation processing in Seurat


```r
# Input: Space Ranger v4.0.1 segmented output, segmented_outputs/filtered_feature_cell_matrix.h5
h5_seg <- lapply(adipose_outs, function(o) {
  file.path(o, "segmented_outputs/filtered_feature_cell_matrix.h5")
})
stopifnot(all(file.exists(unlist(h5_seg))))
 
odir <- "/work/Visium/Adipose_combined_samples/on_disk_mats"
dir.create(odir, recursive = TRUE, showWarnings = FALSE)
 
mat_list <- lapply(names(h5_seg), function(nm) {
  m <- BPCells::open_matrix_10x_hdf5(path = h5_seg[[nm]], feature_type = "Gene Expression")
  out <- file.path(odir, nm)
  BPCells::write_matrix_dir(mat = m, dir = out, overwrite = TRUE)
  BPCells::open_matrix_dir(out)
})
names(mat_list) <- names(h5_seg)
 
obj_list <- lapply(names(mat_list), function(nm) {
  so <- CreateSeuratObject(counts = mat_list[[nm]], project = nm, assay = "RNA",
                           min.cells = 0, min.features = 0)
  so$sample <- nm
  RenameCells(so, add.cell.id = nm)
})
names(obj_list) <- names(mat_list)
 
adipose <- Reduce(function(a, b) merge(a, b), obj_list)
adipose <- subset(adipose, subset = nCount_RNA >= 25)
DefaultAssay(adipose) <- "RNA"
adipose <- NormalizeData(adipose, scale.factor = median(adipose$nCount_RNA), verbose = FALSE)
adipose <- FindVariableFeatures(adipose, verbose = FALSE)
adipose <- SketchData(adipose, ncells = 50000, method = "LeverageScore", sketched.assay = "sketch")
 
DefaultAssay(adipose) <- "sketch"
adipose <- NormalizeData(adipose, verbose = FALSE)
adipose <- FindVariableFeatures(adipose, verbose = FALSE)
adipose <- ScaleData(adipose, verbose = FALSE)
adipose <- RunPCA(adipose, npcs = 50, verbose = FALSE)
adipose <- IntegrateLayers(adipose, method = HarmonyIntegration, orig = "pca",
                           new.reduction = "harmony", dims = 1:30)
adipose <- RunUMAP(adipose, reduction = "harmony", dims = 1:30, reduction.name = "umap")
adipose <- FindNeighbors(adipose, reduction = "harmony", dims = 1:30)
adipose <- FindClusters(adipose, resolution = 0.5)
saveRDS(adipose, "adipose_mergeobj_preProject_sketchHarmony.rds")
 
adipose <- ProjectIntegration(adipose, sketched.assay = "sketch", assay = "RNA", reduction = "harmony")
adipose <- ProjectData(adipose, assay = "RNA", full.reduction = "harmony.full",
                       sketched.assay = "sketch", sketched.reduction = "harmony.full",
                       umap.model = "umap", dims = 1:30,
                       refdata = list(cluster_full = "seurat_clusters"))
saveRDS(adipose, "adipose_mergeobj_preProject_full_sketchHarmony.rds")
```


## 7. Adipose label transfer and curated merged cell labels


```r
# Label transfer from adipose scRNA-seq reference
adipose <- readRDS("adipose_mergeobj_preProject_sketchHarmony.fixed_onDiskPath_joinCounts.rds")
ref_obj <- readRDS("adipose_merged_obj_reference.rds")
DefaultAssay(ref_obj) <- "RNA"
DefaultAssay(adipose) <- "RNA"
 
ref_obj <- FindVariableFeatures(ref_obj, verbose = FALSE)
adipose <- FindVariableFeatures(adipose, verbose = FALSE)
features_use <- intersect(VariableFeatures(ref_obj), VariableFeatures(adipose))
features_use <- head(features_use, 500)
 
anchors <- FindTransferAnchors(reference = ref_obj, query = adipose,
                               features = features_use, dims = 1:30)
pred <- TransferData(anchorset = anchors, refdata = ref_obj$Cellname, dims = 1:30)
adipose <- AddMetaData(adipose, metadata = pred)
adipose$Cellname_final <- ifelse(adipose$prediction.score.max < 0.4,
                                 "Unknown", adipose$predicted.id)
saveRDS(adipose, "adipose_mergeobj_preProject_sketchHarmony_labelTransferred.rds")
 
# Ambiguous labels were reclassified using Harmony SNN reclustering.
DefaultAssay(adipose) <- "SCT"
adipose <- FindNeighbors(adipose, reduction = "harmony", dims = 1:30,
                         k.param = 50, graph.name = "harmony_snn_k50")
adipose <- FindClusters(adipose, graph.name = "harmony_snn_k50", resolution = 0.2,
                        random.seed = 1)
adipose$clusters_k50_r020 <- adipose[["harmony_snn_k50_res.0.2"]][, 1]
 
# Assign broad cluster identities used to resolve mixed/unknown labels.
cl <- as.character(adipose$clusters_k50_r020)
new_name <- rep("Unassigned", length(cl))
new_name[cl == "0"] <- "Adipocyte (mature; lipid-signal mixed)"
new_name[cl == "1"] <- "Endothelial (capillary~general)"
new_name[cl == "2"] <- "Adipose stromal / Fibroblast"
new_name[cl == "3"] <- "Myeloid (Macrophage/DC; MHCII+; incl M2-like)"
new_name[cl == "4"] <- "Vascular mural (VSMC/Pericyte) + Endo mixed (Endo-mural mixed; suspicious)"
new_name[cl == "5"] <- "B cell"
new_name[cl == "8"] <- "Erythroid"
new_name[cl %in% as.character(10:29)] <- "Adipocyte-like / lipid-signal dominated"
adipose$cluster_based_naming <- new_name
 
# Cluster-based curation of mixed/unknown categories.
adipose$Cellname_partial <- as.character(adipose$Cellname_final)
adipose$Cellname_partial[adipose$Cellname_partial == "Adjacent parenchymal cells"] <- "Unknown"
adipose$Cellname_partial[adipose$Cellname_partial == "Endothelial cell"] <- "Endothelial_mixed"
adipose$Cellname_partial[adipose$Cellname_partial == "Macrophage/Monocyte/DC"] <- "Myeloid_mixed"
 
# Final merged categories used downstream.
adipose$Merged_cellname <- as.character(adipose$Cellname_partial)
map_cb_to_merged <- c(
  "Adipocyte (mature; lipid-signal mixed)" = "Adipocyte",
  "Adipocyte-like / lipid-signal dominated" = "Adipocyte",
  "Adipose stromal / Fibroblast" = "Adipose stromal cell",
  "B cell" = "B cell",
  "Endothelial (capillary~general)" = "Endothelial cell",
  "Erythroid" = "Unknown",
  "Myeloid (Macrophage/DC; MHCII+; incl M2-like)" = "Macrophage / Monocyte / DC",
  "Vascular mural (VSMC/Pericyte) + Endo mixed (Endo-mural mixed; suspicious)" = "VSMC/Pericyte"
)
idx_change <- as.character(adipose$Cellname_partial) %in% c("Endothelial_mixed", "Myeloid_mixed", "Unknown")
cb <- as.character(adipose$cluster_based_naming[idx_change])
new_label <- unname(map_cb_to_merged[cb])
new_label[is.na(new_label) | new_label == ""] <- "Unknown"
adipose$Merged_cellname[idx_change] <- new_label
saveRDS(adipose, "adipose_mergeobj_sketchHarmony_SCT_labelTransferred_k50_r020_Naming_260220.rds")
```


## 8. Adipose segmentation geometry and centroid coordinates


```r
suppressPackageStartupMessages({
  library(sf)
  library(jsonlite)
  library(stringr)
  library(tibble)
  library(dplyr)
})
 
add_segmentation_centroids <- function(seu, geojson_path, id_col = "cell_id") {
  cells_sf <- st_read(geojson_path, quiet = TRUE)
  poly_geoms <- st_geometry(cells_sf)
  invalid_idx <- which(!st_is_valid(poly_geoms))
  if (length(invalid_idx) > 0) {
    cells_sf[invalid_idx, ] <- st_make_valid(cells_sf[invalid_idx, ])
    poly_geoms <- st_geometry(cells_sf)
  }
 
  calc_area_px <- function(g) as.numeric(st_area(g))
  cells_sf$seg_area_px <- vapply(poly_geoms, calc_area_px, numeric(1))
  cents <- st_centroid(poly_geoms)
  coords <- st_coordinates(cents)
  cells_sf$centroid_x <- coords[, 1]
  cells_sf$centroid_y <- coords[, 2]
 
  barcodes <- colnames(seu)
  id_num <- as.integer(str_extract(barcodes, "(?<=_)\\d+(?=-)"))
  stopifnot(sum(is.na(id_num)) == 0)
  map_df <- tibble(barcode = barcodes, !!id_col := id_num)
  cell_meta <- cells_sf |>
    st_drop_geometry() |>
    dplyr::select(!!id_col, seg_area_px, centroid_x, centroid_y)
  cell_meta_bar <- inner_join(map_df, cell_meta, by = id_col)
  cell_meta_named <- column_to_rownames(cell_meta_bar, "barcode")
  cell_meta_named <- as.data.frame(cell_meta_named)[barcodes, , drop = FALSE]
  stopifnot(identical(rownames(cell_meta_named), barcodes))
  AddMetaData(seu, metadata = cell_meta_named)
}
 
# Example; repeated for all four adipose samples.
obj_list[["Adipose_Young_1"]] <- add_segmentation_centroids(
  obj_list[["Adipose_Young_1"]],
  "GeoJason/Adipose_Young_1/cell_segmentations.geojson"
)
saveRDS(obj_list[["Adipose_Young_1"]], "Adipose_Young_1_visiumHDsementation_ASSIGNEDonly_withXY.rds")
saveRDS(obj_list[["Adipose_Young_2"]], "Adipose_Young_2_visiumHDsementation_ASSIGNEDonly_withXY.rds")
saveRDS(obj_list[["Adipose_Aged_1"]],  "Adipose_Aged_1_visiumHDsementation_ASSIGNEDonly_withXY.rds")
saveRDS(obj_list[["Adipose_Aged_2"]],  "Adipose_Aged_2_visiumHDsementation_ASSIGNEDonly_withXY.rds")
```


## 9. Aging_1 module score and spatial neighborhood definition


```r
gene_list <- c(
  "Ifit3", "Ifi44", "Isg15", "Ifit1", "Oasl2", "Ifitm3", "Bst2",
  "Ifi204", "Psmb8", "Lgals3bp", "Oasl1", "Irgm2", "Herc6", "Usp18",
  "Irgm1", "Trim30a", "Gbp3", "Mndal", "Irf7", "Xaf1", "Lgals9",
  "Phf11d", "Ifit2", "Stat1", "Gbp7", "Iigp1", "Ifi203", "Cmpk2",
  "Rsad2", "Ifit3b", "Rtp4"
)
 
# Heart: score and 100-um neighborhood in LV ROI
library(RANN)
heart_roi <- readRDS("rds_out/Heart_Aged_1_LV_ROI_norm.rds")
DefaultAssay(heart_roi) <- "Spatial.008um"
gene_use <- intersect(gene_list, rownames(heart_roi))
heart_roi <- AddModuleScore(heart_roi, features = list(gene_use), name = "Aging_1_Score")
 
thr <- 0.3
radius_um <- 100
bin_um <- 8
score_col <- "Aging_1_Score1"
seed_cells <- colnames(heart_roi)[heart_roi@meta.data[[score_col]] >= thr]
tc <- GetTissueCoordinates(heart_roi)[, c("x", "y"), drop = FALSE]
coords_all <- as.matrix(tc)
coords_seed <- coords_all[seed_cells, , drop = FALSE]
spot_px <- heart_roi@images[[1]]@scale.factors$spot
radius_px <- (radius_um / bin_um) * spot_px
nn <- RANN::nn2(data = coords_seed, query = coords_all, k = 1)
d1 <- nn$nn.dists[, 1]
names(d1) <- rownames(tc)
neighbor_cells <- setdiff(names(d1)[d1 <= radius_px], seed_cells)
heart_roi$Aging1_nbhd100um <- "other"
heart_roi$Aging1_nbhd100um[neighbor_cells] <- "neighbor"
heart_roi$Aging1_nbhd100um[seed_cells] <- "seed"
heart_roi$Aging1_nbhd100um <- factor(heart_roi$Aging1_nbhd100um,
                                     levels = c("other", "neighbor", "seed"))
saveRDS(heart_roi, "rds_out/Heart_Aged_1_LV_ROI_norm_Aging1Score.rds")
 
# Adipose: score and 100-um neighborhood by segmentation centroids
adipose_a1 <- readRDS("Adipose_Aged_1_visiumHDsementation_ASSIGNEDonly_withXY.rds")
DefaultAssay(adipose_a1) <- "SCT"
adipose_a1 <- AddModuleScore(adipose_a1, features = list(gene_list), name = "Aging_1_Score")
 
thr <- 0.05
radius_um <- 100
seed_cells <- colnames(adipose_a1)[adipose_a1@meta.data[[score_col]] >= thr]
md <- adipose_a1@meta.data
coords_all <- as.matrix(md[, c("centroid_x", "centroid_y")])
rownames(coords_all) <- rownames(md)
coords_seed <- coords_all[seed_cells, , drop = FALSE]
# Pixel-to-micron conversion estimated from the 6.5-mm capture area.
px_per_um <- mean(c(diff(range(coords_all[,1], na.rm = TRUE)),
                    diff(range(coords_all[,2], na.rm = TRUE)))) / 6500
radius_px <- radius_um * px_per_um
nn <- RANN::nn2(data = coords_seed, query = coords_all, k = 1)
d1 <- nn$nn.dists[, 1]
names(d1) <- rownames(md)
neighbor_cells <- setdiff(names(d1)[d1 <= radius_px], seed_cells)
adipose_a1$Aging1_nbhd100um <- "other"
adipose_a1$Aging1_nbhd100um[neighbor_cells] <- "neighbor"
adipose_a1$Aging1_nbhd100um[seed_cells] <- "seed"
adipose_a1$Aging1_nbhd100um <- factor(adipose_a1$Aging1_nbhd100um,
                                      levels = c("other", "neighbor", "seed"))
saveRDS(adipose_a1, "rds_out/Adipose_Aged_1_norm_Aging1Score.rds")
```


## 10. LIANA export from Seurat objects


```r
# Heart LIANA export: Core, Neighbor_<cell type>, and Distal_<cell type>
export_liana_heart <- function(obj, outdir) {
  DefaultAssay(obj) <- "Spatial.008um"
  obj <- JoinLayers(obj, assay = "Spatial.008um", layers = "data")
  keep_ct <- c("Ventricular cardiac muscle", "B cell lineage", "Cardiac neuron", "Endocardium",
               "Endothelial cell", "Epicardium", "Fibroblast", "Lymphatic endothelial cell",
               "Macrophage / Monocyte / DC", "Mast cell", "Platelet", "T cell lineage",
               "VSMC / Pericyte")
  obj_use <- subset(obj, subset = Merged_cellname %in% keep_ct)
  obj_use$Aging1_role <- dplyr::case_when(
    obj_use$Aging1_nbhd100um == "seed" ~ "Core",
    obj_use$Aging1_nbhd100um == "neighbor" ~ "Neighbor",
    obj_use$Aging1_nbhd100um == "other" ~ "Distal",
    TRUE ~ NA_character_
  )
  obj_use$CCI_group <- dplyr::case_when(
    obj_use$Aging1_role == "Core" ~ "Core",
    obj_use$Aging1_role == "Neighbor" ~ paste0("Neighbor_", obj_use$Merged_cellname),
    obj_use$Aging1_role == "Distal" ~ paste0("Distal_", obj_use$Merged_cellname),
    TRUE ~ NA_character_
  )
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  X <- GetAssayData(obj_use, assay = "Spatial.008um", layer = "data")
  stopifnot(identical(colnames(X), colnames(obj_use)))
  Matrix::writeMM(X, file = file.path(outdir, "X.mtx"))
  writeLines(rownames(X), con = file.path(outdir, "genes.txt"))
  writeLines(colnames(X), con = file.path(outdir, "cells.txt"))
  obs <- data.frame(cell_id = colnames(obj_use), CCI_group = obj_use$CCI_group,
                    stringsAsFactors = FALSE)
  write.csv(obs, file = file.path(outdir, "obs.csv"), row.names = FALSE, quote = FALSE)
}
 
# Adipose LIANA export: Core, Neighbor_<cell type>, and Distal_<cell type>
export_liana_adipose <- function(obj, outdir) {
  assay_use <- "SCT"
  DefaultAssay(obj) <- assay_use
  obj <- tryCatch(JoinLayers(obj, assay = assay_use, layers = "data"), error = function(e) obj)
  exclude_ct <- c("Erythroid cell", "Unknown")
  obj$CCI_group <- NA_character_
  obj$CCI_group[obj$Aging1_nbhd100um == "seed"] <- "Core"
  obj$CCI_group[obj$Aging1_nbhd100um == "neighbor"] <- paste0("Neighbor_", obj$Merged_cellname[obj$Aging1_nbhd100um == "neighbor"])
  obj$CCI_group[obj$Aging1_nbhd100um == "other"] <- paste0("Distal_", obj$Merged_cellname[obj$Aging1_nbhd100um == "other"])
  obj_use <- subset(obj, subset = !is.na(CCI_group) &
                    (CCI_group == "Core" | !(Merged_cellname %in% exclude_ct)))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  X <- GetAssayData(obj_use, assay = assay_use, layer = "data")
  Matrix::writeMM(X, file = file.path(outdir, "X.mtx"))
  writeLines(rownames(X), con = file.path(outdir, "genes.txt"))
  writeLines(colnames(X), con = file.path(outdir, "cells.txt"))
  obs <- obj_use@meta.data[, c("CCI_group"), drop = FALSE]
  obs$cell_id <- rownames(obs)
  obs <- obs[, c("cell_id", "CCI_group")]
  write.csv(obs, file = file.path(outdir, "obs.csv"), row.names = FALSE, quote = FALSE)
}
```


## 11. Python LIANA rank aggregation and Neighbor-versus-Distal filtering


```python
import os
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io as sio
import liana as li
 
 
def load_export(outdir: str) -> sc.AnnData:
    X = sio.mmread(f"{outdir}/X.mtx").tocsr().T.tocsr()  # cells x genes
    genes = pd.read_csv(f"{outdir}/genes.txt", header=None)[0].tolist()
    cells = pd.read_csv(f"{outdir}/cells.txt", header=None)[0].tolist()
    obs = pd.read_csv(f"{outdir}/obs.csv").set_index("cell_id").loc[cells]
    adata = sc.AnnData(X=X)
    adata.var_names = genes
    adata.obs_names = cells
    adata.obs = obs
    return adata
 
 
def neighbor_distal_keep(res_df: pd.DataFrame) -> pd.DataFrame:
    df = res_df.copy()
    df["lr"] = df["ligand_complex"].astype(str) + "→" + df["receptor_complex"].astype(str)
    nb = df[(df["target"] == "Core") & (df["source"].str.startswith("Neighbor_"))].copy()
    ds = df[(df["target"] == "Core") & (df["source"].str.startswith("Distal_"))].copy()
    nb = nb[~nb["source"].isin({"Neighbor_Erythroid cell", "Neighbor_Unknown"})].copy()
    ds = ds[~ds["source"].isin({"Distal_Erythroid cell", "Distal_Unknown"})].copy()
    nb["ct"] = nb["source"].str.replace("^Neighbor_", "", regex=True)
    ds["ct"] = ds["source"].str.replace("^Distal_", "", regex=True)
    nb_mean = nb.groupby(["ct", "lr"], as_index=False)["lrscore"].mean().rename(columns={"lrscore": "mean_neighbor"})
    ds_mean = ds.groupby(["ct", "lr"], as_index=False)["lrscore"].mean().rename(columns={"lrscore": "mean_distal"})
    m = nb_mean.merge(ds_mean, on=["ct", "lr"], how="left")
    m["mean_distal"] = m["mean_distal"].fillna(0.0)
    m["delta"] = m["mean_neighbor"] - m["mean_distal"]
    return m[m["delta"] > 0].copy()
 
 
def run_liana_one(tag: str, outdir: str) -> pd.DataFrame:
    adata = load_export(outdir)
    li.mt.rank_aggregate(
        adata,
        groupby="CCI_group",
        resource_name="mouseconsensus",
        use_raw=False,
        verbose=True,
    )
    res = adata.uns["liana_res"].copy()
    keep = neighbor_distal_keep(res)
    keep["sample"] = tag
    keep.to_csv(f"{tag}_liana_neighbor_gt_distal.csv", index=False)
    return keep
 
exports = {
    "Heart_Aged1":  "Aged1_LIANA_export",
    "Heart_Aged2":  "Aged2_LIANA_export",
    "Heart_Young1": "Young1_LIANA_export",
    "Heart_Young2": "Young2_LIANA_export",
}
 
heart_liana = pd.concat([run_liana_one(tag, outdir) for tag, outdir in exports.items()], ignore_index=True)
heart_liana.to_csv("Heart_LIANA_neighbor_gt_distal_all_samples.csv", index=False)
 
exports = {
    "Adipose_Aged1":  "Adipose_Aged1_LIANA_export_SCT_CoreNeighborDistal",
    "Adipose_Aged2":  "Adipose_Aged2_LIANA_export_SCT_CoreNeighborDistal",
    "Adipose_Young1": "Adipose_Young1_LIANA_export_SCT_CoreNeighborDistal",
    "Adipose_Young2": "Adipose_Young2_LIANA_export_SCT_CoreNeighborDistal",
}
 
adipose_liana = pd.concat([run_liana_one(tag, outdir) for tag, outdir in exports.items()], ignore_index=True)
adipose_liana.to_csv("Adipose_LIANA_neighbor_gt_distal_all_samples.csv", index=False)
```


## 12. Cross-sample LIANA summary for age-group comparison


```python
def combine_two(k1: pd.DataFrame, k2: pd.DataFrame, prefix: str) -> pd.DataFrame:
    x = k1[["ct", "lr", "mean_neighbor"]].rename(columns={"mean_neighbor": f"{prefix}_1"})
    y = k2[["ct", "lr", "mean_neighbor"]].rename(columns={"mean_neighbor": f"{prefix}_2"})
    out = x.merge(y, on=["ct", "lr"], how="outer").fillna(0.0)
    out[f"{prefix}_present_n"] = (out[[f"{prefix}_1", f"{prefix}_2"]] > 0).sum(axis=1)
    out[f"{prefix}_mean"] = out[[f"{prefix}_1", f"{prefix}_2"]].mean(axis=1)
    return out
 
# Example framework after loading per-sample Neighbor>Distal tables.
# Retain ligand-receptor pairs reproducibly observed in either two Aged samples or two Young samples.
aged = combine_two(k_aged1, k_aged2, "Aged")
young = combine_two(k_young1, k_young2, "Young")
combined = aged.merge(young, on=["ct", "lr"], how="outer").fillna(0.0)
combined["max_mean"] = combined[["Aged_mean", "Young_mean"]].max(axis=1)
combined = combined[(combined["Aged_present_n"] >= 2) | (combined["Young_present_n"] >= 2)].copy()
combined.to_csv("LIANA_age_group_summary_neighbor_gt_distal.csv", index=False)
```


# Part III. Other RNA-seq-related and transcriptome analyses


## 1. Regulon ratios, IFN signatures, Reactome/scGSVA, RSS comparisons, and Junb-associated DEGs


```r
############################################################
 
############################
 
# Regulon mean ratios (Aged/Young vs HFD/Normal) per tissue
############################
 
colors_tissue <- c(
  Adipose = "#F8766D",
  Heart   = "#7CAE00",
  Muscle  = "#00BFC4",
  Liver   = "#C77CFF"
)
 
calculate_mean_ratios <- function(aging_obj, diet_obj, regulons, tissue_name) {
  aging_obj$Age <- factor(aging_obj$Age, levels = c("Young", "Aged"))
  diet_obj$Diet <- factor(diet_obj$Diet, levels = c("Normal", "HFD"))
 
  young_cells <- WhichCells(aging_obj, expression = Age == "Young")
  aged_cells  <- WhichCells(aging_obj, expression = Age == "Aged")
 
  normal_cells <- WhichCells(diet_obj, expression = Diet == "Normal")
  hfd_cells    <- WhichCells(diet_obj, expression = Diet == "HFD")
 
  mean_ratios_aging <- data.frame(Gene = regulons, Mean_Ratio = NA_real_)
  mean_ratios_diet  <- data.frame(Gene = regulons, Mean_Ratio = NA_real_)
 
  for (g in regulons) {
    young_mean <- mean(aging_obj@assays$regulonAUC_extended@data[g, young_cells])
    aged_mean  <- mean(aging_obj@assays$regulonAUC_extended@data[g, aged_cells])
    mean_ratios_aging$Mean_Ratio[mean_ratios_aging$Gene == g] <- aged_mean / young_mean
 
    normal_mean <- mean(diet_obj@assays$regulonAUC_extended@data[g, normal_cells])
    hfd_mean    <- mean(diet_obj@assays$regulonAUC_extended@data[g, hfd_cells])
    mean_ratios_diet$Mean_Ratio[mean_ratios_diet$Gene == g]  <- hfd_mean / normal_mean
  }
 
  data.frame(
    Gene             = regulons,
    Mean_Ratio_Aging = mean_ratios_aging$Mean_Ratio,
    Mean_Ratio_Diet  = mean_ratios_diet$Mean_Ratio,
    Tissue           = tissue_name
  )
}
 
all_results <- list(
  calculate_mean_ratios(adipose_aging, adipose_diet, Adipose_Regulons, "Adipose"),
  calculate_mean_ratios(heart_aging,   heart_diet,   Heart_Regulons,   "Heart"),
  calculate_mean_ratios(muscle_aging,  muscle_diet,  Muscle_Regulons,  "Muscle"),
  calculate_mean_ratios(liver_aging,   liver_diet,   Liver_Regulons,   "Liver")
)
 
final_plot_data <- do.call(rbind, all_results)
 
ggplot(final_plot_data,
       aes(x = Mean_Ratio_Aging, y = Mean_Ratio_Diet,
           color = Tissue, label = Gene)) +
  geom_point() +
  geom_text(vjust = -0.5, hjust = 1) +
  scale_color_manual(values = colors_tissue) +
  labs(
    title = "Mean Regulon Activity Ratios: Aging vs Diet",
    x     = "Aged / Young (regulonAUC)",
    y     = "HFD / Normal (regulonAUC)"
  ) +
  theme_minimal()
 
 
############################
 
# Regulon–regulon correlation heatmap (regulonAUC_extended)
############################
 
correlations_DEGs_log <- cor(
  method = "spearman",
  log2(
    t(as.matrix(seurat_all@assays[["regulonAUC_extended"]]@data[Regulons, ])) + 1
  )
)
 
make_spearman_heatmap(
  correlations_DEGs_log,
  palette_low   = "navyblue",
  palette_mid   = "white",
  palette_high  = "red",
  breaks_q_low  = 0.17,
  breaks_q_high = 1.00,
  step          = 0.06,
  cex_row       = 1,
  cex_col       = 1,
  dendrogram    = "both"
)
 
 
############################
 
# IFN signature score & ridge plot (Adipose example)
############################
 
# IFN$Symbol_Merged: vector of IFN-related genes (length 65)
IFN_genes <- head(IFN$Symbol_Merged, n = 65)
 
seurat_adipose_tmp <- seurat_adipose
normalized_matrix_ifn <- as.data.frame(seurat_adipose_tmp[["RNA"]]@data)
 
genes.score <- calculate_signature_score(normalized_matrix_ifn, IFN_genes)
seurat_adipose_tmp <- AddMetaData(seurat_adipose_tmp, genes.score, "IFN_score")
 
plot_data <- as.data.frame(seurat_adipose_tmp@meta.data)
plot_data$IFN_score <- genes.score[rownames(plot_data)]
 
# Reverse Condition order for plotting
unique_condition <- unique(plot_data$Condition)
plot_data$Condition <- factor(plot_data$Condition, levels = rev(unique_condition))
 
options(repr.plot.width = 6, repr.plot.height = 5)
ggplot(plot_data, aes(x = IFN_score, y = Vessel_Type, fill = Condition)) +
  geom_density_ridges(alpha = 0.2, scale = 1) +
  scale_fill_viridis_d(direction = -1) +
  geom_point(aes(color = Condition),
             position = position_dodge(width = 0.75),
             alpha    = 0.5) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border     = element_blank(),
    axis.text        = element_text(size = 12),
    axis.title       = element_text(size = 14),
    legend.text      = element_text(size = 12),
    legend.title     = element_text(size = 14),
    plot.title       = element_text(size = 16),
    strip.text       = element_text(size = 10),
    strip.background = element_rect(fill = "gray90")
  ) +
  labs(y = "Vessel_Type", x = "IFN Signature Score")
 
# t-test per Vessel_Type (Condition groups)
t_test_results <- lapply(unique(seurat_adipose_tmp@meta.data$Vessel_Type), function(category) {
  subset_data <- subset(seurat_adipose_tmp@meta.data, Vessel_Type == category)
  t_test_result <- t.test(genes.score ~ Condition, data = subset_data)
  return(t_test_result$p.value)
})
names(t_test_results) <- unique(seurat_adipose_tmp@meta.data$Vessel_Type)
 
 
############################
 
# Artery vs Capillary RSS: gap from tissue-wise max
############################
 
CapRSS    <- read.csv("/path/to/rss_Aging_Adipose_Capillary.csv")
ArteryRSS <- read.csv("/path/to/rss_Aging_Adipose_Artery.csv")
 
# Example for Aged column; generalized by function
make_plotdf_Aged <- function(CapRSS, ArteryRSS, tissue_label, topn = 100) {
  top_cap <- CapRSS  %>% arrange(desc(Aged)) %>% slice_head(n = topn) %>% pull(X)
  top_art <- ArteryRSS %>% arrange(desc(Aged)) %>% slice_head(n = topn) %>% pull(X)
  genes_to_plot <- union(top_cap, top_art)
 
  cap_vals <- CapRSS  %>%
    filter(X %in% genes_to_plot) %>%
    select(X, Aged) %>%
    rename(Aged_CapRSS = Aged)
 
  art_vals <- ArteryRSS %>%
    filter(X %in% genes_to_plot) %>%
    select(X, Aged) %>%
    rename(Aged_ArteryRSS = Aged)
 
  merge(cap_vals, art_vals, by = "X") %>%
    mutate(Tissue = tissue_label)
}
 
plot_data_all <- dplyr::bind_rows(
  make_plotdf_Aged(AdiposeCapRSS, AdiposeArteryRSS, "Adipose"),
  make_plotdf_Aged(HeartCapRSS,   HeartArteryRSS,   "Heart"),
  make_plotdf_Aged(MuscleCapRSS,  MuscleArteryRSS,  "Muscle")
)
 
plot_gap <- plot_data_all %>%
  group_by(Tissue) %>%
  mutate(
    max_art   = max(Aged_ArteryRSS, na.rm = TRUE),
    max_cap   = max(Aged_CapRSS,    na.rm = TRUE),
    gap_Artery = Aged_ArteryRSS - max_art,
    gap_Cap    = Aged_CapRSS    - max_cap,
    rank_art   = min_rank(desc(Aged_ArteryRSS)),
    rank_cap   = min_rank(desc(Aged_CapRSS)),
    is_top5    = rank_art <= 3 | rank_cap <= 3,
    is_low     = gap_Artery < -0.1 | gap_Cap < -0.1,
    label_gene = ifelse(is_top5 | is_low, X, NA_character_)
  ) %>%
  ungroup() %>%
  mutate(Tissue = factor(Tissue, levels = c("Heart", "Muscle", "Adipose"))) %>%
  arrange(Tissue)
 
min_gap <- min(plot_gap$gap_Artery, plot_gap$gap_Cap, na.rm = TRUE)
 
tissue_cols <- c(
  Adipose = "#F8766D",
  Heart   = "#7CAE00",
  Muscle  = "#00BFC4"
)
 
ggplot(plot_gap, aes(x = gap_Artery, y = gap_Cap)) +
  geom_point(aes(fill = Tissue),
             shape = 21, size = 1.8, stroke = 0.25,
             color = "grey20", alpha = 0.85) +
  geom_text_repel(
    data = subset(plot_gap, !is.na(label_gene)),
    aes(label = label_gene, color = Tissue),
    size = 4,
    box.padding = 0.25,
    point.padding = 0.25,
    max.overlaps = Inf,
    min.segment.length = 0,
    seed = 1,
    show.legend = FALSE
  ) +
  geom_hline(yintercept = 0, color = "black") +
  geom_vline(xintercept = 0, color = "black") +
  scale_fill_manual(values = tissue_cols, name = "Tissue") +
  scale_color_manual(values = tissue_cols, guide = "none") +
  scale_x_continuous(limits = c(min_gap, 0.05), expand = c(0, 0)) +
  scale_y_continuous(limits = c(min_gap, 0.05), expand = c(0, 0)) +
  coord_fixed() +
  labs(
    title = "Gap from Tissue-wise Max RSS (Aged Top100 unions per tissue)",
    x     = "Artery: RSS − max(RSS) within tissue",
    y     = "Capillary: RSS − max(RSS) within tissue"
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))
 
 
############################
 
# Reactome (UCell) × TF regulon correlation & network
############################
 
# Reactome annotation with UCell
# (Run once; result saved and reloaded)
# BiocManager::install("UCell")
library(UCell)
 
msrt <- buildAnnot(species = "mouse", keytype = "SYMBOL", anntype = "Reactome")
# duplicate column 2 to column 3 as in original code (if needed)
msrt@annot[, 3] <- msrt@annot[, 2]
 
seurat_cap_reactome <- readRDS(
  "/path/to/Capillary_Adipose_SCENIC_regulonAUC_extended.rds"
)
 
DefaultAssay(seurat_cap_reactome) <- "SCT"
res <- scgsva(seurat_cap_reactome, msrt, method = "UCell", maxRank = 5000)
seurat_cap_reactome <- AddMetaData(seurat_cap_reactome, res@gsva)
 
saveRDS(
  seurat_cap_reactome,
  "Capillary_Adipose_SCENIC_regulonAUC_extended_Reactome.rds"
)
 
# Reload for correlation analysis
seurat_cap_reactome <- readRDS(
  "Capillary_Adipose_SCENIC_regulonAUC_extended_Reactome.rds"
)
 
# Reactome pathway names
df_reactome <- read.csv("Reactome_annotations.csv", fileEncoding = "UTF-8")
Reactome_path <- df_reactome$x
 
# Example: Irf7 regulon vs Reactome (R²)
reg_data <- GetAssayData(
  seurat_cap_reactome,
  assay = "regulonAUC_extended",
  slot  = "data"
)["Irf7", ]
 
correlations_r2 <- list()
for (path in Reactome_path) {
  pathway_data <- seurat_cap_reactome@meta.data[[path]]
  corr <- cor(reg_data, pathway_data, method = "pearson")
  if (!is.na(corr) && corr > 0) {
    correlations_r2[[path]] <- corr^2
  }
}
r_squared_sorted <- sort(unlist(correlations_r2), decreasing = TRUE)
Irf7_df <- data.frame(
  Pathway   = names(r_squared_sorted),
  R_squared = r_squared_sorted
)
 
top20_df <- head(Irf7_df, 20)
 
ggplot(top20_df,
       aes(x = reorder(Pathway, R_squared), y = R_squared)) +
  geom_bar(stat = "identity", fill = "hotpink3") +
  coord_flip() +
  xlab("scGSVA Reactome Pathway") +
  ylab("R-squared with Irf7 regulon") +
  ggtitle("Irf7 regulon and Reactome Pathways") +
  theme_minimal()
 
# Multi-TF network (Irf7, Irf9, Stat1, Stat2, Junb, Pparg, Cebpd, Jund, Jun, Klf2)
targets <- c("Irf7", "Irf9", "Stat1", "Stat2",
             "Junb", "Pparg", "Cebpd", "Jund", "Jun", "Klf2")
 
result_list <- list()
top20_list  <- list()
 
for (tf in targets) {
  reg_tf <- GetAssayData(
    seurat_cap_reactome,
    assay = "regulonAUC_extended",
    slot  = "data"
  )[tf, ]
 
  corr_list <- list()
  for (path in Reactome_path) {
    pathway_data <- seurat_cap_reactome@meta.data[[path]]
    corr <- cor(reg_tf, pathway_data, method = "spearman")
    if (!is.na(corr)) {
      corr_list[[path]] <- corr
    }
  }
 
  corr_sorted <- sort(unlist(corr_list), decreasing = TRUE)
  df_tf <- data.frame(
    Pathway    = names(corr_sorted),
    Coefficient = corr_sorted
  )
  top20_tf <- head(df_tf, 20)
 
  result_list[[tf]] <- df_tf
  top20_list[[tf]]  <- top20_tf
}
 
pathways_list <- lapply(top20_list, function(df) df$Pathway)
names(pathways_list) <- targets
all_pathways <- unique(unlist(pathways_list))
 
edges <- data.frame(
  from       = character(),
  to         = character(),
  Coefficient = numeric(),
  stringsAsFactors = FALSE
)
 
for (tf in names(top20_list)) {
  df_tf <- top20_list[[tf]]
  for (i in seq_len(nrow(df_tf))) {
    edges <- rbind(edges, data.frame(
      from       = tf,
      to         = df_tf$Pathway[i],
      Coefficient = df_tf$Coefficient[i],
      stringsAsFactors = FALSE
    ))
  }
}
edges <- unique(edges)
 
network_graph <- graph_from_data_frame(edges, directed = FALSE)
 
significant_nodes <- unique(edges$to[edges$Coefficient > 0.3])
 
V(network_graph)$label <- NA
V(network_graph)$label[V(network_graph)$name %in% c(significant_nodes, targets)] <-
  V(network_graph)$name[V(network_graph)$name %in% c(significant_nodes, targets)]
 
V(network_graph)$frame.color <- "grey"
V(network_graph)$frame.width <- 1
V(network_graph)$frame.width[V(network_graph)$name %in% significant_nodes] <- 2
 
dataset_colors <- c(
  Irf7  = "royalblue1", Irf9  = "royalblue1",
  Stat1 = "royalblue1", Stat2 = "royalblue1",
  Junb  = "royalblue1",
  Pparg = "coral", Cebpd = "coral",
  Jund  = "coral", Jun   = "coral",
  Klf2  = "coral"
)
 
V(network_graph)$color <- "pink"
V(network_graph)$color[V(network_graph)$name %in% names(dataset_colors)] <-
  dataset_colors[V(network_graph)$name]
 
common_pathways <- names(table(all_pathways)[table(all_pathways) > 1])
V(network_graph)$color[V(network_graph)$name %in% common_pathways] <- "skyblue"
 
V(network_graph)$size <- 10
V(network_graph)$size[V(network_graph)$name %in% all_pathways] <- 2
 
num_edges <- length(E(network_graph)$Coefficient)
E(network_graph)$width <- 1
E(network_graph)$width[E(network_graph)$Coefficient > 0.3] <- 3
E(network_graph)$color <- gray.colors(
  num_edges, start = 0.9, end = 0.1
)[rank(E(network_graph)$Coefficient)]
 
plot(
  network_graph,
  layout              = layout_with_fr(network_graph),
  vertex.label        = V(network_graph)$label,
  vertex.color        = V(network_graph)$color,
  vertex.size         = V(network_graph)$size,
  edge.width          = E(network_graph)$width,
  edge.color          = E(network_graph)$color,
  edge.arrow.size     = 0.5,
  vertex.frame.color  = V(network_graph)$frame.color,
  vertex.frame.width  = V(network_graph)$frame.width
)
 
 
############################
 
# Junb regulon activity and DEGs in adipose ECs
############################
 
# Junb regulon (regulonAUC_extended assay) in adipose ECs
gene_of_interest <- "Junb"
cutoff <- 0.15
 
DefaultAssay(seurat_adipose) <- "regulonAUC_extended"
junb_reg <- seurat_adipose[["regulonAUC_extended"]]@data[gene_of_interest, ]
 
seurat_adipose$Regulon_Junb <- ifelse(
  junb_reg >= cutoff, "Junb high", "Junb low"
)
 
# Subset by age (3 month vs 1 year) as in original
young_obj <- subset(seurat_adipose, subset = Type5 == "3 month")
aged_obj  <- subset(seurat_adipose, subset = Type5 == "1 year")
 
# DEGs in Aged (Junb high vs Junb low)
DEGs_Junb <- FindMarkers(
  aged_obj,
  group.by        = "Regulon_Junb",
  ident.1         = "Junb high",
  ident.2         = "Junb low",
  logfc.threshold = 0.1
)
 
df_EC <- as.data.frame(DEGs_Junb)
 
plot(df_EC$avg_log2FC, -log10(df_EC$p_val + 1e-308),
     xlab = "log2 Fold change",
     ylab = "-log10(p-value)",
     cex  = 0.5,
     col  = "gray",
     pch  = 16)
```
# Part IV. Xenium spatial analyses


## 1. R setup and organ-specific Xenium objects


```r
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(RANN)
})

# Final organ-specific Xenium objects after manual annotation refinement.
heart_qc <- readRDS("heart_qc_annotated_refined_updated.rds")
muscle_qc <- readRDS("muscle_qc_annotated_refined.rds")
adipose_qc <- readRDS("adipose_qc_annotated_refined_updated.rds")

# Final annotation columns:
#   heart_qc   : heart_annotation_refined
#   muscle_qc  : muscle_annotation_refined
#   adipose_qc : adipose_annotation_refined
#
# FOV names:
#   Heart   : Young_Heart, Aged_Heart
#   Muscle  : Young_Muscle, Aged_Muscle
#   Adipose : Young_Adipose, Aged_Adipose
```


## 2. Cell-centroid coordinates and Bst2/Pira2 expression


```r
# Xenium centroid coordinates are expressed in micrometers.
# Expression values are taken from the normalized data layer of the SCT assay.

make_xenium_spatial_table <- function(obj,
                                      fov,
                                      age,
                                      annotation_col,
                                      assay = "SCT") {
  coords <- GetTissueCoordinates(obj, image = fov)

  if (!"cell" %in% colnames(coords)) {
    coords$cell <- rownames(coords)
  }

  stopifnot(all(c("cell", "x", "y") %in% colnames(coords)))

  coords <- coords[, c("cell", "x", "y"), drop = FALSE]
  coords <- coords[!duplicated(coords$cell), , drop = FALSE]

  md <- obj@meta.data
  stopifnot(annotation_col %in% colnames(md))
  stopifnot(all(coords$cell %in% rownames(md)))

  expr <- GetAssayData(
    obj,
    assay = assay,
    layer = "data"
  )

  stopifnot(all(c("Bst2", "Pira2") %in% rownames(expr)))
  stopifnot(all(coords$cell %in% colnames(expr)))

  expr_use <- as.data.frame(
    t(as.matrix(expr[c("Bst2", "Pira2"), coords$cell, drop = FALSE]))
  )
  expr_use$cell <- rownames(expr_use)

  out <- coords %>%
    left_join(expr_use, by = "cell")

  out$Age <- age
  out$annotation <- as.character(md[out$cell, annotation_col])

  out
}


heart_spatial_df <- bind_rows(
  make_xenium_spatial_table(
    heart_qc,
    fov = "Young_Heart",
    age = "Young",
    annotation_col = "heart_annotation_refined"
  ),
  make_xenium_spatial_table(
    heart_qc,
    fov = "Aged_Heart",
    age = "Aged",
    annotation_col = "heart_annotation_refined"
  )
)

muscle_spatial_df <- bind_rows(
  make_xenium_spatial_table(
    muscle_qc,
    fov = "Young_Muscle",
    age = "Young",
    annotation_col = "muscle_annotation_refined"
  ),
  make_xenium_spatial_table(
    muscle_qc,
    fov = "Aged_Muscle",
    age = "Aged",
    annotation_col = "muscle_annotation_refined"
  )
)

adipose_spatial_df <- bind_rows(
  make_xenium_spatial_table(
    adipose_qc,
    fov = "Young_Adipose",
    age = "Young",
    annotation_col = "adipose_annotation_refined"
  ),
  make_xenium_spatial_table(
    adipose_qc,
    fov = "Aged_Adipose",
    age = "Aged",
    annotation_col = "adipose_annotation_refined"
  )
)
```


## 3. Organ-specific definition of Bst2-high endothelial-containing segments


```r
# Bst2-high cells are selected only from cells annotated as
# "Endothelial-containing segment".
#
# Within each organ, Young and Aged endothelial-containing segments are pooled
# to define one common cutoff. If at least 10% of these cells express Bst2,
# Bst2-high is defined by the organ-specific 90th percentile of Bst2 expression.
# In skeletal muscle, fewer than 10% of endothelial-containing segments express
# Bst2; therefore, all Bst2-positive endothelial-containing segments are retained.

define_bst2_high_ec <- function(spatial_df, organ_name) {
  ec_df <- spatial_df %>%
    filter(annotation == "Endothelial-containing segment")

  stopifnot(nrow(ec_df) > 0)

  positive_fraction <- mean(ec_df$Bst2 > 0, na.rm = TRUE)
  q90 <- unname(quantile(ec_df$Bst2, probs = 0.90, na.rm = TRUE, type = 7))

  if (positive_fraction < 0.10 || q90 <= 0) {
    cutoff <- 0
    rule <- "All Bst2-positive endothelial-containing segments"
    high_df <- ec_df %>%
      filter(Bst2 > 0)
  } else {
    cutoff <- q90
    rule <- "Bst2 expression at or above the organ-specific 90th percentile"
    high_df <- ec_df %>%
      filter(Bst2 > 0, Bst2 >= cutoff)
  }

  summary_df <- data.frame(
    Organ = organ_name,
    Endothelial_n = nrow(ec_df),
    Bst2_positive_n = sum(ec_df$Bst2 > 0, na.rm = TRUE),
    Bst2_positive_fraction = positive_fraction,
    Bst2_high_cutoff = cutoff,
    Bst2_high_n = nrow(high_df),
    Selection_rule = rule,
    stringsAsFactors = FALSE
  )

  list(
    high_cells = high_df,
    summary = summary_df
  )
}


heart_bst2_definition <- define_bst2_high_ec(
  heart_spatial_df,
  organ_name = "Heart"
)

muscle_bst2_definition <- define_bst2_high_ec(
  muscle_spatial_df,
  organ_name = "Muscle"
)

adipose_bst2_definition <- define_bst2_high_ec(
  adipose_spatial_df,
  organ_name = "Adipose"
)

bst2_high_definition_summary <- bind_rows(
  heart_bst2_definition$summary,
  muscle_bst2_definition$summary,
  adipose_bst2_definition$summary
)

bst2_high_definition_summary

# Confirm the muscle-specific rule used in the analysis.
stopifnot(
  muscle_bst2_definition$summary$Bst2_positive_fraction < 0.10
)
```


## 4. Nearest-neighbor distance from monocyte/macrophage-lineage cells to Bst2-high endothelial-containing segments


```r
calculate_nearest_bst2_high_ec_distance <- function(spatial_df,
                                                    bst2_high_ec_df) {
  macrophage_df <- spatial_df %>%
    filter(annotation == "Monocyte/macrophage lineage") %>%
    mutate(
      Pira2_status = ifelse(
        Pira2 > 0,
        "Pira2-positive",
        "Pira2-negative"
      )
    )

  distance_list <- lapply(c("Young", "Aged"), function(age_group) {
    query_df <- macrophage_df %>%
      filter(Age == age_group)

    reference_df <- bst2_high_ec_df %>%
      filter(Age == age_group)

    if (nrow(query_df) == 0 || nrow(reference_df) == 0) {
      return(NULL)
    }

    nn <- RANN::nn2(
      data = as.matrix(reference_df[, c("x", "y")]),
      query = as.matrix(query_df[, c("x", "y")]),
      k = 1
    )

    query_df$nearest_Bst2_high_EC_distance_um <- nn$nn.dists[, 1]
    query_df$nearest_Bst2_high_EC_cell <-
      reference_df$cell[nn$nn.idx[, 1]]

    query_df
  })

  bind_rows(distance_list) %>%
    mutate(
      Age = factor(Age, levels = c("Young", "Aged")),
      Pira2_status = factor(
        Pira2_status,
        levels = c("Pira2-negative", "Pira2-positive")
      )
    )
}
```


## 5. Cumulative proximity curves and bootstrap confidence intervals


```r
summarize_cumulative_proximity <- function(distance_df,
                                           thresholds_um = seq(0, 200, 5)) {
  summary_list <- lapply(
    split(distance_df, list(distance_df$Age, distance_df$Pira2_status),
          drop = TRUE),
    function(df) {
      data.frame(
        Age = as.character(df$Age[1]),
        Pira2_status = as.character(df$Pira2_status[1]),
        Distance_um = thresholds_um,
        Cumulative_percent = vapply(
          thresholds_um,
          function(d) {
            100 * mean(df$nearest_Bst2_high_EC_distance_um <= d)
          },
          numeric(1)
        ),
        Cell_n = nrow(df),
        stringsAsFactors = FALSE
      )
    }
  )

  bind_rows(summary_list) %>%
    mutate(
      Age = factor(Age, levels = c("Young", "Aged")),
      Pira2_status = factor(
        Pira2_status,
        levels = c("Pira2-negative", "Pira2-positive")
      )
    )
}


bootstrap_proximity_difference <- function(distance_df,
                                           thresholds_um = seq(0, 200, 5),
                                           n_boot = 1000,
                                           seed = 123) {
  set.seed(seed)

  age_results <- lapply(c("Young", "Aged"), function(age_group) {
    df_age <- distance_df %>%
      filter(Age == age_group)

    d_pos <- df_age$nearest_Bst2_high_EC_distance_um[
      df_age$Pira2_status == "Pira2-positive"
    ]
    d_neg <- df_age$nearest_Bst2_high_EC_distance_um[
      df_age$Pira2_status == "Pira2-negative"
    ]

    if (length(d_pos) == 0 || length(d_neg) == 0) {
      return(NULL)
    }

    observed_difference <- vapply(
      thresholds_um,
      function(d) {
        100 * (
          mean(d_pos <= d) -
            mean(d_neg <= d)
        )
      },
      numeric(1)
    )

    boot_matrix <- replicate(
      n_boot,
      {
        d_pos_boot <- sample(d_pos, size = length(d_pos), replace = TRUE)
        d_neg_boot <- sample(d_neg, size = length(d_neg), replace = TRUE)

        vapply(
          thresholds_um,
          function(d) {
            100 * (
              mean(d_pos_boot <= d) -
                mean(d_neg_boot <= d)
            )
          },
          numeric(1)
        )
      }
    )

    data.frame(
      Age = age_group,
      Distance_um = thresholds_um,
      Difference_percent = observed_difference,
      CI_lower = apply(
        boot_matrix,
        1,
        quantile,
        probs = 0.025,
        na.rm = TRUE
      ),
      CI_upper = apply(
        boot_matrix,
        1,
        quantile,
        probs = 0.975,
        na.rm = TRUE
      ),
      Pira2_positive_n = length(d_pos),
      Pira2_negative_n = length(d_neg),
      stringsAsFactors = FALSE
    )
  })

  bind_rows(age_results) %>%
    mutate(Age = factor(Age, levels = c("Young", "Aged")))
}


plot_proximity_difference <- function(bootstrap_df, organ_name) {
  age_colors <- c(
    Young = "#F8766D",
    Aged = "#00BFC4"
  )

  ggplot(
    bootstrap_df,
    aes(
      x = Distance_um,
      y = Difference_percent,
      color = Age,
      fill = Age
    )
  ) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.4,
      color = "grey40"
    ) +
    geom_ribbon(
      aes(
        ymin = CI_lower,
        ymax = CI_upper
      ),
      alpha = 0.20,
      color = NA
    ) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = age_colors) +
    scale_fill_manual(values = age_colors) +
    labs(
      title = organ_name,
      x = "Distance to nearest Bst2-high endothelial-containing segment (µm)",
      y = "Pira2-positive minus Pira2-negative monocyte/macrophage-lineage cells (%)"
    ) +
    theme_classic(base_size = 12)
}
```


## 6. Organ-wise Bst2–Pira2 proximity analysis


```r
run_bst2_pira2_proximity <- function(spatial_df,
                                     bst2_definition,
                                     organ_name,
                                     thresholds_um = seq(0, 200, 5),
                                     n_boot = 1000,
                                     seed = 123) {
  distance_df <- calculate_nearest_bst2_high_ec_distance(
    spatial_df = spatial_df,
    bst2_high_ec_df = bst2_definition$high_cells
  )

  cumulative_df <- summarize_cumulative_proximity(
    distance_df = distance_df,
    thresholds_um = thresholds_um
  )

  bootstrap_df <- bootstrap_proximity_difference(
    distance_df = distance_df,
    thresholds_um = thresholds_um,
    n_boot = n_boot,
    seed = seed
  )

  proximity_plot <- plot_proximity_difference(
    bootstrap_df = bootstrap_df,
    organ_name = organ_name
  )

  list(
    bst2_definition = bst2_definition$summary,
    bst2_high_ec = bst2_definition$high_cells,
    macrophage_distance = distance_df,
    cumulative_proximity = cumulative_df,
    bootstrap_difference = bootstrap_df,
    proximity_plot = proximity_plot
  )
}


heart_proximity <- run_bst2_pira2_proximity(
  spatial_df = heart_spatial_df,
  bst2_definition = heart_bst2_definition,
  organ_name = "Heart"
)

muscle_proximity <- run_bst2_pira2_proximity(
  spatial_df = muscle_spatial_df,
  bst2_definition = muscle_bst2_definition,
  organ_name = "Muscle"
)

adipose_proximity <- run_bst2_pira2_proximity(
  spatial_df = adipose_spatial_df,
  bst2_definition = adipose_bst2_definition,
  organ_name = "Adipose"
)

heart_proximity$proximity_plot
muscle_proximity$proximity_plot
adipose_proximity$proximity_plot
```


## 7. Save proximity-analysis results


```r
saveRDS(
  heart_proximity,
  file = "Heart_Bst2_Pira2_proximity_analysis.rds"
)

saveRDS(
  muscle_proximity,
  file = "Muscle_Bst2_Pira2_proximity_analysis.rds"
)

saveRDS(
  adipose_proximity,
  file = "Adipose_Bst2_Pira2_proximity_analysis.rds"
)

write.csv(
  bst2_high_definition_summary,
  file = "Bst2_high_definition_summary.csv",
  row.names = FALSE
)
```
