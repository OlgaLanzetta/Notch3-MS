# scRNAseq Analysis Script — Benamar et al.
# ----------------------------------------
# This script performs single-cell RNA-seq data processing and analysis
# as described in the manuscript "Notch3+ Regulatory T cells Drive Autoimmune Neuroinflammation in Multiple Sclerosis".

# ----------------------
# 1. Load Required Packages
# ----------------------
knitr::opts_chunk$set(echo = TRUE)
library(RColorBrewer)
library(Seurat)
library(DoubletFinder)
library(patchwork)
library(dplyr)
library(ggplot2)
library(tidyr)
library(pheatmap)
library(SingleR)
library(celldex)
library(speckle)
library(DESeq2)
library(xlsx)
library(scales)
library(ComplexHeatmap)
library(circlize)

# Reference dataset for cell classification
ref.se <- celldex::ImmGenData()

# ----------------------
# 2. Data Analysis Setup
# ----------------------
# Specify your sample names and data directory
sample_list <- c("notch3_1", "notch3_2", "notch3_3", "wt_1", "wt_3")
data_dir <- "/path/to/files/"   # TODO: change this to your data path

so_list <- list()
markers_list <- list()

# ----------------------
# 3. Process Each Sample
# ----------------------
for (isample in sample_list) {
  set.seed(123)
  cat("\n------------------------------\n")
  cat(paste0(isample, " - Start Reading!\n"))
  myfile <- file.path(data_dir, isample, "filtered_feature_bc_matrix.h5")
  
  h5_obj <- Read10X_h5(filename = myfile, use.names = TRUE, unique.features = TRUE)
  
  Seurat_obj <- CreateSeuratObject(counts = h5_obj,
                                   project = isample,
                                   min.cells = 0,
                                   min.features = 200)
  
  cat("Original dimensions: ", dim(Seurat_obj), "\n")
  
  # Rename cells to keep sample info
  Seurat_obj <- RenameCells(object = Seurat_obj, add.cell.id = paste0(isample, "_"))
  
  cat("QC and Filtering\n")
  Seurat_obj[["percent.mt"]] <- PercentageFeatureSet(Seurat_obj, pattern = "^mt-")
  
  # Visual QC plots
  VlnPlot(Seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  FeatureScatter(object = Seurat_obj, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA') + geom_smooth(method = "lm")
  
  cat("Number of cells before filtering: ", ncol(Seurat_obj), "\n")
  
  # Filtering low quality cells
  Seurat_obj <- subset(Seurat_obj,
                       subset = nCount_RNA > 500 & nCount_RNA < 30000 &
                                nFeature_RNA > 200 & nFeature_RNA < 3000 &
                                percent.mt < 10)
  
  cat("Number of cells after filtering: ", ncol(Seurat_obj), "\n")
  
  # ----------------------
  # Cell Type Classification
  # ----------------------
  Counts <- GetAssayData(Seurat_obj, slot = "counts")
  pred <- SingleR(test = Counts, ref = ref.se, labels = ref.se$label.main)
  Seurat_obj$main_label <- pred$pruned.labels
  print(table(Seurat_obj$main_label))
  
  Seurat_obj <- subset(Seurat_obj, main_label %in% c("T cells", "NKT"))
  print(table(Seurat_obj$main_label))
  cat("Number of cells after cell classification: ", ncol(Seurat_obj), "\n")
  
  # ----------------------
  # Seurat Standard Workflow
  # ----------------------
  Seurat_obj <- NormalizeData(Seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  Seurat_obj <- FindVariableFeatures(Seurat_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  HVG <- VariableFeatures(Seurat_obj)
  top15 <- head(HVG, 15)
  p2 <- VariableFeaturePlot(Seurat_obj)
  p2 <- LabelPoints(plot = p2, points = top15, repel = TRUE)
  # print(p2)
  
  all.genes <- rownames(Seurat_obj)
  Seurat_obj <- ScaleData(Seurat_obj, features = all.genes)
  Seurat_obj <- RunPCA(Seurat_obj, features = HVG, verbose = FALSE)
  ElbowPlot(Seurat_obj, ndims = 30)
  DimPlot(Seurat_obj, reduction = "pca")
  
  set.seed(123)
  ndims <- 25
  Seurat_obj <- FindNeighbors(Seurat_obj, dims = 1:ndims, verbose = FALSE)
  Seurat_obj <- FindClusters(Seurat_obj, resolution = c(0.2, 0.5, 0.7), verbose = FALSE)
  Seurat_obj <- RunUMAP(Seurat_obj, dims = 1:ndims)
  DimPlot(Seurat_obj, group.by = "RNA_snn_res.0.5", reduction = "umap", label = TRUE)
  Idents(Seurat_obj) <- "RNA_snn_res.0.5"
  
  # ----------------------
  # Doublet Detection with DoubletFinder
  # ----------------------
  sweep.res.list <- paramSweep_v3(Seurat_obj, PCs = 1:ndims, sct = FALSE, num.cores = 8)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- bcmvn %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
  pK <- as.numeric(as.character(pK[[1]]))
  
  annotations <- Seurat_obj@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.075 * nrow(Seurat_obj@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  Seurat_obj <- doubletFinder_v3(Seurat_obj,
                                 PCs = 1:ndims,
                                 pN = 0.25,
                                 pK = pK,
                                 nExp = nExp_poi.adj,
                                 reuse.pANN = FALSE, sct = FALSE)
  
  colnames(Seurat_obj@meta.data) <- gsub("_0..*", "", colnames(Seurat_obj@meta.data))
  print(table(Seurat_obj$DF.classifications))
  
  DimPlot(Seurat_obj, group.by = "DF.classifications", reduction = "umap")
  
  # Filter out doublets
  Seurat_obj <- subset(Seurat_obj, subset = DF.classifications == "Singlet")
  cat("Number of cells after doublet removal: ", ncol(Seurat_obj), "\n")
  cat("Final object dimensions: ", dim(Seurat_obj), "\n")
  
  so_list[[isample]] <- Seurat_obj
  cat(paste0(isample, " .....ended!\n"))
}

# ----------------------
# 4. Save List of Seurat Objects
# ----------------------
saveRDS(so_list, file = file.path(data_dir, "so_list_slim.rds"))


# ----------------------
# 5. Integration with Seurat (CCA)
# ----------------------
all_genes <- unique(unlist(lapply(so_list, rownames)))
cat("Number of unique genes: ", length(all_genes), "\n")

set.seed(123)
features <- SelectIntegrationFeatures(object.list = so_list, nfeatures = 20000)
anchors <- FindIntegrationAnchors(object.list = so_list, anchor.features = features, verbose = FALSE)
seurat.integrated <- IntegrateData(anchorset = anchors, verbose = FALSE)

# Add sample info
seurat.integrated@meta.data$sample <- seurat.integrated$orig.ident
seurat.integrated@meta.data <- separate(seurat.integrated@meta.data, col = sample, into = c("condition", "replicate"), sep = "_")
Idents(seurat.integrated) <- 0
cat("Integrated object dimensions: ", dim(seurat.integrated), "\n")
cat("Integration DONE!\n")

# Save integrated object
saveRDS(seurat.integrated, file = file.path(data_dir, "seurat.integrated.rds"))


# ----------------------
# 6. Downstream Analysis of Integrated Object
# ----------------------
DefaultAssay(seurat.integrated) <- "integrated"
seurat.integrated <- ScaleData(seurat.integrated, features = rownames(seurat.integrated))
seurat.integrated <- RunPCA(seurat.integrated, npcs = 50, verbose = FALSE)
DimPlot(seurat.integrated, reduction = "pca")

set.seed(1234)
ndims <- 30
seurat.integrated <- RunUMAP(seurat.integrated, reduction = "pca", dims = 1:ndims, verbose = FALSE)
DimPlot(seurat.integrated, reduction = "umap")

g1 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'condition')
g2 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'replicate', cols = c('red', 'green', 'blue'))
g1; g2

# Clustering
seurat.integrated@meta.data <- dplyr::select(seurat.integrated@meta.data,
                                             -c("RNA_snn_res.0.2", "RNA_snn_res.0.5", "RNA_snn_res.0.7", "seurat_clusters"))
set.seed(321)
seurat.integrated <- FindNeighbors(seurat.integrated, reduction = "pca", dims = 1:ndims, verbose = FALSE)
seurat.integrated <- FindClusters(seurat.integrated, resolution = 0.15, verbose = FALSE)
seurat.integrated@meta.data$integrated_snn_res.0.15 <- as.factor(paste0("C", seurat.integrated@meta.data$integrated_snn_res.0.15))
Idents(seurat.integrated) <- seurat.integrated@meta.data$integrated_snn_res.0.15

# Define cluster colors
colors <- c("C0" = "#1F77B4", "C1" = "#AEC7E8", "C2" = "darkgreen", "C3" = "#98DF8A",
            "C4" = "#D62728", "C5" = "#FF9896", "C6" = "#9467BD", "C7" = "#C5B0D5",
            "C8" = "#8C564B", "C9" = "burlywood3", "C10" = "#FF7F0E", "C11" = "#FFBB78")

g2a <- DimPlot(seurat.integrated, reduction = 'umap', label = TRUE, repel = TRUE, cols = colors,
               group.by = "integrated_snn_res.0.15", pt.size = 0.8, label.size = 6) +
  guides(color = guide_legend(override.aes = list(size = 8), ncol = 1)) + ggtitle("UMAP")

g3 <- DimPlot(seurat.integrated, reduction = 'umap', split.by = 'condition', label = TRUE,
              repel = TRUE, cols = colors, pt.size = 0.8, group.by = "integrated_snn_res.0.15") +
  guides(color = guide_legend(override.aes = list(size = 12), ncol = 1)) + ggtitle("UMAP")

plots_dir <- file.path(data_dir, "Plots")
if (!dir.exists(plots_dir)) dir.create(plots_dir)
ggsave(filename = file.path(plots_dir, "UMAP_plot_main.pdf"), plot = g2a, device = "pdf", dpi = 600)
ggsave(filename = file.path(plots_dir, "UMAP_plot.pdf"), plot = g3, device = "pdf", width = 20, height = 10, dpi = 1000)

g3.n <- DimPlot(seurat.integrated, reduction = 'umap', split.by = 'condition', group.by = "integrated_snn_res.0.15")
ggsave(filename = file.path(plots_dir, "UMAP_plot2.pdf"), plot = g3.n, device = "pdf", dpi = 600)

table(seurat.integrated@meta.data$integrated_snn_res.0.15)

# ----------------------
# 7. Feature Plots and Marker Identification
# ----------------------
myfeatures <- c("Ikzf2", "Ctla4", "Il2ra", "Foxp3", "Tnfrsf18")
FeaturePlot(seurat.integrated, features = myfeatures, min.cutoff = "q5")

all.markers <- FindAllMarkers(seurat.integrated, logfc.threshold = 0.25, min.pct = 0.15, only.pos = TRUE)
Idents(seurat.integrated) <- seurat.integrated$integrated_snn_res.0.15

features <- all.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

p5 <- DotPlot(seurat.integrated, features = unique(features$gene)) + RotatedAxis()
p5

# ----------------------
# 8. Differential Cell Proportions
# ----------------------
pt <- table(Idents(seurat.integrated), seurat.integrated$orig.ident)
knitr::kable(pt, "simple")

pt_df <- as.data.frame(pt)
prop_plot <- ggplot(pt_df, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = colors) +
  theme(legend.title = element_blank())
prop_plot

ggsave(filename = file.path(plots_dir, "proportion_plot.pdf"), plot = prop_plot, device = "pdf", dpi = 600, height = 5, width = 10)

# Test proportions with speckle::propeller
Test.cellprops <- propeller(clusters = Idents(seurat.integrated),
                            sample = seurat.integrated$orig.ident,
                            group = seurat.integrated$condition)
knitr::kable(Test.cellprops)

seurat.integrated <- ScaleData(seurat.integrated)

# ----------------------
# 9. Heatmap of Marker Genes
# ----------------------


# Number of cells for each cluster
# C0    C1    C2    C3    C4    C5    C6    C7    C8   C9  C10  C11 
# 5464  4725  3587  1138  1098  458   380   275   150   93   71   57 

# To ensure the heatmap has the same number of cells per cluster,
# we randomly downsample each cluster to match the size of the smallest cluster.
cell.list <- WhichCells(seurat.integrated, idents = names(colors), downsample = 57)
DefaultAssay(seurat.integrated) <- "RNA"
subset <- subset(seurat.integrated, cells = cell.list)
mat <- subset[["RNA"]]@data[unique(features$gene), ] %>% as.matrix()
mat <- t(scale(t(mat)))
cluster_anno <- subset@meta.data$integrated_snn_res.0.15

quantile(mat, c(0.1, 0.95))

col_fun <- circlize::colorRamp2(c(-1, 0, 3), c("blue", "black", "#FFFF00"))
d <- ComplexHeatmap::Heatmap(mat, name = "Expression",  
        column_split = factor(cluster_anno),
        cluster_columns = FALSE,
        show_column_dend = FALSE,
        cluster_column_slices = TRUE,
        column_title_gp = gpar(fontsize = 20),
        column_gap = unit(0.5, "mm"),
        cluster_rows = FALSE,
        show_row_dend = FALSE,
        col = col_fun,
        row_names_gp = gpar(fontsize = 20),
        column_title_rot = 0,
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = colors))),
        show_column_names = FALSE,
        use_raster = TRUE,
        raster_quality = 4)
d

# ----------------------
# END OF SCRIPT
# ----------------------
