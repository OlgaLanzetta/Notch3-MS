knitr::opts_chunk$set(echo = TRUE)
library(RColorBrewer)
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(tidyr)
library(pheatmap)
library(DESeq2)
library(xlsx)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)

# ---------------------------------------------------------
# scRNAseq Analysis Script — Benamar et al. (Part II)
# ---------------------------------------------------------
# This script continues the analysis of the 5 scRNA-seq samples.
# It starts from the processed and integrated Seurat object (see Part I)
# and performs differential expression analysis between conditions.

# ----------------------
# 1. Load Integrated Seurat Object
# ----------------------
path <- "/path/to/files/rds/"  # TODO: Update to your path
myfile <- "seurat.integrated_processed_slim.rds"
seurat.integrated <- readRDS(file = paste0(path, myfile))
#head(seurat.integrated@meta.data)

# ----------------------
# 2. Visualization of Clusters and Cell Types
# ----------------------
clust_plot <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'main_label', label = TRUE)
type_plot <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'integrated_snn_res.0.15', label = TRUE)
clust_plot | type_plot

g1 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'integrated_snn_res.0.15',
              label = TRUE, repel = TRUE, split.by = "condition")
ggsave(filename = "UMAP_plot1.pdf", g1, device = "pdf", dpi = 600)

g2 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'integrated_snn_res.0.15',
              split.by = "condition")
ggsave(filename = "UMAP_plot2.pdf", g2, device = "pdf", dpi = 600)

# ----------------------
# 3. Differential Expression Analysis - Pseudobulk with DESeq2
# ----------------------
# Aggregate counts per cluster and sample
DefaultAssay(object = seurat.integrated) <- "RNA"
cts <- AggregateExpression(seurat.integrated, 
                          group.by = c("integrated_snn_res.0.15", "orig.ident"),
                          assays = "RNA",
                          slot = "counts",
                          return.seurat = FALSE)
cts <- cts$RNA
# transpose and convert to data.frame
cts.t <- as.data.frame(t(cts))
# get cluster from rownames
splitRows <- gsub('_.*', '', rownames(cts.t))
# split data.frame per cluster
cts.split <- split.data.frame(cts.t, f = factor(splitRows))
# fix colnames and transpose
cts.split.modified <- lapply(cts.split, function(x){
  rownames(x) <- substring(gsub('*._', '', rownames(x)), first = 2, last = 1000000L)
  t(x)
})

# Automatic DE analysis for all clusters with at least 200 cells (all conditions)
DESeq2_DE_list <- list()  # Differential expression results for notch3 vs wt
cell_in_cluster <- table(Idents(seurat.integrated))
cluster_labels <- levels(Idents(seurat.integrated))

for (iclust in cluster_labels) {
  if (cell_in_cluster[iclust] > 200) {
    print("---------------------------------------")
    print(paste0("Performing DE on cluster: ", iclust))
    
    counts <- cts.split.modified[[iclust]]
    colData <- data.frame(samples = colnames(counts),
                          condition = as.factor(c("notch3", "notch3", "notch3", "wt", "wt")),
                          replicate = c(1, 2, 3, 1, 3))
    # create DESeq2 object
    dds <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = colData,
                                  design = ~ condition)
    dds$condition <- relevel(dds$condition, "wt")
    # filter lowly expressed genes
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep, ]
    print(dim(dds))
    vsd <- vst(dds, blind = FALSE)
    # PCA for QC
    pcaData <- plotPCA(vsd, intgroup = c("condition", "replicate"), returnData = TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    head(pcaData)
    # DESeq2
    dds <- DESeq(dds, quiet = TRUE)
    resultsNames(dds)
    res <- results(dds, name = "condition_notch3_vs_wt", alpha = 0.05)
    print(summary(res))
    res <- res[order(res$padj), ]
    res <- as.data.frame(res)
    res$GeneName <- rownames(res)
    DESeq2_DE_list[[iclust]] <- res[, c(7, 1, 2, 3, 4, 5, 6)]
    head(DESeq2_DE_list[[iclust]], 10)
    #write.xlsx(DESeq2_DE_list[[iclust]], file = paste0("DESeq2_", iclust, "_condition_notch3_vs_wt.xlsx"), 
    #           sheetName = "Sheet1", col.names = TRUE, row.names = FALSE, append = FALSE)
  }
}

# ----------------------
# 4. Volcano Plot Function
# ----------------------
library(ggrepel)
make_volcano <- function(res_df, cluster_label) {
  res_df$Gene <- rownames(res_df)
  # Remove NA
  res_df <- res_df[!is.na(res_df$log2FoldChange) & !is.na(res_df$padj), ]
  res_df$Genes <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1,
                         "Significant", "Not Significant")
  # Top gene labels
  top_genes <- rbind(
    head(res_df[res_df$Genes == "Significant", ], 5),
    tail(res_df[res_df$Genes == "Significant", ], 5)
  )
  max_abs_lfc <- max(abs(res_df$log2FoldChange), na.rm = TRUE)
  xlim_range <- ceiling(max_abs_lfc)
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(data = subset(res_df, Genes == "Not Significant"),
               aes(color = Genes), size = 1.5, alpha = 0.8) +
    geom_point(data = subset(res_df, Genes == "Significant"),
               aes(color = Genes), size = 3, alpha = 0.8) +
    geom_text(
      data = top_genes,
      aes(label = Gene),
      vjust = 1.5,
      size = 5,
      color = "black",
      check_overlap = TRUE
    ) +
    scale_color_manual(values = c("darkgrey", "red")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    scale_x_continuous(limits = c(-xlim_range, xlim_range)) +
    theme_classic(base_size = 14) +
    labs(
      x = "Log2 Fold Change",
      y = "-Log10 adjusted p-value",
      title = paste0("Volcano Plot: ", cluster_label, " (notch3 vs wt)")
    ) +
    guides(color = guide_legend(override.aes = list(size = 3)))
  return(p)
}
# Example volcano plots (for clusters with enough cells)
volcano_C3 <- make_volcano(DESeq2_DE_list[["C3"]], "C3")
volcano_C0 <- make_volcano(DESeq2_DE_list[["C0"]], "C0")
volcano_C1 <- make_volcano(DESeq2_DE_list[["C1"]], "C1")
volcano_C7 <- make_volcano(DESeq2_DE_list[["C7"]], "C7")

print(volcano_C0)
print(volcano_C1)
print(volcano_C7)
print(volcano_C3)

# Define output directory for plots (create if it doesn't exist)
path <- "/path/to/files/"  # TODO: Update to your path
plots_dir <- file.path(path, "Plots")
if (!dir.exists(plots_dir)) dir.create(plots_dir)

# Save volcano plots in plots_dir
ggsave(filename = file.path(plots_dir, "Volcano_C3.pdf"), plot = volcano_C3, width = 10, height = 9)
ggsave(filename = file.path(plots_dir, "Volcano_C0.pdf"), plot = volcano_C0, width = 10, height = 9)
ggsave(filename = file.path(plots_dir, "Volcano_C1.pdf"), plot = volcano_C1, width = 10, height = 9)
ggsave(filename = file.path(plots_dir, "Volcano_C7.pdf"), plot = volcano_C7, width = 10, height = 9)


# ----------------------
# 5. Custom DESeq2 Example: C7 vs C0 in WT samples
# ----------------------
seurat.integrated$cluster_condition <- paste0(Idents(seurat.integrated), "_", seurat.integrated$orig.ident)
cts2 <- AggregateExpression(
  seurat.integrated,
  group.by = "cluster_condition",
  assays = "RNA",
  slot = "counts",
  return.seurat = FALSE
)$RNA

# Select only C0_wt and C7_wt samples
counts_sub <- cts2[, grepl("^C0-wt|^C7-wt", colnames(cts2))]
group_labels <- ifelse(grepl("^C0-", colnames(counts_sub)), "C0", "C7")
colData <- data.frame(
  row.names = colnames(counts_sub),
  group = factor(group_labels, levels = c("C0", "C7"))
)

dds <- DESeqDataSetFromMatrix(countData = counts_sub, colData = colData, design = ~ group)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", "C7", "C0"), alpha = 0.05)
print(summary(res))
res <- res[order(res$padj), ]
res <- as.data.frame(res)
res$GeneName <- rownames(res)
res_df <- dplyr::filter(res, padj < 0.05)

res_df$Genes <- ifelse(
  res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1,
  "Significant",
  "Not Significant"
)
path <- "/path/to/files/rds/"  # TODO: Update to your path
plots_dir <- file.path(path, "Tables")
if (!dir.exists(plots_dir)) dir.create(plots_dir)

write.xlsx(as.data.frame(res), "DESeq2_C7_vs_C0_WT.xlsx", row.names = FALSE)

# Volcano plot for C7 vs C0 WT
library(ggplot2)
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Genes)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(values = c("gray", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_bw(base_size = 14) +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 adjusted p-value",
    title = "Volcano Plot: C7 vs C0 (WT)"
  )

# Clean NA for label
res_df_clean <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), ]
top_genes <- rbind(
  head(res_df_clean[res_df_clean$Genes == "Significant", ], 5),
  tail(res_df_clean[res_df_clean$Genes == "Significant", ], 5)
)

ggplot(res_df_clean, aes(x = log2FoldChange, y = -log10(padj), color = Genes)) +
  geom_point(data = subset(res_df_clean, Genes == "Not Significant"),
             aes(color = Genes), size = 1.5, alpha = 0.8) +
  geom_point(data = subset(res_df_clean, Genes == "Significant"),
             aes(color = Genes), size = 3, alpha = 0.8) +
  geom_text(
    data = top_genes,
    aes(label = GeneName),
    vjust = 1.5,
    size = 5,
    color = "black",
    check_overlap = TRUE
  ) +
  scale_color_manual(values = c("darkgrey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_classic(base_size = 14) +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 adjusted p-value",
    title = "Volcano Plot: C7 vs C0 (WT)"
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))

# ----------------------
# END OF SCRIPT
# ----------------------

