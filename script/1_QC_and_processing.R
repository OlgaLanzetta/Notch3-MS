# load packages 
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
ref.se<-celldex::ImmGenData() ## Database for Cell Classification


## Data Analysis (Part I)

In this document, we report the results of the analysis of  **5** scRNA-seq samples.

- 2 experimental conditions (**notch3**, **wt**)
- 3 replicates per condition **notch3** and 2 replicates for conditions **wt**
We first analyzed the quality of individual samples. Then, we performed the integrated analysis.


### 1. Reading, inspecting and QC filtering individual samples

We used the Seurat pipeline, **SingleR**, and  **DoubletFinder** functions on each sample individually to check its overall quality before proceeding with their integration.

This point includes several steps:

1. QC filtering, where we remove low-quality cells according to mitochondrial percentages, Number of Features, and Number of Reads. This is part of the standard **Seurat** QC workflow.

2. Cell contamination removal. Here, we first applied **SingleR** to classify cell types automatically. Then, according to the experimental protocols,  we select only Cells annotated as "T-Cell" or "NKT".

3. Doublet removal using the R package **DoubletFinder**. To this purpuse, we estimate the DoubletFinder parameters without assuming any ground truth.

sample_list<-c("notch3_1","notch3_2","notch3_3","wt_1","wt_3")

so_list <- list()
markers_list <-list()
for (isample in sample_list) {
  set.seed(123) 
  print("----------------------------------")
  print(paste0(isample, " Start Reading!"))
  myfile=paste0("/path/to/files/",isample,"/filtered_feature_bc_matrix.h5")
  
  h5_obj<-Read10X_h5(filename=myfile, use.names = TRUE, 
                     unique.features = TRUE)
  
  Seurat_obj <- CreateSeuratObject(counts =  h5_obj, 
                                                project = isample,
                                                min.cells = 0,
                                   min.features = 200)
  
  
  print(dim(Seurat_obj))
  ## I rename cells with the sample, to uniquely identify cells since after I will combine the 
  ## available experiments 
  
 Seurat_obj<- RenameCells(object = Seurat_obj, add.cell.id =  paste0(isample,"_"))
  
   print("QC and Filetering")

     Seurat_obj  [["percent.mt"]] <- PercentageFeatureSet(Seurat_obj, pattern = "^mt-")
  
  # Visualize QC metrics as a violin plot
  p1a<-VlnPlot(Seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  print(p1a)

  p1b<-FeatureScatter(object = Seurat_obj, 
                      feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')+
    geom_smooth(method = "lm")
  
#  print(p1b)
  
  ncell_before <- ncol(Seurat_obj)
  print(paste0("Number of cells before filtering:", ncell_before))
  
  ## FILTERING Low quality cell 
  ## note without sample wt2, one can used standard parameters. Instead, with wt2 one must be more stringent
  
   Seurat_obj <- subset(Seurat_obj, 
                        subset = nCount_RNA >500 & nCount_RNA< 30000 &
                          nFeature_RNA > 200 & nFeature_RNA < 3000 &
                          percent.mt < 10)
  
  ncell_after <- ncol(Seurat_obj)
  print(paste0("Number of cells after filtering:", ncell_after))
  
  ###################################################
  ### Classify cell Types using SingleR
  ##################################################
  Counts<-GetAssayData(Seurat_obj,slot="counts")
  
  pred <- SingleR(test = Counts,
                  ref = ref.se,
                  labels = ref.se$label.main)
  
  Seurat_obj$main_label<-pred$pruned.labels
  
  print(table(Seurat_obj$main_label))
  
  Seurat_obj<-subset(Seurat_obj,main_label %in% c("T cells", "NKT"))
  print(table(Seurat_obj$main_label))
  
  ncell_after_cellclass <- ncol(Seurat_obj)
  print(paste0("number of cells after Cell classification:",  ncell_after_cellclass))
  
  
  print("..........Analyzing Individual sample with Seurat.............")
  
### NORMALIZATION  
  
  Seurat_obj <- NormalizeData(Seurat_obj, 
                                   normalization.method = "LogNormalize", 
                                   scale.factor = 10000)
### HVG genes
  
  Seurat_obj<- FindVariableFeatures(Seurat_obj, 
                                          selection.method = "vst", 
                                          nfeatures = 2000, verbose = FALSE)
  
  # Identify the 15 most highly variable genes
  HVG<-VariableFeatures(Seurat_obj)
  top15 <- head(HVG, 15)
  
  #print(top15)
  # plot variable features with and without labels
  p2 <- VariableFeaturePlot(Seurat_obj)
  p2 <- LabelPoints(plot = p2, points = top15, repel = TRUE)
 # print(p2)

  ## SCALE DATA   
  all.genes <- rownames(Seurat_obj)
  Seurat_obj<- ScaleData(Seurat_obj, features = all.genes)
  
## DIMENSION REDUCTION (Da valutare quante dimensioni inserire)
  
  Seurat_obj <- RunPCA(Seurat_obj, 
                            features = HVG,verbose = FALSE)
  
  a1<-ElbowPlot(Seurat_obj,ndims = 30)
  
 # print(a1)
  
  p3<-DimPlot(Seurat_obj, reduction = "pca")
  print(p3)
  
## Preliminary CLUSTERING (this part will be re-done after integration)
  set.seed(123) 
  ndims=25
  Seurat_obj <- FindNeighbors(Seurat_obj, dims = 1:ndims,verbose = FALSE)
 
  Seurat_obj <- FindClusters(Seurat_obj, resolution = c(0.2,0.5,0.7),verbose = FALSE)
  
  Seurat_obj <- RunUMAP(Seurat_obj , dims = 1:ndims)

  p4<-DimPlot(Seurat_obj, group.by = "RNA_snn_res.0.5",
              reduction = "umap",label = TRUE)
  print(p4)
  
  ##We can play with the resolution parameters when choosing the clustering
  Idents(Seurat_obj)<-"RNA_snn_res.0.5" 
#  print(head(Seurat_obj@meta.data))
  
  ###############################################################
  ## Finding and removing potential doublets using doubletFinder 
  ################################################################
## pK Identification (no ground-truth) 

#set.seed(12345)  
sweep.res.list <- paramSweep_v3(Seurat_obj, PCs = 1:ndims, sct = FALSE,num.cores=8)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)


pK <- bcmvn %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

## Homotypic Doublet Proportion Estimate 

annotations <- Seurat_obj@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- 
nExp_poi <- round(0.075*nrow(Seurat_obj@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# run doubletFinder 
Seurat_obj <- doubletFinder_v3(Seurat_obj, 
                                         PCs = 1:ndims, 
                                         pN = 0.25, 
                                         pK = pK, 
                                         nExp = nExp_poi.adj,
                                         reuse.pANN = FALSE, sct = FALSE)


colnames(Seurat_obj@meta.data)<-gsub("_0..*","",colnames(Seurat_obj@meta.data))

print(table(Seurat_obj$DF.classifications))

pd<-DimPlot(Seurat_obj, group.by = "DF.classifications",
              reduction = "umap")
print(pd)
  
#### filtering potential doublets
Seurat_obj <-subset(Seurat_obj, 
                         subset = DF.classifications=="Singlet")

  ncell_after_doublet <- ncol(Seurat_obj)
  print(paste0("number of cells after Doublet Removal:",  ncell_after_doublet))
  
  ####################################
  print(dim(Seurat_obj))
  so_list[[isample]] <-  Seurat_obj 
  

  print(paste0(isample, " .....ended!"))
}

## Save the list of seurat objects
saveRDS(so_list, file = "/path/to/files/so_list_slim.rds")


### 2.  Integrating  samples using **CCA** approach with Seurat Pipeline.

Here, we integrate all the samples and remove the batch effect.
Indeed, preliminary analysis of merged data showed a batch effect primarily associated with the replicates.

There are several ways to proceed with the integration. Here, we have used the standard CCA approach available in Seurat. 

## the so_list already contains Normalized data and HVG useful for selecting anchors
all_genes<-unique(c(rownames(so_list[[1]]),rownames(so_list[[2]]),rownames(so_list[[3]]),rownames(so_list[[4]]),rownames(so_list[[5]])))
                  
length(all_genes) ## Note: Since no gene filtering was applied each dataset should contain the same number of genes (with a fignificant part not expressed)

set.seed(123) 
# select integration features
features <- SelectIntegrationFeatures(object.list = so_list, nfeatures = 20000) ## I set the number of feature to 20.000 which is roughly the number of expressed genes in at least one condition.

# find integration anchors (CCA)
anchors <- FindIntegrationAnchors(object.list = so_list,
                                  anchor.features = features,verbose = FALSE)

# integrate data
seurat.integrated <- IntegrateData(anchorset = anchors, verbose = FALSE)

## adding information about condition and replicates
seurat.integrated@meta.data$sample<-seurat.integrated$orig.ident
seurat.integrated@meta.data <-separate(seurat.integrated@meta.data, 
                               col=sample,into=c("condition","replicate"),sep="_")
Idents(seurat.integrated)<-0
print(dim(seurat.integrated))
print("Integration DONE !")

## Save the CCA integrated Seurat object
saveRDS(seurat.integrated, file = "/path/to/files/seurat.integrated.rds")
```


### 3. Analysis of the integrated dataset: re-run the Seurat pipeline

```{r}
DefaultAssay(seurat.integrated) <- "integrated"

# Scale data, run PCA and UMAP and visualize integrated data
seurat.integrated <- ScaleData(seurat.integrated, features=rownames(seurat.integrated))
seurat.integrated <- RunPCA(seurat.integrated, npcs = 50,verbose = FALSE)

# ElbowPlot(seurat.integrated,ndims = 30)

DimPlot(seurat.integrated, reduction = "pca")

set.seed(1234) 
ndims <- 30
seurat.integrated <- RunUMAP(seurat.integrated, reduction = "pca", dims = 1:ndims,verbose = FALSE)

DimPlot(seurat.integrated, reduction = "umap")
```


```{r}
g1 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'condition')
g2 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'replicate', cols = c('red','green','blue'))
g1

g2
```

### 4. Cluster Cells after integration

seurat.integrated@meta.data<-dplyr::select(seurat.integrated@meta.data,
                                           -c("RNA_snn_res.0.2","RNA_snn_res.0.5","RNA_snn_res.0.7", "seurat_clusters"))
set.seed(321) 
# Cluster the cells
seurat.integrated<- FindNeighbors(seurat.integrated, reduction = "pca", dims = 1:ndims,verbose=FALSE)
seurat.integrated <- FindClusters(seurat.integrated, resolution = 0.15,verbose=FALSE)

## rename clusters as C0,C1,.... for easy reading

seurat.integrated@meta.data$integrated_snn_res.0.15<-as.factor(paste0("C",seurat.integrated@meta.data$integrated_snn_res.0.15))
Idents(seurat.integrated)<-seurat.integrated@meta.data$integrated_snn_res.0.15 

colors <- c("C0" = "#1F77B4",  # Blue
            "C1" = "#AEC7E8",  # Light Blue
            "C2" = "darkgreen",  # Green
            "C3" = "#98DF8A",  # Light Green
            "C4" = "#D62728",  # Red
            "C5" = "#FF9896",  # Light Red
            "C6" = "#9467BD",  # Purple
            "C7" = "#C5B0D5",  # Light Purple
            "C8" = "#8C564B",  # Brown
            "C9" = "burlywood3",  # Pink
            "C10" = "#FF7F0E",  # Orange
            "C11" = "#FFBB78")  # Light Orange




g2 <-DimPlot(seurat.integrated, reduction = 'umap',label = TRUE, repel = TRUE,  cols = colors,
               group.by = "integrated_snn_res.0.15", pt.size = 0.8, label.size =6)

g2 = g2 +  guides(color = guide_legend(override.aes = list(size=8), ncol=1) ) + ggtitle("Umap")

g3 <- DimPlot(seurat.integrated, reduction = 'umap', 
              split.by = 'condition',label = TRUE, repel = TRUE, cols = colors, pt.size = 0.8,
              group.by = "integrated_snn_res.0.15")
g3= g3 +  guides(color = guide_legend(override.aes = list(size=12), ncol=1) ) + ggtitle("Umap")

setwd("/path/to/files/Plots/")
ggsave(filename="/path/to/files/Plots/UMAP_plot_main.pdf",g2,device = "pdf",
       dpi = 600)
ggsave(filename="/path/to/files/Plots/UMAP_plot.pdf",g3,device = "pdf",width = 20, height = 10,
       dpi = 1000)

g3.n <- DimPlot(seurat.integrated, reduction = 'umap', 
              split.by = 'condition', 
              group.by = "integrated_snn_res.0.15")
g3.n

ggsave(filename="/path/to/files/Plots/UMAP_plot2.pdf",g3,device = "pdf",
       dpi = 600)


table(seurat.integrated@meta.data$integrated_snn_res.0.15)

## Save the CCA integrated Seurat object with processing and clustering
#saveRDS(seurat.integrated, file = "/path/to/files/Plots/rds/seurat.integrated_processed_slim.rds")

# read the file
seurat.integrated= readRDS(file = "/path/to/files/Plots/rds/seurat.integrated_processed.rds")


# see the expression of some markers on featureplots

myfeatures<-c("Ikzf2","Ctla4","Il2ra", "Foxp3","Tnfrsf18")
FeaturePlot(seurat.integrated, features = myfeatures, min.cutoff = "q5")

### 5. Identify Cell-Type Markers for each cluster

all.markers<-FindAllMarkers(seurat.integrated,logfc.threshold = 0.25,min.pct = 0.15,
                              only.pos = TRUE)
  Idents(seurat.integrated) <- seurat.integrated$integrated_snn_res.0.15

  features<-all.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
  
  p5<-DotPlot(seurat.integrated, features = unique(features$gene)) + RotatedAxis()
p5  

### 6. Differential Cell Proportions across clusters and conditions

```{r}
pt <- table(Idents(seurat.integrated), seurat.integrated$orig.ident)

knitr::kable(pt, "simple")

pt <- as.data.frame(pt)
pt$Var1 <- pt$Var1

colors

prop_plot <-ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = colors) +
  theme(legend.title = element_blank())

prop_plot

ggsave(filename="/Volumes/LaCie/sc_mehdi_foxp3/Figure/proportion_plot.pdf",prop_plot,device = "pdf",
       dpi = 600, height = 5, width = 10)

## differential proportion via propeller as in the speckle package

Test.cellprops <-propeller(clusters = Idents(seurat.integrated), 
          sample = seurat.integrated$orig.ident, 
          group = seurat.integrated$condition)

knitr::kable(Test.cellprops)

```


## Output relevant files (used for further analysis)

1. **so_list_slim.rds**: This files containg the list of Seurat object for each sample (after quality control and cleaning).
2. **seurat.integrated_slim.rds**: This file contains the integrated cells (after quality control and cleaning).
3. **seurat.integrated_processed_slim.rds**:  This file contains the integrated cells (after quality control and cleaning) with the cluster information.


```{r}
seurat.integrated <- ScaleData(seurat.integrated)

 #C0   C1   C2   C3   C4   C5   C6   C7   C8   C9  C10  C11 
#5464 4725 3587 1138 1098  458  380  275  150   93   71   57 


cell.list3 = WhichCells(seurat.integrated, idents = c("C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7",
                                                    "C8","C9","C10","C11" ), downsample = 57)
DefaultAssay(seurat.integrated) <- "RNA"
mini = subset(seurat.integrated, cells= cell.list3)
mat<- mini[["RNA"]]@data[unique(features$gene), ] %>% as.matrix()

## scale the rows
mat<- t(scale(t(mat)))

cluster_anno<- mini@meta.data$integrated_snn_res.0.15
# what's the value range in the matrix
quantile(mat, c(0.1, 0.95))

col_fun = circlize::colorRamp2(c(-1, 0, 3), c("blue", "black", "#FFFF00"))
library(ComplexHeatmap)
d=ComplexHeatmap::Heatmap(mat, name = "Expression",  
        column_split = factor(cluster_anno),
        cluster_columns = F,
        show_column_dend = FALSE,
        cluster_column_slices = TRUE,
        column_title_gp = gpar(fontsize = 20),
        column_gap = unit(0.5, "mm"),
        cluster_rows = F,
        show_row_dend = FALSE,
        col = col_fun,
        row_names_gp = gpar(fontsize = 20),
        column_title_rot = 0,
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = colors))),
        show_column_names = FALSE,
        use_raster = TRUE,
        raster_quality = 4)
d

```


```{r}
library(Seurat)
library(ggplot2)
library(ggforce) ## version
genes_of_interest <- "Cd3e"
# Define the y-axis limits
y_limits <- c(0, max(FetchData(seurat.integrated, vars = genes_of_interest), na.rm = TRUE))

# Ciclo per creare i violin plot per ogni gene
for (gene in genes_of_interest) {
  # Estrai i dati per il gene corrente
  gene_data <- FetchData(seurat.integrated, vars = c("condition", gene))
  gene_data$cluster <- Idents(seurat.integrated)
  
  # Genera il violin plot per il gene corrente facettato per cluster
  violin_plot <- ggplot(gene_data, aes(x = condition, y = !!sym(gene), fill = condition)) +
    #geom_violin(scale = "width", trim = FALSE) +
    geom_sina(aes(color = condition), size = 1.5) +
    stat_summary(fun = mean, geom = "point", shape = 3, size = 3, color = "black", alpha = 0.6) +  # Add cross for the mean
    facet_wrap(~ cluster, scales = "free_y") +
    scale_fill_manual(values = c("notch3" = "red", "wt" = "black")) +
    scale_color_manual(values = c("notch3" = "red", "wt" = "black")) +  # Add scale_color_manual for geom_sina and crosses
    coord_cartesian(ylim = y_limits) +  # Set the y-axis limits
    theme_classic(base_family = "Arial") +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16, face = "bold"),
      axis.text.y = element_text(size = 16, face = "bold"),
      axis.title.x = element_text(size = 22, face = "bold"),
      axis.title.y = element_text(size = 18, face =  "bold"),
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
      strip.background = element_rect(color = "black", fill = "grey90"),
      strip.text = element_text(size = 20, face = "bold")
    ) +
    labs(title = gene, x = "Condition", y = "Expression Level") +  scale_x_discrete(limits = c("wt", "notch3"))
  
  # Aggiungi il plot alla lista principale
  violin_plots[[gene]] <- violin_plot
}

# Visualizza i pannelli per ciascun gene
for (gene in genes_of_interest) {
  print(violin_plots[[gene]])
}

```


