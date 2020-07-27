library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pbmc.data <- read.delim("/Users/tangxin/Dropbox (Harvard University)/Harvard_research&course/Sig_seq/sig_seq_translation/Enhance_imputation/enhance/denoised-patch-seq-counts.csv", header=TRUE, sep = "\t")
names <- make.unique(as.character(pbmc.data$gene))
rownames(pbmc.data) <- names
pbmc.data <- pbmc.data[,-1] # get rid of old names
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data)
pbmc

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")


# find all markers of cluster 1
#cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
#head(cluster1.markers, n = 5)
# find all markers distinguishing cluster 5 from clusters 0 and 3
#cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, ident.2 = c(0, 3), min.pct = 0.25)
#head(cluster2.markers, n = 5)
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
#cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)



current.cluster.ids <- c(0, 1)
new.cluster.ids <- c('Layer 5','Layer 4')
pbmc@active.ident <- plyr::mapvalues(x =  pbmc@active.ident, from = current.cluster.ids, to = new.cluster.ids)
# cluster re-assignment occurs, which re-assigns clustering in my_levels (assuming you have 12 clusters in total)
my_levels <- c('Layer 5','Layer 4')

# Re-level object@ident
pbmc@active.ident <- factor(x = pbmc@active.ident, levels = my_levels)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 1000, wt = avg_logFC)
p1 <- DoHeatmap(pbmc, features = top10$gene) + NoLegend()
p1_b <- ggplot_build(p1)
p1_b$layout$panel_scales_y[[1]]$range$range
p1 + scale_y_discrete(limits =  rev(p1_b$layout$panel_scales_y[[1]]$range$range))
