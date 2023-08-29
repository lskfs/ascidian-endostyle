
library(future)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

rds.file <- 'zebrafish.inte_1st.rds'
if (!file.exists(rds.file)){
    zf01<-readRDS('../data/zebrafish.develop.rds')
    rownames(zf01@meta.data) <- zf01@meta.data$X
    zf01@meta.data$X <- NULL
    zf01@meta.data$paper <- 'zf01'

    zf04<-readRDS('../data/thyroid.rds')
    zf04@meta.data$paper <- 'zf04'
    zf04@meta.data$Tissue <- 'Thyroid'

    zf06<-readRDS('../data/cochlea.rds')
    zf06@meta.data$paper <- 'zf06'
    zf06@meta.data$Tissue <- 'Cochlea'

    merge.list <- list(zf01, zf04, zf06)
    # normalize and identify variable features for each dataset independently
    merge.list <- lapply(X = merge.list, FUN = function(x) {
                     DefaultAssay(x) <- 'RNA'
                     x <- NormalizeData(x)
                     x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    })
    
    features <- SelectIntegrationFeatures(object.list = merge.list, nfeatures = 2000)
    dataset.anchors <- FindIntegrationAnchors(object.list = merge.list, dims = 1:30, anchor.features = features)
    dataset.merge <- IntegrateData(anchorset = dataset.anchors, dims = 1:30)
    DefaultAssay(dataset.merge) <- "integrated"
    print('integrated finished')

    dataset.merge <- ScaleData(dataset.merge)
    dataset.merge <- RunPCA(dataset.merge)
    dataset.merge <- RunUMAP(dataset.merge, dims = 1:30)
    dataset.merge <- FindNeighbors(dataset.merge, dims = 1:30)
    dataset.merge <- FindClusters(dataset.merge, verbose = FALSE, resolution = 0.8)
    write.table(dataset.merge@meta.data, file = 'zebrafish.inte_1st.metadata.txt', sep = '\t', quote = F, row.names = T)

    saveRDS(dataset.merge, file = rds.file)

    #' save processed seurat object to loom file
    DefaultAssay(dataset.merge) <- 'RNA'
    obj.loom <- as.loom(dataset.merge, filename = paste0('./zebrafish.inte_1st.loom'), 
            overwrite = TRUE)
    obj.loom$close_all()
}else{
    dataset.merge <- readRDS(rds.file)
}

#' create cluster color palette for figure plot function
cluster_number <- length(unique(dataset.merge@meta.data$Tissue))
if (cluster_number <= 25){
    cluster_Palette <- c('dodgerblue2', '#E31A1C', 'green4', '#6A3D9A', '#FF7F00', 'black', 'gold1', 
                         'skyblue2', '#FB9A99', 'palegreen2', '#CAB2D6', '#FDBF6F', 'gray70', 'khaki2', 
                         'maroon', 'orchid1', 'deeppink1', 'blue1', 'steelblue4', 'darkturquoise', 
                         'green1', 'yellow4', 'yellow3','darkorange4', 'brown')
} else if (cluster_number > 25 && cluster_number <= 70){
    qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
    cluster_Palette <- unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
} else if (cluster_number > 70){
    heatmap_Palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
    cluster_Palette <- heatmap_Palette(cluster_number)
}

DimPlot(dataset.merge, reduction = "umap", label = F, repel = F, pt.size = 1.5, cols = cluster_Palette)
ggsave('umap.seurat_clusters.1st.pdf', width = 35, height = 35)

Idents(dataset.merge) <- 'Tissue'
DimPlot(dataset.merge, reduction = "umap", label = F, repel = F, pt.size = 1.5, cols = cluster_Palette)
ggsave('umap.Tissue.1st.pdf', width = 35, height = 35)

DimPlot(dataset.merge, reduction = "umap", group.by = 'Tissue', split.by = 'paper', pt.size = 0.5, cols = cluster_Palette)
ggsave('umap.paper.1st.pdf', width = 30, height = 12)

dataset.merge <- subset(dataset.merge, subset = paper == 'zf01', invert = T)
DimPlot(dataset.merge, reduction = "umap", group.by = 'CellType', split.by = 'paper', pt.size = 0.5, cols = cluster_Palette)
ggsave('umap.CellType.1st.pdf', width = 30, height = 12)


