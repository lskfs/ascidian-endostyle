### Get the parameters
parser = argparse::ArgumentParser(description = 'Script to integrate scRNA and spatial data')
parser$add_argument('-i', dest = 'input', help = 'input directory contains rds files') 
parser$add_argument('-s', dest = 'id', help = 'sample id')
parser$add_argument('-o', dest = 'out', help = 'out directory')
parser$add_argument('-d', dest = 'dims', default = 30, help = 'dims for umap, default 30')
parser$add_argument('-r', dest = 'resolution', default = 0.5, help = 'cluster resolution, default 0.5')
parser$add_argument('-n', dest = 'normalize', default = 'sct', choices = c('sct', 'log'), help = 'normalization method for seurat, sct for sctransform, log for LogNormalize, default sct')
parser$add_argument('-p', dest = 'pointSize', default = 1, help = 'point size in the figure')
opts = parser$parse_args()

###Combine mtx
library(data.table)
library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(patchwork)
library(RColorBrewer)

opts$pointSize <- as.numeric(opts$pointSize)

rds_list <- list.files(opts$input, pattern = '\\.rds$')

obj.list <- list()
for (index in 1:length(rds_list)){
    rds <- rds_list[index]
    rdspath <- paste0(opts$input, '/', rds)
    print(rdspath)
    
    obj <- readRDS(rdspath)

    if (substr(rds, 1, 2) == 'st'){
        obj$Batch <- 'spatial'
        obj@assays$RNA <- obj@assays$Spatial
    }else if (substr(rds, 1, 2) == 'sc'){
        obj$Batch <- 'scRNA'
    }
    DefaultAssay(obj) <- 'RNA'

    colnames(obj@meta.data)[colnames(obj@meta.data) == 'seurat_clusters'] <- 'orig.seurat_clusters'
    
    obj.list[[index]] <- obj
}

#obj.list <- SplitObject(object = merge, split.by = "Batch")

normalize <- function(obj, method = 'sct', nfeatures = 2000){
    if (obj$Batch == 'spatial'){
        #assay = 'Spatial'
        assay = 'RNA'
    }else{
        assay = 'RNA'
    }
    if (method == 'sct'){
        obj <- SCTransform(obj, method = "glmGamPoi", variable.features.n = nfeatures, 
                           assay = assay, verbose = FALSE)
    }else if (method == 'log'){
        obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
        obj <- ScaleData(obj)
        obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures)
    }
    return(obj)
}

for (i in 1:length(x = obj.list)) {
    obj.list[[i]] <- normalize(obj.list[[i]], method = opts$normalize)
    obj.list[[i]] <- RunPCA(obj.list[[i]])
    obj.list[[i]] <- RunUMAP(obj.list[[i]], dims = 1:as.numeric(opts$dims))
}

features <- SelectIntegrationFeatures(object.list = obj.list)
if (opts$normalize == 'log'){
    obj.anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:as.numeric(opts$dims), normalization.method = "LogNormalize", 
                                        anchor.features = features)
    obj.integrated <- IntegrateData(anchorset = obj.anchors, dims = 1:as.numeric(opts$dims), normalization.method = "LogNormalize")
} else if (opts$normalize == 'sct'){
    obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = features)
    obj.anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:as.numeric(opts$dims), normalization.method = "SCT", 
                                        anchor.features = features)
    obj.integrated <- IntegrateData(anchorset = obj.anchors, dims = 1:as.numeric(opts$dims), normalization.method = "SCT")
}

DefaultAssay(object = obj.integrated) <- "integrated"
if (opts$normalize == 'log'){
    obj.integrated <- ScaleData(obj.integrated, verbose = FALSE)
}
obj.integrated <- RunPCA(object = obj.integrated, verbose = FALSE)
obj.integrated <- FindNeighbors(object = obj.integrated)
obj.integrated <- FindClusters(object = obj.integrated, resolution = as.numeric(opts$resolution))
obj.integrated <- RunUMAP(object = obj.integrated, reduction = "pca", dims = 1:as.numeric(opts$dims))

write.table(obj.integrated@meta.data, file=paste0(opts$out, '/', opts$id, '.metadata.txt'), quote=F, sep='\t', row.names=F)

#' create cluster color palette for figure plot function
cluster_number <- length(unique(obj.integrated@meta.data$seurat_clusters))
if (cluster_number <= 25){
    cluster_Palette <- c('dodgerblue2', '#E31A1C', 'green4', '#6A3D9A', '#FF7F00', 'black', 'gold1', 
                         'skyblue2', '#FB9A99', 'palegreen2', '#CAB2D6', '#FDBF6F', 'gray70', 'khaki2', 
                         'maroon', 'orchid1', 'deeppink1', 'blue1', 'steelblue4', 'darkturquoise', 
                         'green1', 'yellow4', 'yellow3','darkorange4', 'brown')
} else if (cluster_number > 25 && cluster_number <= 70){
    qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
    cluster_Palette <- unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
}

DimPlot(obj.integrated, reduction = 'umap', group.by = 'Batch', split.by = 'Batch', cols = cluster_Palette)
ggsave('batch_umap.png', width = 30, height = 12)

DimPlot(obj.integrated, reduction = 'umap', group.by = 'orig.ident', split.by = 'orig.ident', cols = cluster_Palette, ncol = 3)
ggsave('orig.ident_umap.png', width = 30, height = 30)

DimPlot(obj.integrated, reduction = 'umap', label = TRUE, repel = TRUE, cols = cluster_Palette)
ggsave('umap.png', width = 12, height = 12)

saveRDS(obj.integrated, file=paste0(opts$out, '/', opts$id, '.integrated.rds'))

q()
##### find markers in integrated data
DefaultAssay(obj.integrated) <- 'RNA'

markers <- FindAllMarkers(obj.integrated, min.pct = 0.1, logfc.threshold = 0.25)

markers <- markers[with(markers, order(cluster, -avg_log2FC)), ]
write.table(markers, paste0(opts$out, '/', opts$id, '.integrated_AllMarkers.xls'), sep = '\t', quote = FALSE, row.names = F)

topn <- markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
write.table(topn, paste0(opts$out, '/', opts$id, '.integrated_AllMarkers.top30.xls'), sep='\t', quote = FALSE, row.names = F)

##### find markers in scRNA data
sc <- subset(obj.integrated, subset = Batch == 'scRNA')
DefaultAssay(sc) <- 'RNA'

markers <- FindAllMarkers(sc, min.pct = 0.1, logfc.threshold = 0.25)

markers <- markers[with(markers, order(cluster, -avg_log2FC)), ]
write.table(markers, paste0(opts$out, '/', opts$id, '.sc_AllMarkers.xls'), sep = '\t', quote = FALSE, row.names = F)

topn <- markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
write.table(topn, paste0(opts$out, '/', opts$id, '.sc_AllMarkers.top30.xls'), sep='\t', quote = FALSE, row.names = F)

##### find markers in spatial data
st <- subset(obj.integrated, subset = Batch == 'spatial')
DefaultAssay(st) <- 'RNA'

markers <- FindAllMarkers(st, min.pct = 0.1, logfc.threshold = 0.25)

markers <- markers[with(markers, order(cluster, -avg_log2FC)), ]
write.table(markers, paste0(opts$out, '/', opts$id, '.st_AllMarkers.xls'), sep = '\t', quote = FALSE, row.names = F)

topn <- markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
write.table(topn, paste0(opts$out, '/', opts$id, '.st_AllMarkers.top30.xls'), sep='\t', quote = FALSE, row.names = F)

