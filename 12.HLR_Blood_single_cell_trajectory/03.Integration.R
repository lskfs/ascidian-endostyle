library(Seurat)
library(RColorBrewer)
library(patchwork)
library(data.table)
library(dplyr)
library(ggplot2)
set.seed(1)
HLR<-readRDS("HLRdownsp_outputHLRdownsp_seurat.rds")
blood <- readRDS("stcr-blood-2_outputstcr-blood-2_seurat.rds")

colnames(HLR@meta.data)[colnames(HLR@meta.data) == 'seurat_clusters'] <- 'orig.seurat_clusters'
colnames(blood@meta.data)[colnames(blood@meta.data) == 'seurat_clusters'] <- 'orig.seurat_clusters'

obj.list <- list()
obj.list[[1]] <- blood
obj.list[[2]] <- HLR

## merge
obj.integrated <- merge(obj.list[[1]], obj.list[[2]])

set.seed(1)
obj.integrated <- SCTransform(obj.integrated, variable.features.n = 2000, assay = 'RNA', verbose = FALSE)
obj.integrated <- RunPCA(object = obj.integrated, verbose = FALSE)
obj.integrated <- FindNeighbors(object = obj.integrated)
obj.integrated <- FindClusters(object = obj.integrated, resolution = 0.3)
obj.integrated <- RunUMAP(object = obj.integrated, reduction = "pca", dims = 1:15)

write.table(obj.integrated@meta.data, file=paste0('./BloodHLR.metadata.txt'), quote=T, sep='\t', row.names=F)
saveRDS(obj.integrated, file=paste0('BloodAndHLR.integrated.rds'))

obj.integrated <- readRDS("BloodAndHLR.integrated.rds")


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

DimPlot(obj.integrated, reduction = 'umap', split.by = 'orig.ident', group.by = 'seurat_clusters', label = TRUE, label.size = 8, repel=TRUE, cols = cluster_Palette)
#DimPlot(obj.integrated, reduction = 'umap', split.by = 'orig.seurat_clusters', cols = cluster_Palette)
ggsave('p1_umap.png', width = 24, height = 12)

DimPlot(obj.integrated, reduction = 'umap', group.by = 'orig.seurat_clusters', cols = cluster_Palette)
#DimPlot(obj.integrated, reduction = 'umap', group.by = 'orig.ident', split.by = 'orig.ident', cols = cluster_Palette, ncol = 3)
ggsave('p2_umap.png', width = 14, height = 10)

DimPlot(obj.integrated, reduction = 'umap', group.by = 'orig.ident', cols = cluster_Palette)
ggsave('p3_umap.png', width = 10, height = 10)

DimPlot(obj.integrated, reduction = 'umap', label = TRUE, repel = TRUE, label.size = 8, cols = cluster_Palette)
ggsave('p4_umap.png', width = 12, height = 12)


##### find markers in integrated data
DefaultAssay(obj.integrated) <- 'RNA'

markers <- FindAllMarkers(obj.integrated, min.pct = 0.1, logfc.threshold = 0.25)

markers <- markers[with(markers, order(cluster, -avg_log2FC)), ]
write.table(markers, paste0('R0.3_integrated_AllMarkers.xls'), sep = '\t', quote = FALSE, row.names = F)

topn <- markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
write.table(topn, paste0('R0.3_integrated_AllMarkers.top30.xls'), sep='\t', quote = FALSE, row.names = F)

cluster11to2and9.markers <- FindMarkers(obj.integrated, ident.1=11, ident.2=c(2,9), min.pct=0.25, logfc.threshold = 0.25)
cluster11to2and9.markers <- cluster11to2and9.markers[with(cluster11to2and9.markers, order(-avg_log2FC)), ]
head(cluster11to2and9.markers,n=10)
write.table(cluster11to2and9.markers, paste0('Markers-11to2and9.xls'), sep = '\t', quote = FALSE, row.names = T)

cluster2to11and9.markers <- FindMarkers(obj.integrated, ident.1=2, ident.2=c(11,9), min.pct=0.25, logfc.threshold = 0.25)
cluster2to11and9.markers <- cluster2to11and9.markers[with(cluster2to11and9.markers, order(-avg_log2FC)), ]
head(cluster2to11and9.markers,n=10)
write.table(cluster2to11and9.markers, paste0('Markers-2to11and9.xls'), sep = '\t', quote = FALSE, row.names = T)

cluster9to2and11.markers <- FindMarkers(obj.integrated, ident.1=9, ident.2=c(2,11), min.pct=0.25)
cluster9to2and11.markers <- cluster9to2and11.markers[with(cluster9to2and11.markers, order(-avg_log2FC)), ]
head(cluster9to2and11.markers,n=10)
write.table(cluster9to2and11.markers, paste0('Markers-9to2and11.xls'), sep = '\t', quote = FALSE, row.names = T)



