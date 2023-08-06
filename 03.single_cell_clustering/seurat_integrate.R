### Get the parameters
parser = argparse::ArgumentParser(description = "Script to QC and Cluster scRNA data")
parser$add_argument('-i', '--input', help = 'input directory')
parser$add_argument('-o', '--out', help = 'out directory')
parser$add_argument('--dims', dest = 'dims', default = 20, help = 'number of component used')
parser$add_argument('--resolution', dest = 'resolution', default = 0.8, help = 'resolution for clustering, default 0.8')
opts = parser$parse_args()

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(RColorBrewer)

opts$dims <- as.numeric(opts$dims)
opts$resolution <- as.numeric(opts$resolution)

rds_list <- list.files(opts$input, pattern = '\\.rds$')

#vf <- vector()
merge.list <- list()
for (index in 1:length(rds_list)){
    rds <- rds_list[index]
    print(rds)
    
    obj <- readRDS(rds)
    obj$Batch <- obj$orig.ident
    
    colnames(obj@meta.data)[colnames(obj@meta.data) == 'seurat_clusters'] <- 'orig.seurat_clusters'
    
    #vf <- append(vf, VariableFeatures(obj))
    merge.list[[index]] <- obj
}
print('read in finished')

features <- SelectIntegrationFeatures(object.list = merge.list, nfeatures = 2000)
obj.list <- PrepSCTIntegration(object.list = merge.list, anchor.features = features)
dataset.anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:opts$dims, normalization.method = "SCT", 
                                      anchor.features = features)
dataset.merge <- IntegrateData(anchorset = dataset.anchors, dims = 1:opts$dims, normalization.method = "SCT")
DefaultAssay(dataset.merge) <- "integrated"

#dataset.merge <- merge(x = merge.list[[1]], y = merge.list[2:length(merge.list)])
#DefaultAssay(dataset.merge) <- "SCT"
#VariableFeatures(dataset.merge) <- vf

dataset.merge <- RunPCA(dataset.merge)
dataset.merge <- FindNeighbors(dataset.merge, dims = 1:opts$dims)
dataset.merge <- FindClusters(dataset.merge, verbose = FALSE, resolution = opts$resolution)
dataset.merge <- RunUMAP(dataset.merge, dims = 1:opts$dims)

write.table(dataset.merge@meta.data, file = 'dataset.merge.metadata.txt', sep = '\t', quote = F)

#' create heatmap and cluster color palette for figure plot function
cluster_number <- length(levels(dataset.merge@meta.data$seurat_clusters))

if (cluster_number <= 25){
    cluster_Palette <- c('dodgerblue2', '#E31A1C', 'green4', '#6A3D9A', '#FF7F00', 'black', 'gold1', 
                         'skyblue2', '#FB9A99', 'palegreen2', '#CAB2D6', '#FDBF6F', 'gray70', 'khaki2', 
                         'maroon', 'orchid1', 'deeppink1', 'blue1', 'steelblue4', 'darkturquoise', 
                         'green1', 'yellow4', 'yellow3','darkorange4', 'brown')
} else if (cluster_number > 25 && cluster_number <= 70){
    qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
    cluster_Palette <- unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
}

DimPlot(dataset.merge, reduction = "umap", group.by = 'Batch', split.by = 'Batch', ncol = 6, pt.size = 1.5)
ggsave('batch.png', width = 30, height = 12)

DimPlot(dataset.merge, reduction = "umap", label = TRUE, repel = TRUE, split.by = 'Batch', ncol = 6, cols = cluster_Palette, pt.size = 1.5)
ggsave('cluster.png', width = 30, height = 12)

#Idents(dataset.merge) <- 'orig.seurat_clusters'
#DimPlot(dataset.merge, reduction = "umap", label = TRUE, repel = TRUE, split.by = 'Batch', ncol = 6, cols = cluster_Palette, pt.size = 1.5)
#ggsave('orig.ident.png', width = 30, height = 12)

#Idents(dataset.merge) <- 'seurat_clusters'

DefaultAssay(dataset.merge) <- 'RNA'
saveRDS(dataset.merge, file = paste0(opts$out, "/sct.integrated.rds"))

markers <- FindAllMarkers(dataset.merge, min.pct = 0.1, logfc.threshold = 0.25)

#ref <- read.table("/dellfsqd2/ST_OCEAN/USER/hankai/Project/06.Axolotl_spinal_cord/04.stereo-seq/Axolotl_mm_hg_gene_0908.xls",header = T,stringsAsFactors = FALSE)
#markers <- merge(markers, ref, by.x="gene", by.y="Axolotl_ID", sort=F)

markers <- markers[with(markers, order(cluster, -avg_log2FC)), ]
write.table(markers, paste0(opts$out, '/sct.integrated_AllMarkers.xls'), sep = '\t', quote = FALSE)

topn <- markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
write.table(topn, 'sct.integrated_AllMarkers.top30.xls', sep='\t', quote = FALSE, row.names = F)

