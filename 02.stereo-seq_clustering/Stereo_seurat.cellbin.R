########
#
#   Stereo_seurat.R
#   
#   1. load stereo-seq format matrix and convert to seurat object
#   2. add code to do seurat process,
#
########

### Get the parameters
parser = argparse::ArgumentParser(description = 'Script for converting Stereo-seq matrix to seurat format')
parser$add_argument('-i', '--input', dest = 'input', help = 'input tsv filename')
parser$add_argument('-s', '--sample', dest = 'sample', help = 'sample name')
parser$add_argument('-m', '--metadata', dest = 'metadata', help = 'metadata of cell')
parser$add_argument('-o', '--out', dest = 'outdir', help = 'directory where to save the output files, all output files will be indexed by sample ID')

parser$add_argument('--minCount', dest = 'minCount', default = 0, type = 'integer', help = 'minimum UMI number')
parser$add_argument('--maxCount', dest = 'maxCount', type = 'integer', help = 'maximum UMI number')
parser$add_argument('--minFeature', dest = 'minFeature', default = 0, type = 'integer', help = 'minimum Feature number')
parser$add_argument('--maxFeature', dest = 'maxFeature', type = 'integer', help = 'maximum Feature number')
parser$add_argument('--vg', dest = 'vg', default = 3000, type = 'integer', help = 'number of variable genes, default 3000')
parser$add_argument('--pc', dest = 'pc', default = 30, type = 'integer', help = 'number of PC to use, default 30')
parser$add_argument('--resolution', dest = 'resolution', default = 0.8, help = 'cluster resolution, default 0.8')
parser$add_argument('--topn', dest = 'topn', help = 'save the top N markers')

parser$add_argument('--pointSize', dest = 'pointSize', default = 1, help = 'point size of spatial plot, default 0.2')
parser$add_argument('--colors', dest = 'colors', default = 70, type = 'integer', help = 'colors palette, one of c(25, 70), default 70')
opts = parser$parse_args()

library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(dplyr)
library(data.table)
library(Matrix)
library(rjson)
library(RColorBrewer)

opts$pointSize <- as.numeric(opts$pointSize)
opts$resolution <- as.numeric(opts$resolution)
dir.create(opts$outdir, recursive=TRUE)

CellGem2Seurat <- function(cellgem, meta.data = NULL, orig.ident = 'sample', assay = 'Spatial'){

    data <- fread(file = cellgem)

    if ('MIDCounts' %in% colnames(data)) {
        data <- data %>% rename(MIDCount = MIDCounts)
    } else if ('UMICount' %in% colnames(data)) {
        data <- data %>% rename(MIDCount = UMICount)
    } else if ('UMICounts' %in% colnames(data)) {
        data <- data %>% rename(MIDCount = UMICounts)
    }

    if ('label' %in% colnames(data)) {
        data <- data %>% rename(cell = label)
    }

    data <- data[, .(counts=sum(MIDCount)), by = .(geneID, cell)]
    #' create sparse matrix from stereo
    data$geneIdx <- match(data$geneID, unique(data$geneID))
    data$cellIdx <- match(data$cell, unique(data$cell))

    mat <- sparseMatrix(i = data$geneIdx, j = data$cellIdx, x = data$counts, 
                        dimnames = list(unique(data$geneID), unique(data$cell)))

    obj <- CreateSeuratObject(counts = mat, project = orig.ident, assay = assay)
    if (!is.null(meta.data)){
        obj <- AddMetaData(object = obj, metadata = meta.data)
    }

    #' filter out empty cell
    obj <- subset(obj, subset = nCount_Spatial > 0)

    ############# here we finished the creatation of seurat spatial object using pseudo image
    return(obj)
}

#' create heatmap and cluster color palette for figure plot function
heatmap_Palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
#' custom spatial plot on specific feature
iPlot <- function(object, features, pt.size = 0.2){
    plot <- ggplot(object@meta.data, aes_string(x = 'cx', y = 'cy', color = features)) +
            geom_point(shape = 19, size = pt.size) + 
            theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
                  axis.title = element_blank(), axis.line = element_blank(), legend.position = 'right') + 
            coord_fixed()
    if (features %in% c('nCount_Spatial', 'nFeature_Spatial')){
        plot <- plot + scale_color_gradientn(colours = heatmap_Palette(100))
    }else if(features %in% c('seurat_clusters')){
        plot <- plot + scale_color_manual(values = cluster_Palette) +
                guides(colour = guide_legend(override.aes = list(size=3), nrow = 15))
    }
    plot <- plot + theme_void()
    return(plot)
}

#' basic statistics plot function
Preview <- function(object, feature, width = 400, height = 200){
    p1 <- VlnPlot(object, features = feature, pt.size = 0.2) + NoLegend() + theme(axis.text.x=element_blank(), axis.title.x=element_blank())
    #p2 <- SpatialFeaturePlot(obj, pt.size.factor = opts$pointSize, features = feature, stroke = 0) + theme(legend.position = 'right')
    ggsave(paste0(opts$outdir, '/', opts$sample, '_', feature, '_VlnPlot.pdf'), p1)
    ggsave(paste0(opts$outdir, '/', opts$sample, '_', feature, '_VlnPlot.png'), p1)
    p2 <- iPlot(object, features = feature, pt.size = opts$pointSize)
    ggsave(paste0(opts$outdir, '/', opts$sample, '_', feature, '_spatial.pdf'), p2, width = width, height = height, units = 'px')
    ggsave(paste0(opts$outdir, '/', opts$sample, '_', feature, '_spatial.png'), p2, width = width, height = height, units = 'px')
    #patch <- p1 | p2
    #return(patch)
}

#' seurat clustering workflow
Clustering <- function(object, dims = 30){
    object <- RunPCA(object)
    object <- FindNeighbors(object, dims = 1:dims)
    object <- FindClusters(object, verbose = FALSE, resolution = opts$resolution)
    object <- RunUMAP(object, dims = 1:dims)
    return(object)
}

############ start seurat processing

meta.data <- read.table(opts$metadata, sep = '\t', header = TRUE)
meta.data <- meta.data %>% select(-one_of('nCount_Spatial', 'nFeature_Spatial'))
if ('Batch' %in% colnames(meta.data)){
    meta.data <- subset(meta.data, subset = Batch == opts$sample)
}
row.names(meta.data) <- meta.data$label

obj <- CellGem2Seurat(opts$input, meta.data = meta.data, orig.ident = opts$sample, assay = 'Spatial')
print(head(obj@meta.data))

#' filter object based on nCount and nFeature
if (is.null(opts$maxCount)){
    opts$maxCount <- max(obj$nCount_Spatial)
}
if (is.null(opts$maxFeature)){
    opts$maxFeature <- max(obj$nFeature_Spatial)
}
obj <- subset(obj, subset = nCount_Spatial >= opts$minCount & nCount_Spatial <= opts$maxCount & 
              nFeature_Spatial >= opts$minFeature & nFeature_Spatial <= opts$maxFeature)

#' basic statistics plot on nCount and nFeature
Preview(obj, 'nCount_Spatial')
print('nCount_Spatial done')
Preview(obj, 'nFeature_Spatial')
print('nFeature_Spatial done')

#' log normalize for compare
obj <- NormalizeData(obj, verbose = FALSE, assay = 'Spatial')
obj <- ScaleData(obj, verbose = FALSE, assay = 'Spatial')
#' do normalization and clustering
obj <- SCTransform(obj, assay = 'Spatial', variable.features.n = as.numeric(opts$vg), 
                   return.only.var.genes = FALSE, n_genes=NULL, min_cells=5, method='qpoisson')
DefaultAssay(obj) <- 'SCT'
obj <- Clustering(obj, dims = as.numeric(opts$pc))

write.table(obj@meta.data, paste0(opts$outdir, '/', opts$sample, '_metadata.txt'), 
                sep='\t', quote = FALSE, row.names = F)

#' save processed seurat object to rds file
saveRDS(obj, file = paste0(opts$outdir, '/', opts$sample, '_seurat.rds'))
#' save processed seurat object to loom file
obj.loom <- as.loom(obj, filename = paste0(opts$outdir, '/', opts$sample, '_seurat.loom'), 
        overwrite = TRUE)
obj.loom$close_all()

n_clusters <- length(unique(obj@meta.data$seurat_clusters))
if (n_clusters <= 25){
    cluster_Palette <- c('dodgerblue2', '#E31A1C', 'green4', '#6A3D9A', '#FF7F00', 'black', 'gold1', 
                         'skyblue2', '#FB9A99', 'palegreen2', '#CAB2D6', '#FDBF6F', 'gray70', 'khaki2', 
                         'maroon', 'orchid1', 'deeppink1', 'blue1', 'steelblue4', 'darkturquoise', 
                         'green1', 'yellow4', 'yellow3','darkorange4', 'brown')
} else if (n_clusters <= 70){
    qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
    cluster_Palette <- unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
} else if (n_clusters > 70){
    cluster_Palette <- heatmap_Palette(n_clusters)
}
#' plot cluster results on umap and spatial
p3 <- DimPlot(obj, reduction = 'umap', label = TRUE, cols = cluster_Palette) + 
                guides(colour = guide_legend(override.aes = list(size=3), nrow = 15))
ggsave(paste0(opts$outdir, '/', opts$sample, '_UMAP.pdf'), p3)
ggsave(paste0(opts$outdir, '/', opts$sample, '_UMAP.png'), p3)

#' find all markers of each seurat_cluster and save to file
#markers <- FindAllMarkers(obj, assay = 'Spatial', min.pct = 0.1, logfc.threshold = 0.25)
#markers <- markers[with(markers, order(cluster, -avg_log2FC)), ]
#write.table(markers, paste0(opts$outdir, '/', opts$sample, '_AllMarkers.xls'), sep = '\t', quote = FALSE)

#' save the top N markers to file
#if (!is.null(opts$topn)){
#    topn <- markers %>% group_by(cluster) %>% top_n(n = opts$topn, wt = avg_log2FC)
#    write.table(topn, paste0(opts$outdir, '/', opts$sample, '_AllMarkers.top', opts$topn, '.xls'), sep='\t', quote = FALSE, row.names = F)
#}

