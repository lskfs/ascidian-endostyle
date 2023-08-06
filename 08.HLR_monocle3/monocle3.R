
library(monocle3)
library(ggplot2)
library(dplyr)
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(RColorBrewer)

parser = argparse::ArgumentParser(description = 'Script for converting Stereo-seq matrix to seurat format')
parser$add_argument('-i', '--input', dest = 'input', help = 'input h5seurat filename')
parser$add_argument('-m', '--meta', dest = 'meta', help = 'meta filename')
parser$add_argument('-r', '--root', dest = 'root', help = 'root celltype')
parser$add_argument('-s', '--samples', dest = 'sample', help = 'sample ID, will be used as output prefix and seurat object ident')
parser$add_argument('-o', '--out', dest = 'outdir', help = 'directory where to save the output files, all output files will be indexed by sample ID')
opts = parser$parse_args()

if (!file.exists('sinus_selectted.rds')){
    obj <- readRDS('/dellfsqd2/ST_OCEAN/USER/hankai/Project/07.Styela_clava/17.downstream/00.data/Styela_clava.anno.rds')
    DefaultAssay(obj) <- 'RNA'
    Idents(obj) <- 'short_inte_anno'

    selectted_cells <- read.csv('../selectted_cells.obs.txt', sep='\t', header=T)
    obj <- subset(obj, cells = rownames(selectted_cells))
    #obj <- subset(obj, subset = Batch %in% c('sc1', 'sc2', 'sc3'))
    all.cells <- colnames(obj)

    obj.list <- SplitObject(obj, split.by = 'Batch')
    obj.list <- lapply(X = obj.list, FUN = function(x) {
                       x <- NormalizeData(x)
                       x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
                       x <- SCTransform(x, assay = 'RNA', return.only.var.genes = F)
    })

    features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 2000)
    #obj.anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:15, 
    #                                      anchor.features = features)
    #obj <- IntegrateData(anchorset = obj.anchors, dims = 1:15) 
    obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = features)
    obj.anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:30, normalization.method = "SCT", 
                                          anchor.features = features)
    obj <- IntegrateData(anchorset = obj.anchors, dims = 1:30, normalization.method = "SCT")
    DefaultAssay(obj) <- "integrated"

    #obj <- ScaleData(obj)
    obj <- RunPCA(obj, npcs = 30)
    obj <- FindNeighbors(obj, dims = 1:30)
    #obj <- RunUMAP(obj, dims = 1:10, min.dist = 0.01)
    obj <- RunUMAP(obj, dims = 1:30)

    .f <- function(){
    data <- obj@assays$RNA@counts
    cds <- new_cell_data_set(as(data, 'sparseMatrix'),
                         cell_metadata = obj@meta.data, 
                         gene_metadata = data.frame(gene_short_name = row.names(data), 
                                                    row.names = row.names(data)
                                                    )
                         )
    obj <- FindVariableFeatures(obj, nfeatures=2000)
    hvg <- VariableFeatures(obj)
    cds <- preprocess_cds(cds, num_dim = 10, use_genes = hvg)

    cds <- align_cds(cds, alignment_group = 'timepoint', alignment_k = 100)
    cds <- reduce_dimension(cds, preprocess_method = 'PCA', reduction_method = 'UMAP', 
                        umap.min_dist=0.005, umap.n_neighbors=100) #, spread=0.1)
    cds <- cluster_cells(cds)
    cds <- learn_graph(cds, use_partition = FALSE)
    }

    #obj <- SCTransform(obj, assay = 'RNA', verbose = FALSE, variable.features.n = 1000) #, vars.to.regress = 'cell_number')
    #obj <- RunPCA(obj, verbose = FALSE, npcs = 8)
    #obj <- RunUMAP(obj, dims = 1:8, verbose = FALSE, min.dist = 0.01)
    #obj <- FindNeighbors(obj, dims = 1:8, verbose = FALSE)

    cds <- as.cell_data_set(obj)
    cds <- cluster_cells(cds = cds, reduction_method = 'UMAP')
    cds <- learn_graph(cds, use_partition = TRUE)

    batch_Palette <- c('C3-2'='#e64b35', 
                       'C3-3'='#4dbbd5',
                       'C3-4'='#00a087',
                       'C3-5'='#3c5488',
                       'C3-6'='#f39b7f',
                       'C3-8'='#8491b4',
                       'sc1'='#91d1c2',
                       'sc2'='#dc0000',
                       'sc3'='#7e6148'
                      )
    p2 <- plot_cells(cds, color_cells_by = 'Batch', 
                 reduction_method = 'UMAP',
                 label_groups_by_cluster = FALSE, 
                 label_leaves = FALSE, 
                 label_branch_points = FALSE,
                 graph_label_size = 1.5,
                 cell_size = 2,
                 cell_stroke = 0) +
        scale_color_manual(values = batch_Palette) +
        guides(colour = guide_legend(override.aes = list(size=5), nrow = 3))
    ggsave(plot = p2, file = 'batch.png')
    ggsave(plot = p2, file = 'batch.pdf')

    cluster_Palette <- c('BC'='#CAB2D6', 'CPC'='#008B00', 'CSC'='#1C86EE', 'LLF'='#FF83FA', 'M'='#7EC0EE',
                         'MLLC'='#FF7F00', 'NC'='#FFD700', 'NCSC'='#E31A1C', 'PLLC'='#B3B3B3', 
                         'SCA'='#90EE90', 'SCB'='#36648B', 'SLLC'='#B03060'
                         )
    p3 <- plot_cells(cds, color_cells_by = 'short_inte_anno', 
                 reduction_method = 'UMAP',
                 label_groups_by_cluster = FALSE, 
                 label_leaves = FALSE, 
                 label_branch_points = FALSE,
                 graph_label_size = 1.5,
                 cell_size = 2,
                 cell_stroke = 0) +
        scale_color_manual(values = cluster_Palette)
    ggsave(plot = p3, file = 'celltype.png')
    ggsave(plot = p3, file = 'celltype.pdf')

    get_earliest_principal_node <- function(cds, time_bin = '0'){

        cell_ids <- which(colData(cds)[, 'timepoint'] %in% c('36hpa1', '36hpa2'))
    
        closest_vertex <- cds@principal_graph_aux[['UMAP']]$pr_graph_cell_proj_closest_vertex
        closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
        root_pr_nodes <- igraph::V(principal_graph(cds)[['UMAP']])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
        root_pr_nodes

    }
    #cds <- order_cells(cds, root_pr_nodes = get_earliest_principal_node(cds, time_bin = opts$root), reduction_method = 'UMAP')
    #stem.cells <- scan('../anno_stem_cell_in_stereo/stemcell.index.txt', what = character())
    stem.cells <- read.table('../anno_stem_cell_in_stereo/stemcell.index.txt', sep='\t', header=F)
    colnames(stem.cells) <- 'cells'
    stem.cells <- subset(stem.cells, subset = cells %in% all.cells)
    cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = stem.cells$cells)
    p4 <- plot_cells(cds, color_cells_by = 'pseudotime',
                 reduction_method = 'UMAP',
                 label_cell_groups = FALSE,
                 label_leaves = TRUE,
                 label_branch_points = TRUE,
                 graph_label_size = 1.5, 
                 cell_size = 2,
                 cell_stroke = 0) +
        scale_color_viridis_c(option='plasma')
    ggsave(plot = p4, file = 'pseudotime.png')
    ggsave(plot = p4, file = 'pseudotime.pdf')

    rowData(cds)$gene_name <- rownames(cds)
    rowData(cds)$gene_short_name <- rowData(cds)$gene_name

    cds <- cluster_cells(cds, k = 15) # hankai cluster 3
    #cds <- cluster_cells(cds, k = 11) # hankai cluster 4
    cds$clusters <- clusters(cds)
    cds$partitions <- partitions(cds, reduction_method='UMAP')
    p3 <- plot_cells(cds, color_cells_by = 'clusters', 
                 reduction_method = 'UMAP',
                 label_groups_by_cluster = FALSE, 
                 label_leaves = FALSE, 
                 label_branch_points = FALSE,
                 graph_label_size = 1.5,
                 cell_size = 2,
                 cell_stroke = 0)
    ggsave(plot = p3, file = 'clusters.png')
    ggsave(plot = p3, file = 'clusters.pdf')

    saveRDS(cds, file = './monocle_cds.rds')

    meta <- data.frame(cds@colData)
    meta$clusters <- as.character(meta$clusters)
    meta <- meta['clusters']
    print(head(meta))
    obj <- AddMetaData(obj, meta, col.name = 'clusters')
    saveRDS(obj, 'sinus_selectted.rds')
}else{
    obj <- readRDS('sinus_selectted.rds')
    cds <- readRDS('monocle_cds.rds')
}


