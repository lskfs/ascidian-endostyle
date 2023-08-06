
library(monocle)
library(ggplot2)
library(dplyr)
library(Seurat)
library(RColorBrewer)

parser = argparse::ArgumentParser(description = 'Script for converting Stereo-seq matrix to seurat format')
parser$add_argument('-i', '--input', dest = 'input', help = 'input h5seurat filename')
parser$add_argument('-m', '--meta', dest = 'meta', help = 'meta filename')
parser$add_argument('-r', '--root', dest = 'root', help = 'root celltype')
parser$add_argument('-s', '--samples', dest = 'sample', help = 'sample ID, will be used as output prefix and seurat object ident')
parser$add_argument('-o', '--out', dest = 'outdir', help = 'directory where to save the output files, all output files will be indexed by sample ID')

opts = parser$parse_args()

obj <- readRDS('../sct.integrated.rds')
obj@meta.data[obj@meta.data$Tissue == 'Thyriod', ]$Tissue <- 'Thyroid'
obj@meta.data <- obj@meta.data[, which(colnames(obj@meta.data) %in% c('CellType', 'Tissue', 'paper', 'seurat_clusters'))]
obj <- subset(obj, subset = Tissue %in% c('Pharyngeal endoderm') | (Tissue %in% 'Thyroid' & CellType == 'Thyrocyte'))

DefaultAssay(obj) <- 'RNA'
Idents(obj) <- 'Tissue'

data <- obj@assays$RNA@counts
pd <- new("AnnotatedDataFrame", data = obj@meta.data)

feature = data.frame(gene_short_name = row.names(data), 
                         row.names = row.names(data))
fd <- new("AnnotatedDataFrame", data = feature)

cds <- newCellDataSet(as(data, 'sparseMatrix'),
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size()
                      )

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# use cluster DEGs as ordering genes
#diff_test_res <- differentialGeneTest(cds, fullModelFormulaStr = "~Tissue")
#ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
#ordering_genes <- row.names(diff_test_res)[order(diff_test_res$qval)][1:2000]
disp_table <- dispersionTable(cds)
disp_table <- subset(disp_table, mean_expression >= 0.3)
ordering_genes <- disp_table$gene_id

print(length(ordering_genes))

cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds)
ggsave('ordering_genes.png')

plot_pc_variance_explained(cds) + geom_vline(xintercept = 10)
ggsave('pc_variance.png')

cds <- reduceDimension(cds, max_components = 3, method = 'DDRTree')
cds <- orderCells(cds)

GM_state <- function(cds){
    if (length(unique(pData(cds)$State)) > 1){
        T0_counts <- table(pData(cds)$State, pData(cds)$Tissue)[,"Pharyngeal endoderm"]
        return(as.numeric(names(T0_counts)[which
                          (T0_counts == max(T0_counts))]))
    } else {
        return (1)
    }
}

cds <- orderCells(cds, root_state = GM_state(cds))

state_Palette <- c('dodgerblue2', '#E31A1C', 'green4', '#6A3D9A', '#FF7F00', 'black', 'gold1', 
                'skyblue2', '#FB9A99', 'palegreen2', '#CAB2D6', '#FDBF6F', 'gray70', 'khaki2', 
                'maroon', 'orchid1', 'deeppink1', 'blue1', 'steelblue4', 'darkturquoise', 
                'green1', 'yellow4', 'yellow3','darkorange4', 'brown')
plot_cell_trajectory(cds, color_by = 'State', cell_size=2) + scale_color_manual(values = state_Palette) 
ggsave('State.png', width=15, height=15)
ggsave('State.pdf', width=15, height=15)

cluster_Palette <- c('Pharyngeal endoderm'='#e2afaf', 'Thyroid'='#a17569')
plot_cell_trajectory(cds, color_by = 'Tissue', cell_size=2) + scale_color_manual(values = cluster_Palette) 
ggsave('tissue.png', width=15, height=15)
ggsave('tissue.pdf', width=15, height=15)

plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 2) + scale_color_viridis_c(option='plasma')
ggsave('pseudotime.png', width=15, height=15)
ggsave('pseudotime.pdf', width=15, height=15)

saveRDS(cds, file = paste0(opts$outdir, '/monocle_cds.rds'))

data.dir <- '/dellfsqd2/ST_OCEAN/USER/hankai/Project/07.Styela_clava/17.downstream/09.SAMap/00.celltypeDGE_newclass/01.celltypeDGE'
throid.markers <- read.csv(paste0(data.dir, '/Thyroid.txt'), sep='\t', header=T)
throid.markers <- subset(throid.markers, subset = (avg_log2FC > 0) & (p_val < 0.05))

endoderm.markers <- read.csv(paste0(data.dir, '/Pharyngeal_endoderm.txt'), sep='\t', header=T)
endoderm.markers <- subset(endoderm.markers, subset = (avg_log2FC > 0) & (p_val < 0.05))

markers <- c(rownames(endoderm.markers), rownames(throid.markers))
markers <- unique(markers)

# generate pseudotime gene heatmap matrix
newdata <- data.frame(Pseudotime = seq(min(pData(cds)$Pseudotime), 
                                       max(pData(cds)$Pseudotime), length.out = 100))
genSmoothCurves_mat <- genSmoothCurves(cds[markers, ], 
                                       new_data = newdata, 
                                       cores = 10)
write.table(genSmoothCurves_mat, 'genSmoothCurves_mat.txt', sep='\t', quote=F)
pdf('sig_gene.pdf')
pheatmap(log10(genSmoothCurves_mat[b$tree_row$order,]+1), 
         scale = "row",
         cluster_rows = F,
         cluster_cols = F,
         #annotation_col = annocol,
         #annotation_colors = annocolor,
         show_rownames = F,
         show_colnames = F, 
         color = hmcols, 
         breaks = bks)
dev.off()




