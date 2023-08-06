
#library(monocle)
library(ggplot2)
library(dplyr)
library(Seurat)
library(RColorBrewer)
#library(pheatmap)
library(seriation)
library(circlize)
library(ComplexHeatmap)

skip <- function(){
parser = argparse::ArgumentParser(description = 'Script for converting Stereo-seq matrix to seurat format')
parser$add_argument('-i', '--input', dest = 'input', help = 'input h5seurat filename')
parser$add_argument('-m', '--meta', dest = 'meta', help = 'meta filename')
parser$add_argument('-r', '--root', dest = 'root', help = 'root celltype')
parser$add_argument('-s', '--samples', dest = 'sample', help = 'sample ID, will be used as output prefix and seurat object ident')
parser$add_argument('-o', '--out', dest = 'outdir', help = 'directory where to save the output files, all output files will be indexed by sample ID')
opts = parser$parse_args()
}

cds <- readRDS('monocle_cds.rds')

#plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 2) + scale_color_viridis_c(option='plasma')
#ggsave('pseudotime.png', width=15, height=15)
#ggsave('pseudotime.pdf', width=15, height=15)

pt.matrix <- read.table('genSmoothCurves_mat.txt', sep = '\t')
pt.matrix <- log10(pt.matrix+1)
pt.matrix <- t(scale(t(pt.matrix)))
pt.matrix <- pt.matrix[complete.cases(pt.matrix), ]
write.table(pt.matrix, file = 'genSmoothCurves_mat.scaled.txt', sep='\t', quote=F)

#genes.to.anno <- scan('genSmoothCurves_genes_to_anno.txt', character())
genes.to.anno <- read.csv('heatmap_trajectory_genes.txt', sep = '\t')
genes.color <- as.character(genes.to.anno$gene.ptclass.color)
names(genes.color) <- genes.to.anno$gene_name

gene.list <- row.names(pt.matrix)
index <- which(gene.list %in% names(genes.color))
labels <- gene.list[index]
colors <- genes.color[labels]

#cluster_colors <- c('c1'='#00FFFF', 'c2'='#FFFF00', 'c3'='#FF00FF')
row_ha = rowAnnotation(gene_id = anno_mark(at = index, labels = labels, labels_gp = gpar(col=colors)))

#Ward.D2 Hierarchical Clustering
pdf('gene_exp_pseudotime.pdf', width = 5, height = 5)
hthc <- Heatmap(
    pt.matrix,
    name = "z-score",
    col = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 36),
    clustering_method_rows = "ward.D2",
    #clustering_method_rows = "average",
    clustering_method_columns = "ward.D2",
    row_title_rot = 0,
    cluster_rows = TRUE,
    #cluster_rows = FALSE,
    cluster_row_slices = FALSE,
    cluster_columns = FALSE,
    right_annotation = row_ha,
	use_raster = F
	#use_raster = TRUE
    )

draw(hthc)
dev.off()



