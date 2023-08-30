#####################
# 	
#    Monocle2.R
#    http://cole-trapnell-lab.github.io/monocle-release/docs/
#
######## INPUT #######
#
#    Path to rds file containing Seurat object   http://cole-trapnell-lab.github.io/monocle-release/docs/
#
#####################
library(monocle)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(argparse)
library(Seurat)
library(SeuratObject)
library(reshape2)

opts <- list()
opts$input <-  "BloodAndHLR.integrated.2-9-11.rds"
opts$output <- "BloodAndHLRdownsp-Monocle2"
opts$projectname <- "BloodAndHLRdownsp-2911-Monocle2"

seurat_object <- readRDS(opts$input)
seurat_object@meta.data$SampleCluster <- paste("Cluster",seurat_object@meta.data$seurat_clusters,":",seurat_object@meta.data$orig.ident, sep="")

exprs <- as.matrix(seurat_object@assays$RNA@counts) #if you do have UMI data, you should not normalize it yourself prior to creating your cds # Monocle also works "out-of-the-box" with the transcript count matrices produced by CellRanger
phenoData <- data.frame(seurat_object@meta.data$orig.ident, seurat_object@meta.data$seurat_clusters, seurat_object@meta.data$SampleCluster, row.names = row.names(seurat_object@meta.data)) # https://stackoverflow.com/questions/63344149/how-to-fix-duplicate-column-names-preventing-plot-cellscds-from-running
featureData <- data.frame(gene_short_name = row.names(exprs), row.names = row.names(exprs))

#colnames(exprs)<-paste0("X", colnames(exprs))
#rownames(phenoData)<-paste0("X", rownames(phenoData))

pd <- new("AnnotatedDataFrame", data = phenoData)
fd <- new("AnnotatedDataFrame", data = featureData)
cds <- newCellDataSet(as.matrix(exprs),
        phenoData = pd, featureData = fd)

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# Filtering low-quality cells
cds <- detectGenes(cds, min_expr = 0.1)
print(head(fData(cds)))


expressed_genes <- row.names(subset(fData(cds),
    num_cells_expressed >= 10))

print(head(pData(cds)))

pData(cds)$Total_mRNAs <- Matrix::colSums(exprs(cds))

cds <- cds[,pData(cds)$Total_mRNAs < 1e6]

upper_bound <- 10^(mean(log10(pData(cds)$Total_mRNAs)) +
            2*sd(log10(pData(cds)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(cds)$Total_mRNAs)) -
            2*sd(log10(pData(cds)$Total_mRNAs)))

pdf(paste0(opts$output,"01_Total_mRNAs.pdf"))
qplot(Total_mRNAs, data = pData(cds),  geom = "density") +
geom_vline(xintercept = lower_bound) +
geom_vline(xintercept = upper_bound)
dev.off()

cds <- cds[,pData(cds)$Total_mRNAs > lower_bound &
      pData(cds)$Total_mRNAs < upper_bound]
cds <- detectGenes(cds, min_expr = 0.1)

# Log-transform each value in the expression matrix.
L <- log(exprs(cds[expressed_genes,]))

# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))

# Plot the distribution of the standardized gene expression values.
pdf(paste0(opts$output,"02_Stadardized_gene_expression_values.pdf"))
qplot(value, geom = "density", data = melted_dens_df) +
stat_function(fun = dnorm, size = 0.5, color = 'red') +
xlab("Standardized log(FPKM)") +
ylab("Density")
dev.off()

#---Done//--Get started and pre-processing--//Done---#


#-------------Classify and count cells---------------#
# Classifying and Counting Cells - Clustering cells without marker genes
disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
pdf(paste0(opts$output,"03_plot_ordering_genes.pdf"))
plot_ordering_genes(cds)
dev.off()
pdf(paste0(opts$output,"04_plot_pc_variance_explained.pdf"))
plot_pc_variance_explained(cds, return_all = F) # norm_method='log'
dev.off()

cds <- reduceDimension(cds, max_components = 2, num_dim=6, reduction_method = "tSNE", verbose = T)
cds <- clusterCells(cds, num_cluster = 5)
pdf(paste0(opts$output,"05_plot_cell_clusters.pdf"))
plot_cell_clusters(cds, 1, 2)
dev.off()
# plot_cell_clusters(cds, 1, 2, color = "CellType", markers = c("CDC20", "GABPB2"))

#----Done-Classify and count cells-Done--------#

saveRDS(cds, file = paste0(opts$output,"ClassifiedCellImage.rds"))
cds<-readRDS(file = paste0(opts$output,"ClassifiedCellImage.rds"))

#------- Constructing Single Cell Trajectories-------#
# Trajectory step 1: choose genes that define a cell's progress


diff_test_res <- differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr  = "~num_genes_expressed", cores = 4)
ordering_genes <- row.names (subset(diff_test_res, qval < 0.1)) 
cds <- setOrderingFilter(cds, ordering_genes)

pdf(paste0(opts$output,"06_plot_ordering_gene.pdf"))
plot_ordering_genes(cds)
dev.off()

# Trajectory step 2: reduce data dimensionality
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')

saveRDS(cds, file = paste0(opts$output,"ClassifiedCellImage-ForTrajectoryStep3.rds"))

cds<-readRDS(file = paste0(opts$output,"ClassifiedCellImage-ForTrajectoryStep3.rds"))

# Trajectory step 3: order cells along the trajectory
cds <- orderCells(cds)
pdf(paste0(opts$output,"07_plot_cell_trajectory.pdf"))
plot_cell_trajectory(cds,color_by="State")
dev.off()

#pdf(paste0(opts$output,"08_plot_cell_trajectory.pdf"))
#plot_cell_trajectory(cds, color_by = "State") +  facet_wrap(~State, nrow = 1)
#dev.off()
pdf(paste0(opts$output,"08_plot_cell_trajectory.pdf"),width=10,height=8)
plot_cell_trajectory(cds, color_by = "seurat_object.meta.data.orig.ident") +  
	facet_wrap(seurat_object.meta.data.orig.ident~seurat_object.meta.data.seurat_clusters, nrow = 2, ncol= 3) +
	scale_color_manual(breaks = c('HLRdownsp', 'stcr-blood-2'),values=c("#327BC0", "#D5231E")) +
	theme(legend.position = "right")
dev.off()

pdf(paste0(opts$output,"083_plot_cell_trajectory.pdf"),width=10,height=4)
plot_cell_trajectory(cds, color_by = "seurat_object.meta.data.orig.ident") +
        facet_wrap(~seurat_object.meta.data.seurat_clusters, nrow = 2, ncol= 3) +
        scale_color_manual(breaks = c('HLRdownsp', 'stcr-blood-2'),values=c("#327BC0", "#D5231E")) +
        theme(legend.position = "right")
dev.off()

pdf(paste0(opts$output,"082_plot_cell_trajectory.pdf"),width=12,height=5)
plot_cell_trajectory(cds, color_by = "seurat_object.meta.data.seurat_clusters", ) +  facet_wrap(~seurat_object.meta.data.seurat_clusters, nrow = 1)
dev.off()

pdf(paste0(opts$output,"084_plot_cell_trajectory.pdf"),width=9,height=6)
plot_cell_trajectory(cds, color_by="seurat_object.meta.data.SampleCluster") + theme(legend.position = "right")
dev.off()


pdf(paste0(opts$output,"09_plot_cell_trajectory.pdf"))
plot_cell_trajectory(cds, show_cell_names = F, color_by = "seurat_object.meta.data.orig.ident", cell_size = 0.5) +
scale_color_manual(values=c("#327BC0", "#D5231E"))
dev.off()

pdf(paste0(opts$output,"092_plot_cell_trajectory.pdf"),width=12, height=6)
plot_cell_trajectory(cds, show_cell_names = F, color_by = "seurat_object.meta.data.orig.ident") +  facet_wrap(~seurat_object.meta.data.orig.ident, nrow = 1) +
scale_color_manual(values=c("#006b7b", "#ef1828"))
dev.off()

pdf(paste0(opts$output,"10_plot_cell_trajectory.pdf"),width=6,height=6)
plot_cell_trajectory(cds, show_cell_names = F, color_by = "seurat_object.meta.data.seurat_clusters") +
scale_color_manual(values=c("#1A863B", "#96CC88", "#F5BB6F"))
dev.off()

pdf(paste0(opts$output,"102_plot_cell_trajectory.pdf"), width=14, height=5)
plot_cell_trajectory(cds, show_cell_names = F, color_by = "seurat_object.meta.data.seurat_clusters") + facet_wrap(~seurat_object.meta.data.seurat_clusters, nrow = 5, ncol = 6)
dev.off()

# gene choosing  

to_be_tested <- row.names(subset(fData(cds),
	gene_short_name %in% c("evm.model.000026F-arrow-pilon.83.update",  # CLDN7  Claudin-7
				"evm.model.000001F-arrow-pilon.43.update",  # EGR1  Early growth response-1
				"evm.model.000115F-arrow-pilon.27.update",  # GSC2  Homeobox protein goosecoid2
				"evm.model.000098F-arrow-pilon.52.update",  # CTNNB1	catenin beta-1-like [Styela clava]
				"evm.model.000044F_arrow_pilon.72.update"  # FRA2    fos-related antigen-2	
				)
))
cds_subset <- cds[to_be_tested,]

diff_test_res <- differentialGeneTest(cds_subset, fullModelFormulaStr = "~sm.ns(Pseudotime)")

diff_test_res[,c("gene_short_name", "pval", "qval")]

pdf(paste0(opts$output,"11_markergene.pdf"),width=6,height=8)
plot_genes_in_pseudotime(cds_subset, color_by = "seurat_object.meta.data.seurat_clusters") + scale_color_manual(values=c("#1A863B", "#96CC88", "#F5BB6F"))
dev.off()

#----//DONE: Constructing Single Cell Trajectories :DONE\\------#