########
#
#    integration.scdata.R
#    Reference: https://satijalab.org/seurat/articles/integration_introduction.html
#
#    1. Load two or more single-cell-rna-seq datasets(filtered_bc_matrix);
#    2. Identify cross-dataset pairs of cell which are of matched biological state("anchors"),can be used to correct technical differences between datasets(i.e.batch effect) and perform comparative scRNA-seq analysis of acrooss experimental conditions;
#    3. Perform integration and create an 'integrated' data asssy for downstream analysis;
#    4. Identify conserved cell type markers and differential expressed genes across conditions;
#
########
#
#    Please store file path of all folder within a file named read10xpath
#
########

library(Seurat)
#library(SeuratData)
library(patchwork)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

## Get the parameters
parser = argparse::ArgumentParser(description = 'Script for combine Single-cell-rna-seq subject and reference data by Seurat')
parser$add_argument('-d','--dir.path', dest = 'dir.path', help='input single-cell read10xpath directory')
parser$add_argument('-b','--batch', dest = 'batch', help='input batch set name, name your workbatch')
parser$add_argument('-o','--output', dest = 'output', help='output filepath')
parser$add_argument('-r','--resolution', dest = 'resolution',default = 0.8, help='cluster resolution default 0.8')
parser$add_argument('-z','--pointsizeumap', dest = 'pointsizeumap', default=0.2, help='point size for plotting')
parser$add_argument('-D', '--dim', dest = 'dim', default = 10, type = 'integer', help = 'dim in principle component, default 10')

opts = parser$parse_args()
opts$output
opts$resolution

if(FALSE){
opts<-list()
opts$dir.path("/mnt/inspurfs/home/dongb/jiangan/Singlecell/Integrationpipeline/211125-Sc-mito20.20.30-SCT-all/path")
opts$dir.path<-c("/mnt/inspurfs/home/dongb/jiangan/Singlecell/Integrationpipeline/211125-Sc-mito20.20.30-SCT-all/path")
opts$batch <- "Sc-mito20.20.30-SCT-all"
opts$output <- "/mnt/inspurfs/home/dongb/jiangan/Singlecell/Integrationpipeline/211125-Sc-mito20.20.30-SCT-all/"
opts$resolution <- 0.4
opts$pointsizeumap <- 2
opts$dim<-10
}

opts$resolution <- as.numeric(opts$resolution)
opts$pointsizeumap <- as.numeric(opts$pointsizeumap)
opts$dim<- as.numeric(opts$dim)
dir.create(opts$output)
setwd(opts$output)

#-------------------------------#
# Load file
dir <- read.table(opts$dir.path,header=F)[,1]
names(dir) <- 1:length(dir)
counts <- list()
scRNA <- list()
object <- list()
for(i in 1:length(dir)){
	dir[i]
	scRNA[[i]] <- readRDS(dir[i])
	scRNA[[i]]$stim <- i
	object[[i]] <- scRNA[[i]]
	object[[i]] <- SCTransform(object[[i]], n_genes=NULL, min_cells=5, method='qpoisson')
}
object

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = object,nfeatures=3000)
object <- PrepSCTIntegration(object.list=object, anchor.features=features)

# Perform integration
object.anchors <- FindIntegrationAnchors(object.list = object, normalization.method = "SCT",  anchor.features = features)
# Creates an 'integrated' data assay
object.combine <- IntegrateData(anchorset = object.anchors, normalization.method = "SCT")


#-------------------------------#
# Perform an integrated analysis
# we can run a single integrated analysis on all cells

object.combine <- RunPCA(object.combine, verbose = FALSE)

object.combine <- FindNeighbors(object = object.combine)
object.combine <- FindClusters(object = object.combine, resolution = opts$resolution)

table(object.combine@active.ident)
#head(subset(as.data.frame(object.combine@active.ident),object.combine@active.ident=="2"))
#head(WhichCells(object.combine,idents="2"))

object.combine <- BuildClusterTree(object.combine)
Tool(object = object.combine, slot = 'BuildClusterTree')
pdf(paste0(opts$batch,'-ClusterTree.pdf'),width = 8,height = 8)
PlotClusterTree(object.combine)
dev.off()

object.combine <- RunUMAP(object = object.combine, reduction = "pca", dims = 1:opts$dim)
object.combine$stim
opts$output

#Visualization
pdf(paste0(opts$output,opts$batch,"-UMAP",'.pdf'),width = 8,height = 8)
DimPlot(object.combine, reduction = "umap", group.by = "stim",pt.size= opts$pointsizeumap)
#DimPlot(object.combine, reduction = "umap", label=TRUE, repel = TRUE)
dev.off()

pdf(paste0(opts$output,opts$batch,"-UMAP-batch",'.pdf'),width = 16,height = 8)
DimPlot(object.combine, reduction = "umap", split.by = "stim")
dev.off()

saveRDS(object.combine, file = paste0(opts$batch,"-integration-unplot.rds"))
#######Resolution 'for' output loop##########

#inter <- object.combine

for (i in seq(0.1,0.5,0.1)){
#for (i in 0.3){
#  object.combine <- inter
  object.combine <- RunPCA(object.combine, verbose = FALSE)
  object.combine <- FindNeighbors(object = object.combine)
 #object.combine <- FindClusters(object = object.combine, resolution = opts$resolution)
  object.combine <- FindClusters(object.combine, resolution = i)
  object.combine <- RunUMAP(object = object.combine, reduction = "pca", dims = 1:opts$dim)

  filename=paste0(opts$batch,"_Umap_label_resolution_",i,".pdf")
  pdf(filename)
  print(DimPlot(object.combine, reduction = "umap", label=TRUE))
  dev.off()
if(FALSE){
  obj.markers <- FindAllMarkers(object.combine, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.3)
  obj.markers <- obj.markers[with(obj.markers, order(cluster, -avg_log2FC)), ]
  write.table(obj.markers, paste0(opts$batch,'_resolution_',i,'_AllMarkers.xls'), sep ='\t', quote = FALSE, row.names = F)
  topn <- obj.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
  write.table(topn, paste0(opts$batch,'_resolution_',i,'_AllMarkers.top30.xls'), sep='\t', quote = FALSE, row.names = F)

  top10 <- obj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  filename<-paste0(opts$batch,"_resolution_",i,"top10_DOHeatmap.pdf")
  pdf(filename,width=15,height=20)
  plot<-DoHeatmap(object.combine, features = top10$gene,slot = "scale.data", size = 0.5) + NoLegend() 
  print(plot)
  dev.off() 
  print("finished")
}

  plot.gene <- rownames(object.combine)
  loop<-c(1:length(plot.gene))
  dir.create("AllFeaturePlot_0.3")
  setwd("AllFeaturePlot_0.3")
  DefaultAssay(object.combine) <- "integrated"
  for(i in loop){
        plotname <- paste0(plot.gene[i],"-FeaturePlot.jpeg")
        print(plotname)
        jpeg(plotname)
        p <- FeaturePlot(object.combine, features = plot.gene[i])
        print(p)
        dev.off()   }
  setwd(outputpath)
}

#sample.markers <- FindAllMarkers(sample, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#write.table(sample.markers,paste0(opts$batch, '-AllMarkers.csv'))

#top_x <- sample.markers %>% group_by(cluster) %>% top_n(n = 10)
#plotname <- paste0(opts$batch, "-Heatmaptop-", "gene.pdf")
#pdf(plotname,width=100,height=100)
#DoHeatmap(sample, features = top_x$gene) + NoLegend()
#dev.off()

saveRDS(object.combine, file = paste0(opts$batch,"-integration-SCT.rds"))
