#Single-cell-seurat script
#By An Jiang 20211101
#https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
# Version compared to 20210726: Add mitochondial gene filer; Update DimHeatmap for dim choosing
# Version compared to 20211010: Refined structure of pipeline
# SoupX added, please cite: "Young, M.D., Behjati, S. (2020). SoupX removes ambient RNA contamination from droplet-based single-cell RNA sequencing data, GigaScience, Volume 9, Issue 12, December 2020, giaa151bioRxiv, 303727, https://doi.org/10.1093/gigascience/giaa151"
# DoubletFinder added, DoubletFinder is an R package that predicts doublets in single-cell RNA sequencing data.


library(dplyr)
library(Seurat)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(SoupX)
library(argparse)
library(ggplot2)
#library(DoubletFinder)

parser <- ArgumentParser(description = 'Single-cell-seurat script')
parser$add_argument('-i', '--inputpath', dest = 'inputpath', help = 'input data directory, filtered_feature_bc_matrix')
parser$add_argument('-R', '--raw_inputpath',dest = 'raw_inputpath', help = "input data directory, raw_feature_bc_matrix, used for SoupX")
parser$add_argument('-o', '--outputpath', dest = 'outputpath', help = 'output data directory')
parser$add_argument('-P','--projectname', dest = 'projectname', help = 'projectname')
parser$add_argument('-r','--resolution', dest = 'resolution', default = 0.8, type = 'numeric', help = 'resolution for clustering')
parser$add_argument('-c','--nCount_RNAmin', dest = 'nCount_RNAmin', type = 'numeric', default = 0, help = 'minimum Count number')
parser$add_argument('-C','--nCount_RNAmax', dest = 'nCount_RNAmax', type = 'numeric', help = 'maximum Count number')
parser$add_argument('-f','--minfeatures', dest = 'minfeatures', type = 'numeric', default = 0, help = 'minimum feature number')
parser$add_argument('-F','--maxfeatures', dest = 'maxfeatures', type = 'numeric', help = 'maximum feature number')
parser$add_argument('-d', '--dim', dest = 'dim', default = 10, type = 'integer', help = 'dim in principle component, default 10')
parser$add_argument('-t','--top', dest = 'top', type = 'integer', default = 10, help = 'top gene for Heatmap')
parser$add_argument('-m','--percentmt', dest = 'percentmt', type = 'integer', default = 5, help = 'percent of mitochondria filtering')
parser$add_argument('-M','--mtgene',dest = 'mtgene', type = 'character', help = 'mitochondria gene list')
#parser$add_argument('-D','--double', dest = 'double', type = 'numeric', default = 0.075, help = 'DoubletFinder - Assuming 7.5% doublet formation rate - tailor for your dataset')
parser$print_help()
opts = parser$parse_args()
opts

#-----------------------------#
print(paste("Workpath:",opts$inputpath,"\n"))
inputpath<-opts$inputpath
print(paste("Workpath:",opts$raw_inputpath,"\n"))
raw_inputpath<-opts$raw_inputpath
print(paste("Outputpath:",opts$outputpath,"\n"))
outputpath<-opts$outputpath
print(paste("Projectname:",opts$projectname,"\n"))
projectname<-opts$projectname
print(paste("Resolution:",opts$resolution,"\n"))
resolution<-opts$resolution
#print(paste("nFeature:",opts$nfeature,"\n"))
nfeature<-opts$nfeature
print(paste("nCount_RNAmin:",opts$nCount_RNAmin,"\n"))
minncount<-opts$nCount_RNAmin
print(paste("nCount_RNAmax:",opts$nCount_RNAmax,"\n"))
maxncount<-opts$nCount_RNAmax
print(paste("percentmt:",opts$percentmt,"\n"))
percentmt<-opts$percentmt
print(paste("minfeatures:",opts$minfeatures,"\n"))
minfeatures<-opts$minfeatures
print(paste("maxfeatures:",opts$maxfeatures,"\n"))
maxfeatures<-opts$maxfeatures
print(paste("dim:",opts$dim,"\n"))
dim<-opts$dim
print(paste("top:",opts$top,"\n"))
top<-opts$top
print(paste("mito.gene:",opts$mtgene,"\n"))
mtgene<-opts$mtgene
#print(paste("double number of Doublefinder:",opts$double,"\n")); double<-opts$double
#print(paste("cellnumber:",opts$cellnumber,"\n")); #cellnumber<-opts$cellnumber
#***********************************************#

dir.create(outputpath)
setwd(outputpath)
sample.data <- Read10X(data.dir = inputpath)
sample <- CreateSeuratObject(counts = sample.data, project = projectname)
sample.raw <- Read10X(data.dir = raw_inputpath)
#sample.raw <- sample.raw[rownames(sample),]
sample

mt.gene <- read.table(mtgene, header=FALSE)
mt.gene <- mt.gene[,1]
C<-GetAssayData(object = sample, slot = "counts")
mt <-intersect(mt.gene, row.names(C))
percent.mito <- colSums(C[mt,])/colSums(C)*100
sample <- AddMetaData(sample, percent.mito, col.name = "percent.mt")

VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(paste0("0-", projectname,"-QCVlnPlot.png"))

## Filtering
#if (is.null(nCount_RNAmax)){nCount_RNAmax <- max(sample$nCount_RNA)}
#if (is.null(maxfeatures)){maxfeatures <- max(sample$nFeature_RNA)}
sample <- subset(sample, subset = (nFeature_RNA > minfeatures) & (nFeature_RNA < maxfeatures) & (nCount_RNA > minncount) & (nCount_RNA < maxncount) & (sample$percent.mt < percentmt))

VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(paste0("1-", projectname,"-Filtered-QCVlnPlot.png"))

print("*****After filtering:*********")
sample

#########
if(FALSE){Find_doublet <- function(obj, pc = dim){
    sweep.res.list <- parserSweep_v3(obj, PCs = 1:pc, sct = TRUE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)

    nExp_poi <- round(as.numeric(double) * ncol(obj))
    p <- as.numeric(as.vector(bcmvn[bcmvn$MeanBC == max(bcmvn$MeanBC), ]$pK))
    obj <- doubletFinder_v3(obj, PCs = 1:pc, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
    colnames(obj@meta.data)[ncol(obj@meta.data)] = "doublet_info"
    return(obj)
}}
#########

sample <- SCTransform(sample, assay = 'RNA', variable.features.n = 2000, return.only.var.genes = FALSE, n_genes=NULL, min_cells=5, method='qpoisson')
sample <- RunPCA(sample)
#sample <- Find_doublet(sample, pc = dim)
#sample <- subset(sample, subset =doublet_info=='Singlet')
#matx <- sample@meta.data
#sc = SoupChannel(sample.raw, sample.data)
#sc = setClusters(sc, setNames(matx$seurat_clusters, rownames(matx)))
#sc = autoEstCont(sc)

#sample = adjustCounts(sc)

VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(paste0("2-", projectname,"-postdoublet-QCVlnPlot.png"))

plot1 <- FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
ggsave(paste0("3-", projectname,"-FeatureScatter-feature2feature.png"))

sample <- NormalizeData(sample, verbose = FALSE, assay = 'RNA')
sample <- SCTransform(sample, assay = 'RNA', variable.features.n = opts$vg, return.only.var.genes = FALSE, n_genes=NULL, min_cells=5, method='qpoisson')

DefaultAssay(sample) <- 'SCT'
sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(sample)
ggsave(paste0("5-", projectname, "-VFplot.png"))

sample <- RunPCA(sample)
sample <- FindNeighbors(sample, dims = 1:dim)
sample <- FindClusters(sample, resolution = resolution)
ElbowPlot(sample)
ggsave(paste0("6-", projectname,"-ElbowPlot.png"))

sample <- RunUMAP(sample, dims = 1:dim)
png(paste0("7-",projectname,"-UMAP.png"))
DimPlot(sample,reduction="umap")
dev.off()

# find markers for every cluster compared to all remaining cells, report only the positive ones
sample.markers <- FindAllMarkers(sample, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
sample.markers %>% group_by(cluster) %>% top_n(n = 2)
write.table(sample.markers,paste(projectname,'-allmarkers.xls',sep=""),sep="\t")

RDSpath = paste(outputpath,projectname,"_seurat.rds",sep="")
saveRDS(sample, RDSpath)