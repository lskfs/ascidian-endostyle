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
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(SoupX)
library(DoubletFinder)
library(DropletUtils)
library(SeuratDisk)

parser = argparse::ArgumentParser(description = 'Single-cell-seurat script')
parser$add_argument('-i', '--inputpath', dest = 'inputpath', help = 'input data directory, filtered_feature_bc_matrix')
parser$add_argument('-R', '--raw_inputpath',dest = 'raw_inputpath', help = "input data directory, raw_feature_bc_matrix, used for SoupX")
parser$add_argument('-o', '--outputpath', dest = 'outputpath', help = 'output data directory')
parser$add_argument('-P','--projectname', dest = 'projectname', help = 'projectname')
parser$add_argument('-r','--resolution', dest = 'resolution',  type = 'double', default = 0.8, help = 'resolution for clustering')
parser$add_argument('-c','--nCount_RNAmin', dest = 'nCount_RNAmin', type='integer', default = 0, help = 'minimum Count number')
parser$add_argument('-C','--nCount_RNAmax', dest = 'nCount_RNAmax', type='integer', help = 'maximum Count number')
parser$add_argument('-f','--minfeatures', dest = 'minfeatures', type='integer', default = 0, help = 'minimum feature number')
parser$add_argument('-F','--maxfeatures', dest = 'maxfeatures', type='integer', help = 'maximum feature number')
parser$add_argument('-d', '--dim', dest = 'dim', default = 10, type = 'integer', help = 'dim in principle component, default 10')
parser$add_argument('-t','--top', dest = 'top', type = 'integer', default = 10, help = 'top gene for Heatmap')
parser$add_argument('-m','--percentmt', dest = 'percentmt', type = 'integer', default = 5, help = 'percent of mitochondria filtering')
parser$add_argument('-M','--mtgene',dest = 'mtgene', type = 'character', help = 'mitochondria gene list')
parser$add_argument('-D','--double', dest = 'double',  default = 0.075, help = 'DoubletFinder - Assuming 7.5% doublet formation rate - tailor for your dataset')
opts = parser$parse_args()

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
#print(paste("nFeature:",parser$nfeature,"\n"))
 nfeature<-opts$nfeature
print(paste("nCount_RNAmin:",opts$nCount_RNAmin,"\n"))
nCount_RNAmin<-opts$nCount_RNAmin
print(paste("nCount_RNAmax:",opts$nCount_RNAmax,"\n"))
nCount_RNAmax<-opts$nCount_RNAmax
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
print(paste("double number of Doublefinder:",opts$double,"\n"))
double<-opts$double
#print(paste("cellnumber:",parser$cellnumber,"\n")); #cellnumber<-opts$cellnumber
#***********************************************#

dir.create(outputpath)
setwd(outputpath)
getwd()
sample.data <- Read10X(data.dir = inputpath)
sample <- CreateSeuratObject(counts = sample.data, project = projectname)
sample.raw <- Read10X(data.dir = raw_inputpath)
sample.raw <- sample.raw[rownames(sample.data),]
sample

VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
ggsave(paste0("0-", projectname,"-presoupx-QCVlnPlot.png"))

#h5path = paste(outputpath,"/",projectname,"_seuratobject_rawmatrix.h5Seurat",sep="")
#SaveH5Seurat(sample, file = h5path, overwrite=TRUE)

###########SOUPX QC#############
sample <- NormalizeData(sample, normalization.methods = "LogNormalize",scale.factor = 10000)
sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 3000)
sample.genes <- rownames(sample)
sample <- ScaleData(sample, features = sample.genes)

sample <- RunPCA(sample, features = VariableFeatures(sample), npcs = 40, verbose = F)
sample <- FindNeighbors(sample, dims = 1:30)
sample <- FindClusters(sample, resolution = 0.5)
sample <- RunUMAP(sample, dims = 1:30)

matx <- sample@meta.data
sc = SoupChannel(sample.raw, sample.data)
sc = setClusters(sc, setNames(matx$seurat_clusters,rownames(matx)))
sc = autoEstCont(sc)
sample = adjustCounts(sc)

DropletUtils:::write10xCounts(paste0("./soupX_matrix"), sample ,version="3")
############SOUPX DONE###########

sample.data <- Read10X(data.dir =paste0("./soupX_matrix"))
sample <- CreateSeuratObject(counts = sample.data, project = projectname)   

VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(paste0("1-", projectname,"-postsoupx-QCVlnPlot.png"))

######Mitochondria ratio calculation
mt.gene <- read.table(mtgene, header=FALSE)
mt.gene <- mt.gene[,1]
C<-GetAssayData(object = sample, slot = "counts")
mt <-intersect(mt.gene, row.names(C))
percent.mito <- colSums(C[mt,])/colSums(C)*100
sample <- AddMetaData(sample, percent.mito, col.name = "percent.mt")
sample 

#VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#ggsave(paste0("2-", projectname,"-postselfcountmt-QCVlnPlot.png"))

##### Filtering
if (is.null(nCount_RNAmax)){nCount_RNAmax <- max(sample$nCount_RNA)}
if (is.null(maxfeatures)){maxfeatures <- max(sample$nFeature_RNA)}
sample <- subset(sample, subset = (nFeature_RNA > minfeatures) & (nFeature_RNA < maxfeatures) & (nCount_RNA > nCount_RNAmin) & (nCount_RNA < nCount_RNAmax) & (sample$percent.mt < percentmt))

VlnPlot(sample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(paste0("3-", projectname,"-Filtered-QCVlnPlot.png"))

print("*****After filtering:*********")
sample

#########
Find_doublet <- function(obj, pc = dim){
    sweep.res.list <- paramSweep_v3(obj, PCs = 1:pc, sct = TRUE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)

    nExp_poi <- round(as.numeric(double) * ncol(obj))
    p <- as.numeric(as.vector(bcmvn[bcmvn$MeanBC == max(bcmvn$MeanBC), ]$pK))
    obj <- doubletFinder_v3(obj, PCs = 1:pc, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
    colnames(obj@meta.data)[ncol(obj@meta.data)] = "doublet_info"
    return(obj)
}
#########

sample <- Find_doublet(sample, pc = dim)
sample <- subset(sample, subset =doublet_info=='Singlet')
###############Filter done#################################

plot1 <- FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sample, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
ggsave(paste0("4-", projectname,"-FeatureScatter-feature2feature.png"))

sample <- NormalizeData(sample, verbose = FALSE, assay = 'RNA')
sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 3000)
sample.genes <- rownames(sample)
sample <- ScaleData(sample, features = sample.genes)
#sample <- SCTransform(sample, assay = 'RNA', variable.features.n = opts$vg, return.only.var.genes = FALSE, n_genes=NULL, min_cells=5, method='qpoisson')
#DefaultAssay(sample) <- 'SCT'
VariableFeaturePlot(sample)
ggsave(paste0("6-", projectname, "-VFplot.png"))

sample <- RunPCA(sample)
sample <- FindNeighbors(sample, dims = 1:dim)
sample <- FindClusters(sample, resolution = resolution)
ElbowPlot(sample)
ggsave(paste0("7-", projectname,"-ElbowPlot.png"))

if(FALSE){
#' create heatmap and cluster color palette for figure plot function
cluster_number <- length(unique(obj@meta.data$seurat_clusters))
if (cluster_number <= 25){
    cluster_Palette <- c('dodgerblue2', '#E31A1C', 'green4', '#6A3D9A', '#FF7F00', 'black', 'gold1', 'skyblue2', '#FB9A99', 'palegreen2', '#CAB2D6', '#FDBF6F', 'gray70', 'khaki2',                 'maroon', 'orchid1', 'deeppink1', 'blue1', 'steelblue4','darkturquoise','green1', 'yellow4', 'yellow3','darkorange4', 'brown')
} else if (cluster_number > 25 && cluster_number <= 70){
    qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
    cluster_Palette <- unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
}
DimPlot(obj, reduction = 'umap', cols = cluster_Palette)
}

sample <- RunUMAP(sample, dims = 1:dim)
DimPlot(sample, reduction = 'umap')
ggsave(paste0("8-",projectname,"-UMAP.png"))

sample <- RunTSNE(sample, features = VariableFeatures(object = sample))
DimPlot(sample, reduction = 'tsne')
ggsave(paste0("9-",projectname,"-tSNE.png"),width=20,height=16,dpi =300,units = 'in')

#sample <- JackStraw(sample, num.replicate = 100)
#sample <- ScoreJackStraw(sample, dims = 1:20)
#JackStrawPlot(sample, dims = 1:dim)
#ggsave(paste0("10-", projectname,"-JackStrawPlot.pdf"))

if(FALSE){
# Finding differentially expressed features (cluster biomarkers)
# find all markers of cluster 2
cluster2.markers <- FindMarkers(sample, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster2.markers <- FindMarkers(sample, ident.1 = 2, ident.2 = c(5,8), min.pct = 0.25)
print("**************find all markers distinguishing cluster 2 from clusters 5 and 8********")
plot(head(cluster2.markers, n = 10))
ggsave(paste0("10.1-",projectname,"cluster2to5&8.png"))

cluster5.markers <- FindMarkers(sample, ident.1 = 5, ident.2 = c(2,8), min.pct = 0.25)
print("**************find all markers distinguishing cluster 5 from clusters 2 and 8********")
plot(head(cluster5.markers, n = 10))
ggsave(paste0("10.1-",projectname,"cluster5to2&8.png"))


cluster8.markers <- FindMarkers(sample, ident.1 = 8, ident.2 = c(2,5), min.pct = 0.25)
print("**************find all markers distinguishing cluster 8 from clusters 2 and 5********")
plot(head(cluster8.markers, n = 10))
ggsave(paste0("10.1-",projectname,"cluster8to2&5.png"))
}

# find markers for every cluster compared to all remaining cells, report only the positive ones
sample.markers <- FindAllMarkers(sample, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
sample.markers %>% group_by(cluster) %>% top_n(n = 2)
write.table(sample.markers,paste(projectname,'allmarkers.csv',sep=""))

cluster0.markers <- FindMarkers(sample, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

#plotname10 <- paste("10", projectname,"VlnPlot-Gene1.pdf",sep="")
#pdf(plotname10)
#VlnPlot(sample, features = c("MSTRG.10271"))
#dev.off()

# you can plot raw counts as well
#plotname11 <- paste("11",projectname,"VlnPlot-Gene2.pdf",sep="")
#pdf(plotname11)
#VlnPlot(sample, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
#dev.off()

if(FALSE){
plot.gene <- rownames(sample)
loop<-c(1:length(plot.gene))
dir.create("12-FeaturePlot")
setwd("12-FeaturePlot")
for(i in loop){
	plotname12 <- paste0(plot.gene[i],"-", projectname,"-FeaturePlot.jpeg")
	print(plotname12)
	jpeg(plotname12)
	p <- FeaturePlot(sample, features = plot.gene[i])
	print(p)
	dev.off()
}
setwd(outputpath)
}
top_x <- sample.markers %>% group_by(cluster) %>% top_n(n = top)
DoHeatmap(sample, features = top_x$gene) + NoLegend()
ggsave(paste0("10-", projectname,"Heatmaptop",top,"gene.png"))

#new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", 
#                     "NK", "DC", "Platelet")
#names(new.cluster.ids) <- levels(sample)
#DimPlot(sample, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
#dev.off()

#dir.create(paste0(outputpath,'/final_matrix'))
if(FALSE){
write.table(data.frame(rownames(sample),rownames(sample)),file = paste0(outputpath,'/final_matrix/genes.tsv'),quote = F,sep = '\t',col.names = F,row.names = F)
write.table(colnames(sample),file = paste0(outputpath,'/final_matrix/barcodes.tsv'),quote = F,col.names = F,row.names = F)

file="matrix.mtx"
sink(file)
cat("%%MatrixMarket matrix coordinate integer general\n")
cat("%\n")
cat(paste(nrow(sample),ncol(sample),sum(sample>0),"\n")) 
sink()
tmp=sample[1:5,1:4]
tmp
tmp=do.call(rbind,lapply(1:ncol(sample),function(i){
  return(data.frame(row=1:nrow(sample),
                    col=i,
                    exp=sample[,i]))
}) )
tmp=tmp[tmp$exp>0,]
head(tmp)
write.table(tmp, file = paste0(outputpath,'/final_matrix/matrix.mtx'), quote = F, col.names = F,row.names = F,append = T )
}

#DropletUtils:::write10xCounts(paste0("./final_matrix"), sample ,version="3")

RDSpath = paste(outputpath,"/",projectname,"_seuratobject_final.rds",sep="")
saveRDS(sample, file = RDSpath)
