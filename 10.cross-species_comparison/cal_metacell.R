library(Seurat)
library(SeuratDisk)
library(dplyr)
library(tidyr)

### function defination for meta cell aggretation
pseudocell <- function(obj, group="orig.ident", group_n=10){
    obj@meta.data <- obj@meta.data %>% mutate(cellN = colnames(obj)) %>% group_by(get(group)) %>%
        mutate(bins.no=rep(1:ceiling(n()/group_n), each=group_n, length.out=n())) %>%
        mutate(pseudo_label=paste0(get(group), "_N.", bins.no)) %>% tibble::column_to_rownames('cellN') %>%
        data.frame()
    obj.mini <- AverageExpression(object=obj, group.by="pseudo_label", slot='data', return.seurat=TRUE)
    return(object.mini)
}

### read in zebrafish rds file
zf.orig <- readRDS('./zebrafish.final.rds')
DefaultAssay(zf.orig)<-'RNA'

### extract data for zebrafish thyroid cells
tissue <- c('Pharyngeal endoderm', 'Thymus', 'Intestine', 'Liver', 'Pancreas')
clusters <- c(7, 10, 16, 18, 22, 23, 25, 28, 29, 31, 32, 34, 35, 39, 42, 45)
zf <- subset(zf.orig, subset = Tissue %in% tissue | (Tissue %in% 'Thyroid' & seurat_clusters %in% clusters))

zf <- pseudocell(zf, group='Tissue', group_n=10)
zf@meta.data['Tissue'] <- as.character(zf@meta.data$orig.ident)

SaveH5Seurat(zf, "zf.thyroid.h5Seurat")
Convert('zf.thyroid.h5Seurat', dest='h5ad')

### extract data for zebrafish hair cells
celltypes <- c('neuromast hair cell 1', 'macula hair cell', 'hair cell progenitor', 'neuromast hair cell 2', 'crista hair cell')
zf <- subset(zf.orig, subset = CellType %in% celltypes)

zf <- pseudocell(zf, group_n=10)
zf@meta.data['CellType'] <- as.character(zf@meta.data$orig.ident)

SaveH5Seurat(zf, "zf.cochlea.h5Seurat")
Convert('zf.cochlea.h5Seurat', dest='h5ad')


### read in endostyle rds file
sc.orig <- readRDS('../data/Styela_clava.anno.rds')
DefaultAssay(sc.orig) <- 'RNA'

### extract data for endostyle thyroid
sc <- subset(sc.orig, subset = inte_anno == 'Thyroid like cells')

sc <- pseudocell(sc, group='inte_anno', group_n=10)
sc@meta.data['inte_anno'] <- as.character(sc@meta.data$orig.ident)

SaveH5Seurat(sc, "sc.TLC.h5Seurat")
Convert('sc.TLC.h5Seurat', dest='h5ad')

### extract data for endostyle cochlea
sc <- subset(sc.orig, subset = inte_anno == 'Hair cells')

sc <- pseudocell(sc, group='inte_anno', group_n=10)
sc@meta.data['inte_anno'] <- as.character(sc@meta.data$orig.ident)

SaveH5Seurat(sc, "sc.HCLC.h5Seurat")
Convert('sc.HCLC.h5Seurat', dest='h5ad')
