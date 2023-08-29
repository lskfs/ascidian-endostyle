
library(future)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

rds.file <- 'zebrafish.final.rds'
if (!file.exists(rds.file)){
    obj <- readRDS('zebrafish.inte_1st.rds')

    early.lineage <- c("Hypoblast","Hypochord","Lateral plate mesoderm","Lens placode","Mesenchyme-related_arch or fin bud","Neural crest","Notochord","Paraxial mesoderm","Pectoral fin bud","Pharyngeal arch","Pharyngeal endoderm","Placode","Tailbud")
    early.timepoints <- c('1a_1dpf', '1b_1dpf', '2a_2dpf', '2b_2dpf')

    later.tissue <- c("Blood","Blood vessel","Central Nervous System","Connective tissue","Integument","Intestine","Kidney","Liver","Muscle","Myeloid lineage","Pancreas","Pineal gland","Retina","Spleen","Thymus", "Neuron", "Visceral mesoderm")
    later.timepoints <- c('5a_5dpf', '5b_5dpf')

    obj <- subset(obj, subset = (Tissue %in% early.lineage & Sample %in% early.timepoints) | 
                            (Tissue %in% later.tissue & Sample %in% later.timepoints) | 
                            (Tissue %in% c('Cochlea', 'Thyroid'))
              )
    obj@meta.data[obj@meta.data$Tissue %in% c('Central Nervous System', 'Neuron'), ]$Tissue <- 'Nervous system'

    saveRDS(obj, file = rds.file)
    write.table(obj@meta.data, file = 'zebrafish.final.metadata.txt', sep='\t', quote=F, row.names = T)

    #' save processed seurat object to loom file
    DefaultAssay(dataset.merge) <- 'RNA'
    obj.loom <- as.loom(dataset.merge, filename = paste0('./zebrafish.final.loom'), 
            overwrite = TRUE)
    obj.loom$close_all()
}
