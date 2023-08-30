library(Seurat)

obj <- readRDS("BloodAndHLR.integrated.rds")
tail(Idents(obj))

objsub <- subset(obj,idents=c(2,9,11),invert = FALSE)
saveRDS(objsub, file="BloodAndHLR.integrated.2-9-11.rds")
write.table(objsub@meta.data, file="BloodAndHLR.integrated.2-9-11.metadata.txt")

