library(SCENIC)
library(Seurat)
library(SeuratObject)


AddSCENICAssay = function(obj){ 
ls = readRDS("int/3.4_regulonAUC.Rds")
ls = S4ToList(ls)
regulonAUC = ls$assays$data$listData$AUC
rownames(regulonAUC) = unlist(lapply(rownames(regulonAUC), function(X) strsplit(X, " ")[[1]][1] ))
regulonAUC = as.data.frame(regulonAUC)
regulonAUC$Names = rownames(regulonAUC)
regulonAUC = regulonAUC[!grepl("extended", regulonAUC$Names), ]
regulonAUC$Names = NULL
obj[['SCENIC']] = CreateAssayObject(regulonAUC)
return(obj)
save(obj, file = 'Seurat_SCENIC_assayObject.RData')
}


CalculateUMAP = function(obj, dims, res) {
DefaultAssay(obj) = "SCENIC"
obj <- NormalizeData(object = obj)
obj <- FindVariableFeatures(object = obj)
obj <- ScaleData(object = obj)
obj <- RunPCA(object = obj)
obj <- FindNeighbors(object = obj, dims = 1:dims)
obj <- FindClusters(object = obj, resolution=res)
obj <- RunUMAP(object = obj, dims = 1:dims)
save(obj, file = paste0("Seurat_SCENIC_assayObject", "res: ", res, "dims: ", dims, ".RData"))
return(obj)
}

CalculateOverlayClusters  = function(obj, SCENIC_clusters, RNA_clusters) {
library(dplyr)
s_t = obj@meta.data
by_regulonAUC_res0.1 = s_t %>% group_by_at(c(SCENIC_clusters, RNA_clusters))
return(by_regulonAUC_res0.1 %>% tally())
}


