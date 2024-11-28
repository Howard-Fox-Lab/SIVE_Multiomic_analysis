library(SCENIC)

#https://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Running.html
createSeurat = function(filename_h5, filename_metadata, mm) {
obj = Read10X_h5(filename_h5)
newMetaData = read.delim(filename_metadata, sep=mm)
colnames(obj) = paste0(newMetaData$Cell.name, "_", newMetaData$Sample.name)
obj = CreateSeuratObject(obj)
obj = AddMetaData(obj, newMetaData)
save(obj, file=paste0(filename_h5, ".RData"))
return(obj)
}

export_heatmap = function(filepath, cluster_on, filename) {
setwd(filepath)

load('pt_2.RData')
load("pt_1_.RData")

Idents(object = obj) <- cluster_on

cellInfo <- data.frame(seuratCluster=Idents(obj))
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seuratCluster),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

save(cellInfo, file = "cellInfo.saved.RData")
write.csv(cellInfo, file = "cellInfo_graph.csv")
print(paste0(getwd(), '/',filename))
png(file = paste0(getwd(), '/',filename), width = 1000, height = 1000)
a = ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity")
return(regulonActivity_byCellType_Scaled)
}

export_image = function(filepath, cluster_on, filename) {
	png(filename, width = 1000, height = 1000)
	a = export_heatmap(filepath, cluster_on, filename)
	print(a)
	ComplexHeatmap::Heatmap(a, name="Regulon activity")
	invisible(dev.off())
	setwd("../../")
	# set to run the Caudatev3_SCENIC directoru
	return(ComplexHeatmap::Heatmap(a, name="Regulon activity"))
}


seuratGrnBoost = function(linkList_tsv = "int/1.4_GENIE3_linkList-9.tsv") {
	GRNBoost_LinkList <- importArboreto(linkList_tsv)
	colnames(GRNBoost_LinkList) = c("TF", "Target", "weight")
	if(length(colnames(GRNBoost_LinkList)) == 4) {
		GRNBoost_LinkList[[4]] = NULL
	}	
	saveRDS(GRNBoost_LinkList, file = "int/1.4_GENIE3_linkList.Rds")
}

seurat_scenicV2 =  function(obj_filename, cluster_on, session_name, cores, topWhat) {

## first part calculate the pt_1 specs
load(obj_filename)
exprMat = as.matrix(obj[['RNA']]$counts)
save(obj, exprMat, file = "pt_1_.RData")
library(SCENIC)
## solve this relation here
scenicOptions <- initializeScenic(org="mgi", dbDir="cisTarget_databases", nCores=cores)
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

### Co-expression network
genesKept <- geneFiltering(as.matrix(exprMat), scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]

runCorrelation(exprMat_filtered, scenicOptions)
exportsForArboreto(exprMat_filtered, scenicOptions, dir = "int")
exprMat_filtered_log <- log2(exprMat_filtered+1)
save(exprMat_filtered_log, scenicOptions, exprMat_filtered, file = "pt_2.RData")

# add more info here in the future .... 

}



seurat_scenic = function(obj_filename, cluster_on, session_name, cores, topWhat) {
# first, load the object you want

load(obj_filename)
exprMat = as.matrix(obj[['RNA']]$counts) #grab the counts

# Step 0 done
save(obj, exprMat, file = "pt_1_.RData")

library(SCENIC)
scenicOptions <- initializeScenic(org="hgnc", dbDir="cisTarget_databases", nCores=cores)

scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

### Co-expression network
genesKept <- geneFiltering(as.matrix(exprMat), scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)

save(scenicOptions, exprMat_filtered_log, file = 'pt_2.RData')
### Build and score the GRN
exprMat_log <- log2(exprMat+1)
#scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod=c(topWhat))# Toy run settings
runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

# first run
save.image(file=session_name)

Idents(object = obj) <- cluster_on
cellInfo <- data.frame(seuratCluster=Idents(obj))
regulons <- loadInt(scenicOptions, "aucell_regulons")
runSCENIC_4_aucell_binarize(scenicOptions)

save.image(file=session_name)
}
