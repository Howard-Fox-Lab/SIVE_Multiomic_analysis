library(Seurat)
library(SCENIC)
library(AUCell)
library(ComplexHeatmap)
library(Seurat)

# https://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Running.html
# usage example
# a = export_image("AstroMicroCaudate/SalineSIV", "AstroMicro.0.1", "SalineSIV.png", true)
# a
# dev.off()
# https://stackoverflow.com/questions/7334644/sort-columns-of-a-dataframe-by-column-name
# https://stackoverflow.com/questions/20295787/how-can-i-use-the-row-names-attribute-to-order-the-rows-of-my-dataframe-in-r
# https://stackoverflow.com/questions/38466276/why-is-row-names-preferred-over-rownames 

export_heatmap = function(filepath, cluster_on, filename, nonduplicate) {
setwd(filepath)

load('pt_2.RData')
load("pt_1_.RData")
load('regulons.RData')
# this is how you sort the Idents in Seurat
Idents(obj) <- factor(x = Idents(obj), levels = sort(levels(obj)))
Idents(object = obj) <- cluster_on

cellInfo <- data.frame(seuratCluster=Idents(obj))
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
# rownames(regulonAUC) = sort(rownames(regulonAUC))
if (nonduplicate) {
	regulonAUC <- regulonAUC[rownames(regulonAUC),]
}
else {
	regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
}

regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seuratCluster),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))


print(regulonActivity_byCellType_Scaled)

# ordering the matrix here
regulonActivity_byCellType_Scaled = regulonActivity_byCellType_Scaled[,order(colnames(regulonActivity_byCellType_Scaled))]
regulonActivity_byCellType_Scaled = regulonActivity_byCellType_Scaled[order(rownames(regulonActivity_byCellType_Scaled)), ]

save(cellInfo, file = "cellInfo.saved.RData")

print(paste0(getwd(), '/',filename))
png(file = paste0(getwd(), '/',filename), width = 1000, height = 1000)
return(regulonActivity_byCellType_Scaled)
}

calc_rss = function(setName, cluster_on, filepath, filename, ret_Dir) {
setwd(filepath)

load('pt_2.RData')
load("pt_1_.RData")
load('regulons.RData')
# this is how you sort the Idents in Seurat
Idents(obj) <- factor(x = Idents(obj), levels = sort(levels(obj)))
Idents(object = obj) <- cluster_on
png(filename, width = 1000, height = 1000)
cellInfo <- data.frame(seuratCluster=Idents(obj))

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "seuratCluster"])
rssPlot <- plotRSS(rss)
setwd(ret_Dir)
return(plotRSS_oneSet(rss, setName))
}

export_bin_heatmap = function(minPerc, filepath, filename, cluster_on, ret_Dir) {
setwd(filepath)

load("pt_2.RData")
load("pt_1_.RData")
load("regulons.RData")
png(filename, width = 1000, height = 1000)

Idents(obj) <- factor(x = Idents(obj), levels = sort(levels(obj)))
Idents(object = obj) <- cluster_on

cellInfo <- data.frame(seuratCluster=Idents(obj))

binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$seuratCluster), 
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]
write.csv(binaryActPerc_subset, file = paste0("binaryActPerc_subset_aucMTX", ".csv"))  
setwd(ret_Dir)
return(ComplexHeatmap::Heatmap(binaryActPerc_subset, 
			       name="Regulon activity (%)", col = c("white","pink","red"), 
			       row_order = order(rownames(binaryActPerc_subset)), 
			       column_order=order(colnames(binaryActPerc_subset))))

}

get_regulon = function(filepath, ret_Dir) {
	setwd(filepath)
	load("pt_2.RData")
	load("pt_1_.RData")
	regulons <- loadInt(scenicOptions, "regulons")
	save(regulons, file="regulons.RData")
	setwd(ret_Dir)
}

export_image = function(filepath, cluster_on, filename, csv_exp, nonduplicate, ret_Dir) {
	png(filename, width = 1000, height = 1000)
	a = export_heatmap(filepath, cluster_on, filename, nonduplicate)
	write.csv(a, csv_exp)
	ComplexHeatmap::Heatmap(a, name="Regulon activity")
	dev.off()
	setwd(ret_Dir)
	# set to run the Caudatev3_SCENIC directoru
	return(ComplexHeatmap::Heatmap(a, name="Regulon activity", row_order=order(rownames(a)),column_order=order(colnames(a))))
}
