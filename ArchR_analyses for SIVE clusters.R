library(ArchR)
library(tidyr)
library(dplyr)
library(ggplot2)
library(GenomicFeatures)
set.seed(1234)
addArchRThreads(threads = 16)
setwd("/work/foxlab/xiaoke/scATAC/atac_enceph_new/")

proj_4 <- readRDS("/work/foxlab/xiaoke/scATAC/atac_enceph_new/Myeloid_cells/Save-ArchR-Project.rds")
proj_4 <- addImputeWeights(proj_4, reducedDims = 'Harmony_combined')

###Compare the expression of MAMU and homeostatic genes
markerGenes  <- c(
  "P2RY12",'CX3CR1','GPR34','SALL1', #Microglia
  "MAMU-DRA",'MAMU-DRB1','MAMU-DRB5','CD74' # CAM
)

p1 <- plotGroups(
  ArchRProj =  proj_4,
  groupBy = "Label_2",
  colorBy = "GeneScoreMatrix",
  name = markerGenes,
  plotAs = "violin",
  pal = c("darkred","darkblue","darkgreen","darkviolet")
)

p2 <- lapply(p1, function(x){
  x + guides(color = "none", fill = "none") + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.text.y=element_text(size = 20), 
      axis.title = element_blank()
    )+
    labs(title = NULL)
})

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/ATAC.png", width = 4000, height = 2000, units = 'px',
    res = 300)
do.call(cowplot::plot_grid, c(list(ncol = 4),p2))
dev.off()


#### Analyses for CAM-like cluster and Microglia-like cluster
### subset the cell clusters
id <- BiocGenerics::which(proj_4$Label_2 %in% c("CAM_like","Micro_like"))
cells <- proj_4$cellNames[id]
proj_5 <- subsetArchRProject(proj_4, cells = cells, outputDirectory = "SIVE_cluster", force = T,dropCells =FALSE)
proj_5 <- readRDS("/work/foxlab/xiaoke/scATAC/atac_enceph_new/SIVE_cluster/Save-ArchR-Project.rds")
proj_5 <- addImputeWeights(proj_5, reducedDims = 'Harmony_combined')

# Plot UMAP
p1 <- plotEmbedding(ArchRProj = proj_5, colorBy = "cellColData", name = "Label_2", embedding = "UMAP_combined", baseSize = 5)
png(filename = "/work/foxlab/xiaoke/scATAC/Figures/UMAP.png", width = 4000, height = 4000, units = 'px',
    res = 600)
p1
dev.off()

### Find markers for atac data
markerGS <- getMarkerFeatures(
  ArchRProj = proj_5,
  useMatrix = 'GeneScoreMatrix',
  groupBy = 'Label_2',
  bias = c('TSSEnrichment','log10(nFrags)'),
  testMethod = 'wilcoxon'
)
markerList <- getMarkers(markerGS, cutOff = 'FDR <= 0.01 & Log2FC >= 1.25')
CAM_like_atac <- markerList$CAM_like[c(5,7)]
Micro_like_atac <- markerList$Micro_like[c(5,7)]
write.csv(markerList, file = '/work/foxlab/xiaoke/scATAC/excel/markers_atac.csv')

markerRNA <- getMarkerFeatures(
  ArchRProj = proj_5,
  useMatrix = 'GeneExpressionMatrix',
  groupBy = 'Label_2',
  bias = c('TSSEnrichment','log10(nFrags)'),
  testMethod = 'wilcoxon'
)
markerList <- getMarkers(markerRNA, cutOff = 'FDR <= 0.01 & Log2FC >= 1.25')
CAM_like_RNA <- markerList$CAM_like[c(5,7)]
Micro_like_RNA <- markerList$Micro_like[c(5,7)]
write.csv(markerList, file = '/work/foxlab/xiaoke/scATAC/excel/markers_RNA.csv')

# Find the same markers for RNA and ATAC
CAM_like <- intersect(CAM_like_RNA$name,CAM_like_atac$name)
Micro_like <- intersect(Micro_like_RNA$name,Micro_like_atac$name)
write.csv(CAM_like, file = '/work/foxlab/xiaoke/scATAC/excel/CAMlike_both.csv')
write.csv(Micro_like, file = '/work/foxlab/xiaoke/scATAC/excel/Microike_both.csv')


marker <- c("IGF2R","VCAN","S100A6","ANTXR2","EREG","FN1","S100A4",
            "S100A10","S100A2","CD38","ILR1","SELL","CXCL14","S100A9","S100A8","CXCL3","MPP1",
            "SPP1","C1QA","C1QC","C1QB","CCL5","ITGB5","IGF1","CD8A","TLR3")
marker_RNA <- c("CCR2","VCAN","PECAM1","FGR","LYZ","S100A6","CD44","S100A9",
            "VIM","S100A10","S100A8","ITGAX","IFI30","ITGB2","CCL5","IL6","IL18",
            "CCL8","GPR34","SPP1","P2RY12","CCL2","MERTK","CXCL8","TLR3","TLR7","TLR1")

# Plot Heatmap
library(ComplexHeatmap)
Marker_Complx_htmap <- function(matrix, col, label, type){
  mtx <- matrix
  lb <- label
  myCol <- col
  type <- type
  ha <- ComplexHeatmap::rowAnnotation(foo = anno_mark(at = which(rownames(mtx) %in% lb),
                                                      labels = rownames(mtx)[rownames(mtx)%in%lb],
                                                      labels_gp = gpar(fontsize = 15)),
                                      width =  unit(50,"mm"))
  myCol <- colorRampPalette(myCol)(100)
  myBreaks <- seq(-2, 2, length.out = 100)
  col_fun <- circlize::colorRamp2(myBreaks, myCol)
  column_ha <- ComplexHeatmap::HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 1:ncol(mtx)), labels =colnames(mtx),
                                                                  labels_gp = gpar(col = "white", fontsize = 20, fontface = "bold")),
                                                 height = unit(2.5,"cm"))
  if(type == "RNA"){
    p <- ComplexHeatmap::Heatmap(mtx, right_annotation = ha, show_row_names = F, col = col_fun,
                                 cluster_rows = F, cluster_columns = F,column_split = colnames(mtx),
                                 column_title = NULL,top_annotation = column_ha,show_column_names = F,
                                 heatmap_legend_param = list(title = "RNA Expression", direction = "horizontal",
                                                             at = c(-1.5, 0, 1.5), 
                                                             labels = c("Low", "Median", "High"),
                                                             labels_gp = gpar(fontsize = 15),
                                                             title_gp = gpar(fontsize = 15),
                                                             legend_height = unit(30, "mm"),
                                                             legend_width = unit(40, "mm")))
  }else if (type == "ATAC"){
    p <- ComplexHeatmap::Heatmap(mtx, right_annotation = ha, show_row_names = F, col = col_fun,
                                 cluster_rows = F, cluster_columns = F,column_split = colnames(mtx),column_title = NULL,top_annotation = column_ha,
                                 show_column_names = F, heatmap_legend_param = list(title = "Predicted Gene Activity", direction = "horizontal",
                                                                                    at = c(-1.5, 0, 1.5), 
                                                                                    labels = c("Low", "Median", "High"),
                                                                                    labels_gp = gpar(fontsize = 15),
                                                                                    title_gp = gpar(fontsize = 15),
                                                                                    legend_height = unit(30, "mm"),
                                                                                    legend_width = unit(40, "mm")))
  }else{
    stop("Please choose between RNA and ATAC as type")
  }
  p <- ComplexHeatmap::draw(p, heatmap_legend_side = "bottom")
  return(p)
}

heatmapGS <- plotMarkerHeatmap(
  seMarker = markerGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
  returnMatrix = T,
  plotLog2FC = T
)
allmarkers <- as.data.frame(markerList)
allmarkers %>%
  arrange(group_name,desc(Log2FC)) -> allmarkers
#Set the order of heatmappk matrix same with allmarkers
heatmapGS <- heatmapGS[match(allmarkers$name,rownames(heatmapGS)), ]
heatmapGS <- heatmapGS[unique(rownames(heatmapGS)),]

p1 <- Marker_Complx_htmap(matrix = heatmapGS, label = marker, col = ArchRPalettes$blueYellow, type = "ATAC")
png(filename = "/work/foxlab/xiaoke/scATAC/Figures/Heatmap_atac.png", width = 4000, height = 4000, units = 'px',
    res = 300)
p1
dev.off()


heatmapRNA <- plotMarkerHeatmap(
  seMarker = markerRNA, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  returnMatrix = T,
  plotLog2FC = T
)
allmarkers <- as.data.frame(markerList)
allmarkers %>%
  arrange(group_name,desc(Log2FC)) -> allmarkers
#Set the order of heatmappk matrix same with allmarkers
heatmapRNA <- heatmapRNA[match(allmarkers$name,rownames(heatmapRNA)), ]
heatmapRNA <- heatmapRNA[unique(rownames(heatmapRNA)),]

#Use functiont to plot
p2 <- Marker_Complx_htmap(matrix = heatmapRNA, label = marker_RNA, col = ArchRPalettes$solarExtra, type = "RNA")
png(filename = "/work/foxlab/xiaoke/scATAC/Figures/Heatmap_rna.png", width = 4000, height = 4000, units = 'px',
    res = 300)
p2
dev.off()

# Violin plots
markerGenes <- c("S100A8","S100A9","S100A6","CD36","IL1R1","CD38",
                 "IL6","IL18","SPP1","CCL5","CCL3","CCL2")
p1 <- plotGroups(
  ArchRProj =  proj_5,
  groupBy = "Label_2",
  colorBy = "GeneExpressionMatrix",
  name = markerGenes,
  plotAs = "violin"
)

p2 <- lapply(p1, function(x){
  x + guides(color = "none", fill = "none") + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.text.y=element_text(size = 20), 
      axis.title = element_blank()
    )+
    labs(title = NULL)
})

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/RNA_expression.png", width = 6000, height = 5000, units = 'px',
    res = 300)
do.call(cowplot::plot_grid, c(list(ncol = 6),p2))
dev.off()

# chromatin accessibility
pkmarkers <- getMarkerFeatures(
  ArchRProj = proj_5,
  useMatrix = 'PeakMatrix',
  groupBy = 'Label_2',
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = 'wilcoxon'
)

p <- plotBrowserTrack(
  ArchRProj = proj_5,
  groupBy = "Label_2", 
  geneSymbol = c("S100A4","S100A8"),
  loops = getPeak2GeneLinks(proj_5),
  features =  getMarkers(pkmarkers,cutOff = "FDR <= 0.01 & Log2FC >= 1",returnGR = TRUE),
  baseSize = 7,
  sizes = c(8,1.5,2,2),
  facetbaseSize = 10,
  upstream = 2000,
  downstream = 40000,
)
png(filename = "/work/foxlab/xiaoke/scATAC/Figures/s100A4_6.png", width = 3500, height = 3500, units = 'px',
    res = 600)
grid::grid.draw(p$S100A4)
dev.off()

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/s100A8_9.png", width = 3500, height = 3500, units = 'px',
    res = 600)
grid::grid.draw(p$S100A8)
dev.off()

p <- plotBrowserTrack(
  ArchRProj = proj_5,
  groupBy = "Label_2", 
  geneSymbol = c("SPP1","CCL5"),
  loops = getPeak2GeneLinks(proj_5),
  features =  getMarkers(pkmarkers,cutOff = "FDR <= 0.01 & Log2FC >= 1",returnGR = TRUE),
  baseSize = 7,
  sizes = c(8,1.5,2,2),
  facetbaseSize = 10,
  upstream = 10000,
  downstream = 20000 
)

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/SPP1.png", width = 3500, height = 3500, units = 'px',
    res = 600)
grid::grid.draw(p$SPP1)
dev.off()

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/CCL5.png", width = 3500, height = 3500, units = 'px',
    res = 600)
grid::grid.draw(p$CCL5)
dev.off()


# Plot embedding
markerGenes  <- c("SPP1","CCL5","CCL2","CCL3","S100A8","S100A9","S100A4","S100A6")
p1 <- plotEmbedding(
  ArchRProj = proj_4,
  colorBy = "GeneExpressionMatrix",
  name = markerGenes[1:4],
  continuousSet = 'horizonExtra',
  embedding = "UMAP_combined",
  imputeWeights = getImputeWeights(proj_4)
)
p2 <- lapply(p1, function(x){
  x + guides(color = "none", fill = "none") + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
    theme(
      axis.text=element_blank(), 
      axis.title = element_blank(),
      axis.ticks = element_blank(),
    )+
    labs(title = NULL)
})

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/CCL.png", width = 6000, height = 1500, units = 'px',
    res = 300)
do.call(cowplot::plot_grid, c(list(ncol = 4),p2))
dev.off()

#Plot SIV
p <- plotBrowserTrack(
  ArchRProj = proj_5,
  groupBy = "Label", 
  geneSymbol = "SIV",
  upstream = 5000,
  downstream = 15000,
  loops = getPeak2GeneLinks(proj_5),
  baseSize = 7,
  facetbaseSize = 10
)
png(filename = "/work/foxlab/xiaoke/scATAC/Figures/SIV.png", width = 3500, height = 3500, units = 'px',
    res = 600)
grid::grid.draw(p$SIV)
dev.off()
