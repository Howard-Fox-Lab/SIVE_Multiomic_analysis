library(ArchR)
library(tidyr)
library(dplyr)
library(ggplot2)
library(GenomicFeatures)
set.seed(1234)
addArchRThreads(threads = 16)
setwd("/work/foxlab/xiaoke/scATAC/atac_enceph_new/")

### set up for genome and gene annotation
library("BSgenome.rheMac10.plusSIV.refMmul10")
Txdb <- loadDb("/work/foxlab/xiaoke/TxDb.Mmulatta.refseq.rheMac10.SIV.sqlite")
genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.rheMac10.plusSIV.refMmul10, filter = F)
id <- c(grep("NC_",seqnames(genomeAnnotation$chromSizes)), grep("SIV239",seqnames(genomeAnnotation$chromSizes)))
genomeAnnotation$chromSizes <- genomeAnnotation$chromSizes[id]
geneAnnotation <- createGeneAnnotation(
  TSS = transcripts(Txdb),
  exons = exons(Txdb),
  genes = GRanges(symbol=genes(Txdb)$gene_id,genes(Txdb))
)

### The path for fragment files
files <- list.files(path = "/work/foxlab/xiaoke/scATAC/atac_enceph/fragments/", pattern = ".gz",
                   recursive = F, full.names = T)
name <- list.files(path = "/work/foxlab/xiaoke/scATAC/atac_enceph/fragments/", pattern = ".gz",
                   recursive = F, full.names = F)
name <- gsub('_atac_fragments.tsv.gz','',name)

for (i in 1:length(files)) {
  names(files)[i] <- name[i]
}

addArchRChrPrefix(chrPrefix = FALSE)

### Create Arrowfiles
Arrowfile <- createArrowFiles(
  inputFiles = files,
  sampleNames = names(files),
  minTSS = 4,
  minFrags = 1000, ##Less than 1000 unique nuclear fragments will be filtered
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  subThreading = FALSE,
  excludeChr = c('NC_005943.1') #Exclude mitocondria genes
)

# Doublet identification
doubScores <- addDoubletScores(
  input = Arrowfile,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)

# Create ArchR project
proj_1 <- ArchRProject(
  ArrowFiles = Arrowfile, 
  outputDirectory = "proj_1",
  copyArrows = TRUE,
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation
)

## Separate columns in metadata
samplename <- data.frame(proj_1$Sample)
samplename <- separate(samplename, col = 'proj_1.Sample', into = c('species','name','condition'), sep = '_')
proj_1$name <- samplename$name
proj_1$condition <- samplename$condition
proj_1@cellColData

## save the ArchR project
saveArchRProject(ArchRProj = proj_1, outputDirectory = "proj_1", load = FALSE)

## number of cells in each sample
table(proj_1$name)

## Plotting QC metrics
proj_1 <- readRDS("/work/foxlab/xiaoke/scATAC/atac_enceph/proj_1/Save-ArchR-Project.rds")
df <- getCellColData(proj_3, select = c("log10(nFrags)", "TSSEnrichment"))
p <- ggPoint(
  x = df[,1], 
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")+
  theme(axis.title = element_text(size = 15), axis.text = element_text(size = 10), legend.key.size = unit(1,"cm"))
png(filename = "/work/foxlab/xiaoke/scATAC/Figures/frag_vs_TSS_1.png", width = 4000, height = 4000, units = 'px',
    res = 300)
p
dev.off()

# Remove doublets
proj_2 <- filterDoublets(proj_1)

### Integration with scRNA-seq data
files <- list.files(path = "/work/foxlab/xiaoke/scATAC/atac_enceph/GEX/",
                    recursive = F, full.names = T)
files <- files[c(1,2,5,6)]
name <- list.files(path = "/work/foxlab/xiaoke/scATAC/atac_enceph/GEX/",
                   recursive = F, full.names = F)
name <- name[c(1,2,5,6)]
name <- gsub('_filtered_feature_bc_matrix.h5','',name)

for (i in 1:length(files)) {
  names(files)[i] <- name[i]
}

seRNA <- import10xFeatureMatrix(
  input = files,
  names = names(files)
)

proj_2 <- addGeneExpressionMatrix(input = proj_2, seRNA = seRNA, force = TRUE)
proj_2 <- proj_2[!is.na(proj_2$Gex_nUMI)]
table(proj_2$name)

saveArchRProject(ArchRProj = proj_2, outputDirectory = "proj_2_wRNA", load = FALSE)

######Add mitochondria percentage into the arch R dataset ############
# Calculate %mitocondria in seurat
library(Seurat)
library(SeuratDisk)
files <- list.files(path = "/work/foxlab/xiaoke/scATAC/atac_enceph/GEX", 
                    recursive = F, full.names = F)
files <- files[c(1,2,5,6)]

for (x in files){
  name <- gsub('_filtered_feature_bc_matrix.h5','',x)
  x <- Read10X_h5(paste0("/work/foxlab/xiaoke/scATAC/atac_enceph/GEX/",x))
  assign(name,CreateSeuratObject(counts = x$`Gene Expression`, min.cell =10))
}

merged_seurat <- merge(RM_104T_Uninfected, y = c(RM_106T_Uninfected,RM_21T_Infected, RM_34T_Infected),
                       add.cell.ids = ls()[12:15])

mt <- read.csv('/work/foxlab/xiaoke/seurat/mt.genes.txt')
mt <- gsub("_","-",mt$mt)
merged_seurat$mtUMI <- Matrix::colSums(merged_seurat[which(rownames(merged_seurat) %in% mt),], na.rm = T)
merged_seurat$mitoPercent <- merged_seurat$mtUMI*100/merged_seurat$nCount_RNA
SaveH5Seurat(merged_seurat, file = "/work/foxlab/xiaoke/scATAC/atac_enceph_new/merged_seurat.h5seurat")

# Fetch the cell names and mitochondria percentage from seurat and archR
seurat_cells <- FetchData(merged_seurat, vars = c('mtUMI','mitoPercent'))
archR_cells <- gsub("#","_",getCellNames(proj_2))

#Filter the seurat cells based on cells from archR
seurat_cells <- seurat_cells[which(rownames(seurat_cells) %in% archR_cells),]
#Change the format of cell barcode back to archR format
seurat_cells$cell_id <-rownames(seurat_cells)
seurat_cells$barcode <- sub("([^_]+)_([^_]+)_([^_]+)_","#", seurat_cells$cell_id)
seurat_cells$name <- sub("([^_]+)_([^_]+)_([^_]+).*","\\1_\\2_\\3", seurat_cells$cell_id)
seurat_cells$cell_id <- paste0(seurat_cells$name, seurat_cells$barcode)
#Add to the metadata of archR project 
proj_2$Gex_MitoRatio <- seurat_cells[match(row.names(proj_2@cellColData),seurat_cells$cell_id),"mitoPercent"]
proj_2$Gex_MitoUMI <- seurat_cells[match(row.names(proj_2@cellColData),seurat_cells$cell_id),"mtUMI"]

#Save ArchR project
saveArchRProject(ArchRProj = proj_2, outputDirectory = "proj_2_wRNA", load = FALSE)

#########################################
# Filter based on RNA data
proj_3 <- proj_2[proj_2$Gex_nUMI>400 & 
                   proj_2$Gex_nGenes>400 &  proj_2$Gex_MitoRatio <15] 
saveArchRProject(ArchRProj = proj_3, outputDirectory = "proj_3_wRNA(filter)", load = FALSE)
table(proj_3$name)

# QC plots for individual sample
# Plot TSS score
proj_3 <- readRDS("/work/foxlab/xiaoke/scATAC/atac_enceph_new/proj_3_wRNA(filter)/Save-ArchR-Project.rds")
p1 <- plotGroups(
  proj_3,
  groupBy = 'name',
  name = 'TSSEnrichment',
  plotAs = 'volin'
)+theme(axis.title = element_blank(), axis.text = element_text(size = 15), 
        plot.title = element_text(hjust = 0.5, size = 18),plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"))+
  ggtitle("TSS Encihemnt Score")
# Plot Fragment size
p2 <- plotGroups(
  proj_3,
  groupBy = 'name',
  name = 'log10(nFrags)',
  plotAs = 'violin'
)+theme(axis.title = element_blank(), axis.text = element_text(size = 15), 
        plot.title = element_text(hjust = 0.5, size = 18),plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"))+
  ggtitle("Unique Fragments")

# Plot for GEX_UMI
p3 <- plotGroups(
  ArchRProj =  proj_3,
  groupBy = "name",
  name = "log10(Gex_nUMI)",
  plotAs = "violin"
)+theme(axis.title = element_blank(), axis.text = element_text(size = 15), 
        plot.title = element_text(hjust = 0.5, size = 18),plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"))+
  ggtitle("UMI Counts")

# Plot for GEX_Genes
p4 <- plotGroups(
  ArchRProj =  proj_3,
  groupBy = "name",
  name = "log10(Gex_nGenes)",
  plotAs = "violin"
)+theme(axis.title = element_blank(), axis.text = element_text(size = 15), 
        plot.title = element_text(hjust = 0.5, size = 18),plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"))+
  ggtitle("Gene Counts")

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/QC_persample.png", width = 3000, height = 3000, units = 'px',
    res = 300)
cowplot::plot_grid(p1,p2,p3, p4, ncol = 2, nrow = 2)
dev.off()

### Reduce dimensionality
#LSI for ATAC data
proj_3 <- addIterativeLSI(
  ArchRProj = proj_3, 
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "TileMatrix", 
  depthCol = "nFrags",
  name = "LSI_ATAC",
  force = TRUE
)

#LSI for RNA data
proj_3 <- addIterativeLSI(
  ArchRProj = proj_3, 
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "GeneExpressionMatrix", 
  depthCol = "Gex_nUMI",
  varFeatures = 2500,
  firstSelection = "variable",
  binarize = FALSE,
  name = "LSI_RNA",
  force = TRUE
)

# Combine two Dims
proj_3 <- addCombinedDims(proj_3, reducedDims = c("LSI_ATAC", "LSI_RNA"), name =  "LSI_Combined")
proj_3@reducedDims

### Batch Effect correction
proj_3 <- addHarmony(
  ArchRProj = proj_3,
  reducedDims = "LSI_Combined",
  name = "Harmony_combined",
  groupBy = "name"
)

### Add clustering
proj_3 <- addUMAP(proj_3, reducedDims = "Harmony_combined", name = "UMAP_combined", minDist = 0.8, force = TRUE)
proj_3 <- addClusters(proj_3, reducedDims = 'Harmony_combined', name = "Clusters_combined", resolution = 0.2, force = TRUE)

#Plot Embeding
p <- plotEmbedding(ArchRProj = proj_3, colorBy = "cellColData", name = "Clusters_combined", embedding = "UMAP_combined",
                    labelSize = 8, labelAsFactors = F)
p1 <- p+theme(legend.text = element_text(size = 15), legend.position = "none", plot.margin = unit(c(0,0,0,0),"cm"))+
  labs(title = NULL)

p2 <- plotEmbedding(ArchRProj = proj_3, colorBy = "cellColData", name = "name", embedding = "UMAP_combined", labelSize = 5)+
  theme(legend.text = element_text(size = 15),legend.position = "none", plot.margin = unit(c(0,0,0,0),"cm"))+
  labs(title = NULL)

p <- plotEmbedding(ArchRProj = proj_3, colorBy = "cellColData", name = "condition", embedding = "UMAP_combined", 
                   labelSize = 8, labelAsFactors = F)
p3 <- p+ theme(legend.text = element_text(size = 15),legend.position = "none",plot.margin = unit(c(0,0,0,0),"cm"))+
  labs(title = NULL)

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/UMAP_cluster.png", width = 3000, height = 3000, units = 'px',
    res = 600)
p1
dev.off()

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/UMAP_name.png", width = 3000, height = 3000, units = 'px',
    res = 600)
p2
dev.off()

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/UMAP_condi.png", width = 3000, height = 3000, units = 'px',
    res = 600)
p3
dev.off()

saveArchRProject(ArchRProj = proj_3, outputDirectory = "proj_3_wRNA(filter)", load = FALSE)

### Characterization
# FeaturePlot for RNA expression
proj_3 <- addImputeWeights(proj_3, reducedDims = 'Harmony_combined')
markerGenes  <- c(
  "P2RY12",'CX3CR1','GPR34', #Microglia
  "MAMU-DRA",'MAMU-DRB1','MAMU-DRB5', # CAM
  'CD3D','CD3E','GZMB' # T cells
)

#### Gene expression/Score from scRNA data
p1 <- plotEmbedding(
  ArchRProj = proj_3,
  colorBy = "GeneExpressionMatrix",
  name = markerGenes,
  continuousSet = 'horizonExtra',
  embedding = "UMAP_combined",
  imputeWeights = getImputeWeights(proj_3)
)

p2 <- lapply(p1, function(x){
  x + guides(color = "none", fill = "none") + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm")) +
    theme(
      axis.text=element_blank(), 
      axis.title = element_blank()
    )+
    labs(title = NULL)
})

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/markers_RNA.png", width = 5000, height = 5000, units = 'px',
    res = 300)
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
dev.off()

p1 <- plotEmbedding(
  ArchRProj = proj_3,
  colorBy = "GeneScoreMatrix",
  name = markerGenes,
  continuousSet = 'horizonExtra',
  embedding = "UMAP_combined",
  imputeWeights = getImputeWeights(proj_3)
)

p2 <- lapply(p1, function(x){
  x + guides(color = "none", fill = "none") + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm")) +
    theme(
      axis.text=element_blank(), 
      axis.title = element_blank()
    )+
    labs(title = NULL)
})

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/markers_ATAC.png", width = 5000, height = 5000, units = 'px',
    res = 300)
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
dev.off()


# Violin Plots
p1 <- plotGroups(
  ArchRProj =  proj_3,
  groupBy = "Clusters_combined",
  colorBy = "GeneExpressionMatrix",
  name = markerGenes,
  plotAs = "violin"
)
p2 <- lapply(p1, function(x){
  x + guides(color = "none", fill = "none") + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
    theme(
      axis.text=element_text(size = 18, angle = 90), 
      axis.title = element_blank()
    )+
    labs(title = NULL)
})

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/RNA_violin.png", width = 4000, height = 4000, units = 'px',
    res = 300)
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
dev.off()


p1 <- plotGroups(
  ArchRProj =  proj_3,
  groupBy = "Clusters_combined",
  colorBy = "GeneScoreMatrix",
  name = markerGenes,
  plotAs = "violin"
)
p2 <- lapply(p1, function(x){
  x + guides(color = "none", fill = "none") + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
    theme(
      axis.text=element_text(size = 18, angle = 90), 
      axis.title = element_blank()
    )+
    labs(title = NULL)
})

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/ATAC_violin.png", width = 4000, height = 4000, units = 'px',
    res = 300)
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))
dev.off()

saveArchRProject(ArchRProj = proj_3, outputDirectory = "proj_3_wRNA(filter)", load = FALSE)

#Exclude C1 and C11
id <- BiocGenerics::which(proj_3$Clusters_combined %in% c('C2','C3','C4','C5','C6','C7','C8','C9','C10','C12'))
cells <- proj_3$cellNames[id]
proj_4 <- subsetArchRProject(proj_3, cells = cells, outputDirectory = "Myeloid_cells", force = T,dropCells =FALSE)

#### Rename the clusters
label_old <- unique(proj_4$Clusters_combined)
label_new <- c("CAM_like_1","Micro_like","CAM_like_2","CAM_like_3","Micro_1",'CAM_1',"Micro_2",
               "Micro_3","Micro_4","CAM_2")
proj_4$Label <- mapLabels(proj_4$Clusters_combined, newLabels = label_new, oldLabels = label_old)
saveArchRProject(ArchRProj = proj_4, outputDirectory = "Myeloid_cells", load = FALSE)

label_old <- unique(proj_4$Label)
label_new <- c("CAM_like","Micro_like","CAM_like","CAM_like","Micro",'CAM',"Micro",
               "Micro","Micro","CAM")
proj_4$Label_2 <- mapLabels(proj_4$Label, newLabels = label_new, oldLabels = label_old)
saveArchRProject(ArchRProj = proj_4, outputDirectory = "Myeloid_cells", load = FALSE)

proj_4$Label_3 <- paste0(paste0(proj_4$Label, "("), paste0(proj_4$Clusters_combined,")"))

p1 <- plotEmbedding(proj_4, colorBy = "cellColData", name = 'Label_3', embedding = "UMAP_combined",
                    labelAsFactors = F, labelSize = 6)+
  theme(legend.position = "none", plot.title = element_blank(), axis.title = element_text(size = 15))

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/UMAP_label.png", width = 5500, height = 5500, units = 'px',
    res = 600)
p1
dev.off()

### Find markers for ATAC data
proj_4 <- readRDS("/work/foxlab/xiaoke/scATAC/atac_enceph_new/Myeloid_cells/Save-ArchR-Project.rds")
markerGenes  <- c(
  "P2RY12",'CX3CR1','GPR34','SALL1', #Microglia
  "MAMU-DRA",'MAMU-DRB1','MAMU-DRB5','CD74' # CAM
)

#### Gene expression/Score from scRNA data
p1 <- plotEmbedding(
  ArchRProj = proj_4,
  colorBy = "GeneExpressionMatrix",
  name = markerGenes,
  continuousSet = 'horizonExtra',
  embedding = "UMAP_combined",
  imputeWeights = getImputeWeights(proj_4)
)

p2 <- lapply(p1, function(x){
  x + guides(color = "none", fill = "none") + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm")) +
    theme(
      axis.text=element_blank(), 
      axis.title = element_blank()
    )+
    labs(title = NULL)
})

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/marekers_RNA.png", width = 6000, height = 3000, units = 'px',
    res = 300)
do.call(cowplot::plot_grid, c(list(ncol = 4),p2))
dev.off()

p1 <- plotEmbedding(
  ArchRProj = proj_4,
  colorBy = "GeneScoreMatrix",
  name = markerGenes,
  continuousSet = 'horizonExtra',
  embedding = "UMAP_combined",
  imputeWeights = getImputeWeights(proj_4)
)

p2 <- lapply(p1, function(x){
  x + guides(color = "none", fill = "none") + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm")) +
    theme(
      axis.text=element_blank(), 
      axis.title = element_blank()
    )+
    labs(title = NULL)
})

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/marekers_ATAC.png", width = 6000, height = 3000, units = 'px',
    res = 300)
do.call(cowplot::plot_grid, c(list(ncol = 4),p2))
dev.off()

markerGS <- getMarkerFeatures(
  ArchRProj = proj_4,
  useMatrix = 'GeneScoreMatrix',
  groupBy = 'Label',
  bias = c('TSSEnrichment','log10(nFrags)'),
  testMethod = 'wilcoxon'
)
markerList <- getMarkers(markerGS, cutOff = 'FDR <= 0.01 & Log2FC >= 1.25')
write.csv(markerList, file = '/work/foxlab/xiaoke/scATAC/excel/markers_atac.csv')

markerRNA <- getMarkerFeatures(
  ArchRProj = proj_4,
  useMatrix = 'GeneExpressionMatrix',
  groupBy = 'Label',
  bias = c('TSSEnrichment','log10(nFrags)'),
  testMethod = 'wilcoxon'
)
markerList <- getMarkers(markerRNA, cutOff = 'FDR <= 0.01 & Log2FC >= 1.25')
write.csv(markerList, file = '/work/foxlab/xiaoke/scATAC/excel/markers_RNA.csv')

# Plot Gene browser
pkmarkers <- getMarkerFeatures(
  ArchRProj = proj_4,
  useMatrix = 'PeakMatrix',
  groupBy = 'Label_2',
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = 'wilcoxon'
)

p <- plotBrowserTrack(
  ArchRProj = proj_4,
  groupBy = "Label_2",
  plotSummary = c("bulkTrack","featureTrack","geneTrack"),
  geneSymbol = c("P2RY12","GPR34","CX3CR1","SALL1"),
  features =  getMarkers(pkmarkers,cutOff = "FDR <= 0.05 & Log2FC >= 1",returnGR = TRUE),
  baseSize = 7,
  sizes = c(8,1.5,2),
  facetbaseSize = 10,
  upstream = 20000,
  downstream = 3000 
)
png(filename = "/work/foxlab/xiaoke/scATAC/Figures/P2RY12.png", width = 3000, height = 3000, units = 'px',
    res = 600)
grid::grid.draw(p$P2RY12)
dev.off()

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/GPR34.png", width = 3000, height = 3000, units = 'px',
    res = 600)
grid::grid.draw(p$GPR34)
dev.off()

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/CX3CR1.png", width = 3000, height = 3000, units = 'px',
    res = 600)
grid::grid.draw(p$CX3CR1)
dev.off()

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/SALL1.png", width = 3000, height = 3000, units = 'px',
    res = 600)
grid::grid.draw(p$SALL1)
dev.off()


#Plotting the heatmap
library(ComplexHeatmap)
#Define a new function for plotting
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
                                                                  labels_gp = gpar(col = "white", fontsize = 15, fontface = "bold")))
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

marker_ATAC <- c("NCAM1","CD109","CD69","CD82","IL1R1","SELL","CSTA","GZMA","NLRP2","S100P","VCAN",
                 "S100A12","S100A8","S100A9","CD38","MCF2L","SALL1","TLR3","CCL20","NCR1","OLFM1","SPP1","JAK3","CCL5")

marker_RNA <- c("MAMU-DRA","MAMU-DRB1","MAMU-DQA1","CD74","CD36",
                "ITGA4","FGR","ITGAX","ITGAL","PECAM1","CD47","ITGB2",
                "VIM","IFI6","CD44","B2M","IFI30","CCR2","LYZ",
                "VCAN","S100A8","S100A9","S100A10","S100A6",
                "P2RY12","GPR34","CX3CR1","MERTK","MRC1","SELPLG",
                "SPP1","APOE","CCL5","CXCL8","CD163","CCL2","CCL8","CXCL10")

heatmapGS <- plotMarkerHeatmap(
  seMarker = markerGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  returnMatrix = T
)
allmarkers <- as.data.frame(markerList)
allmarkers %>%
  arrange(group_name,desc(Log2FC)) -> allmarkers
#Set the order of heatmappk matrix same with allmarkers
heatmapGS <- heatmapGS[match(allmarkers$name,rownames(heatmapGS)), ]
heatmapGS <- heatmapGS[unique(rownames(heatmapGS)),]

#Use functiont to plot
p1 <- Marker_Complx_htmap(matrix = heatmapGS, label = marker_ATAC, col = ArchRPalettes$blueYellow, type = "ATAC")
png(filename = "/work/foxlab/xiaoke/scATAC/Figures/Heatmap_atac.png", width = 5000, height = 4000, units = 'px',
    res = 300)
p1
dev.off()

heatmapRNA <- plotMarkerHeatmap(
  seMarker = markerRNA, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  returnMatrix = T
)
allmarkers <- as.data.frame(markerList)
allmarkers %>%
  arrange(group_name,desc(Log2FC)) -> allmarkers
#Set the order of heatmappk matrix same with allmarkers
heatmapRNA <- heatmapRNA[match(allmarkers$name,rownames(heatmapRNA)), ]
heatmapRNA <- heatmapRNA[unique(rownames(heatmapRNA)),]

#Use functiont to plot
p2 <- Marker_Complx_htmap(matrix = heatmapRNA, label = marker_RNA, col = ArchRPalettes$solarExtra, type = "RNA")
png(filename = "/work/foxlab/xiaoke/scATAC/Figures/Heatmap_rna.png", width = 5000, height = 4000, units = 'px',
    res = 300)
p2
dev.off()

#########SIV positive cells---------------------------
# FeaturePlot for gene expression (RNA)
# FeaturePlot for gene expression (atac)
proj_4 <- addImputeWeights(proj_4, reducedDims = 'Harmony_combined')
p <- plotEmbedding(
  ArchRProj = proj_4, 
  colorBy = "GeneExpressionMatrix", 
  name = 'SIV', 
  embedding = "UMAP_combined",
  imputeWeights = getImputeWeights(proj_4)
)

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/SIV_geneexp.png", width = 2000, height = 2000, units = 'px',
    res = 300)
p
dev.off()

# FeaturePlot for gene expression (atac)
p <- plotEmbedding(
  ArchRProj = proj_4, 
  colorBy = "GeneScoreMatrix", 
  name = 'SIV', 
  embedding = "UMAP_combined",
  imputeWeights = getImputeWeights(proj_4)
)

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/SIV_genescore.png", width = 2000, height = 2000, units = 'px',
    res = 300)
p
dev.off()


## calculate percentage of SIV positive cells
# Grab certain gene expression from gene matrix
getAvailableMatrices(proj_4)
GenMtrix <- getMatrixFromProject(proj_4, useMatrix = 'GeneScoreMatrix')
RNAMtrix <- getMatrixFromProject(proj_4, useMatrix = "GeneExpressionMatrix")
rownames(assay(GenMtrix, withDimnames=F)) <- rowData(GenMtrix)$name
rownames(assay(RNAMtrix, withDimnames=F)) <- rowData(RNAMtrix)$name
SIV_atac <- data.frame(SIV = assay(GenMtrix, withDimnames = F)['SIV',], 
                  Cluster = colData(GenMtrix)$Label,
                  Sample = colData(GenMtrix)$name)
SIV_RNA <- data.frame(SIV = assay(RNAMtrix, withDimnames = F)['SIV',], 
                      Cluster = colData(RNAMtrix)$Label,
                      Sample = colData(RNAMtrix)$name)

SIV.pos_atac <- filter(SIV_atac, SIV>0)
SIV.pos_RNA <- filter(SIV_RNA, SIV>0)

table(SIV.pos_RNA$Sample)
table(SIV.pos_atac$Sample)

# Add SIV counts into archR metadata
library(Seurat)
library(SeuratDisk)
merged_seurat <- LoadH5Seurat("/work/foxlab/xiaoke/scATAC/atac_enceph_new/merged_seurat.h5seurat")

siv_cells <- data.frame(FetchData(merged_seurat, vars = c('SIV')))
merged_seurat$SIV_count <- siv_cells$SIV
seurat_cells <- data.frame(FetchData(merged_seurat, vars = c('SIV','nCount_RNA')))
archR_cells <- gsub("#","_",getCellNames(proj_4))
#Filter seurat_cells
seurat_cells <- seurat_cells[which(rownames(seurat_cells) %in% archR_cells),]
#Change the format
seurat_cells$cell_id <-rownames(seurat_cells)
seurat_cells$barcode <- sub("([^_]+)_([^_]+)_([^_]+)_","#", seurat_cells$cell_id)
seurat_cells$name <- sub("([^_]+)_([^_]+)_([^_]+).*","\\1_\\2_\\3", seurat_cells$cell_id)
seurat_cells$cell_id <- paste0(seurat_cells$name, seurat_cells$barcode)
#Add to the metadata of archR project 
proj_4$SIV_count <- seurat_cells[match(row.names(proj_4@cellColData),seurat_cells$cell_id),"SIV"]
saveArchRProject(ArchRProj = proj_4, outputDirectory = "Myeloid_cells", load = FALSE)

#Calculate the co-expression
SIV <- data.frame(SIV_DNA = assay(GenMtrix, withDimnames = F)['SIV',], 
                  SIV_RNA = colData(GenMtrix)$SIV_count,
                  Cluster = colData(GenMtrix)$Label,
                  Sample = colData(GenMtrix)$name)
SIV.pos <- filter(SIV, SIV_RNA>0 & SIV_DNA>0)
table(SIV.pos$Cluster)
