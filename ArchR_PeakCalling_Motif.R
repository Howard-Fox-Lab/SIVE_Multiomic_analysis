library(ArchR)
library(tidyr)
library(dplyr)
library(ggplot2)
library(GenomicFeatures)
set.seed(1234)
addArchRThreads(threads = 16)
setwd("/work/foxlab/xiaoke/scATAC/atac_enceph_new/")
library("BSgenome.rheMac10.plusSIV.refMmul10")

proj_4 <- readRDS("/work/foxlab/xiaoke/scATAC/atac_enceph_new/Myeloid_cells/Save-ArchR-Project.rds")

#####################Peak Calling------------------------------
# Add Pseudo-bulk replicates
library("BSgenome.rheMac10.plusSIV.refMmul10")
proj_4 <- addGroupCoverages(ArchRProj = proj_4, groupBy = 'Label')
saveArchRProject(ArchRProj = proj_4, outputDirectory = "Myeloid_cells", load = FALSE)

##### Peak Calling using MACS2------------------------------------------
# Change the conda environment
# conda activate macs
# Load gene and genome annotation
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
pathToMacs2 <- findMacs2()
proj_4 <- readRDS("/work/foxlab/xiaoke/scATAC/atac_enceph_new/Myeloid_cells/Save-ArchR-Project.rds")

proj_4 <- addReproduciblePeakSet(
  ArchRProj = proj_4, 
  groupBy = "Label", 
  pathToMacs2 = pathToMacs2,
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  excludeChr = c('NC_005943.1'),
  genomeSize = '2.2e09' 
  # 78% of total genome size of rhesus macaques
)

# Add peak matrix
proj_4 <- addPeakMatrix(proj_4)
saveArchRProject(ArchRProj = proj_4, outputDirectory = "Myeloid_cells", load = FALSE)

# Calculate peak regions
pk <- getPeakSet(proj_4)
write.csv(pk, file = "/work/foxlab/xiaoke/scATAC/excel/pkcalling.csv", row.names = T)
df <- data.frame(table(pk$peakType), Var2 = rep("all",4))

files <- list.files("/work/foxlab/xiaoke/scATAC/atac_enceph_new/Myeloid_cells/PeakCalls/",
                    pattern = ".rds")
for (i in files) {
  names <- gsub("-reproduciblePeaks.gr.rds","",i)
  i <- readRDS(paste0("/work/foxlab/xiaoke/scATAC/atac_enceph_new/Myeloid_cells/PeakCalls/",i))
  x <- data.frame(table(i$peakType), Var2 = rep(names, 4))
  df <- rbind(df,x)
}

# ggplot
colors <- c("#403990","#00405b", "#008dca","#80a6e2", "#c0beb8","#f46f43", "#d70000", "#7d0000")
colors <- setNames(colors, levels(df$Var1))
plot1 <- ggplot(df, aes(x = Freq/1000, y = Var2, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_bar(stat = "identity") +
  theme(legend.title = element_blank())+
  scale_fill_manual(values = colors)+
  theme(axis.text.y = element_text(size = 15, face = 'bold'),
        axis.title.y = element_blank())+
  xlab("Number of peaks (x10^3)")
png(filename = "/work/foxlab/xiaoke/scATAC/Figures/pk_region.png", width = 2000, height = 3000, units = 'px',
    res = 300)
plot1
dev.off()


#obtain peak matrix for downstream analyses
getAvailableMatrices(proj_4)
pkMtrix <- getMatrixFromProject(proj_4, useMatrix = 'PeakMatrix')
pkdf <- as.data.frame(assay(pkMtrix, withDimnames = F))

##### Identify marker Peaks
pkmarkers <- getMarkerFeatures(
  ArchRProj = proj_4,
  useMatrix = 'PeakMatrix',
  groupBy = 'Label',
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = 'wilcoxon'
)
markerList <- getMarkers(pkmarkers, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
write.csv(markerList, file = '/work/foxlab/xiaoke/scATAC/excel/peakmarkers.csv')

# Plotting Peak Heatmap
heatmappk <- plotMarkerHeatmap(
  seMarker = pkmarkers,
  cutOff = "FDR <= 0.05 & Log2FC >= 1",
  transpose = TRUE,
  nLabel = 0
)
p1 <- draw(heatmappk, heatmap_legend_side = "bot")
png(filename = "/work/foxlab/xiaoke/scATAC/Figures/pk_heatmap.png", width = 2000, height = 2000, units = 'px',
    res = 300)
p1
dev.off()

######### Motif Analyses---------------------------------
proj_4 <- addMotifAnnotations(ArchRProj = proj_4, motifSet = "JASPAR2020", force = TRUE,
                              annoName = "Motif", species = "Homo sapiens")
proj_4 <- addMotifAnnotations(ArchRProj = proj_4, motifSet = "vierstra", collection = "archetype", force = TRUE,
                              annoName = "Motif_vier", species = "Homo Sapiens") #v2 archetype clustered pwms

saveArchRProject(ArchRProj = proj_4, outputDirectory = "Myeloid_cells", load = FALSE)

pkmarkers <- getMarkerFeatures(
  ArchRProj = proj_4,
  useMatrix = 'PeakMatrix',
  groupBy = 'Label',
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = 'wilcoxon'
)

enrichMotifs <- peakAnnoEnrichment(
  seMarker = pkmarkers,
  ArchRProj = proj_4,
  peakAnnotation = "Motif_vier",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

enrichMotifs <- peakAnnoEnrichment(
  seMarker = pkmarkers,
  ArchRProj = proj_4,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

heatmapEM <- plotEnrichHeatmap(enrichMotifs, returnMatrix = T)

#For Vierstra:
rownames(heatmapEM) <- ifelse(
  nchar(rownames(heatmapEM)) > 20,
  paste0(str_sub(rownames(heatmapEM),1,20),"..."),
  rownames(heatmapEM))

#For JASPAR:
rownames(heatmapEM) <- gsub("_.*","",rownames(heatmapEM))

myCol <- ArchRPalettes$comet
myCol <- colorRampPalette(myCol)(4)
myBreaks <- seq(0, 100, length.out = 4)
col_fun <- circlize::colorRamp2(myBreaks, myCol)

column_ha <- ComplexHeatmap::HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 1:10), labels =colnames(heatmapEM),
                                                                labels_gp = gpar(col = "white", fontsize = 14, fontface = "bold"),
                                                                height = unit(1, "cm")))

p <- ComplexHeatmap::Heatmap(heatmapEM,row_names_gp = gpar(fontsize = 13), col = col_fun, cluster_rows = F, cluster_columns = F,column_split = colnames(mtx),
                             column_title = NULL,top_annotation = column_ha,show_column_names = F,width = ncol(mtx)*unit(35,"mm"), height = nrow(mtx)*unit(6,"mm"),
                             heatmap_legend_param = list(title = "Norm.Enrichment -log(P-adj)[0-Max]", direction = "horizontal",at = c(0, 50, 100), 
                                                         labels = c("Low", "Median", "High"),
                                                         labels_gp = gpar(fontsize = 15),
                                                         title_gp = gpar(fontsize = 15), 
                                                         legend_height = unit(30, "mm"),
                                                         legend_width = unit(40, "mm")))
p <- ComplexHeatmap::draw(p, heatmap_legend_side = "bottom")
png(filename = "/work/foxlab/xiaoke/scATAC/Figures/motifjaspr_heatmap.png", width = 5000, height = 2500, units = 'px',
    res = 300)
p
dev.off()

# Export the motif enrichment for each cluster
df1 <- assays(enrichMotifs)$mlog10Padj
write.csv(df1, file = "/work/foxlab/xiaoke/scATAC/excel/motif.enrchment.csv")

## Correct for deviations for motif enrichment
proj_4 <- addBgdPeaks(proj_4, method = "ArchR",force = TRUE)
saveArchRProject(ArchRProj = proj_4, outputDirectory = "Myeloid_cells", load = FALSE)
proj_4 <- addDeviationsMatrix(proj_4, peakAnnotation = "Motif", force = TRUE)
proj_4 <- addDeviationsMatrix(proj_4, peakAnnotation = "Motif_vier", force = TRUE)

saveArchRProject(ArchRProj = proj_4, outputDirectory = "Myeloid_cells", load = FALSE)
getAvailableMatrices(proj_4)

## Plot the variability of the motif
options(ggrepel.max.overlaps = Inf)
library(ggrepel)
df <- getVarDeviations(proj_4, name = "MotifMatrix", plot = F)
df <- getVarDeviations(proj_4, name = "Motif_vierMatrix", plot = F)
#For Jaspar
df$name <- sub("_[^_]+$", "",df$name)
#For vierstra
df <- as.data.frame(df)
df <- df %>%
  separate(col=name, into=c("B","C","D"), sep="\\s*\\|\\s*", fill="right")


p <- ggplot(df, aes(x = rank, y = combinedVars, color = combinedVars))+
  geom_point(size = 3)+
  geom_label_repel(data = df[rev(seq_len(40)),], aes(x = rank, y = combinedVars, label = name),
                   size = 5, nudge_x = 2, color = "black")+
  theme_ArchR()+
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))+
  ylab("Variability")+
  xlab("Rank Sorted Annotation")+
  theme(legend.position = "none", plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"),
        axis.title = element_text(size = 15), axis.text = element_text(size = 12))

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/motif_varibility.png", width = 3000, height = 3000, units = 'px',
    res = 600)
p
dev.off()

# Heatmap for motif deviation
getAvailableMatrices(proj_4)
mtfmkr <- getMarkerFeatures(
  ArchRProj = proj_4,
  useMatrix = "Motif_vierMatrix",
  groupBy = 'Label',
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = 'wilcoxon',
  useSeqnames = "z"
)

markerList <- getMarkers(mtfmkr, cutOff = 'FDR <= 0.001 & MeanDiff >=1')
write.csv(markerList, file = '/work/foxlab/xiaoke/scATAC/excel/markermotif_jaspr.csv')

heatmapmotif <- plotMarkerHeatmap(
  seMarker = mtfmkr, 
  cutOff = "FDR <0.001 & MeanDiff >= 1", 
  log2Norm = FALSE,
  returnMatrix = T
)

#The prediction from SCENIC
df <- read.csv("/work/foxlab/xiaoke/scATAC/Figures/regulon activity.csv")
rownames(df) <- df$X
df <- df[,-1]
df <- df[,c(6,2,9,8,7,1,4,3,5,10)]
same <- intersect(rownames(df),rownames(heatmapmotif))
diff <- setdiff(rownames(df), rownames(heatmapmotif))
df.col <- data.frame(TF = c(same,diff), color = c(rep("red",10),rep("black",15)))
df.col <- df.col[match(rownames(df),df.col$TF),]

#Arrange the allmarkers by log2FC
allmarkers <- as.data.frame(markerList)
allmarkers %>%
  arrange(group_name,desc(MeanDiff))-> allmarkers

allmarkers$name_2 <- ifelse(
  nchar(allmarkers$name) > 20,
  paste0(str_sub(allmarkers$name,1,20),"..."),
  allmarkers$name)

allmarkers %>%
  group_by(group_name)%>%
  slice_head(n=10)%>%
  ungroup()-> top10

label <- data.frame(TF = c(same,top10$name), 
          color = c(rep("red",length(same)),rep("black",length(top10$name))))
label <- unique(label)

#Manipulate with heatmapmotif matrix
rownames(heatmapmotif) <- gsub("_.*","",rownames(heatmapmotif))

#Set the order of heatmapmotif matrix same with allmarkers
heatmapmotif <- heatmapmotif[match(allmarkers$name,rownames(heatmapmotif)), ]
# For vierstra only: rownames(heatmapmotif) <- allmarkers$name_2
heatmapmotif <- heatmapmotif[unique(rownames(heatmapmotif)),]

#MATCH the label
label <- label[match(rownames(heatmapmotif), label$TF),]
label <- na.omit(label)

library(ComplexHeatmap)
myCol <- ArchRPalettes$solarExtra
myCol <- colorRampPalette(myCol)(100)
myBreaks <- seq(-2, 2, length.out = 100)
col_fun <- circlize::colorRamp2(myBreaks, myCol)

# Heatmap for different Treatment
ha <- ComplexHeatmap::rowAnnotation(foo = anno_mark(at = which(rownames(heatmapmotif) %in% top10$name_2),
                                                    labels = rownames(heatmapmotif)[rownames(heatmapmotif)%in%top10$name_2],
                                                    labels_gp = gpar(fontsize = 15)),
                                    width =  unit(50,"mm"))

column_ha <- ComplexHeatmap::HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 1:10), labels =colnames(heatmapmotif),
                                                                labels_gp = gpar(col = "white", fontsize = 14, fontface = "bold"),
                                                                height = unit(1, "cm")))

p <- ComplexHeatmap::Heatmap(heatmapmotif, show_row_names = F, col = col_fun, cluster_rows = F, cluster_columns = F,column_split = colnames(heatmapmotif),
                             column_title = NULL,top_annotation = column_ha,show_column_names = F, right_annotation = ha, width = ncol(heatmapmotif)*unit(35,"mm"),
                             heatmap_legend_param = list(title = "Deviation z score", direction = "horizontal",at = c(-1.5, 0, 1.5), 
                                                         labels = c("Low", "Median", "High"),
                                                         labels_gp = gpar(fontsize = 15),
                                                         title_gp = gpar(fontsize = 15), 
                                                         legend_height = unit(30, "mm"),
                                                         legend_width = unit(40, "mm")))
p <- ComplexHeatmap::draw(p, heatmap_legend_side = "bottom")
png(filename = "/work/foxlab/xiaoke/scATAC/Figures/motif_heatmap_vier.png", width = 6000, height = 7000, units = 'px',
    res = 300)
p
dev.off()

# Heatmap for SCENIC prediction
myCol <- ArchRPalettes$solarExtra
myCol <- colorRampPalette(myCol)(100)
myBreaks <- seq(0, 1, length.out = 100)
col_fun <- circlize::colorRamp2(myBreaks, myCol)

mtx <- as.matrix(df)
column_ha <- ComplexHeatmap::HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 1:10), labels =colnames(mtx),
                                                                labels_gp = gpar(col = "white", fontsize = 14, fontface = "bold"),
                                                                height = unit(1, "cm")))

p <- ComplexHeatmap::Heatmap(mtx,row_names_gp = gpar(col = df.col$color, fontsize = 20), col = col_fun, cluster_rows = T, cluster_columns = F,column_split = colnames(mtx),
                             column_title = NULL,top_annotation = column_ha,show_column_names = F,width = ncol(mtx)*unit(35,"mm"), height = nrow(mtx)*unit(8,"mm"),
                             heatmap_legend_param = list(title = "Regulon activity score", direction = "horizontal",at = c(0, 0.5, 1), 
                                                         labels = c("Low", "Median", "High"),
                                                         labels_gp = gpar(fontsize = 15),
                                                         title_gp = gpar(fontsize = 15), 
                                                         legend_height = unit(30, "mm"),
                                                         legend_width = unit(40, "mm")))
p <- ComplexHeatmap::draw(p, heatmap_legend_side = "bottom")
png(filename = "/work/foxlab/xiaoke/scATAC/Figures/motif_heatmap_scenic.png", width = 5000, height = 4000, units = 'px',
    res = 300)
p
dev.off()

#####Plot marker motifs in embedding
motif <- c("FOSL1","JUND","BACH2","NFE2","CEBPB","CEBPE","CTCF","SPIB",
           "IRF1","STAT1..STAT2","NFKB1","NFKB2","REL","RELA")
motif <- c("IRF7","STAT1","NFKB1","IKZF1","FLI1")

motif <- c("FOSL/JUND","BACH/NFE","NFKB/RELA","SPI/BCL11A","CEBPB/CEBPE","ELF/SPIB","CTCF/CTCFL","IRF/STAT")
motif <- c("IRF1","IRF3","IRF7","IRF9")

markermotifs <-getFeatures(proj_4, select = paste(motif, collapse = "|"), useMatrix = "Motif_vierMatrix")
markermotifs <-getFeatures(proj_4, select = paste(motif, collapse = "|"), useMatrix = "MotifMatrix")

markermotifs <- grep("z:", markermotifs, value = TRUE)

motif <- paste0("z:",motif)
motif_1 <- gsub("_.*","",markermotifs)
markermotifs <- markermotifs[which(motif_1 %in% motif)]

p1 <- plotEmbedding(
  proj_4,
  colorBy = "MotifMatrix",
  name = markermotifs[c(3,5,2,1,4)],
  embedding = "UMAP_combined",
  imputeWeights = getImputeWeights(proj_4)
)

p1 <- plotGroups(
  ArchRProj =  proj_4,
  groupBy = "Label",
  colorBy = "MotifMatrix",
  name = markermotifs[c(3,5,2,1,4)],
  plotAs = "violin"
)

p2 <- plotEmbedding(
  proj_4,
  colorBy = "GeneScoreMatrix",
  name = motif,
  embedding = "UMAP_combined",
  imputeWeights = getImputeWeights(proj_4)
)

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/ELF.png", width = 4000, height = 4000, units = 'px',
    res = 600)
p1[8]
dev.off()

p3 <- lapply(p1, function(x){
  x + guides(color = "none", fill = "none") + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x = element_text(angle = 45),
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank(),
      axis.title = element_blank()
    )
})

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/archR.png", width = 7500, height = 1500, units = 'px',
    res = 600)
do.call(cowplot::plot_grid, c(list(ncol = 5),p3))
dev.off()

#####FootPrinting----------------
motif <- c("NFKB1","NFKB2","STAT1..STAT2","IRF1","REL","RELA")
motif <- c("FOSL/JUND","BACH/NFE","NFKB/RELA","SPI/BCL11A","CEBPB/CEBPE","ELF/SPIB","CTCF/CTCFL","IRF/STAT")

motifpositions <- getPositions(proj_4, name ="Motif_vier")

markerMotifs <- unlist(lapply(motif, function(x) grep(x, names(motifpositions), value = TRUE)))
markerMotifs <- markerMotifs[c(1:6)]
motifpositions <- motifpositions[!seqnames(motifpositions) %in% c("NC_041761.1")]

seFoot <- getFootprints(
  ArchRProj = proj_4,
  positions = motifpositions[markerMotifs],
  groupBy = "Label_2"
)

p1 <- plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj_4, 
  normMethod = "Subtract",
  addDOC = FALSE,
  smoothWindow = 6,
  plot = FALSE,
  force = T
)

names(p1) <- gsub("_.*","", names(p1))


png(filename = "/work/foxlab/xiaoke/scATAC/Figures/FOSL.png", width = 2000, height = 2000, units = 'px',
    res = 600)
cowplot::plot_grid(p1[[1]])
dev.off()


#####################Co-accessibility analyses--------------
#Store co-accessibility information in the ArchR project
proj_4 <- addCoAccessibility(
  proj_4,
  reducedDims = "Harmony_combined"
)
proj_4 <- addPeak2GeneLinks(
  proj_4,
  reducedDims = "Harmony_combined",
  useMatrix = "GeneExpressionMatrix"
)
saveArchRProject(ArchRProj = proj_4, outputDirectory = "Myeloid_cells", load = FALSE)

#Retrieve peak co-accessibility information
CA <- getCoAccessibility(
  proj_4,
  corCutOff = 0.5,
  resolution = 1000,
  returnLoops = TRUE
)


p <- plotBrowserTrack(
  ArchRProj = proj_4, 
  groupBy = "Label", 
  geneSymbol = "CD4", 
  upstream = 50000,
  downstream = 50000,
  loops = getCoAccessibility(proj_4)
)

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/coacess_CD4.png", width = 2000, height = 2000, units = 'px',
    res = 300)
grid::grid.draw(p$`CD4`)
dev.off()

# Retrieve peak-to-gene accessibility
p2g <- getPeak2GeneLinks(
  proj_4,
  corCutOff = 0.45,
  resolution = 1000,
  returnLoops = TRUE
)
p <- plotBrowserTrack(
  ArchRProj = proj_4, 
  groupBy = "Label", 
  geneSymbol = "SIV", 
  upstream = 500,
  downstream = 10000,
  loops = getPeak2GeneLinks(proj_4)
)
png(filename = "/work/foxlab/xiaoke/scATAC/Figures/coacessG_SIV.png", width = 2000, height = 2000, units = 'px',
    res = 300)
grid::grid.draw(p$`SIV`)
dev.off()


# Plot a heatmap for peak-to-gene links
p <- plotPeak2GeneHeatmap(
  proj_4,
  groupBy = "Label"
)
png(filename = "/work/foxlab/xiaoke/scATAC/Figures/pk2gen_heatmap.png", width = 4000, height = 6000, units = 'px',
    res = 300)
p
dev.off()

############ Identify positive regulators-----------
seGroupMotif <- getGroupSE(proj_4, useMatrix = "Motif_Matrix",
                           groupBy = "Label")
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames == "z"]

rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs
# Genescore matrix correlates with motif matrix
corGSM_MM <- correlateMatrices(
  proj_4,
  useMatrix1 = "GeneScoreMatrix",
  useMatrix2 = "Motif_vierMatrix",
  reducedDims = "Harmony_combined"
)

# Gene expression matrix correlates with motif matrix
corGIM_MM <- correlateMatrices(
  proj_4,
  useMatrix1 = "GeneExpressionMatrix",
  useMatrix2 = "MotifMatrix",
  reducedDims = "Harmony_combined"
)

# Add Maxium Delta Deviation to the correlation dataframe
corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]

# Identify positive TF regulators (motif ~ chromatin accessibility)
corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
corGSM_MM$TFRegulator <- "NO"
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"

df1 <- data.frame(corGSM_MM)

library(ggrepel)
options(ggrepel.max.overlaps = Inf)
p1 <- ggplot(df1, aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  guides(color = guide_legend(override.aes = list(shape = 16)))+
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGSM_MM$maxDelta)*1.05)
  )
p <- p1 +  geom_label_repel(aes(label = ifelse(TFRegulator == "YES", GeneScoreMatrix_name,'')), show.legend = F)

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/pos_TF_regulator.png", width = 2000, height = 2000, units = 'px',
    res = 300)
p
dev.off()

# Identify positive TF regulators (motif ~ RNA expression)
corGIM_MM <- corGIM_MM[order(abs(corGIM_MM$cor), decreasing = TRUE), ]
corGIM_MM <- corGIM_MM[which(!duplicated(gsub("\\-.*","",corGIM_MM[,"MotifMatrix_name"]))), ]
corGIM_MM$TFRegulator <- "NO"
corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.5 & corGIM_MM$padj < 0.01 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.75))] <- "YES"

df2 <- data.frame(corGIM_MM)

p1 <- ggplot(df2, aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Expression") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGIM_MM$maxDelta)*1.05)
  )
p <- p1 +  geom_label_repel(aes(label = ifelse(TFRegulator == "YES", GeneExpressionMatrix_name,'')), show.legend = F)


png(filename = "/work/foxlab/xiaoke/scATAC/Figures/pos_TF_regulator_1.png", width = 2000, height = 2000, units = 'px',
    res = 300)
p
dev.off()

# The positive TFs found in both correlations
library(ggvenn)
pos_GSM = sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])
pos_GIM = sort(corGIM_MM[corGIM_MM$TFRegulator=="YES",1])
lis <- list(pos_GSM=pos_GSM, pos_GIM=pos_GIM)

plot1 <- ggvenn(lis,show_elements = T, label_sep = "\n", 
                fill_color = c("lightblue","lightyellow","pink"), text_size = 4,
                stroke_size =2)

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/Venn.png", width = 3000, height = 2000, units = 'px',
    res = 300)
plot1
dev.off()


