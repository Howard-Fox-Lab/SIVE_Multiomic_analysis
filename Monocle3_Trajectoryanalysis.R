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
unique(proj_4$Label)


# Trajectory Analyses
Trajectory <- c("Micro_like","CAM_like_1","CAM_like_3","CAM_like_2")
proj_4 <- addTrajectory(
  proj_4,
  name = "InfectedU",
  groupBy = "Label",
  trajectory = Trajectory,
  embedding = "UMAP_combined",
  force = TRUE
)

####-----Unsupervised Trajactory analyses-------------------
# change conda env to monocle3
library(monocle3)
library(SingleCellExperiment)
sce <- SingleCellExperiment(
  assays = SimpleList(
    counts = as(matrix(rnorm(nCells(proj_4) * 3), ncol = nCells(proj_4), nrow = 3), "dgCMatrix")
  ),
  colData = getCellColData(proj_4)
)

cds <- methods::new(
  "cell_data_set", 
  assays = SummarizedExperiment::Assays(list(counts = methods::as(assay(sce), "dgCMatrix"))), 
  colData = colData(sce), 
  int_elementMetadata = int_elementMetadata(sce), 
  int_colData = int_colData(sce), 
  int_metadata = int_metadata(sce), 
  metadata = metadata(sce), 
  NAMES = NULL, 
  elementMetadata = elementMetadata(sce)[, 0], 
  rowRanges = rowRanges(sce)
)
metadata(cds)$cds_version <- Biobase::package.version("monocle3")
reducedDims(cds)$UMAP <- getEmbedding(proj_4, embedding = "UMAP_combined")

# Run clustering
clusterParams <- NULL
clusterParams$cds <- cds
cds <- do.call(monocle3::cluster_cells, clusterParams)
rm(clusterParams)

# Learning Graph
graphParams <- NULL
graphParams$cds <- cds
cds <- do.call(monocle3::learn_graph, graphParams)
rm(graphParams)

# Get Principal Node
pCells <- which(colData(cds)[, "Label"] == "CAM_2")
closestVertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closestVertex <- as.matrix(closestVertex[colnames(cds), ])
rootNodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closestVertex[pCells,]))))]

# order cells
cds <- order_cells(cds, root_pr_nodes = rootNodes)

# Get pseudotime
cds@principal_graph_aux[[1]]$pseudotime <- ArchR:::.getQuantiles(cds@principal_graph_aux[[1]]$pseudotime) * 100


proj_4 <- addMonocleTrajectory(
  proj_4,
  name = "Mono_CAM2",
  groupBy = "Label",
  useGroups = c("Micro_like","CAM_like_1","CAM_like_2","CAM_like_3","Micro_1","Micro_2",
                "Micro_3","Micro_4","CAM_1","CAM_2"),
  monocleCDS = cds
)

saveArchRProject(ArchRProj = proj_4, outputDirectory = "Myeloid_cells", load = FALSE)

p1 <- plotTrajectory(proj_4, trajectory = "Mono_CAM2", colorBy = "cellColData",embedding = "UMAP_combined",name = "Mono_CAM2")
p2 <- plotTrajectory(proj_4, trajectory = "Mono_CAM2",colorBy = "GeneScoreMatrix",embedding = "UMAP_combined",name = "CEBPB", continuousSet = "horizonExtra")
p3 <- plotTrajectory(proj_4, trajectory = "Mono_CAM2",colorBy = "GeneExpressionMatrix",embedding = "UMAP_combined",name = "CEBPB", continuousSet = "blueYellow")

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/Traj_all_CAM2.png", width = 6000, height = 2000, units = 'px',
    res = 300)
ggAlignPlots(p1[[1]],p2[[1]],p3[[1]], type = "h")
dev.off()

#---------------------------------------------------
p1 <- plotTrajectory(proj_4, trajectory = "InfectedU", colorBy = "cellColData",embedding = "UMAP_combined",name = "InfectedU")
p2 <- plotTrajectory(proj_4, trajectory = "InfectedU",colorBy = "GeneScoreMatrix",embedding = "UMAP_combined",name = "CEBPB", continuousSet = "horizonExtra")
p3 <- plotTrajectory(proj_4, trajectory = "InfectedU",colorBy = "GeneExpressionMatrix",embedding = "UMAP_combined",name = "CEBPB", continuousSet = "blueYellow")

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/Traj_Infected.png", width = 6000, height = 2000, units = 'px',
    res = 300)
ggAlignPlots(p1[[1]],p2[[1]],p3[[1]], type = "h")
dev.off()

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/Traj_Infected_1.png", width = 6000, height = 2000, units = 'px',
    res = 300)
ggAlignPlots(p1[[2]],p2[[2]],p3[[2]], type = "h")
dev.off()


# Pseudotime Heatmap for Transcrition factor
trajMM  <- getTrajectory(ArchRProj = proj_4, name = "InfectedU", useMatrix = "MotifMatrix", log2Norm = FALSE)
trajpk <- getTrajectory(ArchRProj = proj_4, name = "InfectedU", useMatrix = "PeakMatrix", log2Norm = FALSE)
trajGSM <- getTrajectory(ArchRProj = proj_4, name = "InfectedU", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
trajGIM <- getTrajectory(ArchRProj = proj_4, name = "InfectedU", useMatrix = "GeneExpressionMatrix", log2Norm = FALSE)

p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "blueYellow"))
p2 <- plotTrajectoryHeatmap(trajpk, pal = paletteContinuous(set = "solarExtra"))
p3 <-  plotTrajectoryHeatmap(trajGSM, pal = paletteContinuous(set = "horizonExtra"))
p4 <-  plotTrajectoryHeatmap(trajGIM, pal = paletteContinuous(set = "horizonExtra"))

P1 <- ComplexHeatmap::draw(p1, heatmap_legend_side = "bot", annotation_legend_side = "bot")
P2 <- ComplexHeatmap::draw(p2, heatmap_legend_side = "bot", annotation_legend_side = "bot")
P3 <- ComplexHeatmap::draw(p3, heatmap_legend_side = "bot", annotation_legend_side = "bot")
P4 <- ComplexHeatmap::draw(p4, heatmap_legend_side = "bot", annotation_legend_side = "bot")

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/Traj_heatmap_geneexpression.png", width = 2500, height = 3000, units = 'px',
    res = 300)
P4
dev.off()


# Integrative pseudo-time analyses Identify positive TF regulators (GSM)
traj_corGSM_MM <- correlateTrajectories(trajGSM,trajMM)
trajGSM2 <- trajGSM[traj_corGSM_MM[[1]]$name1, ]
trajMM2 <- trajMM[traj_corGSM_MM[[1]]$name2, ]
trajCombined <- trajGSM2

assay(trajCombined, withDimnames =F) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))
ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)

P1 <- ComplexHeatmap::draw(ht1, heatmap_legend_side = "bot", annotation_legend_side = "bot")
P2 <- ComplexHeatmap::draw(ht2, heatmap_legend_side = "bot", annotation_legend_side = "bot")


png(filename = "/work/foxlab/xiaoke/scATAC/Figures/GSM_MM_posTF_M.png", width = 2000, height = 3000, units = 'px',
    res = 300)
P2
dev.off()

# Integrative pseudo-time analyses Identify positive TF regulators (GIM)
traj_corGIM_MM <- correlateTrajectories(trajGIM, trajMM)

trajGIM2 <- trajGIM[traj_corGIM_MM[[1]]$name1, ]
trajMM2 <- trajMM[traj_corGIM_MM[[1]]$name2, ]
trajCombined <- trajGIM2

assay(trajCombined, withDimnames =F) <- t(apply(assay(trajGIM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)

rowOrder <- match(rownames(combinedMat), rownames(trajGIM2))
ht1 <- plotTrajectoryHeatmap(trajGIM2,  pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)

P1 <- ComplexHeatmap::draw(ht1, heatmap_legend_side = "bot", annotation_legend_side = "bot")
P2 <- ComplexHeatmap::draw(ht2, heatmap_legend_side = "bot", annotation_legend_side = "bot")

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/GIM_MM_posTF_M.png", width = 2000, height = 3000, units = 'px',
    res = 300)
P2
dev.off()

# VennPlot
library(ggvenn)
pos_Traj_GSM <- rownames(trajGSM2)
pos_Traj_GIM <- rownames(trajGIM2)

pos_Traj_GSM <- gsub(".*:","",pos_Traj_GSM)
pos_Traj_GIM <- gsub(".*:","",pos_Traj_GIM)

pos_GSM <- sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])
pos_GIM <- sort(corGIM_MM[corGIM_MM$TFRegulator=="YES",1])



lis <- list(pos_Traj_GSM = pos_Traj_GSM, pos_Traj_GIM = pos_Traj_GIM)

plot1 <- ggvenn(lis,show_elements = T, label_sep = "\n", 
                fill_color = c("lightblue","lightyellow"), text_size = 6,
                stroke_size =2)

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/Venn.png", width = 4000, height = 4000, units = 'px',
    res = 300)
plot1
dev.off()


