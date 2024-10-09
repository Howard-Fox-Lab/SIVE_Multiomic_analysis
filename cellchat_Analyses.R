proj_4 <- readRDS("/work/foxlab/xiaoke/scATAC/atac_enceph_new/Myeloid_cells/Save-ArchR-Project.rds")
# Get gene expression matrix
RNAMtrix <- getMatrixFromProject(proj_4, useMatrix = "GeneExpressionMatrix")
saveRDS(RNAMtrix, file = "/work/foxlab/xiaoke/scATAC/atac_enceph_new/RNAMtrix.rds")
# Change conda environment
# conda activate cellchat
# Load packages for cellchat
library(CellChat)
library(patchwork)
library(ggplot2)
library(SummarizedExperiment)
library(Seurat)
library(SeuratDisk)
options(stringsAsFactors = FALSE)

# Read in RNAMatrix
merged_seurat <- LoadH5Seurat("/work/foxlab/xiaoke/scATAC/atac_enceph_new/merged_seurat.h5seurat")
RNAMtrix <- readRDS("/work/foxlab/xiaoke/scATAC/atac_enceph_new/RNAMtrix.rds")
meta <- as.data.frame(colData(RNAMtrix))
archR_cells <- gsub("#","_",rownames(meta))
merged_seurat$Label <- meta[match(rownames(merged_seurat@meta.data),archR_cells),"Label"]
Idents(merged_seurat) <- merged_seurat$Label
merged_seurat <- merged_seurat[,!is.na(merged_seurat$Label)]
SaveH5Seurat(merged_seurat, file = "/work/foxlab/xiaoke/scATAC/atac_enceph_new/merged_seurat.h5seurat")

cells <- subset(merged_seurat, idents = c("CAM_like_1","CAM_like_2","CAM_like_3","Micro_like"))
cells <- NormalizeData(cells)

cells <- subset(merged_seurat, idents = c("Micro_1","Micro_2","Micro_3","Micro_4","CAM_1","CAM_2"))
cells <- NormalizeData(cells)

# Get cellchat inputs from seuratobject
data.input <- cells[["RNA"]]@data
labels <- Idents(cells)
meta <- data.frame(labels = labels, row.names = names(labels))

# Create cellchat object
cellchat <- createCellChat(object = cells, group.by = "ident", assay = "RNA")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels")

CellChatDB <- CellChatDB.human #Get database

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB)
cellchat@DB <- CellChatDB.use

# Subset cellchat
cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Compute the communication probability
cellchat <- computeCommunProb(cellchat, type = "triMean")

# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, file = "/work/foxlab/xiaoke/scATAC/atac_enceph_new/cellchat/cellchat_SIVE.rds")
saveRDS(cellchat, file = "/work/foxlab/xiaoke/scATAC/atac_enceph_new/cellchat/cellchat_Uninfect.rds")

# Visualization of aggregated cell-cell communication network
cellchat <- readRDS("/work/foxlab/xiaoke/scATAC/atac_enceph_new/cellchat/cellchat_SIVE.rds")
cellchat <- readRDS("/work/foxlab/xiaoke/scATAC/atac_enceph_new/cellchat/cellchat_Uninfect.rds")

groupSize <- as.numeric(table(cellchat@idents))
p1 <- netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                       label.edge= F, title.name = "Number of interactions")
p2 <- netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                       label.edge= F, title.name = "Interaction weights/strength")

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/cellchat_weight.png", width = 2000, height = 2000, units = 'px',
    res = 300)
p2
dev.off()


# access all signaling pathway
cellchat@netP$pathway
pathways.show <- c("CCL")
# Chord Plot
group.cellType <- c(rep("CAM_like",3),"Micro_like")
group.cellType <- c("Micro","CAM",rep("Micro",3),"CAM")
names(group.cellType) <- levels(cellchat@idents)

p1 <- netVisual_chord_cell(cellchat, signaling = pathways.show, lab.cex =1, 
                           group = group.cellType, title.name = "")
png(filename = "/work/foxlab/xiaoke/scATAC/Figures/CCL_signaling.png", width = 2500, height = 2500, units = 'px',
    res = 600)
p1
dev.off()

# The contribution of ligand-receptor pair
pariLR.CCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pariLR.CCL[1,]
p1 <- netVisual_individual(cellchat,group = group.cellType ,signaling = pathways.show, pairLR.use = LR.show, 
                           layout = "chord", big.gap = 5)
png(filename = "/work/foxlab/xiaoke/scATAC/Figures/CCL3_CCR1.png", width = 2000, height = 2000, units = 'px',
    res = 600)
p1
dev.off()

# All the significant interactions from micro_like/CAM-like/CAM/Microglia
levels(cellchat@idents)
p1 <- netVisual_chord_gene(cellchat,sources.use = 10, targets.use = c(7:9), 
                           lab.cex = 0.8,show.legend = F,small.gap = 3, annotationTrackHeight = c(0.05),
                           signaling = c("APP","SPP1","FN1","COMPLEMENT","CD45","ApoE","CCL","PTPRM"))


p2 <- netVisual_chord_gene(cellchat, sources.use = c(7:9), targets.use = 10, 
                           lab.cex = 0.8,show.legend = F,small.gap = 3, annotationTrackHeight = c(0.05),
                           signaling = c("APP","SPP1","FN1","COMPLEMENT","CD45","ApoE","CCL","PTPRM"))

p3 <- netVisual_chord_gene(cellchat, sources.use = c(1,3,4,5), targets.use = c(2,6), 
                           lab.cex = 0.8,show.legend = F,small.gap = 3, annotationTrackHeight = c(0.05),
                           signaling = c("APP","SPP1","FN1","COMPLEMENT","CD45","ApoE","CCL","PTPRM"))

p4 <- netVisual_chord_gene(cellchat, sources.use = c(2,6), targets.use = c(1,3,4,5), 
                           lab.cex = 0.8,show.legend = F,small.gap = 3, annotationTrackHeight = c(0.05),
                           signaling = c("APP","SPP1","FN1","COMPLEMENT","CD45","ApoE","CCL","PTPRM"))

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/MicroLike.png", width = 4000, height = 4000, units = 'px',
    res = 600)
p1
dev.off()
png(filename = "/work/foxlab/xiaoke/scATAC/Figures/CAMLike.png", width = 4000, height = 4000, units = 'px',
    res = 600)
p2
dev.off()
png(filename = "/work/foxlab/xiaoke/scATAC/Figures/Microglia.png", width = 4000, height = 4000, units = 'px',
    res = 600)
p3
dev.off()
png(filename = "/work/foxlab/xiaoke/scATAC/Figures/CAM.png", width = 4000, height = 4000, units = 'px',
    res = 600)
p4
dev.off()

# Identify dominant senders and receivers
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
saveRDS(cellchat, file = "/work/foxlab/xiaoke/scATAC/atac_enceph_new/cellchat/cellchat_SIVE.rds")
saveRDS(cellchat, file = "/work/foxlab/xiaoke/scATAC/atac_enceph_new/cellchat/cellchat_Uninfect.rds")

gg1 <- netAnalysis_signalingRole_scatter(cellchat, label.size = 5)

gg2 <- gg1+xlab("ability of sending signals")+ylab("ability of responding signals")+
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 15))

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/main_contributor.png", width = 4000, height = 3000, units = 'px',
    res = 600)
gg2
dev.off()

# Identify signals contributing the most to outgoing 
# or incoming signaling of certain cell groups
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", height =10, width = 6, font.size = 10,font.size.title = 10)

ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", height =10, width = 6, font.size =10,font.size.title = 10)
png(filename = "/work/foxlab/xiaoke/scATAC/Figures/main_interactions.png", width = 4500, height = 4000, units = 'px',
    res = 600)
ht1+ht2
dev.off()


