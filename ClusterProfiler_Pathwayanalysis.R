#-----------------------Gene ontology (GO) analysis------------------------------
library(clusterProfiler)
library("org.Mmu.eg.db")
library(enrichplot)

# Pathway analyses for RNA data
# Read in csv
data <- read.csv("/work/foxlab/xiaoke/scATAC/excel/markers_RNA.csv")
CAM_like <- data[which(data$group_name %in% "CAM_like"),c(8,10)]
Micro_like <- data[which(data$group_name %in% "Micro_like"),c(8,10)]


# Creat a list
sample <- list("Micro_like" = Micro_like$name,
               "CAM_like" = CAM_like$name)

#Change the format of the gene
sample$`Micro_like` <- bitr(sample$`Micro_like`, fromType = "SYMBOL", toType = 'ENTREZID', OrgDb = "org.Mmu.eg.db")
sample$`CAM_like` <- bitr(sample$`CAM_like`, fromType = "SYMBOL", toType = 'ENTREZID', OrgDb = "org.Mmu.eg.db")

genelist <- list('CAM_like'= sample$`CAM_like`$ENTREZID,
                 'Micro_like'= sample$`Micro_like`$ENTREZID)

genelist$Micro_like <- sort(genelist$Micro_like, decreasing = T)
genelist$CAM_like <- sort(genelist$CAM_like, decreasing = T)

Gocluster <- compareCluster(geneCluster = genelist, fun = 'enrichGO', OrgDb = "org.Mmu.eg.db",ont = "BP")
edox <- setReadable(Gocluster,"org.Mmu.eg.db", 'ENTREZID' )

KEGGcluster <- compareCluster(geneCluster = genelist, fun = 'enrichKEGG', organism = 'mcc')

plot2 <- dotplot(Gocluster, showCategory=10)+
  theme(axis.title = element_blank())

plot1 <- dotplot(KEGGcluster, showCategory=10)+
  theme(axis.title = element_blank())

plot3 <- cnetplot(edox, cex_category = 1.5)

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/GO_compare.png", width = 2500, height = 3500, units = 'px',
    res = 300)
plot2
dev.off()

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/KEGG_compare.png", width = 2500, height = 3500, units = 'px',
    res = 300)
plot1
dev.off()


# Pathway analyses for atac data
# Read in csv
data <- read.csv("/work/foxlab/xiaoke/scATAC/excel/markers_atac.csv")
CAM_like <- data[which(data$group_name %in% "CAM_like"),c(8,10)]
Micro_like <- data[which(data$group_name %in% "Micro_like"),c(8,10)]


# Creat a list
sample <- list("Micro_like" = Micro_like$name,
               "CAM_like" = CAM_like$name)

#Change the format of the gene
sample$`Micro_like` <- bitr(sample$`Micro_like`, fromType = "SYMBOL", toType = 'ENTREZID', OrgDb = "org.Mmu.eg.db")
sample$`CAM_like` <- bitr(sample$`CAM_like`, fromType = "SYMBOL", toType = 'ENTREZID', OrgDb = "org.Mmu.eg.db")

genelist <- list('CAM_like'= sample$`CAM_like`$ENTREZID,
                 'Micro_like'= sample$`Micro_like`$ENTREZID)

genelist$Micro_like <- sort(genelist$Micro_like, decreasing = T)
genelist$CAM_like <- sort(genelist$CAM_like, decreasing = T)

Gocluster <- compareCluster(geneCluster = genelist, fun = 'enrichGO', OrgDb = "org.Mmu.eg.db",ont = "BP")

KEGGcluster <- compareCluster(geneCluster = genelist, fun = 'enrichKEGG', organism = 'mcc')

plot2 <- dotplot(Gocluster, showCategory=10)+
  theme(axis.title = element_blank())

plot1 <- dotplot(KEGGcluster, showCategory=10)+
  theme(axis.title = element_blank())


png(filename = "/work/foxlab/xiaoke/scATAC/Figures/GO_compare.png", width = 2000, height = 4000, units = 'px',
    res = 300)
plot2
dev.off()

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/KEGG_compare.png", width = 2000, height = 4000, units = 'px',
    res = 300)
plot1
dev.off()

