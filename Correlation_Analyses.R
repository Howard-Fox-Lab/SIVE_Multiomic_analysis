library(ArchR)
library(ggplot2)
proj_4 <- readRDS("/work/foxlab/xiaoke/scATAC/atac_enceph_new/Myeloid_cells/Save-ArchR-Project.rds")
#Get motif matrix
motifmtrx <- getMatrixFromProject(proj_4, useMatrix = "Motif_vierMatrix")
motifmtrx <- getMatrixFromProject(proj_4, useMatrix = "MotifMatrix")
RNAmtrx <- getMatrixFromProject(proj_4, useMatrix = "GeneExpressionMatrix")
Genemtrx <- getMatrixFromProject(proj_4, useMatrix = "GeneScoreMatrix")

rownames(assay(RNAmtrx, withDimnames=F)) <- rowData(RNAmtrx)$name
rownames(assay(Genemtrx, withDimnames=F)) <- rowData(Genemtrx)$name

mtfname <- rowData(motifmtrx)
mtfname <- as.data.frame(mtfname)
#For vierstra
mtfname <- mtfname %>%
  separate(col=name, into=c("B","C","D"), sep="\\s*\\|\\s*", fill="right")
rownames(assay(motifmtrx, withDimnames = F)) <- mtfname$C
#For JASPAR
mtfname$name <- sub("_[^_]+$", "",mtfname$name)
rownames(assay(motifmtrx, withDimnames = F)) <- mtfname$name

df_RNA <- data.frame(IL1b = assay(RNAmtrx, withDimnames = F)['IL1B',],
                     CXCL1 = assay(RNAmtrx, withDimnames = F)['CXCL1',],
                     IL10 = assay(RNAmtrx, withDimnames = F)['IL10',])

df_motif <- data.frame(NFKB = assay(motifmtrx, withDimnames = F)['NFKB/RELA',],
                       STAT = assay(motifmtrx, withDimnames = F)['STAT1..STAT2',])

df_Gene <- data.frame(IFNA1 = assay(Genemtrx, withDimnames = F)['IFNA13',],
                      IFNA8 = assay(Genemtrx, withDimnames = F)['IFNA8',],
                      IFNE = assay(Genemtrx, withDimnames = F)['IFNE',])

df1 <- cbind(df_RNA, df_motif)
df2 <- cbind(df_Gene, df_motif)
df3 <- as.data.frame(proj_4@cellColData)
df1$cluster <- df3[match(rownames(df1), rownames(df3)),"Label_2"]
df2$cluster <- df3[match(rownames(df2), rownames(df3)),"Label_2"]

df1 %>%
  group_by(cluster)%>%
  summarise(IL1b = mean(IL1b),
            CXCL1 = mean(CXCL1),
            IL10 = mean(IL10),
            NFKB = mean(NFKB),
            STAT = mean(STAT)) -> mean_value

df2 %>%
  group_by(cluster)%>%
  summarise(IFNA1 = mean(IFNA1),
            IFNA8 = mean(IFNA8),
            IFNE = mean(IFNE),
            STAT = mean(STAT),
            NFKB = mean(NFKB)) -> mean_value_1

round(cor(mean_value$IL1b,mean_value$NFKB), 2)
round(cor(mean_value$IL10,mean_value$NFKB), 2)
round(cor(mean_value$CXCL1,mean_value$NFKB), 2)
round(cor(mean_value_1$IFNA1,mean_value_1$STAT), 2)
round(cor(mean_value_1$IFNA8,mean_value_1$STAT), 2)
round(cor(mean_value_1$IFNE,mean_value_1$STAT), 2)

cor.test(mean_value$IL1b,mean_value$NFKB)
cor.test(mean_value$IL10,mean_value$NFKB)
cor.test(mean_value$CXCL1,mean_value$NFKB)
cor.test(mean_value_1$IFNA1,mean_value_1$STAT)
cor.test(mean_value_1$IFNA8,mean_value_1$STAT)
cor.test(mean_value_1$IFNE,mean_value_1$STAT)


plot1 <- ggplot(mean_value_1, #mean_value
                aes(y=STAT, x=IFNE, color = cluster)) + #y = NFKB x = IFNA1/IFNA8/IFNE/IL1b/IL10/CXCL1
  geom_point(size = 5) + 
  geom_smooth(method="lm", col="black") + 
  theme_ArchR()+
  theme(axis.title = element_text(size = 15),legend.text = element_text(size = 15),
        legend.title = element_blank(),plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"))+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  xlab("Average predicted gene activity of IFNE")+
  ylab("Average Deviation Z-score of STAT1/STAT2")

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/cor_ifne.png", width = 2000, height = 2000, units = 'px',
    res = 300)
plot1
dev.off()


######--------Correlation between cell cluster (RNA Expression)-----------
id1 <- BiocGenerics::which(proj_4$Label_2 %in% "CAM") #Microglia
id2 <- BiocGenerics::which(proj_4$Label_2 %in% "CAM_like") #Microglia-like

Uninf <- proj_4$cellNames[id1]
Infec <- proj_4$cellNames[id2]

RNA_Uninf <- getMatrixFromProject(proj_4[Uninf,], useMatrix = "GeneExpressionMatrix")
RNA_Infec <- getMatrixFromProject(proj_4[Infec,], useMatrix = "GeneExpressionMatrix")
RNA_all <- getMatrixFromProject(proj_4, useMatrix = "GeneExpressionMatrix")

rownames(assay(RNA_Uninf, withDimnames=F)) <- rowData(RNA_Uninf)$name
rownames(assay(RNA_Infec, withDimnames=F)) <- rowData(RNA_Infec)$name
rownames(assay(RNA_all, withDimnames=F)) <- rowData(RNA_all)$name

df_Uninf <- assay(RNA_Uninf, withDimnames = F)
df_Infec <- assay(RNA_Infec, withDimnames = F)
df_RNAall <- assay(RNA_all, withDimnames = F)

RNAall <- data.frame(RNA = apply(df_RNAall,1,mean))
RNAall <- RNAall[order(RNAall$RNA, decreasing = T),,drop = F]
gene <- rownames(RNAall)[1:4000]

Uninf <- data.frame(RNA = apply(df_Uninf,1,mean))
Uninf <- Uninf[order(Uninf$RNA, decreasing = T),,drop = F]
Uninf_sub <- Uninf[which(rownames(Uninf) %in% gene),,drop = F]

Infec <- data.frame(RNA = apply(df_Infec,1,mean))
Infec <- Infec[order(Infec$RNA, decreasing = T),,drop = F]
Infec_sub <- Infec[which(rownames(Infec) %in% gene),,drop = F]
Infec_sub <- Infec_sub[match(rownames(Uninf_sub), rownames(Infec_sub)),,drop = F]

df <- data.frame(Uninfect = Uninf_sub$RNA, Infect = Infec_sub$RNA, gene = rownames(Uninf_sub))

round(cor(df$Uninfect,df$Infect), 2)
round(coef(lm(df$Uninfect ~ df$Infect))[1],2) #intercept
round(coef(lm(df$Uninfect ~ df$Infect))[2],2) #slope

library(ggforce)
plot1 <- ggplot(df, aes(y=Uninfect, x=Infect)) + 
  geom_point(size = 3) +
  geom_smooth(method="lm", col="black") + 
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 2) +
  theme_ArchR()+
  theme(axis.title = element_text(size = 15),legend.text = element_text(size = 15),
        legend.title = element_blank(),plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"))+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  xlab("Average Expression in CAM cluters")+
  ylab("Average Expression CAM-like clusters")+
  facet_zoom(xlim = c(0,10), ylim = c(0,10))

png(filename = "/work/foxlab/xiaoke/scATAC/Figures/cor_CAMvsCAMLike.png", width = 2000, height = 2000, units = 'px',
    res = 300)
plot1
dev.off()
