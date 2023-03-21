library(ggplot2)
library(scales)
library(viridis)
library(reshape)
library(RColorBrewer)
library(pheatmap)
library(Rtsne)
library(dplyr)
library(tidyr)
library(cluster)
library(ggfortify)
library(plotly)

# Setting the working directory to the directory which contains this script
if (exists("RStudio.Version")){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(getSrcDirectory()[1])
}
#Theme for plotting
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(),    #strip major gridlines
    panel.grid.minor = element_blank(),    #strip minor gridlines
    plot.title=element_text(size=26, face="bold"),
    legend.title = element_text(size=24,face="bold"),
    legend.text = element_text(size=18)
  )
dir.create('../results/plots/chassis_strain')
source('plotVennDiagram.R')
#Get heatmap for targets matrix
allTargetsMat <- read.csv('../results/production_targets/targetsMatrix_compatible.txt',sep='\t',stringsAsFactors = TRUE)
allTargetsMat <- read.csv('../results/production_targets/targetsMatrix_mech_validated.txt',sep='\t',stringsAsFactors = TRUE)
#Discard down-regulatio

#Top modifications
OEs<- which(allTargetsMat$shortNames=='TPI1')
KDs<- which(allTargetsMat$shortNames=='URA7')
KOs<- which(allTargetsMat$shortNames=='SER3')
df.genes<- allTargetsMat[c(OEs,KDs,KOs),]

OEs <- which(colSums(allTargetsMat[OEs,5:ncol(allTargetsMat)])==3)
KDs <- which(colSums(allTargetsMat[KDs,5:ncol(allTargetsMat)])==2)
KOs <- which(colSums(allTargetsMat[KOs,5:ncol(allTargetsMat)])==1)
OEprod <- colnames(df.genes)[OEs]
KDprod <- colnames(df.genes)[KDs]
KOprod <- colnames(df.genes)[KOs]
overlap  <- intersect(intersect(OEs,KDs),KOs)
prods_1 <- colnames(df.genes)[overlap]

lista <- list(OEprod,KDprod,KOprod)
colores <- cividis(11)
pdf(paste('../results/plots/chassis_strain/topMod1_strain.pdf',sep=''),width = 6.5, height = 6.5)
overlap <- plotVennDiagram(lista,c('TPI1','URA7','SER3'),c(colores[11],colores[6],colores[1]),rep(2,7),3)
dev.off()

#Top modifications 1 and 2
OEs<- which(allTargetsMat$shortNames=='TPI1' | allTargetsMat$shortNames=='PRO1')
KDs<- which(allTargetsMat$shortNames=='URA7'| allTargetsMat$shortNames=='ODC2')
KOs<- which(allTargetsMat$shortNames=='SER3'| allTargetsMat$shortNames=='GLK1')
df.genes<- allTargetsMat[c(OEs,KDs,KOs),]

OEs <- which(colSums(allTargetsMat[OEs,5:ncol(allTargetsMat)])==6)
KDs <- which(colSums(allTargetsMat[KDs,5:ncol(allTargetsMat)])==4)
KOs<- which(colSums(allTargetsMat[KOs,5:ncol(allTargetsMat)])==2)

overlap  <- intersect(intersect(OEs,KDs),KOs)
prods_2 <- colnames(df.genes)[overlap]
OEprod <- colnames(df.genes)[OEs]
KDprod <- colnames(df.genes)[KDs]
KOprod <- colnames(df.genes)[KOs]

lista <- list(OEprod,KDprod,KOprod)
colores <- cividis(11)
pdf(paste('../results/plots/chassis_strain/topMod2_strain.pdf',sep=''),width = 6.5, height = 6.5)
overlap <- plotVennDiagram(lista,c('TPI1 & PRO1','URA7 & ODC2','SER3 & GLK1'),c(colores[11],colores[6],colores[1]),rep(2,7),3)
dev.off()