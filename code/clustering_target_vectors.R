library(ggplot2)
library(scales)
library(viridis)
library(reshape)
library(RColorBrewer)
library(pheatmap)
library(Rtsne)
library(tibble)
library(dplyr)
library(tidyr)
library(cluster)
library(ggfortify)
library(plotly)
library(matrixStats)
library(htmlwidgets)
#library(qgraph)
library(reshape2)
library(umap)

if (!require("processx")) install.packages("processx")


# Setting the working directory to the directory which contains this script
if (exists("RStudio.Version")){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(getSrcDirectory()[1])
}
families <- c('amino acid','alkaloid','organic acid','protein','alcohol','terpene','fatty acids and lipids','flavonoid','aromatic','bioamine','stilbenoids')
codes    <- c('_AA','_alk','_oAc','_pro','_alc','_ter','_FAL','_fla','_aro','_bio','_stb')
#Load targets summary
filename        <- paste('../results/production_targets/targets_summary.txt',sep='')
targets_summary <- read.csv(filename,sep='\t',stringsAsFactors = FALSE)
classes <- unique(targets_summary$chemClass)
counts  <- c()
for (i in 1:length(classes))
{
  counts <- c(counts,sum(targets_summary$chemClass==classes[i]))
}
#Get data frame for plotting
df <- data.frame(classes,counts)
df <- df[order(-counts),]
df$classes <- factor(df$classes,levels<-df$classes) 
df$perc <- round(df$counts*100/sum(df$counts),digits=1)
#Get pallete of colors
colourCount <- length(unique(df$classes))
getPalette <- colorRampPalette(brewer.pal(colourCount,"Paired"))

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

#Analyse targets matrix
allTargetsMat <- read.csv('../results/production_targets/targetsMatrix_L3.txt',sep='\t',stringsAsFactors = TRUE)
#allTargetsMat <- read.csv('../results/production_targets/targetsMatrix_mech_validated.txt',sep='\t',stringsAsFactors = TRUE)
targetsMat <- allTargetsMat
targetsMat <- targetsMat[rowSums(targetsMat[,5:ncol(targetsMat)])>0,]
cases      <- ncol(targetsMat) - 4

#Dimensionality reduction
newDF <- targetsMat[,5:ncol(targetsMat)]
genes <- targetsMat$shortNames
idxs  <- which(rowSums(targetsMat[,5:ncol(targetsMat)])!=cases)
extra <- colnames(newDF)
rownames(newDF) <- genes
newDF <- as.data.frame(t(newDF))
newDF$extra <- extra
newDF<- newDF %>% separate(extra, c("chemical", "family"), "_fam_")
newDF <- newDF[sort(newDF$chemical),]
##Remove duplicate vectors 
#(sum random small quantities to avoid duplicatio)
M1<-matrix(rnorm(nrow(newDF)*(ncol(newDF)-2),mean = 0.001,sd=0.0001),nrow=nrow(newDF))
newDF[,1:(ncol(newDF)-2)] <- newDF[,1:(ncol(newDF)-2)]  + M1
idxs <- which(!duplicated(newDF[,1:(ncol(newDF)-3)]))
duplicates <- which(duplicated(newDF[,1:(ncol(newDF)-3)]))
duplicates <- newDF[duplicates,]
newDF <- newDF[idxs,]
#load origin data
chem_origin <- read.csv('../results/production_targets/chemicals_origin.txt',sep='\t',stringsAsFactors = FALSE)
chem_origin <- chem_origin[order(chem_origin$b1),]
origin <- chem_origin$b3 == 'native'
origin <- as.numeric(origin)

prodCap <- read.csv('../results/production_capabilities/prodCapabilities_allChemicals.txt',sep='\t',stringsAsFactors = FALSE)
prodCap <- prodCap[order(prodCap$compound),]
idxs <- which(prodCap$cFlux_l>= 0.95)
Plim <- rep(1,nrow(prodCap))
Plim[idxs] <- 0

#Plim <- prodCap$Pburden/max(prodCap$Pburden)
#Plim[Plim!=1] <- 0

origin <- as.numeric(origin)

idxs <- c()
idxM <- c()
for (i in 1:nrow(newDF)){
  chemical <- rownames(newDF)[i]
  chemical <- substr(chemical, 1, (nchar(chemical)-8))
  index <- which(chem_origin$b1==chemical)
  index2 <- which(prodCap$compound==chemical)
  idxs <- c(idxs,index)
  idxM <- c(idxM,index2)
}
origin <- origin[idxs]
Plim <- Plim[idxM]
#Add product family info
#origin<-origin[order(newDF$family)]
#newDF <- newDF[order(newDF$family),]
famLvls <- as.numeric(unique(factor(newDF$family)))
famLvls <- (unique(factor(newDF$family)))
famLvls <- famLvls[order(famLvls)]
newDF$family <- factor(newDF$family,levels = famLvls)
orgStr <- as.character(origin)
orgStr[orgStr=='1']<-'N'
orgStr[orgStr=='0']<-'H'
Plim[Plim==0] <- 8
Plim[Plim==1] <- 12
newDF$Plim <- factor(Plim,levels = unique(Plim))
newDF$origin <- factor(orgStr,levels = unique(orgStr))
#PCA
PCAdata   <- prcomp(newDF[,1:(ncol(newDF)-4)], center = TRUE,scale = FALSE,retx=TRUE)
p         <- autoplot(PCAdata,data = newDF,colour = 'family',size = 4)
plotTitle <- paste('../results/plots/PCA_allTargets.png',sep='')
png(plotTitle,width = 600, height = 600)
plot(p)
dev.off()
#tSNE

for (i in 1:25)
{
set.seed(18) # Set a seed if you want reproducible results
perplxty <- i
tsne_out <- Rtsne(newDF[,1:(ncol(newDF)-4)],dims=2,perplexity=perplxty, max_iter = 10000,theta=0.1) # Run TSNE
tsne_plot <- data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2],newDF$chemical,newDF$family,newDF$origin,newDF$origin)
#Define color pallete
colourCount2 <- length(unique(tsne_plot$family))
getPalette2  <- colorRampPalette(brewer.pal(colourCount, "Paired"))
colnames(tsne_plot)[(ncol(tsne_plot)-2)]<- 'family'
colnames(tsne_plot)[(ncol(tsne_plot)-1)]<- 'Plim'
colnames(tsne_plot)[(ncol(tsne_plot))]<- 'origin'
p <- plot_ly(tsne_plot,x=~x, y=~y, text =~newDF.chemical, type="scatter", mode="markers",symbol=~origin,symbols=c(8,20), color=~family,colors = getPalette2(colourCount),size=Plim)%>% 
layout(title= list(text = paste('tsne_',i,'_perplexity')))
#plotTitle <- paste('../results/plots/cluster_analysis/tSNE_allTargets_',i,'.pdf',sep='')
#orca(p, plotTitle,width=900,height=600)
name <-paste("../results/plots/tsne/tsne_",i,".html",sep = '')
saveWidget(p, name, selfcontained = F, libdir = "lib")
}



