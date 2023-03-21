#stoich_distance analysis
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
library(htmlwidgets)
#library(qgraph)
library(reshape2)
library(readr)

if (exists("RStudio.Version")){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(getSrcDirectory()[1])
}

families <- c('amino acid','alkaloid','organic acid','protein','alcohol','terpene','fatty acids and lipids','flavonoid','aromatic','bioamine')
codes    <- c('_AA','_alk','_oAc','_pro','_alc','_ter','_FA','_fla','_aro','_bioAm')
#Load chemicals info
filename        <- '../data/chemicals_info.txt'
chemicals_info <- read_delim(filename,"\t", escape_double = FALSE, na = "NA",trim_ws = TRUE)
#reformat
chemicals_info$ecModel<- gsub('-','_',chemicals_info$ecModel)
#chemicals_info$ecModel<- gsub('(S)_reticuline','S_reticuline',chemicals_info$ecModel)
chemicals_info$ecModel<- gsub(',','_',chemicals_info$ecModel)
#chemicals_info$ecModel<- gsub('ecGlutamate.mat','ecGlutamate',chemicals_info$ecModel)
#chemicals_info$ecModel<- gsub('.mat','',chemicals_info$ecModel)

#load targets matrix
allTargetsMat <- read.csv('../results/production_targets/targetsMatrix_L3.txt',sep='\t',stringsAsFactors = TRUE)
targetsMat    <- allTargetsMat
cases <- ncol(targetsMat) - 4

#Distance matrix
x <- targetsMat[,5:ncol(targetsMat)]
#x[x==0] <- 1E-12
#x <- log10(x)

colnames(x) <- colnames( targetsMat[,5:ncol(targetsMat)])
x <- t(x)
distMat <- dist(x, method = "euclidean", diag = FALSE, upper = FALSE)
distMat <- as.data.frame(as.matrix(distMat),stringsAsFactors = FALSE)
rownames(distMat) <- rownames(x)
#var1 <- rownames(distMat)
dd <- as.dist(distMat)
hc <- hclust(dd)
distMat <-distMat[hc$order, hc$order]
var1 <- rownames(distMat)
# Melt the distance matrix
melted_distMat <- melt(distMat)
melted_distMat <- cbind(var1,melted_distMat)
#Rearrange for plotting
copyMelted <- melted_distMat
melted_distMat$var1 <- gsub('_fam_','_',melted_distMat$var1)
melted_distMat$variable <- gsub('_fam_','_',melted_distMat$variable)
melted_distMat$var1 <- factor(melted_distMat$var1,levels = unique(melted_distMat$var1))
melted_distMat$variable <- factor(melted_distMat$variable,levels = unique(melted_distMat$variable))
#melted_distMat$value <- melted_distMat$value/max(melted_distMat$value)
#melted_distMat <- as.data.frame(apply(melted_distMat, 2, rev))
temp <- melted_distMat
temp$value <- temp$value/max(temp$value)
temp$value <- 1-temp$value
p <- ggplot(data = temp, aes(variable, var1, fill = value))+
  geom_tile(color = "white") + scale_fill_viridis(discrete = F, option = "E")
p <- p+theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 28, hjust = 1),
        axis.text.y = element_text(size = 28))
plotTitle <- paste('../results/plots/targets_ALL_distMat_L3_discrete.png',sep='')
png(plotTitle,width = 6000, height = 6000)
plot(p)
dev.off()
#Plot histogram of distances
p <- ggplot(temp, aes(x=value)) + geom_histogram(binwidth=0.1,fill='grey')
p <-  p + theme_bw(base_size = 2*12)+xlab('Gene expression distance') + ylab('Frequency')
plotTitle <- paste('../results/plots/histogram_targetsSimilarity_L3_discrete.png',sep='')
png(plotTitle,width = 600, height = 600)
plot(p)
dev.off()
#simplify distMat
temp$value[temp$value<0.75] <- 0
temp$value[temp$value>=0.75]  <- 1 

p <- ggplot(data = temp, aes(variable, var1, fill = value))+
  geom_tile(color = "white") + scale_fill_viridis(discrete = F, option = "E")
p <- p+theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 28, hjust = 1),
        axis.text.y = element_text(size = 28))
plotTitle <- paste('../results/plots/targets_ALL_distMat_simp_L3.png',sep='')
png(plotTitle,width = 6000, height = 6000)
plot(p)
dev.off()
#regenerate targets dist. tabldis
targets.dist<- copyMelted
#get rid of product family info
targets.dist<- targets.dist %>% separate(var1, c("var1", "fam1"), "_fam_")
targets.dist<- targets.dist %>% separate(variable, c("variable", "fam2"), "_fam_")
#normalize distance
targets.dist$value <- targets.dist$value/max(targets.dist$value)
#targets.dist$value <- 1-targets.dist$value
#open fluxDist_dist matrix
filename  <- paste('../results/fluxDist_distance_allChemicals_noBio.txt',sep='')
flux.dist <- read.csv(filename,sep='\t',stringsAsFactors = FALSE)
rownames(flux.dist) <- flux.dist$Row
flux.dist <- flux.dist[,2:ncol(flux.dist )]
var1 <- rownames(flux.dist)
# Melt the distance matrix
flux.dist <- melt(flux.dist)
flux.dist <- cbind(var1,flux.dist)
flux.dist$var1 <- as.character(flux.dist$var1)
flux.dist$variable <- as.character(flux.dist$variable)

#Create a new df with both distance metrics for each chemicals pair
targets.dist$var1 <- gsub('_','',targets.dist$var1)
targets.dist$variable <- gsub('_','',targets.dist$variable)
newMatrix <- data.frame(chem1=as.character(targets.dist$var1),chem2=as.character(targets.dist$variable),distTargets=targets.dist$value,distStoich=rep(NA,nrow(targets.dist)),fam1=as.character(targets.dist$fam1),fam2=as.character(targets.dist$fam2),stringsAsFactors = FALSE)
newMatrix <- newMatrix %>% unite("pair", chem1:chem2, sep= '_and_',na.rm = FALSE, remove = TRUE)
#newMatrix<-unite(newMatrix, var1, variable, sep = "_", remove = TRUE, na.rm = FALSE)
for (i in 1:nrow(flux.dist)){
   str<- paste(flux.dist$var1[i],'_and_',flux.dist$variable[i],sep='')
  indexes <- match(str,newMatrix$pair)
  if (length(indexes)==1 & !is.na(indexes[1])){
      #print(indexes)
      newMatrix$distStoich[indexes]<- flux.dist$value[i]
  }else{
    print(str)
  }
}
#Check which chemical pairs correspond to the same family
indexes <- which(newMatrix$fam1==newMatrix$fam2)
newMatrix$related <- rep('No',nrow(newMatrix))
newMatrix$related[indexes] <- 'Yes'
newMatrix$distStoich <- newMatrix$distStoich/max(newMatrix$distStoich)
#newMatrix$distStoich <- 1-newMatrix$distStoich 
#newMatrix <- newMatrix[newMatrix$distStoich<=0.25,]

p <- ggplot(newMatrix, aes(x=distStoich, y=distTargets,color=related)) +
  geom_point(size=1) +scale_color_manual(values=c(rgb(0.65,0.65,0.65,0.25),rgb(0.7,0,0.2,0.6))) + 
  theme_bw(base_size = 32) +
  xlab('Flux distance') + ylab('Gene expression distance') +geom_smooth(aes(group = 1),se=TRUE, fullrange=TRUE)
plotName <- '../results/plots/gDist_vs_fDist_L3.pdf'
pdf(plotName,width = 10, height = 7.5)
plot(p)
dev.off()
distance.lm = lm(distStoich ~ distTargets, data=newMatrix)
summary(distance.lm)$r.squared 
dev.off()




