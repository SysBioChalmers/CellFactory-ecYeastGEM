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
#load chemicals info
filename <- paste('../data/chemicals_info.txt',sep='')
chemicals <- read.csv(filename,sep = '\t',stringsAsFactors = FALSE)
families <- unique(chemicals$class)
codes    <- c('_AA','_alk','_oAc','_pro','_alc','_ter','_FA','_fla','_aro','_bioAm')
pathWays <- c('Oxidative phosphorylation','Glycolysis','TCA cycle','Pentose phosphate pathway')
#load enzyme info
filename <- paste('../data/enzymeTable.txt',sep='')
enzTable <- read.csv(filename,sep = '\t',stringsAsFactors = FALSE)
#load KEGG pathways info
filename <- paste('../data/keggPathways.txt',sep='')
keggDF <- read.csv(filename,sep = '\t',stringsAsFactors = FALSE)
#Load targets summary
filename        <- paste('../results/production_targets/targets_summary.txt',sep='')
targets_summary <- read.csv(filename,sep='\t',stringsAsFactors = FALSE)
#Get pie chart for families of chemicals
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
getPalette  <- colorRampPalette(brewer.pal(5, "Set1"))
colourCount <- length(unique(df$classes))
#Create blank theme
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
#Get pie chart for chemical classes
p <- ggplot(df, aes(x='',y=perc,fill=classes))+
  geom_bar(width = 1, stat = 'identity')
pie <- p + coord_polar("y", start=0,direction=1)
pie <- pie + blank_theme 
pie <- pie +scale_fill_manual(values = getPalette(colourCount))#scale_fill_viridis(discrete = T, option = "E")
pie <- pie + labs(fill = 'Chemical families')
png('../results/plots/chemicalFamilies.png',width = 800, height = 800)
plot(pie)
dev.off()
#Get heatmap for targets matrix (binary values KOs and OEs)
targetsMat                  <- allTargetsMat
targetsMat[targetsMat == 2] <- 0
targetsMat                  <- targetsMat[rowSums(targetsMat[,5:ncol(targetsMat)])>0,]
targetsMat[targetsMat == 1] <- -1
targetsMat[targetsMat == 3] <- 1
plotTitle                   <- paste('../results/plots/targetsMatrix.png',sep='')
png(plotTitle,width = 5000, height = 13000)
rownames(targetsMat) <- targetsMat$shortNames
pheatmap(targetsMat[,5:ncol(targetsMat)],cluster_cols = T,cluster_rows = T, show_rownames = TRUE,scale='none',fontsize = 28)
dev.off()
#Dimensionality reduction (tSNE)
targetType <- c('KO','KD','OE')
for (i in 1:3){
  #Isolate targets group 
  allTargetsMat <- read.csv('../results/production_targets/targetsMatrix_mech_validated.txt',sep='\t',stringsAsFactors = FALSE)
  targetsMat <-(allTargetsMat[,5:ncol(allTargetsMat)])
  targetsMat <- as.matrix(targetsMat)
  targetsMat[targetsMat != i] <- 0
  targetsMat[targetsMat == i] <- 1
  #Extract numerical values
  newDF <- targetsMat[rowSums(targetsMat)>0,5:ncol(targetsMat)]
  #newDF <- newDF[rowSums(newDF)>0,]
  newDF <- as.data.frame(t(newDF))
  idxs  <- which(!duplicated(newDF))
  newDF <- newDF[idxs,]
  
  newDF$extra <- rownames(newDF)
  newDF <- newDF %>% separate(extra, c("chemical", "family"), "_del_")
  newDF$family <- gsub('_',' ',newDF$family)
  #t <- data.frame(newDF$chemical,newDF$family)
  set.seed(18) # Set a seed if you want reproducible results
  tsne_out  <- Rtsne(newDF[,1:(ncol(newDF)-3)],dims=3,perplexity=nrow(newDF)/length(famLvls),check_duplicates=FALSE) # Run TSNE
  tsne_plot <- data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2],z = tsne_out$Y[,3],newDF$chemical,newDF$family)
  colnames(tsne_plot)[ncol(tsne_plot)]<- 'family'
  plot_ly(x=tsne_plot$x, y=tsne_plot$y, z=tsne_plot$z, type="scatter3d", mode="markers", color=tsne_plot$family,text = tsne_plot$newDF.chemical)
  #str <- paste('../results/plots/tsne_',targetType[i],'.html',sep='')
  #htmlwidgets::saveWidget(as_widget(p),str)
}
