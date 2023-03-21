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
library(stringr)

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
data <- read.csv('../results/met_cofactors_turnovers_allChemicals.txt',sep='\t',stringsAsFactors = TRUE)
data <- read.csv('../results/met_precursors_turnovers_allChemicals.txt',sep='\t',stringsAsFactors = TRUE)

#data <- t(data)
#data <- as.data.frame(data)
#data <- write.table(data,'../results/met_precursors_turnovers_allChemicals_T.txt',sep='\t')

chemicals <- data$Var1
data <- data[,2:ncol(data)]
dataWT <- c(1, 0.7776,0.0604,0.0308,1.5030,1.5030,0.7928,0.7379,0.1074,1.2212,6.9708E-7,0.2457)
#for (i in 1:nrow(data)){
#  data[i,] <- data[i,]*dataWT[i]
#}
data[data>100]=100
data[data<0.01]=0.01
#data <- data+1E-3
data <- log(data,10)
vector <- data[1,]
#data <- scale(data,center=F,scale =t(vector))
chemicals <- c('G6P','R5P','G3P','PEP','CoA','SUC','F6P','E4P','3PG','PYR','OXO','OAC')
#chemicals <- c('ATP','NADH','NADPH','coenzyme A')

rownames(data) <- chemicals
#Isolate numeric values
t.data <- t(data)
#get a heatmap
colores <- cividis(100)
plotName <- '../results/met_turnover/cof_turnover_heatmap.pdf'
plotName <- '../results/met_turnover/met_turnover_heatmap.pdf'

pdf(plotName,width = 6, height = 20)

pheatmap(t.data,
                  method = c("pearson"),
                  clustering_method = "complete",
                  cluster_row = T,
                  cluster_col = T,
                  show_rownames = T,
                  show_colnames = T,
                  legend = T,
                  fontsize = 12,
                  color = colores)

dev.off()
df <- t(data)
#t.data <- t.data[,2:ncol(t.data)]
PCAdata   <- prcomp(t.data, center = FALSE,scale = FALSE,retx=TRUE)
t.data <- as.data.frame(t.data)


p         <- autoplot(PCAdata,data = t.data)
plotTitle <- paste('../results/met_turnover/PCA_cof_TO.png',sep='')
png(plotTitle,width = 600, height = 600)
plot(p)
dev.off()
#Manual table
data2 <- read.csv('../results/met_precursors_turnovers_allChemicals_wBio.txt',sep='\t',stringsAsFactors = TRUE)
data2[,1] <- c('G6P','F6P','R5P','E4P','G3P','3PG','PEP','PYR','ACC','AKG','SUC','OXO')
data2 <- t(data2)
colnames(data2) <- data2[1,]
data2 <- as.data.frame(data2[2:nrow(data2),])
data2 <- cbind(rownames(data2),data2)
colnames(data2)[1] <- 'Chemical'
data2 <- data2[order(data2$Chemical),] 
#data2 <- data2[3:nrow(data2),]
clusters <- read.csv('../results/cluster_strains/product_clusters.txt',sep='\t',stringsAsFactors = TRUE)
data2$cluster <- clusters$cluster
write.table(data2,'../results/prec_mTO_clusters.txt',sep='\t')


chem_origin <- read.csv('../results/production_targets/chemicals_origin.txt',sep='\t',stringsAsFactors = FALSE)
chem_origin <- chem_origin[order(chem_origin$b1),]
chem_origin$b1 <- gsub('_','',chem_origin$b1)
origin <- chem_origin$b3 == 'native'
origin <- as.numeric(origin)

prodCap <- read.csv('../results/production_capabilities/prodCapabilities_allChemicals_wBio.txt',sep='\t',stringsAsFactors = FALSE)
prodCap <- prodCap[order(prodCap$compound),]
prodCap$compound <- gsub('_','',prodCap$compound )

#Plim <- prodCap$Pburden/max(prodCap$Pburden)
idxs <- which(prodCap$cFlux_l>= 0.99)
Plim <- rep(1,nrow(prodCap))
Plim[idxs] <- 0
#Plim[Plim!=1] <- 0

origin <- as.numeric(origin)
results <- data.frame(chem_origin$b1,chem_origin$b2,chem_origin$b3,Plim) 
colnames(results) <- c('chem','family','origin','Plim')
results$cluster <- clusters$cluster
#results <- results[order(results$precGroup,results$cluster),] 
#results <- results[order(-results$Plim),] 
# for (i in 1:nrow(clusters)){
#   comps <- str_split(clusters$internal_ids[i], ' // ', n = Inf, simplify = FALSE)
#   comps <- as.data.frame(comps)
#   colnames(comps) <- 'names'
#   #comps$names <- gsub('_','',comps$names)
#   idxs <- match(comps$names,results$chem)
#   results$cluster[idxs] <- i
#   #print(comps)
#   print(idxs)
# }
#results$cluster[results$cluster==0] <- NA
results$cluster <- as.character(results$cluster)

metDF <- data2[,2:ncol(data2)]
metTurnover <- cbind(results,metDF)
metTurnover$Plim[metTurnover$Plim==0]<-5
metTurnover$Plim[metTurnover$Plim==1]<-10
#metTurnover$precGroup[is.na(metTurnover$precGroup )]<-0
#metTurnover$precGroup <- as.character(metTurnover$precGroup)
metTurnover$cluster <- as.character(metTurnover$cluster)

famLvls <- as.numeric(unique(factor(metTurnover$family)))
famLvls <- (unique(factor(metTurnover$family)))
famLvls <- famLvls[order(famLvls)]
metTurnover$family <- factor(metTurnover$family,levels = famLvls)
write.table(metTurnover,'../results/met_turnover/metTO_clusters.txt',sep='\t',row.names = FALSE)

for (i in 1:25)
{
  set.seed(18) # Set a seed if you want reproducible results
  perplxty <- i
  tsne_out  <- Rtsne(metTurnover,dims=2,perplexity=perplxty, max_iter = 10000,theta=0.1) # Run TSNE
  tsne_plot <- data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2],compound=metTurnover$chem,family=metTurnover$family,origin=metTurnover$origin,Plim=metTurnover$Plim,cluster =metTurnover$cluster)
  #Define color pallete
  colourCount <- length(unique(tsne_plot$family))
  getPalette2  <- colorRampPalette(brewer.pal(colourCount, "Paired"))
  colnames(tsne_plot)[4] <- 'family'
  colnames(tsne_plot)[3] <- 'compound'
  p <- plot_ly(tsne_plot,x=~x, y=~y, text =~compound, type="scatter", mode="markers", color=~cluster,colors = getPalette2(colourCount),symbol=~metTurnover$origin,symbols=c(20,8),size=metTurnover$Plim)%>% 
    layout(title= list(text = paste('tsne_',i,'_perplexity')))
  plotTitle <- paste('../results/plots/metTO/tSNE_metTO_cluster_',i,'.html',sep='')
  #orca(p, plotTitle,width=900,height=600)
  name <-paste(plotTitle,sep = '')
  saveWidget(p, plotTitle, selfcontained = F, libdir = "lib")
}
