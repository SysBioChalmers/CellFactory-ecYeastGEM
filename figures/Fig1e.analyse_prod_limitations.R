#production_capabilities
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


# Setting the working directory to the directory which contains this script
if (exists("RStudio.Version")){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(getSrcDirectory()[1])
}

filename        <- paste('../results/production_capabilities/prodCapabilities_allChemicals.txt',sep='')
df <- read.csv(filename,sep='\t',stringsAsFactors = FALSE)
df$yield2 <- (df$prodRate_ec/df$cFlux_h)*(df$MW/180)
df$prodYield_ec <- df$prodYield_ec*(df$MW/180)
#df[, "maxYield"] <- apply(df[, c(9,16)], 1, max)
df$maxYield <- pmax(df$yield2, df$prodYield_ec, na.rm = TRUE)
df$maxYield    <- 1/df$maxYield
df$protYield <- df$Pburden/(df$prodRate_ec*(df$MW/1000))
df$prodYield_ec <- 1/df$prodYield_ec
#df$MW <- df$MW/max(df$MW)
df$CCMratio <- 1-df$CCMratio

families <- c('amino acid','alcohol','alkaloid','aromatic','bioamine','fatty acids and lipids','flavonoid','organic acid','stilbenoids','terpene')
codes    <- c('AAs','alc','alk','aro','bio','FAL','fla','oAc','stb','ter')
colores  <- c('#A5CEE3','#1F78B4','#B2DF8A','#33A02B','#FB9A99','#E3211C','#FDBF6F','#FF7F00','#CAB2D6','#6A3D9A')
for (i in 1:length(families)){
  idxs <- grep(families[i],df$family)
  df$family[idxs]<-codes[i]
}

#df$family <- factor(df$family,levels = unique(df$family))
df$family <- factor(df$family,levels = codes)
colourCount2 <- length(codes)
getPalette2  <- colorRampPalette(brewer.pal(colourCount2, "Paired"))
df$limitation <- FALSE
df$limitation[df$Pburden==max(df$Pburden)] <- TRUE
#df$MW <-df$MW/180
#df$limitation[df$Pburden==max(df$Pburden)] <- 'protein'
p <- ggplot(df, aes(x=maxYield,y=protYield,shape=type,size=MW,color=family)) +
     geom_point(alpha=0.7) + theme_bw(base_size = 2*12)+ scale_x_log10()  + scale_y_log10() +
     scale_color_brewer(palette="Paired")+xlab('Substrate cost [g glucose/g product]') + ylab('Protein cost [g protein/g product]') +
     scale_size_continuous(range = c(1, 12)) + scale_size_continuous(range = c(1, 12), breaks = c(100, 200, 300))+guides(
                                    shape = guide_legend(order = 1,override.aes = list(size = 4)),
                                    color = guide_legend(order = 2,override.aes = list(size = 5)),  
                                    size = guide_legend(order = 3))  #+ scale_shape(aes(solid = limitation))
pdf('../results/production_capabilities/plots/protCost_vs_subsCost.pdf',width = 9, height = 8)
plot(p)
dev.off()

newDF <- data.frame(df$compound,df$type,df$family,df$limitation)
het_df <- newDF[newDF$df.type=='heterologous',]
het_Pconst <- (sum(het_df$df.limitation==TRUE)-1)/(nrow(het_df)-1)

nat_df <- newDF[newDF$df.type=='native',]
nat_Pconst <- sum(nat_df$df.limitation==TRUE)/nrow(nat_df)

Plim_perc <- c()
nRowFam  <- c()
for (code in codes){
  temp <- newDF[newDF$df.family==code,]
  percentage <- (sum(temp$df.limitation==TRUE))/(nrow(temp))
  Plim_perc <- c(Plim_perc,percentage)
  nRowFam  <- c(nRowFam,nrow(temp))
}
pLims <- data.frame(codes,nRowFam,Plim_perc)


