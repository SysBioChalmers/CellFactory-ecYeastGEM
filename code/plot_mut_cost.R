#production_capabilities
#library(ggplot)
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

filename <- paste('../results/production_capabilities/mutant_costs.txt',sep='')
df       <- read.csv(filename,sep='\t',stringsAsFactors = FALSE)
# plot for psilocybin
df_psil  <- df[df$product=='psilocybin',]
df_psil$O2_in <- df_psil$O2_in/max(df_psil$O2_in)
df_psil$Pcost <- df_psil$Pcost
p <- ggplot(df_psil, aes(x=Ccost,y=Pcost,size=OEf,color = O2_in)) +
     geom_point(alpha=0.6) + theme_bw(base_size = 2*12)+ xlab('Substrate cost [g glucose/g product]') + ylab('Protein cost [g protein/g product]') +
     scale_color_continuous(name='Norm.  O2 in')  +
     scale_size_continuous(breaks = c(1,2,5,10,20,50,100),range = c(5,20),name='OE factor')
pdf('../results/production_capabilities/plots/mutCosts_psilocybin.pdf',width =9, height = 8)
plot(p)
dev.off()

# plot for valencene
df_valen <- df[df$product=='valencene',]
df_valen$O2_in <- df_valen$O2_in/max(df_valen$O2_in)
df_valen$Pcost <- df_valen$Pcost
p <- ggplot(df_valen, aes(x=Ccost,y=Pcost,size=OEf,color = O2_in)) +
  geom_point(alpha=0.6) + theme_bw(base_size = 2*12)+ xlab('Substrate cost [g glucose/g product]') + ylab('Protein cost [g protein/g product]') +
  scale_color_continuous(name='Norm.  O2 in')  +
  scale_size_continuous(breaks = c(1,2,5,10,20,50,100),range = c(5,20),name='OE factor')
pdf('../results/production_capabilities/plots/mutCosts_valencene.pdf',width =9, height = 8)
plot(p)
dev.off()