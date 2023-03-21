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
library(htmlwidgets)
library(reshape2)


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

filename  <- paste('../results/production_targets/targets_summary.txt',sep='')
targetsDF <- read.csv(filename,sep='\t',stringsAsFactors = FALSE)
medianTargets <- c()
for (j in 1:length(families)){
  chemClass <- families[j] 
  famCode   <- codes[j]
  targets_summary <- targetsDF
  if (nchar(chemClass)>1){
   targets_summary <- targets_summary[targets_summary$chemClass==chemClass,]
  }
  nCompounds <- nrow(targets_summary)
  if (nrow(targets_summary)>=1){
    
    vector1 <- c(targets_summary$cand_1_OE,targets_summary$cand_2_OE,targets_summary$cand_3_OE)
    vector2 <- c(rep('Lvl1',nCompounds),rep('Lvl2',nCompounds),rep('Lvl3',nCompounds))
    vector3 <- c(rep('OE',length(vector1)))
    df      <- data.frame(vector1,vector2,vector3)
    values  <- cbind(mean(targets_summary$cand_1_OE),mean(targets_summary$cand_2_OE),mean(targets_summary$cand_3_OE))
    values  <- round(values,2)
    str     <- gsub('_','',famCode)
    type    <- 'OE'
    OErow   <- data.frame(str,values,type)
    
    vector1 <- c(targets_summary$cand_1_dR,targets_summary$cand_2_dR,targets_summary$cand_3_dR)
    vector2 <- c(rep('Lvl1',nCompounds),rep('Lvl2',nCompounds),rep('Lvl3',nCompounds))
    vector3 <- c(rep('KD',length(vector1)))
    df2     <- data.frame(vector1,vector2,vector3)
    values  <- cbind(mean(targets_summary$cand_1_dR),mean(targets_summary$cand_2_dR),mean(targets_summary$cand_3_dR))
    values  <- round(values,2)
    str     <- gsub('_','',famCode)
    type <- 'KD'
    KDrow   <- data.frame(str,values,type)
    
    vector1 <- c(targets_summary$cand_1_del,targets_summary$cand_2_del,targets_summary$cand_3_del)
    vector2 <- c(rep('Lvl1',nCompounds),rep('Lvl2',nCompounds),rep('Lvl3',nCompounds))
    vector3 <- c(rep('KO',length(vector1)))
    df3     <- data.frame(vector1,vector2,vector3)
    df      <- rbind(df,df2,df3)
    values  <- cbind(mean(targets_summary$cand_1_del),mean(targets_summary$cand_2_del),mean(targets_summary$cand_3_del))
    values  <- round(values,2)
    str     <- gsub('_','',famCode)
    type <- 'KO'
    KOrow   <- data.frame(str,values,type)
    
    colors <- cividis(11)
    #Generate box plots with all targets (KOs, KDs and OEs)
    p <-ggplot(df, aes(x=vector2, y=vector1, fill=vector3)) + geom_boxplot()
    p <- p + theme_bw(base_size = 2*12) + xlab('') +
    ylab('# of targets')+ylim(c(0,100))+labs(fill = 'Modification') +
    scale_fill_manual(values = c(colors[11],colors[6],colors[3]))
    plotTitle <- paste('../results/plots/targetsNumber/targetDistributions',famCode,'.pdf',sep='')
    pdf(plotTitle,width = 8, height = 6)
    plot(p)
    dev.off()
    medianTargets <- rbind(medianTargets,OErow,KDrow,KOrow)
  } 
}
colnames(medianTargets) <- c('family','Lvl1','Lvl2','LVl3','type')

actions <- c('OE','KD','KO')

for (action in actions){
temp <- medianTargets[which(medianTargets$type==action),1:(ncol(medianTargets)-1)]
temp<-melt(temp,id.vars="family")

p <- ggplot(temp,aes(x=family,y=value,fill=factor(variable)))+
  geom_bar(stat="identity",position="dodge")+theme_bw(base_size = 2*12) +
  scale_fill_manual(name='',values = c(colors[10],colors[6],colors[3]))+
  xlab("Family")+ylab("Mean number of targets per product")
plotTitle <- paste('../results/plots/targetsNumber/meanTargets_AllFams',action,'.pdf',sep='')
pdf(plotTitle,width = 8.5, height = 7)
plot(p)
dev.off()
}
write.table(medianTargets,'../results/processed_results/average_targets_number.txt',quote=FALSE,sep='\t',row.names = FALSE)


