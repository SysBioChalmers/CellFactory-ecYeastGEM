library(ggplot2)
library(scales)
library(viridis)

# Setting the working directory to the directory which contains this script
if (exists("RStudio.Version")){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(getSrcDirectory()[1])
}
families <- c('amino acid','alkaloid','organic acid','protein','alcohol','terpene','fatty acid','natural pigment','flavonoid','aromatic','bioamine')
codes    <- c('_AA','_alk','_oAc','_pro','_alc','_ter','_FA','_nPg','_fla','_aro','_bioAm')
#Get pie chart for families of chemicals
filename <- paste('../results/targets_summary.txt',sep='')
targets_summary <- read.csv(filename,sep='\t',stringsAsFactors = FALSE)

classes <- unique(targets_summary$chemClass)
counts  <- c()
for (i in 1:length(classes))
{
  counts <- c(counts,sum(targets_summary$chemClass==classes[i]))
}

df <- data.frame(classes,counts)
df <- df[order(-counts),]
df$classes <- factor(df$classes,levels<-df$classes) 
df$perc <- round(df$counts*100/sum(df$counts),digits=1)

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=26, face="bold")
  )
  
p <- ggplot(df, aes(x='',y=perc,fill=classes))+
  geom_bar(width = 1, stat = 'identity')
pie <- p + coord_polar("y", start=0,direction=1)
pie <- pie + blank_theme + theme(axis.text.x=element_blank())#+
      #geom_text(aes(y = perc + c(0, cumsum(perc)[-length(perc)]), 
      #          label = perc), size=5)
pie <- pie +scale_fill_viridis(discrete = T, option = "E")
pie <- pie + labs(fill = 'Chemical families')
pie


filename    <- paste('../results/targets_summary.txt',sep='')
targetsSummary <- read.csv(filename,sep='\t',stringsAsFactors = FALSE)
for (j in 1:length(families)){
  chemClass <- families[j] 
  famCode   <- codes[j]
  targets_summary <- targetsSummary
if (nchar(chemClass)>1){
  targets_summary <- targets_summary[targets_summary$chemClass==chemClass,]
}
nCompounds <- nrow(targets_summary)
vector1 <- c()
vector2 <- c()
vector3 <- c()
vector1 <- c(targets_summary$cand_1_OE,targets_summary$cand_2_OE,targets_summary$cand_3_OE)
vector2 <- c(rep('step1',nrow(targets_summary)),rep('step2',nrow(targets_summary)),rep('step3',nrow(targets_summary)))
vector3 <- c(rep('OE',length(vector1)))
df <- data.frame(vector1,vector2,vector3)
vector1 <- c(targets_summary$cand_1_del,targets_summary$cand_2_del,targets_summary$cand_3_del)
vector2 <- c(rep('step1',nrow(targets_summary)),rep('step2',nrow(targets_summary)),rep('step3',nrow(targets_summary)))
vector3 <- c(rep('deletions',length(vector1)))
df2 <- data.frame(vector1,vector2,vector3)
df <- rbind(df,df2)
#colnames(df) <- c('chemicals','step1','step2','step3')
p<-ggplot(df, aes(x=vector2, y=vector1, fill=vector3)) +
  geom_boxplot()
p <- p+theme_bw(base_size = 2*12) + xlab('') +
  ylab('# of targets')+ylim(c(0,150))+labs(fill = 'Modification')
plotTitle <- paste('../results/plots/allTargets',famCode,'.png',sep='')
png(plotTitle,width = 600, height = 600)
plot(p)
dev.off()

filename    <- paste('../results/all_deletions.txt',sep='')
all_deletions <- read.csv(filename,sep='\t',stringsAsFactors = FALSE)
if (nchar(chemClass)>1){
  all_deletions <- all_deletions[all_deletions$chem_class_del==chemClass,]
}
uniqueDels <- unique(all_deletions$del_targets)
counts <- c()
for (i in 1:length(uniqueDels))
{
  counts <- c(counts,sum(all_deletions$del_targets==uniqueDels[i]))
}
df <- data.frame(uniqueDels,counts,stringsAsFactors = FALSE)
df$counts <- df$counts/nCompounds
dt <- df[order(-counts),]
dt$uniqueDels <- factor(dt$uniqueDel, levels<- dt$uniqueDel)
dt <- dt[1:10,]

p <- ggplot(data=dt, aes(x=uniqueDels,y=(counts))) + geom_bar(stat='identity',fill='#00bfc4')+ 
  theme_bw(base_size = 2*12)+xlab('Gene targets') +
  ylab('Relative frequency')#+ scale_y_continuous(breaks = (seq(0,5,by = 0.5)))
plotTitle <- paste('../results/plots/topDeletions',famCode,'.png',sep='')
png(plotTitle,width = 950, height = 600)
plot(p)
dev.off()

filename    <- paste('../results/all_OEs.txt',sep='')
all_OEs <- read.csv(filename,sep='\t',stringsAsFactors = FALSE)
if (nchar(chemClass)>1){
  all_OEs <- all_OEs[all_OEs$chem_class_OE==chemClass,]
}
uniqueOEs <- unique(all_OEs$OE_targets)
counts <- c()
for (i in 1:length(uniqueOEs))
{
  counts <- c(counts,sum(all_OEs$OE_targets==uniqueOEs[i]))
}
df <- data.frame(uniqueOEs,counts,stringsAsFactors = FALSE)
df$counts <- df$counts/nCompounds
dt <- df[order(-counts),]
dt$uniqueOEs <- factor(dt$uniqueOEs, levels<- dt$uniqueOEs)
dt <- dt[1:10,]

p <- ggplot(data=dt, aes(x=uniqueOEs,y=(counts))) + geom_bar(stat='identity',fill='#f8766d')+ 
  theme_bw(base_size = 2*12)+xlab('Gene targets') +
  ylab('Relative frequency')#+ scale_y_continuous(breaks = (seq(0,5,by = 0.5)))
plotTitle <-paste('../results/plots/topOEs',famCode,'.png',sep='')
png(plotTitle,width = 950, height = 600)
plot(p)
dev.off()
}
