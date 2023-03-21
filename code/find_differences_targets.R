library(ggplot2)
library(scales)
library(viridis)
library(reshape)
library(RColorBrewer)
library(tibble)
library(dplyr)
library(tidyr)
library(matrixStats)
library(reshape2)

if (!require("processx")) install.packages("processx")


# Setting the working directory to the directory which contains this script
if (exists("RStudio.Version")){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(getSrcDirectory()[1])
}

#colourCount <- length(unique(df$classes))
#getPalette <- colorRampPalette(brewer.pal(colourCount,"Paired"))

#Load targets summary
filename        <- paste('../results/production_targets/targets_summary.txt',sep='')
targets_summary <- read.csv(filename,sep='\t',stringsAsFactors = FALSE)
#classes <- unique(targets_summary$chemClass)
#counts  <- c()
#for (i in 1:length(classes))
#{
#  counts <- c(counts,sum(targets_summary$chemClass==classes[i]))
#}


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

targetAverage <- read.csv('../results/production_targets/targets_summary.txt',sep='\t',stringsAsFactors = FALSE)
x <- as.matrix(targetAverage[,((ncol(targetAverage)-2)):ncol(targetAverage)])
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
Plim <- as.character(Plim)
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

newDF$Plim <- factor(Plim,levels = unique(Plim))
newDF$origin <- factor(orgStr,levels = unique(orgStr))
values <- newDF[,1:(ncol(newDF)-4)]
allData <- newDF
modF <- c(4,0.25,0)
mods <- c('OE','KD','KO')
modifications <- matrix(0,nrow(newDF),3)
modifications[,1:3] <- as.matrix(targetAverage[,((ncol(targetAverage)-2)):ncol(targetAverage)])

newDF <- allData[,((ncol(allData)-3)):ncol(allData)]
total <- rowSums(modifications)
newDF <- cbind(newDF,modifications,total)

newDF <- as.data.frame(newDF)
colnames(newDF) <- c('chemical','family','Plim','origin','OE','KD','KO','total')
mods <- c('OE','KD','KO','total')
variables <- c('Plim','origin')
#newDF$Plim <- newDF$Plim+1
#newDF$origin <- newDF$origin+1
#newDF$Plim[newDF$Plim==0] <- 'No'
#newDF$Plim[newDF$Plim==1] <- 'Yes'

colors <- cividis(11)

for (i in 1:4){
  for (j in 1:2){
    values <- newDF[,4+i]
    type <- newDF[,2+j]
    if (j==1){
      dist1 <- values[which(type=='0')]
      dist2 <- values[which(type=='1')]
      #type[type=='0'] <- 'No'
      #type[type=='1'] <- 'Yes'
    }
    if (j==2){
      dist1 <- values[which(type=='N')]
      dist2 <- values[which(type=='H')]
    }
    plot_data <- as.data.frame(type,values)
    p <- ks.test(dist2,dist1,alternative='less')
    print(mods[i])
    print(variables[j])
    print(p)
    #colnames(plot_data) <- c('type','values')
    #plot_data$type <- as.factor(plot_data$type)
    
    #plot_data$type <- as.factor(plot_data$type,levels=c('No','Yes'))
    
    #Generate box plots with all targets (KOs, KDs and OEs)
    p <-ggplot(plot_data, aes(x=type, y=values,fill=as.factor(plot_data$type))) + geom_boxplot()
    p <- p + theme_bw(base_size = 2*12) + xlab('') +
    ylab('# of targets')+ylim(c(0,50))+labs(fill = variables[j]) +
    scale_fill_manual(values = c(colors[3],colors[9]))
    plotTitle <- paste('../results/plots/targetsNumber/distributions_targets_',mods[i],'_',variables[j],'.pdf',sep='')
    pdf(plotTitle,width = 8, height = 6)
    plot(p)
    dev.off()
  }
}
