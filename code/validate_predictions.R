library(matrixStats)
library(tibble)
library(dplyr)
library(tidyr)
library(reshape2)
library(stringr)
# Setting the working directory to the directory which contains this script
if (exists("RStudio.Version")){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(getSrcDirectory()[1])
}
#families <- c('amino acid','alkaloid','organic acid','protein','alcohol','terpene','fatty acids and lipids','flavonoid','aromatic','bioamine','stilbenoids')
#codes    <- c('_AA','_alk','_oAc','_pro','_alc','_ter','_FAL','_fla','_aro','_bio','_stb')
#Load targets summary
filename        <- paste('../data/genetic_background.txt',sep='')
genetic_background <- read.csv(filename,sep='\t',stringsAsFactors = FALSE)
#classes <- unique(targets_summary$chemClass)
#Analyse targets matrix
allTargetsMat <- read.csv('../results/production_targets/targetsMatrix_L2.txt',sep='\t',stringsAsFactors = FALSE)
targetsMat <- allTargetsMat
cases      <- ncol(targetsMat) - 4
targetsMat <- targetsMat[rowSds(as.matrix(targetsMat[,5:ncol(targetsMat)]))!=0,]

newDF <- targetsMat[,5:ncol(targetsMat)]
genes <- targetsMat$genes
extra <- colnames(newDF)
rownames(newDF) <- genes
newDF <- as.data.frame(t(newDF))
newDF$extra <- extra
newDF<- newDF %>% separate(extra, c("chemical", "family"), "_fam_")
newDF <- newDF[sort(newDF$chemical),]
chemIDs <- newDF$chemical
expTargets <- newDF[,1:(ncol(newDF)-2)]
expTargets[,] <- 1
rownames(expTargets) <- chemIDs
results <- c()

genes <- colnames(newDF)
genes <-genes[1:(length(genes)-2)]
idxs <- match(genes,allTargetsMat$genes)
genes2 <- allTargetsMat$genes[idxs]
shortNames <- allTargetsMat$shortNames[idxs]
#colnames(newDF) <-shortNames
allTargets <- c()
genetic_background$down <- paste(genetic_background$downregulation,',',genetic_background$deletion)
for (chem in newDF$chemical){
  idx <- which(genetic_background$internal_ids==chem)
  row <- which(chemIDs==chem)
  if (length(idx)>1){
    print(chem)
  }else{
    countOE <- 0
    countKD <- 0
    #countKO <- 0
    
    OEs <- str_split(genetic_background$overexpression[idx],",")
    if (nchar(OEs[[1]][1])>=1){
      idxs <- match(OEs[[1]],colnames(newDF))
      idxs <- idxs[!is.na(idxs)]
      print(colnames(newDF)[idxs])
      position <- which(newDF[row,idxs]==4)
      position <- idxs[position]
      gTargetsOE <- shortNames[position]
      countOE <- length(position)
      allTargets <- c(allTargets,gTargetsOE)
      gTargetsOE <- str_c(gTargetsOE, collapse = ",")
      #expTargets[row,idxs] <- 4
    }

    KDs <- str_split(genetic_background$down[idx],",")
    if (nchar(KDs[[1]][1])>=1){
      idxs <- match(KDs[[1]],colnames(newDF))
      idxs <- idxs[!is.na(idxs)]
      print(colnames(newDF)[idxs])
      position <- which(newDF[row,idxs]<=0.25)
      position <- idxs[position]
      gTargetsKD <- shortNames[position]
      countKD <- length(position)
      allTargets <- c(allTargets,gTargetsKD)
      gTargetsKD <- str_c(gTargetsKD, collapse = ",")
      #expTargets[row,idxs] <- 0.25
    }
    
    #KOs <- str_split(genetic_background$deletion[idx],",")
    #if (nchar(KOs[[1]][1])>=1){
    #  idxs <- match(KOs[[1]],colnames(newDF))
    #  idxs <- idxs[!is.na(idxs)]
    #  countKO <- length(which(newDF[row,idxs]==0))
    #}
    if (sum(c(countOE,countKD))>0){
    results <- rbind(results,cbind(chem,gTargetsOE,gTargetsKD))
    }
  }
}
allTargets <- unique(allTargets)
#sumas <- rowSums(as.numeric(results[,2:4]))
#results <- results[,]
write.table(results,'../results/processed_results/validated_targets.txt',sep='\t',row.names = FALSE)



