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
#Theme for plotting
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
dir.create('../results/plots/chassis_strain')

#Get heatmap for targets matrix
allTargetsMat <- read.csv('../results/production_targets/targetsMatrix_compatible.txt',sep='\t',stringsAsFactors = TRUE)
allTargetsMat <- read.csv('../results/production_targets/targetsMatrix_mech_validated.txt',sep='\t',stringsAsFactors = TRUE)
#Discard down-regulations (as they're not generalizable)
allTargetsMat[allTargetsMat == 2] <- 0

chass_strain_chem <- list()
actions <- c('OE','KO')
colors    <- cividis(11) 
counter <- 0
#Get target matrix for OEs and KOs for all genes that are predicted to work for more than one product
for (action in actions){
  matrix <- allTargetsMat
  if (action=='OE'){
    matrix[matrix == 1] <- 0
    matrix[matrix == 3] <- 1
  }else{
    matrix[matrix == 1] <- 1
    matrix[matrix == 3] <- 0
  }
  matrix <- matrix[,5:ncol(matrix)]
  rownames(matrix) <- allTargetsMat$shortNames
  colnames(matrix) <- colnames(allTargetsMat)[5:ncol(allTargetsMat)]
  matrix <- matrix[rowSums(matrix)>0,]
  idxs   <- sort(rowSums(matrix),decreasing= TRUE,index.return=TRUE)
  matrix <- matrix[idxs$ix,]
  matrix <- matrix[rowSums(matrix)>1,]
  if (action=='OE'){OE_mat <- matrix}
  else{KO_mat <- matrix}
}
chemicals   <- data.frame(product=colnames(matrix))
chemicals   <- chemicals %>% separate(product, c("product", "family"), "_fam_")
lastElement <- min(nrow(OE_mat),nrow(KO_mat))
modNumber   <- 18

OEs_idx   <- 1
KOs_idx   <- 1
remaining <- 1:nrow(chemicals)
geneTargets   <- c()
modifications <- c()
prodNumber    <- c()
remChems      <- list()
for (i in 1:modNumber){
  idxs_OE_i <- which(OE_mat[OEs_idx,]>0)
  idxs_KO_i <- which(KO_mat[KOs_idx,]>0)
  
  temp_OE <- intersect(remaining,idxs_OE_i)
  temp_KO <- intersect(remaining,idxs_KO_i)
  if (length(temp_OE)>=length(temp_KO)){
    remaining     <- temp_OE
    geneTargets   <- c(geneTargets,rownames(OE_mat)[OEs_idx])
    modifications <- c(modifications,'OE')
    OEs_idx   <- OEs_idx+1
  }else{
    remaining <- temp_KO
    geneTargets   <- c(geneTargets,rownames(KO_mat)[KOs_idx])
    modifications <- c(modifications,'KO')
    KOs_idx   <- KOs_idx+1
  }
  prodNumber <- c(prodNumber,length(remaining))
  counter <- counter +1
  chass_strain_chem[[counter]] <- colnames(matrix)[remaining]
  #Get a pie chart with the product families distribution for the remaining 
  #product targets aimed by the chassis_strain
  if (length(remaining)>=1){
    newTable   <- chemicals[remaining,]
    remChems[[i]] <- (chemicals[remaining,])
    #Get pie chart for families of chemicals
    classes <- unique(newTable$family)
    counts  <- c()
    for (j in 1:length(classes))
    {
      counts <- c(counts,sum(newTable$family==classes[j]))
    }
  #Get data frame for plotting
  df <- data.frame(classes,counts)
  df <- df[order(-counts),]
  df$classes <- factor(df$classes,levels<-df$classes) 
  df$perc <- round(df$counts*100/sum(df$counts),digits=1)
  #Get pallete of colors
  getPalette  <- colorRampPalette(brewer.pal(5, "Set1"))
  colourCount <- length(unique(df$classes))
  p <- ggplot(df, aes(x='',y=perc,fill=classes))+
    geom_bar(width = 1, stat = 'identity')
  pie <- p + coord_polar("y", start=0,direction=1)
  pie <- pie + blank_theme 
  pie <- pie +scale_fill_manual(values = getPalette(colourCount))#scale_fill_viridis(discrete = T, option = "E")
  pie <- pie + labs(fill = 'Chemical families')
  titleStr <- paste('Chassis strain (',OEs_idx-1,' OEs, ',KOs_idx-1,' KOs) for ', length(remaining),' products',sep='')
  pie <- pie + ggtitle(titleStr)
  png(paste('../results/plots/chassis_strain/products_chassis_strain_',i,'.png',sep=''),width = 650, height = 600)
  plot(pie)
  dev.off()
  write.table(remChems[[i]],paste('../results/chassis_strain_',i,'_chemicals.txt',sep=''),sep='\t',quote = F,row.names = F)
  }
}
chassis_strain <- data.frame(geneTargets,modifications,prodNumber)
write.table(chassis_strain,'../results/chassis_strains_modifications.txt',sep='\t',quote = F,row.names = F)

