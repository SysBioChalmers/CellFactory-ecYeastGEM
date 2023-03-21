#Gene-centric analysis of metabolic targets
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
library(fmsb)

# Setting the working directory to the directory which contains this script
if (exists("RStudio.Version")){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(getSrcDirectory()[1])
}
dir.create('../results/gene_centric')

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

#load matrix indicating the relation between gene targets and chemical products
allTargetsMat <- read.csv('../results/production_targets/targetsMatrix_L3.txt',sep='\t',stringsAsFactors = TRUE)
prot_lims <- read.csv('../results/production_capabilities/proteinLimitations_allChemicals.txt',sep='\t',stringsAsFactors = TRUE)
#scatter plot for prot limitations
p <- ggplot(prot_lims, aes(x=prot_lims$Prod_FC, y=prot_lims$Prot_cost)) +
  geom_point(size=2, shape=23)
plot(p)
#allTargetsMat <- read.csv('../results/production_targets/targetsMatrix_ecFSEOF.txt',sep='\t',stringsAsFactors = TRUE)
targetsMat <- allTargetsMat
direction  <- c('OE','KD','KO')
colors     <- cividis(11) 
threshold <- 10

productsNmr <- c()
chemicals <- data.frame(product = colnames(allTargetsMat)[5:(ncol(allTargetsMat))])
chemicals <- chemicals %>% separate(product, c("product", "family"), "_fam_")
families  <- unique(chemicals$family)
families  <- c(families,'ALL')
#families  <- families[1]
commonTargets <- c()
columns <- c()
for (fam in families){
  newCol <- c()
  for (action in direction){
    matrix <- targetsMat[,5:ncol(targetsMat)]
    
    if (action=='OE'){
      matrix[matrix != 4] <- 0
      matrix[matrix == 4] <- 1
      barColor <- colors[10]
      j <- 1
    }  
    if (action=='KD'){
      matrix[matrix != 0.25] <- 0
      matrix[matrix == 0.25] <- 1
      barColor <- colors[6]
      j <- 2
    }
    if (action=='KO'){
      matrix[matrix != 0] <- 100
      matrix[matrix == 0] <- 1
      matrix[matrix == 100] <- 0
      barColor <- colors[2]
      j <- 3
    }
    targetsMat <- allTargetsMat
    rownames(matrix) <- targetsMat$shortNames
    colnames(matrix) <- colnames(targetsMat)[5:ncol(targetsMat)]
    if (fam!='ALL'){
      fam_str <- paste('_fam_',fam,sep='')
      matrix  <- matrix[,grep(fam_str,colnames(matrix))]
    }
    if (is.numeric(matrix)){
      matrix <- matrix[matrix>0] 
      matrix <- data.frame(matrix,suma=matrix)
    } else{
      matrix <- matrix[rowSums(matrix)>0,]
      #Keep those genes that are targets
      matrix$suma <- rowSums(matrix)
    }
    matrix <- matrix[order(matrix$suma,decreasing=TRUE),]
    top10  <- matrix[1:threshold,]
    sumas  <- colSums(top10)
    newDF  <- data.frame(genes = rownames(matrix),occurrence = matrix$suma)
    nChems <- ncol(matrix)-1
    panGenes_num <- length(which(newDF$occurrence>=(nChems)))
    print(fam)
    print(action)
    print(newDF$genes[which(newDF$occurrence>=(nChems))])
    newCol <- c(newCol,panGenes_num)
    newDF  <- newDF[1:threshold,]
    if (fam!='ALL'){
      newDF$occurrence <- newDF$occurrence/(ncol(matrix)-1)
    }
    #generate bar plots for top genes
    p <- ggplot(newDF,aes(x= reorder(genes,-occurrence),occurrence))+geom_bar(stat ="identity",fill=barColor)+  
      theme_bw(base_size = 2*12)+xlab('Gene targets') +
      ylab('Occurences')#+ scale_y_continuous(breaks = (seq(0,5,by = 0.5)))
    if (fam!='ALL'){
      p <- p + ylim(c(0,1))
    }
    plotTitle <- paste('../results/gene_centric/top', action,'s_',fam,'.png',sep='')
    png(plotTitle,width = 850, height = 600)
    plot(p)
    dev.off()
    #Generate stacked bar plots for top represented genes, indicating products family
    if (fam=='ALL'){
      chemProds <- data.frame(product = colnames(matrix)[1:(ncol(matrix)-1)])
      chemProds <- chemProds %>% separate(product, c("product", "family"), "_fam_")
      geneNames <- c()
      gene_fam  <- c()
      for (k in 1:nrow(top10)){
        indexes <- which(top10[k,1:nrow(chemProds)]>0)
        if (length(indexes)>=1){
          for (l in indexes){
            geneNames <- c(geneNames,rownames(top10)[k])
            gene_fam  <- c(gene_fam,chemProds$family[l])
          }
        }
      }
      topGenesDF <- data.frame(geneNames,Family=gene_fam,stringsAsFactors = TRUE)
      topGenesDF$geneNames <- factor(topGenesDF$geneNames,levels = rownames(top10),ordered = TRUE)
      #plot
      colourCount <- length(unique(topGenesDF$Family))
      getPalette  <- colorRampPalette(brewer.pal(colourCount, "Paired"))
      p <- ggplot(topGenesDF, aes(fill=Family, x=(geneNames))) + 
        geom_bar(colour = 'black',position="stack", stat="count") + theme_bw(base_size = 2*12)+
        xlab('') + ylab('Number of chemicals') + scale_fill_manual(values = getPalette(colourCount))
      pdf(paste('../results/gene_centric/topGenes_chemFam_',action,'_',fam,'.pdf',sep=''),width = 12, height = 6)
      plot(p)
      dev.off()
    }
    #Get histogram for product-gene specificity
    if (fam=='ALL'){
      #matrix  <- targetsMat[,5:ncol(targetsMat)
      nonZ <- rowSums(matrix[,1:(ncol(matrix)-1)])>0
      sumas   <- rowSums(matrix[,1:(ncol(matrix)-1)])[nonZ]
    #  sumas   <- rowSums(matrix)#[rowSums(matrix)>0]
      geneNmr <- c()
      newDF <- data.frame(chemicals = rownames(matrix)[nonZ],productNmr =sumas)
      p <- ggplot(newDF, aes(x=productNmr)) + geom_histogram(binwidth=1,fill=barColor)
      p <-  p + theme_bw(base_size = 2*12)+xlab('Number of products') + ylab('Number of gene targets')
      #p <-  p + xlim(c(1,100))
      pdf(paste('../results/gene_centric/target_specificity_',action,'.pdf',sep=''),width = 6.5, height = 6)
      plot(p)
      dev.off()
    }
  }
  if (nChems>2 & fam!='ALL'){
    commonTargets <- cbind(commonTargets,(newCol))
    columns <- c(columns,fam)
  }
}
rows <- c('OE','KD','KO')
rownames(commonTargets) <- rows
colnames(commonTargets) <- columns
maxLim <- 20
minLim <- 0
commonTargets <- rbind(rep(maxLim,ncol(commonTargets)),rep(0,ncol(commonTargets)),commonTargets)
commonTargets <- as.data.frame(commonTargets,stringsAsFactors = FALSE)

#plot spider plot
colors_border = c(rgb(0.8,0.6,0,0.8), rgb(0.4,0.4,0.40,0.8),rgb(0.1,0,0.8,0.8))
colors_in     = c(rgb(0.8,0.6,0,0.2), rgb(0.4,0.4,0.40,0.2),rgb(0.1,0,0.8,0.2))
plotName <- '../results/plots/panGenes_by_chemFam.pdf'
pdf(plotName,width = 5.5, height = 5)
radarchart( commonTargets , axistype=1 , 
            #custom polygon
            pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="black", caxislabels=seq(minLim,maxLim*100,(maxLim-minLim)/4), cglwd=1.5,
            #custom labels
            vlcex=2, calcex = 1.5)
#plot(p)
legend(x=1.1, y=1.1, legend = rows, bty = "n", pch=20 , col=colors_in , text.col = "black", cex=1.5, pt.cex=3)
dev.off()



