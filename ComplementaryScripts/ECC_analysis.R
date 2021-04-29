library(ggplot2)
library(tidyverse)
#library(hongR)
library(readxl)
library(viridis)

getSingleReactionFormula <- function(description, reaction, ko) {
  index <- vector()
  result <- vector()
  tt <- vector()
  for (i in 1:length(ko)){
    if(length(match(ko[i],reaction))){
      index <- match(ko[i],reaction)
      tt <- description[index]
      result[i] <- paste0(tt, collapse = ";")
    } else{
      
      result[i] <- NA
    }
  }
  return(result)
}




# Setting the working directory to the directory which contains this script
if (exists("RStudio.Version")){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(getSrcDirectory()[1])
}

# check the input!!!!!!!!
# merge all the ECC result together
# datafile
ECC_dir <- "../results/ECC/"
all_product <- list.files(ECC_dir)
All_gene <- vector()

for (x in all_product){
  print(x)
  #x <- "serine_ECCs.txt"
  ss0 <- paste(ECC_dir, x, sep = "")
  print(ss0)
  ECC0 <- read.table(ss0, header = TRUE, stringsAsFactors = FALSE)
  # remove the rows with zero ECC in low and high conditions
  ECC1 <- filter(ECC0, CC_lowGlc !=0 | CC_highGlc !=0)
  gene0 <- ECC1$genes
  All_gene <- c(All_gene, gene0)
}

All_gene_unique <- unique(All_gene)



# creat a dataframe
# low ECC
ECC_low_df <- data.frame(gene=All_gene_unique, stringsAsFactors = FALSE)
for (x in all_product){
  print(x)
  #x <- "serine_ECCs.txt"
  ss0 <- paste(ECC_dir, x, sep = "")
  print(ss0)
  ECC0 <- read.table(ss0, header = TRUE, stringsAsFactors = FALSE)
  # remove the rows with zero ECC in low and high conditions
  ECC1 <- filter(ECC0, CC_lowGlc !=0 | CC_highGlc !=0)
  # product name
  pp <- str_replace(x, "_ECCs.txt", "")
  ECC_low_df[[pp]] <- getSingleReactionFormula(ECC1$CC_lowGlc, ECC1$genes, ECC_low_df$gene)
  ECC_low_df[[pp]] <- as.numeric(ECC_low_df[[pp]])
}



# High ECC
ECC_high_df <- data.frame(gene=All_gene_unique, stringsAsFactors = FALSE)
for (x in all_product){
  print(x)
  #x <- "serine_ECCs.txt"
  ss0 <- paste(ECC_dir, x, sep = "")
  print(ss0)
  ECC0 <- read.table(ss0, header = TRUE, stringsAsFactors = FALSE)
  # remove the rows with zero ECC in low and high conditions
  ECC1 <- filter(ECC0, CC_lowGlc !=0 | CC_highGlc !=0)
  # product name
  pp <- str_replace(x, "_ECCs.txt", "")
  ECC_high_df[[pp]] <- getSingleReactionFormula(ECC1$CC_highGlc, ECC1$genes, ECC_high_df$gene)
  ECC_high_df[[pp]] <- as.numeric(ECC_high_df[[pp]])
}


## choose the glucose uptake rate used for the ECC analysis
ECC_input_df <- ECC_high_df
cut_off0 <- 24 #for high glucose uptake rate
# or choose:
# ECC_input_df <- ECC_low_df
# cut_off0 <- 52 #for low glucose uptake rate;  
## start analysis
ECC_input_df1 <- as.data.frame(t(ECC_input_df), stringsAsFactors = FALSE)
colnames(ECC_input_df1) <- ECC_input_df1[1,]
ECC_input_df1 <-  ECC_input_df1[-c(1),]
ECC_input_df1[] <- sapply(ECC_input_df1, as.numeric)
ECC_input_df1[is.na(ECC_input_df1)] <- 0
# add product annotation information
product_df <- data.frame(product=rownames(ECC_input_df1), stringsAsFactors = FALSE)
chemicals_info <- read_excel("../ComplementaryData/chemicals_info.xlsx")
chemicals_info$Name0 <- str_replace_all(chemicals_info$ecModel, ".mat", "")
chemicals_info$Name0 <- str_replace_all(chemicals_info$Name0, "^ec", "")
chemicals_info$Name0 <- str_to_lower(chemicals_info$Name0)
# product classification by Iven
filename <- paste('../results/targets_summary.txt',sep='')
targets_summary <- read.csv(filename,sep='\t',stringsAsFactors = FALSE)
targets_summary$name0 <- str_replace_all(targets_summary$models, "^ec", "")
#product_df$class <- getSingleReactionFormula(chemicals_info$class,chemicals_info$Name0,product_df$product)
product_df$class <- getSingleReactionFormula(targets_summary$chemClass,targets_summary$name0,product_df$product)
# it is found some products are not grouped
# also need a unique name of product
ECC_input_df2 <- ECC_input_df1
ECC_input_df2$product <- product_df$product
ECC_input_df2$class <- product_df$class
# PCA analysis
# it seems there is no good classification based on product families?
library(ggfortify)
# firstly remove genes with no ECCs in all products
sum_result <- mapply(sum,ECC_input_df1[,])
sum_result1 <- sum_result[which(sum_result >0)]
ECC_input_df10 <- ECC_input_df1[, names(sum_result1)]

autoplot(prcomp(ECC_input_df10), data = ECC_input_df2, colour = 'class') +
  theme(axis.text=element_text(size=20, family="Arial"),
        axis.title=element_text(size=20, family="Arial") ) +
  ggtitle('') +
  theme(panel.background = element_rect(fill = "white", color="black", size = 1)) +
  theme(legend.text=element_text(size=15),
        legend.title =element_text(size=15))
# other statistical analysis
InfluencedProducts <- c()
for ( gene in All_gene_unique ){
  print(gene)
  product_info <- ECC_input_df2[, gene]
  product_info <- product_info[product_info > 0 | product_info < 0]
  product_num <- length(product_info)
  InfluencedProducts <- c(InfluencedProducts, product_num)
}

product_num_gene <- data.frame(gene=All_gene_unique, products=InfluencedProducts, stringsAsFactors = FALSE)
# plot
# choose the top 10 genes which could affect most products in ECC analysis
product_num_gene_top_15 <- product_num_gene[product_num_gene$products > 24,]  # cut_off0 <- 52 for low glucose uptake rate;  cut_off0 <- 24 for high glucose uptake rate

Factor <- product_num_gene_top_15
Factor <- Factor[order(Factor$products,decreasing = TRUE),]
product_num_gene_top_15$gene <-factor(product_num_gene_top_15$gene, levels=Factor$gene)

fileName <- paste('../results/plots/top15_controlGenes_histogram.png',sep='')
png(fileName,width=700, height=600)
colores <- cividis(11)
colour <- colores[6]
p <- ggplot(data=product_num_gene_top_15, aes(x=gene, y=products)) +
  geom_bar(stat="identity",fill=colour) +
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
  theme(panel.background = element_rect(fill = "NA")) +
  theme(panel.grid.major = element_line(colour = "grey90")) +
  theme(axis.line = element_line(size = 1, colour = "black")) +
  theme(legend.position = c(0.85, 0.2)) +
  theme(axis.text.y=element_text(size=18, family="Arial"),
        axis.text.x=element_text(size=12, family="Arial"),
        axis.title=element_blank(),
        legend.text = element_text(size=10, family="Arial"),
        plot.title = element_text(size=18, family="Arial")) +
  ggtitle('') #+
plot(p)
dev.off()
#theme(panel.background = element_rect(fill = "white", color="black", size = 1)) 

# build a dataframe of top 15 genes which affect the most products with product families information
gene_product_family <- data.frame()
for (x in product_num_gene_top_15$gene){
  print(x)
  ss0 <- ECC_input_df2[,c(x, 'class')]
  # remove the zero ECCs
  ss0 <- ss0[ss0[[x]] > 0,]
  new_df <- as.data.frame(table(ss0$class), stringsAsFactors = FALSE)
  new_df$gene <- x
  gene_product_family <- rbind.data.frame(gene_product_family, new_df)
}
colnames(gene_product_family) <- c("class","num","gene")
gene_product_family$gene <- factor(gene_product_family$gene, levels=Factor$gene)
# Stacked
ggplot(gene_product_family, aes(fill=class, y=num, x=gene)) + 
  geom_bar(position="stack", stat="identity")+
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
  theme(axis.text=element_text(size=10, family="Arial"),
        axis.title=element_text(size=12,family="Arial"),
        legend.text = element_text(size=10, family="Arial")) +
  ggtitle('') #+

# Stacked + percent
fileName <- paste('../results/plots/top15_controlGenes_pathway.png',sep='')
png(fileName,width=900, height=600)
p <- ggplot(gene_product_family, aes(fill=class, y=num, x=gene)) + 
  geom_bar(position="fill", stat="identity")+
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
  theme(panel.background = element_rect(fill = "NA")) +
  theme(panel.grid.major = element_line(colour = "grey90")) +
  theme(axis.line = element_line(size = 1, colour = "black")) +
  theme(axis.text.y=element_text(size=18, family="Arial"),
        axis.text.x=element_text(size=12, family="Arial"),
        axis.title=element_blank(),
        legend.text = element_text(size=12, family="Arial"),
        plot.title = element_text(size=18, family="Arial")) +
  ggtitle('') #+
plot(p)
dev.off()




# check the ECC distribution of each gene across all products that this gene has effect
all_ECC <- vector()
all_ECC_gene <- vector()

for (x in colnames(ECC_input_df1)){
  print(x)
  ecc0 <- ECC_input_df1[[x]]
  if(sum(ecc0)>0){
    all_ECC <- c(all_ECC, ecc0)
    ecc_gene0 <- rep(x,each=length(ecc0))
    all_ECC_gene <- c(all_ECC_gene, ecc_gene0)
  }
}

ECC_gene_df_input <- data.frame(ECC=all_ECC, gene = all_ECC_gene, stringsAsFactors = FALSE)
# still choose the top 15 genes
ECC_gene_df_input1 <- ECC_gene_df_input[which(ECC_gene_df_input$gene %in% product_num_gene_top_15$gene),]
# plot the box plot
ECC_gene_df_input1$gene <- factor(ECC_gene_df_input1$gene, levels=Factor$gene)
#ECC_gene_df_input1$gene <- as.factor(ECC_gene_df_input1$gene)
fileName <- paste('../results/plots/top15_controlGenes_FCC_boxPlots.png',sep='')
png(fileName,width=700, height=600)
p <- ggplot(ECC_gene_df_input1, aes(x=gene, y=ECC)) + 
  geom_boxplot(fill='light blue') +
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
  theme(panel.background = element_rect(fill = "NA")) +
  theme(panel.grid.major = element_line(colour = "grey90")) +
  theme(axis.text.y=element_text(size=18, family="Arial"),
        axis.text.x=element_text(size=12, family="Arial"),
        axis.title.x =element_blank(),
        axis.title.y =element_text(size=18, family="Arial"),
        legend.text = element_text(size=12, family="Arial"),
        plot.title = element_text(size=18, family="Arial")) +
  ylab('Flux control coefficient') +
  ggtitle('')
plot(p)
dev.off()


# general analysis
product_num_gene$group <- NA
product_num_gene$group[product_num_gene$products==1] <- "a. 1 product"
product_num_gene$group[product_num_gene$products > 1 &  product_num_gene$products <=3] <- "b. 2-3 products"
product_num_gene$group[product_num_gene$products > 3 &  product_num_gene$products <=5] <- "c. 4-5 products"
product_num_gene$group[product_num_gene$products > 5 &  product_num_gene$products <=11] <- "d. 6-11 products"
product_num_gene$group[product_num_gene$products > 11 &  product_num_gene$products <= 19] <- "e. 12-19 products"
product_num_gene$group[product_num_gene$products > 19] <- "f. over 20 products"

product_num_gene1 <- product_num_gene[!is.na(product_num_gene$group),]
product_num_df <- as.data.frame(table(product_num_gene1$group), stringsAsFactors = FALSE)
# add a new row
product_num_df[nrow(product_num_df) + 1,] = c("g. over 60 products", length(which(product_num_gene$products > 60)))
product_num_df$Freq <- as.numeric(product_num_df$Freq)
colnames(product_num_df) <- c("Group", "Number")
# plot
fileName <- paste('../results/plots/FCC_product_specificity.png',sep='')
png(fileName,width=700, height=600)
p <- ggplot(data=product_num_df, aes(x=Group, y=Number)) +
  geom_bar(stat="identity",fill = colour) +
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
  theme(panel.background = element_rect(fill = "NA")) +
  theme(panel.grid.major = element_line(colour = "grey90")) +
  theme(axis.text.y=element_text(size=18, family="Arial"),
        axis.text.x=element_text(size=18, family="Arial"),
        axis.title.x =element_blank(),
        axis.title.y =element_text(size=18, family="Arial"),
        legend.text = element_text(size=12, family="Arial"),
        plot.title = element_text(size=18, family="Arial")) +
  ylab('Enzymes with FCC>0') +
  ggtitle('') #+
plot(p)
dev.off()
#theme(panel.background = element_rect(fill = "white", color="black", size = 1)) 

# heatmap
# choose the top 15 genes and prepare a heatmap
top_15_genes <- as.character(product_num_gene_top_15$gene)
ECC_input_df1_top15 <- ECC_input_df1[, top_15_genes]
# randome choose 40 products
ECC_input_df1_top15_20products <- ECC_input_df1_top15[sample(nrow(ECC_input_df1_top15), 40), ]
library(pheatmap)
pheatmap(ECC_input_df1_top15_20products,
         method = c("pearson"),
         clustering_method = "complete",
         treeheight_row = 40,
         treeheight_col = 40,
         cluster_row = FALSE,
         cluster_col = FALSE,
         show_rownames = T,
         show_colnames = T,
         legend = T,
         fontsize = 12,
         color = colorRampPalette(c("white", "SandyBrown", "firebrick3"))(100))




# tSNE plot

library(Rtsne) # Load package
library(plotly)
ECC_input_tsne <- ECC_input_df10[!duplicated(ECC_input_df10[,1:ncol(ECC_input_df10)]),]
metadata <- data.frame(sample_id = rownames(ECC_input_tsne),
                       colour = NA)
metadata$colour <- getSingleReactionFormula(targets_summary$chemClass,targets_summary$name0, metadata$sample_id)


data <- as.matrix(ECC_input_tsne)
set.seed(1)
tsne <- Rtsne(data)
df <- data.frame(x = tsne$Y[,1],
                 y = tsne$Y[,2],
                 colour = metadata$colour)
# 2D
ggplot(df, aes(x, y, colour = colour)) +
  geom_point()

# 3D
set.seed(8)
tsne_out <- Rtsne(data, dims=3, perplexity=10,check_duplicates=FALSE) # perplexity <- 10 seems good!
tsne_plot <- data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2],z = tsne_out$Y[,3],metadata$sample_id,metadata$colour)
colnames(tsne_plot)[ncol(tsne_plot)]<- 'family'
plot_ly(x=tsne_plot$x, y=tsne_plot$y, z=tsne_plot$z, type="scatter3d", mode="markers", color=tsne_plot$family,text = tsne_plot$metadata.sample_id)
