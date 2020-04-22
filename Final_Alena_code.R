#Loading RNA Seq counts data

monkey <- read.csv(file="monkey.csv")
colnames(monkey) <- c("gene", "EXMC_m", "PE_m", "TE_m", "EPI_m", "VE_m")

human_EPI_PE <- read.delim("human_EPI_PE.txt")
human_TE <- read.delim("human_TE.txt")

hESC_BMP4 <- read.delim("hESC_BMP4.tabular")
colnames(hESC_BMP4) <- c("geneID", "BMP4_1", "BMP4_2", "gene", "an")
hESC_BMP4 <- hESC_BMP4[ , c(4,2,3)] #get rid of unneeded cols


#Building merge on a example
#Data1 <- data.frame(gene = c("Gene1" , "Gene2"), v1 = c(0, 2), v2 = c(9,9))
#Data2 <- data.frame(geneX = c("Gene2" , "Gene3"), v1 = c(99, 2), v2 = c(99,3))
#merge.df <- merge(Data1, Data2, by.x = "gene", by.y = "geneX")

#Merge all the data based on the gene column
df <- merge(human_EPI_PE, human_TE, by.x = "gene", by.y = "gene")
df <- merge(df, hESC_BMP4, by.x = "gene", by.y = "gene")
df <- merge(df, monkey, by.x = "gene", by.y = "gene")

rm(list=setdiff(ls(), "df")) #remove everything except the merged df


#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

dim(df)[1] #we have 9220 genes for analysis

#To make data from different sets comparible, I will normalize 
#them and also multiply by 1e6 so that numbers look better
normalization <- function(x) {
  y <- floor(((x - min(x))/(max(x) - min(x)))*1e6)
  return(y)
}

#x <- c(0,5,10,15)
#normalization(x) #check how it works
df[ ,2:17] <- apply(df[ ,2:17], 2, FUN = normalization)

#install.packages("data.table")
library(data.table)
t.df <- transpose(df) #transpose
colnames(t.df) <- t.df[1, ] #gene names to the colnames
t.df <- t.df[-1, ] #delete a row with gene names

dm <- data.matrix(t.df, rownames.force = NA)

#put all info in one list
dl <- list(Genes = dm, Cells = list(Study = c(rep("Human", 9), rep("hESCs", 2), rep("Monkey", 5)), 
                                    Sample = colnames(df)[-1], 
                              Type = c("EPI", "PE", "EPI", "PE", "EPI", "PE", "TE", 
                                       "TE", "TE", "BMP4", "BMP4", "EXCM", "PE", "TE", "EPI", "VE")))

rm(list=setdiff(ls(), "dl")) #remove everything except the list dl

#Normalize libraries on ~Ref genes with edgeR
#BiocManager::install("edgeR")
library(edgeR)

dge <- DGEList(counts=t(dl$Genes), group = dl$Cells$Sample)
dim(dge) #9220 genes

#First, I want to keep only genes with sufficients reads
#keep <- rowSums(cpm(dge)>500) >= 2
#dge <- dge[keep,]
#dim(dge)  #785 genes
#head(dge)

#dge$samples$lib.size <- colSums(dge$counts) #record new col size
#dge$samples

dge <- calcNormFactors(dge)
dge$samples
#These norm factors seem insufficient for me in taking account the lib sizes, so I 
#will try to use DeSeq2 method instead of edgeR 
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("DESeq")
library("DESeq")

cds <- newCountDataSet(data.frame(dge$counts), dge$samples$group )
cds <- estimateSizeFactors( cds )
ds.factors <- sizeFactors( cds )
ds.factors #DESeq seems much better in considering the differences in lib sizes of samples

dge$samples[ ,2]/ds.factors #however, samples 11 and 12 seem to be too big in lib size though


x <- dge$counts %*% diag(1/ds.factors) #divide by scaling factors
x <- floor(t(x))
subset(x, select = "ACTB") #try to normalize on  housekeeping gene #work better

actB <- subset(x, select = "ACTB")/max(subset(x, select = "ACTB"))
actB <- as.numeric(actB)
y <- t(x) %*% diag(1/actB) #divide by ActB scaling factor
y <- floor(t(y))

#actB / x[ , 5] *1e6

dl <- list(Genes = y, Cells = list(Study = c(rep("Human", 9), rep("hESCs", 2), rep("Monkey", 5)), 
                                    Sample = dl$Cells$Sample, 
                                    Type = c("EPI", "PE", "EPI", "PE", "EPI", "PE", "TE", 
                                             "TE", "TE", "BMP4", "BMP4", "EXCM", "PE", "TE", "EPI", "VE")))


#Plot gene dispersion
cds <- estimateDispersions( cds , method="blind")
plotDispEsts(cds)


#Install Bioconductor (http://bioconductor.org/install/)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.10")

library(BiocManager) #load it

#Install mixOmics package (http://bioconductor.org/packages/release/bioc/html/mixOmics.html)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("mixOmics")

library(mixOmics) #load


genes <- dl$Genes/1e6
genes <- genes[,-5]
trans.pca <- pca(genes, ncomp = 10, center = TRUE, scale = FALSE)
trans.pca

# Load ggplot2
library(ggplot2)

# Create PCs-barplot
data <- data.frame(
  PC= colnames(trans.pca$x),  
  Variation=trans.pca$explained_variance)
positions <- rev(colnames(trans.pca$x))

# Barplot
ggplot(data, aes(x=PC, y=Variation, fill =PC)) + 
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette="Spectral") +
  coord_flip() +
  scale_x_discrete(limits = positions) + 
  theme(legend.position = "none") +
  ggtitle("PCA components") +
  theme(plot.title = element_text(hjust = 0.5))



plot(trans.pca)
plotIndiv(trans.pca, comp = c(1,2), ind.names = FALSE, 
          group = dl$Cells$Study, ellipse = FALSE, pch.levels = dl$Cells$Type,
          pch = c(16, 2, 16, 2, 16, 2, rep(8,3), rep(17,2), 15, 2, 8, 16, 0),
          legend.title = "Study", legend.title.pch = "Sample",
          legend = TRUE, title = 'Overall gene expression (785 best genes), PCA comp 1-2')

#install.packages("factoextra")
library("factoextra")
fviz_pca_ind(trans.pca, geom.ind = "point", pointshape = 21, addEllipses = FALSE,
             invisible="quali",
             pointsize = 5, 
             fill.ind = dl$Cells$Type, 
             palette = "jco", 
             label = "var",
             repel = TRUE,
             legend.title = "Sample") +
  ggtitle("Overall gene expression (9220 genes), PCA comp 1-2") +
  labs(x = "PC1 (59% expl. var)", y = "PC2 (15% expl. var)") +
  theme(plot.title = element_text(hjust = 0.5)
)

fviz_pca_ind(trans.pca, geom.ind = "point", pointshape = 21,
             invisible="quali",
             pointsize = 5, 
             fill.ind = dl$Cells$Study, 
             palette = "jco", 
             addEllipses = TRUE,
             label = "var",
             repel = TRUE,
             legend.title = "Study") +
  ggtitle("Overall gene expression (9220 genes), PCA comp 1-2") +
  labs(x = "PC1 (59% expl. var)", y = "PC2 (15% expl. var)") +
  theme(plot.title = element_text(hjust = 0.5)
)


load_scores <- trans.pca$rotation[,1]
#since we are insterested in both genes that push data
#in the left and in the right
sorted <- sort(abs(load_scores), decreasing = TRUE)
top_100 <- names(sorted[1:150])
top_100  #top_100 different genes

#find TB genes and see their loading scores
y <- subset(genes, select = top_100)
yy <- data.frame(dl$Cells$Sample, y)
rownames(yy) <- yy[ ,1] #sample names to row names
yy <- yy[ ,-1]

trans.pca$rotation[top_100, 1] #abs loading scores 

#Waterfalls
plotLoadings(trans.pca, comp = 1, method = 'mean', contrib = 'max', border = FALSE,
             xlim = c(-0.4, 0.4), title = "Genes influence on PC1 (Loading scores)",
             size.title = 1, ndisplay = 100, size.name = 0.5, size.legend = 0.3)

plotLoadings(trans.pca, comp = 2, method = 'mean', contrib = 'max', border = FALSE,
             xlim = c(-0.4, 0.4), title = "Genes influence on PC2 (Loading scores)",
             size.title = 1, ndisplay = 100, size.name = 0.5, size.legend = 0.3)

#https://www.datanovia.com/en/lessons/heatmap-in-r-static-and-interactive-visualization/
#Scale data
#install.packages("pheatmap")
library("pheatmap")

#install.packages("RColorBrewer")
library(RColorBrewer)
#display.brewer.all()
#display.brewer.pal(n = 3, name = 'PuBuGn')


yy <- scale(yy)
pheatmap(yy, cutree_rows = 8, color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(50),
         cluster_cols = FALSE, fontsize_col = 6, clustering_method = "median")





TE <- subset(genes, select = c("CDX2", "TFAP2A", "GATA3",
                                "TFAP2C", "KRT7", "SPINT1"))
rownames(TE) <- rownames(yy)
TE <- scale(TE)

#Annotation for general heatmap 

annotation_col = data.frame(Study = factor(dl$Cells$Study))
rownames(annotation_col) <- rownames(yy)
annotation_colors = list(Study=c(Human="#F0F0F0",hESCs= "#BDBDBD",Monkey="#636363"))[1]

pheatmap(TE, cluster_cols = FALSE, cutree_rows = 3, clustering_method = "median", fontsize = 11,
         annotation_row =  annotation_col , color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(50),
         annotation_colors = annotation_colors, angle_col = 0, fontsize_col = 12, fontsize_row = 12,
          main = "Expression of TE-specific genes in different samples and studies")

#kmeans_k = 6      to force 6 clusters

#color = colorRampPalette(brewer.pal(n = 3, name = 'Set1'))(50),

#brewer.pal(n = 3, name = 'Greys')






