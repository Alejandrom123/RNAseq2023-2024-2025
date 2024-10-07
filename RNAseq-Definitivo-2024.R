#Script definitivo con: 
  #Muestras sin modificaciones
  #Muestras L4 y L9 eliminadas 
  #Muestras P3 eliminada 

library("ggplot2")
library("gridExtra")
library("edgeR")
library("UpSetR")
library("dplyr")
library(gplots)
library(ggplot2)
library(RColorBrewer)
library(dendextend)

####1. Distancias Poisson Lambdacialotrina####

countTable = read.table("C:/Users/alejandro/Desktop/Segundo artículo Transcriptómica Aedes/After_bam-R/1_Conteos/Counts_RNAseq_subsample_2024_V2.txt", sep="\t", header=T, row.names = 1)

#Remover columnas de counTable
countTable$Chr <- NULL
countTable$Start <- NULL
countTable$End <- NULL
countTable$Strand <- NULL
countTable$Length <- NULL

#Mantener columnas de counTable de Columna 4 a 12
countTable <-countTable[,4:12]


colData = read.table("C:/Users/Alejandro/Desktop/Segundo artículo Transcriptómica Aedes/After_bam-R/Metadata/Factors_Metadata_especificos_subsample_V2.txt", sep="\t", header=T, row.names = 1)

#Mantener filas lambdacialotrina 4 a 12 
colData = colData[(4:12), , drop = F]


#Confirmar nombres de columnas y filas
all(colnames(countTable) %in% rownames(colData))

all(colnames(countTable) == rownames(colData))

countTable[, "max"]= apply(countTable[, 1:ncol(countTable)], 1, max)
countTable=countTable[countTable[,ncol(countTable)]>32,] 
countTable= countTable [,-ncol(countTable)]



library(DESeq2)
dds<-DESeqDataSetFromMatrix(countData= countTable,colData= colData,design= ~ estado_susceptibilidad)
dds<-DESeq(dds)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

# For samples
correlation= cor(normalized_counts, method= "spearman")
colors = colorRampPalette(brewer.pal(9, "RdBu"))(100)
heatmap.2(correlation, main= "Heatmap correlación Spearman entre muestras", trace="none", col = colors , margin=c(10, 10), scale = "none")
heatmap(correlation, scale= "column", main= "Heatmap de correlación entre muestras DEG", margins = c(15,10))

#Distancia poisson con todos los genes 


#install.packages("pheatmap")
library("pheatmap")
install.packages("PoiClaClu")
library("PoiClaClu")
library("RColorBrewer")

# Heatmap, Poisson
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
colnames(samplePoisDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors,
         main = "Poisson distances")

#Distancia poission con solo DEG´s 

####2. Distancias Poisson Permetrina####

countTable = read.table("C:/Users/alejandro/Desktop/Segundo artículo Transcriptómica Aedes/After_bam-R/1_Conteos/Counts_RNAseq_subsample_2024_V2.txt", sep="\t", header=T, row.names = 1)

#Remover columnas de counTable
countTable <-countTable[,18:26]

colData = read.table("C:/Users/Alejandro/Desktop/Segundo artículo Transcriptómica Aedes/After_bam-R/Metadata/Factors_Metadata_especificos_subsample_V2.txt", sep="\t", header=T, row.names = 1)

colData = colData[(13:21), , drop = F]

all(colnames(countTable) %in% rownames(colData))

all(colnames(countTable) == rownames(colData))

countTable[, "max"]= apply(countTable[, 1:ncol(countTable)], 1, max)
countTable=countTable[countTable[,ncol(countTable)]>32,] 
countTable= countTable [,-ncol(countTable)]

library(DESeq2)
dds<-DESeqDataSetFromMatrix(countData= countTable,colData= colData,design= ~ estado_susceptibilidad)
dds<-DESeq(dds)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
colnames(samplePoisDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors,
         main = "Poisson distances")
###3. Correlación muestras Lambdacialotrina exclusivamente - aquí con L4 Y L9####

countTable = read.table("C:/Users/alejandro/Desktop/Segundo artículo Transcriptómica Aedes/After_bam-R/1_Conteos/Counts_RNAseq_subsample_2024_V2.txt", sep="\t", header=T, row.names = 1)

#Remover columnas de counTable
countTable$Chr <- NULL
countTable$Start <- NULL
countTable$End <- NULL
countTable$Strand <- NULL
countTable$Length <- NULL
countTable$L4_filtered.bam <- NULL
countTable$L9_filtered.bam <- NULL
#Mantener columnas de counTable de Columna 4 a 12
countTable <-countTable[,4:12]


colData = read.table("C:/Users/Alejandro/Desktop/Segundo artículo Transcriptómica Aedes/After_bam-R/Metadata/Factors_Metadata_especificos_subsample_V2.txt", sep="\t", header=T, row.names = 1)

#Mantener filas lambdacialotrina 4 a 12 
colData = colData[(4:12), , drop = F]

#Remover la fila 4 y la 9
colData <- colData[-4, , drop = FALSE]
colData <- colData[-8, , drop = FALSE]

#Confirmar nombres de columnas y filas
all(colnames(countTable) %in% rownames(colData))

all(colnames(countTable) == rownames(colData))

countTable[, "max"]= apply(countTable[, 1:ncol(countTable)], 1, max)
countTable=countTable[countTable[,ncol(countTable)]>32,] 
countTable= countTable [,-ncol(countTable)]

library(DESeq2)
dds<-DESeqDataSetFromMatrix(countData= countTable,colData= colData,design= ~ estado_susceptibilidad)
dds<-DESeq(dds)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

#Esta siguiente no tiene sentido con todas las muestras y los conteos normalizados
clusters <- hclust(dist(t(normalized_counts), method = "correlation"))
plot(clusters, hang = -1)
rect.hclust(clusters , k = 5, border = 4:4)
Dendogram <- as.dendrogram(clusters)
Dendogram <- color_branches(clusters, k = 3)
plot(Dendogram)
grid.draw(dendrogram(Dendogram, horiz = TRUE, leaflab = "none"))

###Para todos los genes
# For samples
correlation= cor(normalized_counts, method= "spearman")
colors = colorRampPalette(brewer.pal(9, "RdBu"))(100)
heatmap.2(correlation, main= "Heatmap correlación Spearman entre muestras", trace="none", col = colors , margin=c(15, 10), scale = "none")
heatmap(correlation, scale= "column", main= "Heatmap de correlación entre muestras DEG", margins = c(15,10))


#Correlación muestras Lambdacialotrina exclusivamente - aquí SIN L4 Y L9
###4. Correlación muestras Lambdacialotrina sin L4 y L9####


countTable = read.table("C:/Users/alejandro/Desktop/Segundo artículo Transcriptómica Aedes/After_bam-R/1_Conteos/Counts_RNAseq_subsample_2024_V2.txt", sep="\t", header=T, row.names = 1)

#Remover columnas de counTable
countTable$Chr <- NULL
countTable$Start <- NULL
countTable$End <- NULL
countTable$Strand <- NULL
countTable$Length <- NULL

#Mantener columnas de counTable de Columna 4 a 12
countTable <-countTable[,4:12]

countTable$L4_filtered.bam <- NULL
countTable$L9_filtered.bam <- NULL

colData = read.table("C:/Users/Alejandro/Desktop/Segundo artículo Transcriptómica Aedes/After_bam-R/Metadata/Factors_Metadata_especificos_subsample_V2.txt", sep="\t", header=T, row.names = 1)

#Mantener filas lambdacialotrina 4 a 12 
colData = colData[(4:12), , drop = F]

#Remover la fila 4 y la 9
colData <- colData[-4, , drop = FALSE]
colData <- colData[-8, , drop = FALSE]

#Confirmar nombres de columnas y filas
all(colnames(countTable) %in% rownames(colData))

all(colnames(countTable) == rownames(colData))

countTable[, "max"]= apply(countTable[, 1:ncol(countTable)], 1, max)
countTable=countTable[countTable[,ncol(countTable)]>32,] 
countTable= countTable [,-ncol(countTable)]

library(DESeq2)
dds<-DESeqDataSetFromMatrix(countData= countTable,colData= colData,design= ~ estado_susceptibilidad)
dds<-DESeq(dds)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

#Esta siguiente no tiene sentido con todas las muestras y los conteos normalizados
clusters <- hclust(dist(t(normalized_counts), method = "correlation"))
plot(clusters, hang = -1)
rect.hclust(clusters , k = 5, border = 4:4)
Dendogram <- as.dendrogram(clusters)
Dendogram <- color_branches(clusters, k = 3)
plot(Dendogram)
grid.draw(dendrogram(Dendogram, horiz = TRUE, leaflab = "none"))

###Para todos los genes
# For samples
correlation= cor(normalized_counts, method= "spearman")
colors = colorRampPalette(brewer.pal(9, "RdBu"))(100)
heatmap.2(correlation, main= "Heatmap correlación Spearman entre muestras", trace="none", col = colors , margin=c(15, 10), scale = "none")
heatmap(correlation, scale= "column", main= "Heatmap de correlación entre muestras", margins = c(15,10))




#

####5. Gráficas y demás con DEG´s con todos los genes sin modificaciones#####

countTable = read.table("C:/Users/alejandro/Desktop/Segundo artículo Transcriptómica Aedes/After_bam-R/1_Conteos/Counts_RNAseq_subsample_2024_V2.txt", sep="\t", header=T, row.names = 1)

colData = read.table("C:/Users/Alejandro/Desktop/Segundo artículo Transcriptómica Aedes/After_bam-R/Metadata/Factors_Metadata_especificos_subsample_V2.txt", sep="\t", header=T, row.names = 1)

#Remover columnas de counTable
countTable$Chr <- NULL
countTable$Start <- NULL
countTable$End <- NULL
countTable$Strand <- NULL
countTable$Length <- NULL

all(colnames(countTable) %in% rownames(colData))

all(colnames(countTable) == rownames(colData))

countTable[, "max"]= apply(countTable[, 1:ncol(countTable)], 1, max)
countTable=countTable[countTable[,ncol(countTable)]>32,] 
countTable= countTable [,-ncol(countTable)]

library(DESeq2)
dds<-DESeqDataSetFromMatrix(countData= countTable,colData= colData,design= ~ estado_susceptibilidad)
dds<-DESeq(dds)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

library(dplyr)
##Comparaciones con R25 y S75
res1<-results(dds,contrast=c("estado_susceptibilidad","L_resistente75","L_resistente25"))
Lambda_Res75vRes25= as.data.frame(res1)
Lambda_Res75vRes25= subset.data.frame(Lambda_Res75vRes25, padj<0.05 & abs(log2FoldChange)>1)
#write.table(as.data.frame(Lambda_Res75vRes25),file="Lambda_Resistente75vResistente25.csv")

res2<-results(dds,contrast=c("estado_susceptibilidad","L_resistente75","L_susceptible75"))
Lambda_Res75vSusc75= as.data.frame(res2)
Lambda_Res75vSusc75= subset.data.frame(Lambda_Res75vSusc75, padj<0.05 & abs(log2FoldChange)>1)
#write.table(as.data.frame(Lambda_Res75vSusc75),file="Lambda_Resistente75vSuseptible75.csv")

res3<-results(dds,contrast=c("estado_susceptibilidad","P_resistente75","P_resistente25"))
Perm_Res75vRes25= as.data.frame(res3)
Perm_Res75vRes25= subset.data.frame(Perm_Res75vRes25, padj<0.05 & abs(log2FoldChange)>1)
#write.table(as.data.frame(Perm_Res75vRes25),file="Perm_Resistente75vResistente25.csv")

res4<-results(dds,contrast=c("estado_susceptibilidad","P_resistente75","P_susceptible75"))
Perm_Res75vSusc75= as.data.frame(res4)
Perm_Res75vSusc75= subset.data.frame(Perm_Res75vSusc75, padj<0.05 & abs(log2FoldChange)>1)
#write.table(as.data.frame(Perm_Res75vSusc75),file="Perm_Resistente75vSusceptible75.csv")

###Tablas faltantes###
res5<-results(dds,contrast=c("estado_susceptibilidad","P_susceptible75","P_resistente25"))
Perm_Susc75vRes25= as.data.frame(res5)
Perm_Susc75vRes25= subset.data.frame(Perm_Susc75vRes25, padj<0.05 & abs(log2FoldChange)>1)
#write.table(as.data.frame(Perm_Susc75vRes25),file="Perm_Susceptible75vResistente25.csv")

res6<-results(dds,contrast=c("estado_susceptibilidad","L_susceptible75","L_resistente25"))
Lambda_Susc75vRes25= as.data.frame(res6)
Lambda_Susc75vRes25= subset.data.frame(Lambda_Susc75vRes25, padj<0.05 & abs(log2FoldChange)>1)
#write.table(as.data.frame(Lambda_Susc75vRes25),file="Lambda_Susceptible75vResistente25.csv")

###
genes_nameres1 = rownames(Lambda_Res75vRes25)
genes_nameres2 = rownames(Lambda_Res75vSusc75)
genes_nameres3 = rownames(Perm_Res75vRes25)
genes_nameres4 = rownames(Perm_Res75vSusc75)
genes_nameres5 = rownames(Perm_Susc75vRes25)
genes_nameres6 = rownames(Lambda_Susc75vRes25)

names_genes= c(genes_nameres1,genes_nameres2,genes_nameres3,genes_nameres4,genes_nameres5,genes_nameres6)
normalized_subset <- normalized_counts[rownames(normalized_counts) %in% names_genes, ]

#write.table(normalized_subset, file="normalized_counts_DEG.txt", sep="\t", quote=F)


###Correlación spearman con DEGs de comparaciones con solo tratados solo genes DEG

correlation= cor(normalized_subset, method= "spearman")
colors = colorRampPalette(brewer.pal(9, "RdBu"))(100)
heatmap.2(correlation, main= "Heatmap correlación spearman entre TODAS las muestras", trace="none", col = colors , margin=c(10, 10), scale = "none")

###Distancia Euclidiana###
clusters <- hclust(dist(t(normalized_subset), method = "euclidian"))
plot(clusters, hang = -1)
rect.hclust(clusters , k = 5, border = 4:4)
Dendogram <- as.dendrogram(clusters)
Dendogram <- color_branches(clusters, k = 7)
plot(Dendogram)
grid.draw(dendrogram(Dendogram, horiz = TRUE, leaflab = "none"))

# ================================== #
# == 5.Principal component analysis  #
# ================================== #

#Conteos normalizados 
all_Z=t(scale(t(normalized_counts)))

pca <- prcomp(t(all_Z))
pca.var1a <- pca$sdev^2
pca.var.per <- round(pca.var1a/sum(pca.var1a)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2],
                       Group=c("Susceptible", "Susceptible","Susceptible", "L_resistente25", "L_resistente25","L_resistente25", "L_susceptible75", "L_susceptible75", "L_susceptible75", "L_resistente75", "L_resistente75", "L_resistente75", "P_resistente25", "P_resistente25", "P_resistente25", "P_susceptible75", "P_susceptible75", "P_susceptible75", "P_resistente75", "P_resistente75", "P_resistente75"))
pca.data

ggplot(data=pca.data, aes(x=X, y=Y)) +
  geom_point(aes(colour= Group), size= 5 ) +
  geom_text(aes(label=Sample)) +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  ggtitle("PCA resistentes y susceptibles TODOS Norm. Counts.")

#Conteos CON SOLO DEG´S  
DEG_Z=t(scale(t(normalized_subset)))

pca <- prcomp(t(DEG_Z))
pca.var1a <- pca$sdev^2
pca.var.per <- round(pca.var1a/sum(pca.var1a)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2],
                       Group=c("Susceptible", "Susceptible","Susceptible", "L_resistente25", "L_resistente25","L_resistente25", "L_susceptible75", "L_susceptible75", "L_susceptible75", "L_resistente75", "L_resistente75", "L_resistente75", "P_resistente25", "P_resistente25", "P_resistente25", "P_susceptible75", "P_susceptible75", "P_susceptible75", "P_resistente75", "P_resistente75", "P_resistente75"))
pca.data

ggplot(data=pca.data, aes(x=X, y=Y)) +
  geom_point(aes(colour= Group), size= 5 ) +
  geom_text(aes(label=Sample)) +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  ggtitle("PCA resistentes y susceptibles TODAS las muestras DEG´S.")

###6. Gráficas y demás con DEG´s con todos los genes SIN CONTROLES susceptibles####

countTable = read.table("C:/Users/alejandro/Desktop/Segundo artículo Transcriptómica Aedes/After_bam-R/1_Conteos/Counts_RNAseq_subsample_2024_V2.txt", sep="\t", header=T, row.names = 1)

#Remover columnas de counTable
countTable <-countTable[,9:26]

colData = read.table("C:/Users/Alejandro/Desktop/Segundo artículo Transcriptómica Aedes/After_bam-R/Metadata/Factors_Metadata_especificos_subsample_V2.txt", sep="\t", header=T, row.names = 1)

colData = colData[(4:21), , drop = F]

all(colnames(countTable) %in% rownames(colData))

all(colnames(countTable) == rownames(colData))

countTable[, "max"]= apply(countTable[, 1:ncol(countTable)], 1, max)
countTable=countTable[countTable[,ncol(countTable)]>32,] 
countTable= countTable [,-ncol(countTable)]

library(DESeq2)
dds<-DESeqDataSetFromMatrix(countData= countTable,colData= colData,design= ~ estado_susceptibilidad)
dds<-DESeq(dds)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

library(dplyr)
##Comparaciones con R25 y S75
res1<-results(dds,contrast=c("estado_susceptibilidad","L_resistente75","L_resistente25"))
Lambda_Res75vRes25= as.data.frame(res1)
Lambda_Res75vRes25= subset.data.frame(Lambda_Res75vRes25, padj<0.05 & abs(log2FoldChange)>1)
#write.table(as.data.frame(Lambda_Res75vRes25),file="Lambda_Resistente75vResistente25.csv")

res2<-results(dds,contrast=c("estado_susceptibilidad","L_resistente75","L_susceptible75"))
Lambda_Res75vSusc75= as.data.frame(res2)
Lambda_Res75vSusc75= subset.data.frame(Lambda_Res75vSusc75, padj<0.05 & abs(log2FoldChange)>1)
#write.table(as.data.frame(Lambda_Res75vSusc75),file="Lambda_Resistente75vSuseptible75.csv")

res3<-results(dds,contrast=c("estado_susceptibilidad","P_resistente75","P_resistente25"))
Perm_Res75vRes25= as.data.frame(res3)
Perm_Res75vRes25= subset.data.frame(Perm_Res75vRes25, padj<0.05 & abs(log2FoldChange)>1)
#write.table(as.data.frame(Perm_Res75vRes25),file="Perm_Resistente75vResistente25.csv")

res4<-results(dds,contrast=c("estado_susceptibilidad","P_resistente75","P_susceptible75"))
Perm_Res75vSusc75= as.data.frame(res4)
Perm_Res75vSusc75= subset.data.frame(Perm_Res75vSusc75, padj<0.05 & abs(log2FoldChange)>1)
#write.table(as.data.frame(Perm_Res75vSusc75),file="Perm_Resistente75vSusceptible75.csv")

###Tablas faltantes###
res5<-results(dds,contrast=c("estado_susceptibilidad","P_susceptible75","P_resistente25"))
Perm_Susc75vRes25= as.data.frame(res5)
Perm_Susc75vRes25= subset.data.frame(Perm_Susc75vRes25, padj<0.05 & abs(log2FoldChange)>1)
#write.table(as.data.frame(Perm_Susc75vRes25),file="Perm_Susceptible75vResistente25.csv")

res6<-results(dds,contrast=c("estado_susceptibilidad","L_susceptible75","L_resistente25"))
Lambda_Susc75vRes25= as.data.frame(res6)
Lambda_Susc75vRes25= subset.data.frame(Lambda_Susc75vRes25, padj<0.05 & abs(log2FoldChange)>1)
#write.table(as.data.frame(Lambda_Susc75vRes25),file="Lambda_Susceptible75vResistente25.csv")

###
genes_nameres1 = rownames(Lambda_Res75vRes25)
genes_nameres2 = rownames(Lambda_Res75vSusc75)
genes_nameres3 = rownames(Perm_Res75vRes25)
genes_nameres4 = rownames(Perm_Res75vSusc75)
genes_nameres5 = rownames(Perm_Susc75vRes25)
genes_nameres6 = rownames(Lambda_Susc75vRes25)

names_genes= c(genes_nameres1,genes_nameres2,genes_nameres3,genes_nameres4,genes_nameres5,genes_nameres6)
normalized_subset <- normalized_counts[rownames(normalized_counts) %in% names_genes, ]

#write.table(normalized_subset, file="normalized_counts_DEG.txt", sep="\t", quote=F)


###Correlación spearman con DEGs de comparaciones con solo tratados solo genes DEG

correlation= cor(normalized_subset, method= "spearman")
colors = colorRampPalette(brewer.pal(9, "RdBu"))(100)
heatmap.2(correlation, main= "correlación spearman entre DEG´S muestras TRATADAS", trace="none", col = colors , margin=c(10, 10), scale = "none")

###Distancia Euclidiana###
clusters <- hclust(dist(t(normalized_subset), method = "euclidian"))
plot(clusters, hang = -1)
rect.hclust(clusters , k = 5, border = 4:4)
Dendogram <- as.dendrogram(clusters)
Dendogram <- color_branches(clusters, k = 6)
plot(Dendogram)
grid.draw(dendrogram(Dendogram, horiz = TRUE, leaflab = "none"))


# ================================== #
# == 5.Principal component analysis  #
# ================================== #

#Conteos normalizados 
DEG_Z=t(scale(t(normalized_subset)))

pca <- prcomp(t(DEG_Z))
pca.var1a <- pca$sdev^2
pca.var.per <- round(pca.var1a/sum(pca.var1a)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2],
                       Group=c("L_resistente25", "L_resistente25","L_resistente25", "L_susceptible75", "L_susceptible75", "L_susceptible75", "L_resistente75", "L_resistente75", "L_resistente75", "P_resistente25", "P_resistente25", "P_resistente25", "P_susceptible75", "P_susceptible75", "P_susceptible75", "P_resistente75", "P_resistente75", "P_resistente75"))
pca.data

ggplot(data=pca.data, aes(x=X, y=Y)) +
  geom_point(aes(colour= Group), size= 5 ) +
  geom_text(aes(label=Sample)) +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  ggtitle("PCA resistentes y susceptibles DEG´s solo Tratados.")
###7. Gráficas y demás con DEG´s con SOLO GENES LAMBDACIALOTRINA SIN MODIFICACIONES####
countTable = read.table("C:/Users/alejandro/Desktop/Segundo artículo Transcriptómica Aedes/After_bam-R/1_Conteos/Counts_RNAseq_subsample_2024_V2.txt", sep="\t", header=T, row.names = 1)

#Remover columnas de counTable
countTable <-countTable[,9:17]

colData = read.table("C:/Users/Alejandro/Desktop/Segundo artículo Transcriptómica Aedes/After_bam-R/Metadata/Factors_Metadata_especificos_subsample_V2.txt", sep="\t", header=T, row.names = 1)

colData = colData[(4:12), , drop = F]

all(colnames(countTable) %in% rownames(colData))

all(colnames(countTable) == rownames(colData))

countTable[, "max"]= apply(countTable[, 1:ncol(countTable)], 1, max)
countTable=countTable[countTable[,ncol(countTable)]>32,] 
countTable= countTable [,-ncol(countTable)]

library(DESeq2)
dds<-DESeqDataSetFromMatrix(countData= countTable,colData= colData,design= ~ estado_susceptibilidad)
dds<-DESeq(dds)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

library(dplyr)
##Comparaciones con R25 y S75
res1<-results(dds,contrast=c("estado_susceptibilidad","L_resistente75","L_resistente25"))
Lambda_Res75vRes25= as.data.frame(res1)
Lambda_Res75vRes25= subset.data.frame(Lambda_Res75vRes25, padj<0.05 & abs(log2FoldChange)>1)
#write.table(as.data.frame(Lambda_Res75vRes25),file="Lambda_Resistente75vResistente25.csv")

res2<-results(dds,contrast=c("estado_susceptibilidad","L_resistente75","L_susceptible75"))
Lambda_Res75vSusc75= as.data.frame(res2)
Lambda_Res75vSusc75= subset.data.frame(Lambda_Res75vSusc75, padj<0.05 & abs(log2FoldChange)>1)
#write.table(as.data.frame(Lambda_Res75vSusc75),file="Lambda_Resistente75vSuseptible75.csv")

res6<-results(dds,contrast=c("estado_susceptibilidad","L_susceptible75","L_resistente25"))
Lambda_Susc75vRes25= as.data.frame(res6)
Lambda_Susc75vRes25= subset.data.frame(Lambda_Susc75vRes25, padj<0.05 & abs(log2FoldChange)>1)
#write.table(as.data.frame(Lambda_Susc75vRes25),file="Lambda_Susceptible75vResistente25.csv")

###
genes_nameres1 = rownames(Lambda_Res75vRes25)
genes_nameres2 = rownames(Lambda_Res75vSusc75)
genes_nameres6 = rownames(Lambda_Susc75vRes25)

names_genes= c(genes_nameres1,genes_nameres2,genes_nameres6)
normalized_subset <- normalized_counts[rownames(normalized_counts) %in% names_genes, ]

#write.table(normalized_subset, file="normalized_counts_DEG.txt", sep="\t", quote=F)


###Correlación spearman con DEGs de comparaciones con solo tratados solo genes DEG

correlation= cor(normalized_subset, method= "spearman")
colors = colorRampPalette(brewer.pal(9, "RdBu"))(100)
heatmap.2(correlation, main= "correlación spearman entre DEG´S muestras TRATADAS", trace="none", col = colors , margin=c(10, 10), scale = "none")

###Distancia Euclidiana### ESTO NO DIO
clusters <- hclust(dist(t(normalized_subset), method = "euclidian"))
plot(clusters, hang = -1)
rect.hclust(clusters , k = 3, border = 4:6)
Dendogram <- as.dendrogram(clusters)
Dendogram <- color_branches(clusters, k = 3)
plot(Dendogram)

# ================================== #
# == 5.Principal component analysis  #
# ================================== #

#Conteos normalizados 
DEG_Z=t(scale(t(normalized_subset)))

pca <- prcomp(t(DEG_Z))
pca.var1a <- pca$sdev^2
pca.var.per <- round(pca.var1a/sum(pca.var1a)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2],
                       Group=c("L_resistente25", "L_resistente25","L_resistente25", "L_susceptible75", "L_susceptible75", "L_susceptible75", "L_resistente75", "L_resistente75", "L_resistente75"))
pca.data

ggplot(data=pca.data, aes(x=X, y=Y)) +
  geom_point(aes(colour= Group), size= 5 ) +
  geom_text(aes(label=Sample)) +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  ggtitle("PCA resistentes y susceptibles DEG´s solo lambdacialotrina")
###8. Gráficas y demás con DEG´s con SOLO GENES PERMETRINA SIN MODIFICACIONES####
countTable = read.table("C:/Users/alejandro/Desktop/Segundo artículo Transcriptómica Aedes/After_bam-R/1_Conteos/Counts_RNAseq_subsample_2024_V2.txt", sep="\t", header=T, row.names = 1)

#Remover columnas de counTable
countTable <-countTable[,18:26]

colData = read.table("C:/Users/Alejandro/Desktop/Segundo artículo Transcriptómica Aedes/After_bam-R/Metadata/Factors_Metadata_especificos_subsample_V2.txt", sep="\t", header=T, row.names = 1)

colData = colData[(13:21), , drop = F]

all(colnames(countTable) %in% rownames(colData))

all(colnames(countTable) == rownames(colData))

countTable[, "max"]= apply(countTable[, 1:ncol(countTable)], 1, max)
countTable=countTable[countTable[,ncol(countTable)]>32,] 
countTable= countTable [,-ncol(countTable)]

library(DESeq2)
dds<-DESeqDataSetFromMatrix(countData= countTable,colData= colData,design= ~ estado_susceptibilidad)
dds<-DESeq(dds)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

library(dplyr)
##Comparaciones con R25 y S75
res3<-results(dds,contrast=c("estado_susceptibilidad","P_resistente75","P_resistente25"))
Perm_Res75vRes25= as.data.frame(res3)
Perm_Res75vRes25= subset.data.frame(Perm_Res75vRes25, padj<0.05 & abs(log2FoldChange)>1)
#write.table(as.data.frame(Perm_Res75vRes25),file="Perm_Resistente75vResistente25.csv")

res4<-results(dds,contrast=c("estado_susceptibilidad","P_resistente75","P_susceptible75"))
Perm_Res75vSusc75= as.data.frame(res4)
Perm_Res75vSusc75= subset.data.frame(Perm_Res75vSusc75, padj<0.05 & abs(log2FoldChange)>1)
#write.table(as.data.frame(Perm_Res75vSusc75),file="Perm_Resistente75vSusceptible75.csv")

###Tablas faltantes###
res5<-results(dds,contrast=c("estado_susceptibilidad","P_susceptible75","P_resistente25"))
Perm_Susc75vRes25= as.data.frame(res5)
Perm_Susc75vRes25= subset.data.frame(Perm_Susc75vRes25, padj<0.05 & abs(log2FoldChange)>1)
#write.table(as.data.frame(Perm_Susc75vRes25),file="Perm_Susceptible75vResistente25.csv")

genes_nameres3 = rownames(Perm_Res75vRes25)
genes_nameres4 = rownames(Perm_Res75vSusc75)
genes_nameres5 = rownames(Perm_Susc75vRes25)

names_genes= c(genes_nameres3,genes_nameres4,genes_nameres5)
normalized_subset <- normalized_counts[rownames(normalized_counts) %in% names_genes, ]

#write.table(normalized_subset, file="normalized_counts_DEG.txt", sep="\t", quote=F)


###Correlación spearman con DEGs de comparaciones con solo tratados solo genes DEG

correlation= cor(normalized_subset, method= "spearman")
colors = colorRampPalette(brewer.pal(9, "RdBu"))(100)
heatmap.2(correlation, main= "correlación spearman entre DEG´S muestras TRATADAS", trace="none", col = colors , margin=c(10, 10), scale = "none")

###Distancia Euclidiana###
clusters <- hclust(dist(t(normalized_subset), method = "euclidian"))
plot(clusters, hang = -1)
rect.hclust(clusters , k = 5, border = 4:4)
Dendogram <- as.dendrogram(clusters)
Dendogram <- color_branches(clusters, k = 6)
plot(Dendogram)
grid.draw(dendrogram(Dendogram, horiz = TRUE, leaflab = "none"))

#Conteos normalizados 
DEG_Z=t(scale(t(normalized_subset)))

pca <- prcomp(t(DEG_Z))
pca.var1a <- pca$sdev^2
pca.var.per <- round(pca.var1a/sum(pca.var1a)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2],
                       Group=c("P_resistente25", "P_resistente25","P_resistente25", "P_susceptible75", "P_susceptible75", "P_susceptible75", "P_resistente75", "P_resistente75", "P_resistente75"))
pca.data

ggplot(data=pca.data, aes(x=X, y=Y)) +
  geom_point(aes(colour= Group), size= 5 ) +
  geom_text(aes(label=Sample)) +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  ggtitle("PCA resistentes y susceptibles DEG´s solo Permetrina")

###8.1 Gráficas y demás con DEG´s con SOLO GENES PERMETRINA SIN MUESTRA P3####
countTable = read.table("C:/Users/alejandro/Desktop/Segundo artículo Transcriptómica Aedes/After_bam-R/1_Conteos/Counts_RNAseq_subsample_2024_V2.txt", sep="\t", header=T, row.names = 1)

#Remover columnas de counTable
countTable <-countTable[,18:26]

#Remover columnas de counTable

countTable$P3_filtered.bam <- NULL

colData = read.table("C:/Users/Alejandro/Desktop/Segundo artículo Transcriptómica Aedes/After_bam-R/Metadata/Factors_Metadata_especificos_subsample_V2.txt", sep="\t", header=T, row.names = 1)

colData = colData[(13:21), , drop = F]
colData <- colData[-3, , drop = FALSE]

all(colnames(countTable) %in% rownames(colData))

all(colnames(countTable) == rownames(colData))

countTable[, "max"]= apply(countTable[, 1:ncol(countTable)], 1, max)
countTable=countTable[countTable[,ncol(countTable)]>32,] 
countTable= countTable [,-ncol(countTable)]

library(DESeq2)
dds<-DESeqDataSetFromMatrix(countData= countTable,colData= colData,design= ~ estado_susceptibilidad)
dds<-DESeq(dds)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

library(dplyr)
##Comparaciones con R25 y S75
res3<-results(dds,contrast=c("estado_susceptibilidad","P_resistente75","P_resistente25"))
Perm_Res75vRes25= as.data.frame(res3)
Perm_Res75vRes25= subset.data.frame(Perm_Res75vRes25, padj<0.05 & abs(log2FoldChange)>1)
#write.table(as.data.frame(Perm_Res75vRes25),file="Perm_Resistente75vResistente25.csv")

res4<-results(dds,contrast=c("estado_susceptibilidad","P_resistente75","P_susceptible75"))
Perm_Res75vSusc75= as.data.frame(res4)
Perm_Res75vSusc75= subset.data.frame(Perm_Res75vSusc75, padj<0.05 & abs(log2FoldChange)>1)
#write.table(as.data.frame(Perm_Res75vSusc75),file="Perm_Resistente75vSusceptible75.csv")

###Tablas faltantes###
res5<-results(dds,contrast=c("estado_susceptibilidad","P_susceptible75","P_resistente25"))
Perm_Susc75vRes25= as.data.frame(res5)
Perm_Susc75vRes25= subset.data.frame(Perm_Susc75vRes25, padj<0.05 & abs(log2FoldChange)>1)
#write.table(as.data.frame(Perm_Susc75vRes25),file="Perm_Susceptible75vResistente25.csv")

genes_nameres3 = rownames(Perm_Res75vRes25)
genes_nameres4 = rownames(Perm_Res75vSusc75)
genes_nameres5 = rownames(Perm_Susc75vRes25)

names_genes= c(genes_nameres3,genes_nameres4,genes_nameres5)
normalized_subset <- normalized_counts[rownames(normalized_counts) %in% names_genes, ]

#write.table(normalized_subset, file="normalized_counts_DEG.txt", sep="\t", quote=F)


###Correlación spearman con DEGs de comparaciones con solo tratados solo genes DEG

correlation= cor(normalized_subset, method= "spearman")
colors = colorRampPalette(brewer.pal(9, "RdBu"))(100)
heatmap.2(correlation, main= "correlación spearman entre DEG´S muestras TRATADAS", trace="none", col = colors , margin=c(10, 10), scale = "none")

###Distancia Euclidiana###
clusters <- hclust(dist(t(normalized_subset), method = "euclidian"))
plot(clusters, hang = -1)
rect.hclust(clusters , k = 5, border = 4:4)
Dendogram <- as.dendrogram(clusters)
Dendogram <- color_branches(clusters, k = 6)
plot(Dendogram)
grid.draw(dendrogram(Dendogram, horiz = TRUE, leaflab = "none"))

#Conteos normalizados 
DEG_Z=t(scale(t(normalized_subset)))

pca <- prcomp(t(DEG_Z))
pca.var1a <- pca$sdev^2
pca.var.per <- round(pca.var1a/sum(pca.var1a)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2],
                       Group=c("P_resistente25", "P_resistente25", "P_susceptible75", "P_susceptible75", "P_susceptible75", "P_resistente75", "P_resistente75", "P_resistente75"))
pca.data

ggplot(data=pca.data, aes(x=X, y=Y)) +
  geom_point(aes(colour= Group), size= 5 ) +
  geom_text(aes(label=Sample)) +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  ggtitle("PCA resistentes y susceptibles DEG´s solo Permetrina")




####9. Análisis de Transcriptoma definitivo con tratados con muestra P3 eliminada del análisis####
#Se hará el siguiente análisis utilizando los tratados y para efectos de mayor resolución el PCA se hará a partir de los datos individualizados de Permetrina con los DEG´s Permetrina y de Lambdacialotrina con los DEG´s Lambdacialotrina. 

countTable = read.table("C:/Users/alejandro/Desktop/Segundo artículo Transcriptómica Aedes/After_bam-R/1_Conteos/Counts_RNAseq_subsample_2024_V2.txt", sep="\t", header=T, row.names = 1)

#Remover columnas de counTable
countTable <-countTable[,9:26]
countTable$P3_filtered.bam <- NULL


colData = read.table("C:/Users/Alejandro/Desktop/Segundo artículo Transcriptómica Aedes/After_bam-R/Metadata/Factors_Metadata_especificos_subsample_V2.txt", sep="\t", header=T, row.names = 1)

#Mantener filas lambdacialotrina 4 a 12 
colData = colData[(4:21), , drop = F]
colData <- colData[-12, , drop = FALSE]


all(colnames(countTable) %in% rownames(colData))

all(colnames(countTable) == rownames(colData))

countTable[, "max"]= apply(countTable[, 1:ncol(countTable)], 1, max)
countTable=countTable[countTable[,ncol(countTable)]>32,] 
countTable= countTable [,-ncol(countTable)]

### guardar el universo de genes para GO
#gene_universe <-rownames(countTable)
#write.table(as.data.frame(gene_universe),file="gene_universe_countTable.csv")

library(DESeq2)
dds<-DESeqDataSetFromMatrix(countData= countTable,colData= colData,design= ~ estado_susceptibilidad)
dds<-DESeq(dds)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
plotDispEsts(dds)


library(dplyr)
##Comparaciones con R25 y S75
res1<-results(dds,contrast=c("estado_susceptibilidad","L_resistente75","L_resistente25"))
Lambda_Res75vRes25= as.data.frame(res1)
#write.table(as.data.frame(res1), file = "res1_HLRvsLLR_1.csv")
Lambda_Res75vRes25= subset.data.frame(Lambda_Res75vRes25, padj<0.05 & abs(log2FoldChange)>1)
#write.table(as.data.frame(Lambda_Res75vRes25),file="Lambda_Resistente75vResistente25.csv")

res2<-results(dds,contrast=c("estado_susceptibilidad","L_resistente75","L_susceptible75"))
Lambda_Res75vSusc75= as.data.frame(res2)
#write.table(as.data.frame(res2), file = "res2_HLRvsHLS.csv")
Lambda_Res75vSusc75= subset.data.frame(Lambda_Res75vSusc75, padj<0.05 & abs(log2FoldChange)>1)
#write.table(as.data.frame(Lambda_Res75vSusc75),file="Lambda_Resistente75vSuseptible75.csv")

res3<-results(dds,contrast=c("estado_susceptibilidad","P_resistente75","P_resistente25"))
Perm_Res75vRes25= as.data.frame(res3)
#write.table(as.data.frame(res3),file="res3_HPRvsLPR.csv")
Perm_Res75vRes25= subset.data.frame(Perm_Res75vRes25, padj<0.05 & abs(log2FoldChange)>1)
#write.table(as.data.frame(Perm_Res75vRes25),file="Perm_Resistente75vResistente25.csv")

res4<-results(dds,contrast=c("estado_susceptibilidad","P_resistente75","P_susceptible75"))
Perm_Res75vSusc75= as.data.frame(res4)
#write.table(as.data.frame(res4),file="res4_HPRvsHPS.csv")
Perm_Res75vSusc75= subset.data.frame(Perm_Res75vSusc75, padj<0.05 & abs(log2FoldChange)>1)
#write.table(as.data.frame(Perm_Res75vSusc75),file="Perm_Resistente75vSusceptible75.csv")

###Tablas faltantes###
res5<-results(dds,contrast=c("estado_susceptibilidad","P_susceptible75","P_resistente25"))
Perm_Susc75vRes25= as.data.frame(res5)
#write.table(as.data.frame(res5),file="res5_HPSvsLPR.csv")
Perm_Susc75vRes25= subset.data.frame(Perm_Susc75vRes25, padj<0.05 & abs(log2FoldChange)>1)
#write.table(as.data.frame(Perm_Susc75vRes25),file="Perm_Susceptible75vResistente25.csv")

res6<-results(dds,contrast=c("estado_susceptibilidad","L_susceptible75","L_resistente25"))
Lambda_Susc75vRes25= as.data.frame(res6)
#write.table(as.data.frame(res6),file="res6_HLSvsLPS.csv")
Lambda_Susc75vRes25= subset.data.frame(Lambda_Susc75vRes25, padj<0.05 & abs(log2FoldChange)>1)
#write.table(as.data.frame(Lambda_Susc75vRes25),file="Lambda_Susceptible75vResistente25.csv")



###
genes_nameres1 = rownames(Lambda_Res75vRes25)
genes_nameres2 = rownames(Lambda_Res75vSusc75)
genes_nameres3 = rownames(Perm_Res75vRes25)
genes_nameres4 = rownames(Perm_Res75vSusc75)
genes_nameres5 = rownames(Perm_Susc75vRes25)
genes_nameres6 = rownames(Lambda_Susc75vRes25)

names_genes= c(genes_nameres1,genes_nameres2,genes_nameres3,genes_nameres4,genes_nameres5,genes_nameres6)
normalized_subset <- normalized_counts[rownames(normalized_counts) %in% names_genes, ]
normalized_subset_Permetrina <- normalized_subset[,10:17]
normalized_subset_Lambdacialotrina <- normalized_subset[,1:9]

#write.table(normalized_subset, file="normalized_counts_DEG.txt", sep="\t", quote=F)

##Correlación spearman y PCA con DEGs de comoparaciones con solo Permetrina##

#Conteos normalizados 
DEG_Z_Permetrina=t(scale(t(normalized_subset_Permetrina)))

pca <- prcomp(t(DEG_Z_Permetrina))
pca.var1a <- pca$sdev^2
pca.var.per <- round(pca.var1a/sum(pca.var1a)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2],
                       Group=c("P_resistente25", "P_resistente25", "P_susceptible75", "P_susceptible75", "P_susceptible75", "P_resistente75", "P_resistente75", "P_resistente75"))
pca.data

ggplot(data=pca.data, aes(x=X, y=Y)) +
  geom_point(aes(colour= Group), size= 5 ) +
  geom_text(aes(label=Sample)) +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  ggtitle("PCA  DEG´s Permetrina")

###Correlación spearman con DEGs de comparaciones con solo tratados solo genes DEG

correlation= cor(normalized_subset_Permetrina, method= "spearman")
colors = colorRampPalette(brewer.pal(9, "RdBu"))(100)
heatmap.2(correlation, main= "Correlación spearman entre DEG´S Permetrina", trace="none", col = colors , margin=c(10, 10), scale = "none")

###Distancia Euclidiana### Da horrible ni la hice

##Correlación spearman y PCA con DEGs de comoparaciones con solo Lambdacialotrina##

#Conteos normalizados 
DEG_Z_Lambdacialotrina=t(scale(t(normalized_subset_Lambdacialotrina)))

pca <- prcomp(t(DEG_Z_Lambdacialotrina))
pca.var1a <- pca$sdev^2
pca.var.per <- round(pca.var1a/sum(pca.var1a)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2],
                       Group=c("L_resistente25", "L_resistente25","L_resistente25", "L_susceptible75", "L_susceptible75", "L_susceptible75", "L_resistente75", "L_resistente75", "L_resistente75"))
pca.data

ggplot(data=pca.data, aes(x=X, y=Y)) +
  geom_point(aes(colour= Group), size= 5 ) +
  geom_text(aes(label=Sample)) +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  ggtitle("PCA resistentes y susceptibles DEG´s solo Tratados.")

###
DEG_Z_PCA=t(scale(t(normalized_subset)))

pca <- prcomp(t(DEG_Z_PCA))
pca.var1a <- pca$sdev^2
pca.var.per <- round(pca.var1a/sum(pca.var1a)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2],
                       Group=c("L_resistente25", "L_resistente25","L_resistente25", "L_susceptible75", "L_susceptible75", "L_susceptible75", "L_resistente75", "L_resistente75", "L_resistente75","P_resistente25", "P_resistente25", "P_susceptible75", "P_susceptible75", "P_susceptible75", "P_resistente75", "P_resistente75", "P_resistente75"))
pca.data

ggplot(data=pca.data, aes(x=X, y=Y)) +
  geom_point(aes(colour= Group), size= 5 ) +
  geom_text(aes(label=Sample)) +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  ggtitle("PCA resistentes y susceptibles DEG´s solo Tratados.")

#write.table(as.data.frame(pca.data), file = "PCA.data_DEGS_Tratamientos.csv")

# ===================#
# == HEATMAP CORREGULACION#
# ===================#

DE_Z= as.matrix(DE_Z)

hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
hmcol2 = colorRampPalette(brewer.pal(11, "RdBu"))(100)
heatmap.2(DE_Z, col = hmcol2, trace="none", margin=c(10, 6))
heatmap.2(DE_Z, trace = "none",  col = hmcol2, labRow = F )


all_Z2= as.matrix(all_Z)

hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
hmcol2 = colorRampPalette(brewer.pal(11, "RdBu"))(100)
heatmap.2(all_Z2, col = hmcol2, trace="none", margin=c(10, 6))
heatmap.2(all_Z2, trace = "none",  col = hmcol2, labRow = F )

####10. Análisis de Transcriptoma definitivo con TRATADOS Y SUSCEPTIBLES con muestra P3 eliminada del análisis####
#Se hará el siguiente análisis utilizando los tratados (Sin la muestra P3) y los susceptibles el valor p ajustado será <0.01 y log2FC >2 

countTable = read.table("C:/Users/alejandro/Desktop/Segundo artículo Transcriptómica Aedes/After_bam-R/1_Conteos/Counts_RNAseq_subsample_2024_V2.txt", sep="\t", header=T, row.names = 1)

#Remover columnas de counTable
countTable <-countTable[,6:26]
countTable$P3_filtered.bam <- NULL

colData = read.table("C:/Users/Alejandro/Desktop/Segundo artículo Transcriptómica Aedes/After_bam-R/Metadata/Factors_Metadata_especificos_subsample_V2.txt", sep="\t", header=T, row.names = 1)

#Mantener filas lambdacialotrina 4 a 12 
colData <- colData[-15, , drop = FALSE]


all(colnames(countTable) %in% rownames(colData))

all(colnames(countTable) == rownames(colData))

countTable[, "max"]= apply(countTable[, 1:ncol(countTable)], 1, max)
countTable=countTable[countTable[,ncol(countTable)]>32,] 
countTable= countTable [,-ncol(countTable)]

library(DESeq2)
dds<-DESeqDataSetFromMatrix(countData= countTable,colData= colData,design= ~ estado_susceptibilidad)
dds<-DESeq(dds)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

colnames(normalized_counts) <- c("S1", "S2", "S3","LLR-1", "LLR-2", "LLR-3", "HLS-1", "HLS-2", "HLS-3", "HLR-1", "HLR-2", "HLR-3", "LPR-1", "LPR-2", "HPS-1", "HPS-2", "HPS-3", "HPR-1", "HPR-2", "HPR-3") 

write.table(as.data.frame(normalized_counts),file="normalized_counts.csv")
getwd()

library(dplyr)
##Comparaciones con los susceptibles##
##Lambdacialotrina
res5<-results(dds,contrast=c("estado_susceptibilidad","L_resistente75","Susceptible"))
Lambda_Res75v_Susc= as.data.frame(res5)
Lambda_Res75v_Susc= subset.data.frame(Lambda_Res75v_Susc, padj<0.01 & abs(log2FoldChange)>2)
#write.table(as.data.frame(Lambda_Res75v_Susc),file="L_Resistente75v_Susceptible.csv")

res6<-results(dds,contrast=c("estado_susceptibilidad","L_susceptible75","Susceptible"))
Lambda_Susc75v_Susc= as.data.frame(res6)
Lambda_Susc75v_Susc= subset.data.frame(Lambda_Susc75v_Susc, padj<0.01 & abs(log2FoldChange)>2)
#write.table(as.data.frame(Lambda_Susc75v_Susc),file="L_Susceptible75v_Susceptible.csv")

res7<-results(dds,contrast=c("estado_susceptibilidad","L_resistente25","Susceptible"))
Lambda_Res25v_Susc= as.data.frame(res7)
Lambda_Res25v_Susc= subset.data.frame(Lambda_Res25v_Susc, padj<0.01 & abs(log2FoldChange)>2)
#write.table(as.data.frame(Lambda_Res25v_Susc),file="L_Resistente25v_Susceptible.csv")

##Permetrina
res8<-results(dds,contrast=c("estado_susceptibilidad","P_resistente75","Susceptible"))
Perm_Res75v_Susc= as.data.frame(res8)
write.table(as.data.frame(res8), file= "res8_")
Perm_Res75v_Susc= subset.data.frame(Perm_Res75v_Susc, padj<0.01 & abs(log2FoldChange)>2)
#write.table(as.data.frame(Perm_Res75v_Susc),file="P_Resistente75v_Susceptible.csv")

res9<-results(dds,contrast=c("estado_susceptibilidad","P_susceptible75","Susceptible"))
Perm_Susc75v_Susc= as.data.frame(res9)
Perm_Susc75v_Susc= subset.data.frame(Perm_Susc75v_Susc, padj<0.01 & abs(log2FoldChange)>2)
#write.table(as.data.frame(Perm_Susc75v_Susc),file="P_Susceptible75v_Susceptible.csv")

res10<-results(dds,contrast=c("estado_susceptibilidad","P_resistente25","Susceptible"))
Perm_Res25v_Susc= as.data.frame(res10)
Perm_Res25v_Susc= subset.data.frame(Perm_Res25v_Susc, padj<0.01 & abs(log2FoldChange)>2)
#write.table(as.data.frame(Perm_Res25v_Susc),file="P_Resistente25v_Susceptible.csv")

res12<-results(dds,contrast=c("estado_susceptibilidad","P_resistente75","P_resistente25"))
Perm_Res75v_Res25= as.data.frame(res12)
write.table(as.data.frame(res12), file= "res12_HPRvLPR.csv")
Perm_Res75v_Susc= subset.data.frame(Perm_Res75v_Susc, padj<0.01 & abs(log2FoldChange)>2)
#write.table(as.data.frame(Perm_Res75v_Susc),file="P_Resistente75v_Susceptible.csv")

###
genes_nameres5 = rownames(Lambda_Res75v_Susc)
genes_nameres6 = rownames(Lambda_Susc75v_Susc)
genes_nameres7 = rownames(Lambda_Res25v_Susc)
genes_nameres8 = rownames(Perm_Res75v_Susc)
genes_nameres9 = rownames(Perm_Susc75v_Susc)
genes_nameres10 = rownames(Perm_Res25v_Susc)

names_genes= c(genes_nameres5,genes_nameres6,genes_nameres7,genes_nameres8,genes_nameres9,genes_nameres10)

normalized_subset <- normalized_counts[rownames(normalized_counts) %in% names_genes, ]

#write.table(normalized_subset, file="normalized_counts_DEG.txt", sep="\t", quote=F)


#Conteos normalizados PCA

all_Z=t(scale(t(normalized_counts)))

pca <- prcomp(t(all_Z))
pca.var1a <- pca$sdev^2
pca.var.per <- round(pca.var1a/sum(pca.var1a)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2],
                       Group=c("Susceptible", "Susceptible","Susceptible", "L_resistente25", "L_resistente25","L_resistente25", "L_susceptible75", "L_susceptible75", "L_susceptible75", "L_resistente75", "L_resistente75", "L_resistente75", "P_resistente25", "P_resistente25", "P_susceptible75", "P_susceptible75", "P_susceptible75", "P_resistente75", "P_resistente75", "P_resistente75"))
pca.data

ggplot(data=pca.data, aes(x=X, y=Y)) +
  geom_point(aes(colour= Group), size= 5 ) +
  geom_text(aes(label=Sample)) +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  ggtitle("PCA resistentes y susceptibles Todos los GENES")

###Conteos PCA DEGS
DEG_Z=t(scale(t(normalized_subset)))

pca <- prcomp(t(DEG_Z))
pca.var1a <- pca$sdev^2
pca.var.per <- round(pca.var1a/sum(pca.var1a)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2],
                       Group=c("Susceptible", "Susceptible","Susceptible", "L_resistente25", "L_resistente25","L_resistente25", "L_susceptible75", "L_susceptible75", "L_susceptible75", "L_resistente75", "L_resistente75", "L_resistente75", "P_resistente25", "P_resistente25", "P_susceptible75", "P_susceptible75", "P_susceptible75", "P_resistente75", "P_resistente75", "P_resistente75"))
pca.data

ggplot(data=pca.data, aes(x=X, y=Y)) +
  geom_point(aes(colour= Group), size= 5 ) +
  geom_text(aes(label=Sample)) +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  ggtitle("PCA resistentes y susceptibles DEG´s")

###Correlación spearman con DEGs de comparaciones TRATADOS Y SUSCEPTIBLES genes DEG 

correlation= cor(normalized_subset, method= "spearman")
colors = colorRampPalette(brewer.pal(9, "RdBu"))(100)
heatmap.2(correlation, main= "Correlación spearman entre DEG´S Permetrina", trace="none", col = colors , margin=c(10, 10), scale = "none")

####11. Union de Anotaciones con tabla con log2FC para MA plot en graphpad####

###                                                                     ###
### Genes HLR (High Lambda Resistant) v Suscebtible(SOLO LOS DEG´S q tienen anotacion )
###                                                                     ###

###Se crea una tabla con los identificadores de los genes y se guarda como TSV y con el siguiente código se lee. 
genes_HLRvS = read.table("C:/Users/Alejandro/Desktop/Segundo artículo Transcriptómica Aedes/After_bam-R/2. Codigo R/RNAseq_definitivo_Ttm-S/2. High-R Lambda v Susceptible/Genes_anotados_Resistenve75vSusceptible.txt", sep="\t", header=T)

### Genes totales de Genes_HPRvS_anotados_log2FC_meanbase_PARAR: Esta corresponde a otra tabla con todos los genes de la comparación HPRvS y en una segunda columna se ponen los log2FC
ALL_genes_HLRvS = read.table("C:/Users/Alejandro/Desktop/Segundo artículo Transcriptómica Aedes/After_bam-R/2. Codigo R/RNAseq_definitivo_Ttm-S/2. High-R Lambda v Susceptible/Genes_HLR_paraR.txt", sep="\t", header=T)

###Combinar ambas listas
merged_genes_HLRvS = merge (genes_HPRvS, ALL_genes_HPRvS, by ="genes_HPRvS", all.x = T)
###Guardar la lista
#install.packages("openxlsx")
library(openxlsx)

write.xlsx(merged_genes_HPRvS, "C:/Users/Alejandro/Desktop/Segundo artículo Transcriptómica Aedes/After_bam-R/2. Codigo R/RNAseq_definitivo_Ttm-S/2. High-R Lambda v Susceptible", row.names = FALSE)

###---------------------------------------------------------------------###
### Genes HPR (High Permet Resistant) v Suscebtible(SOLO LOS DEG´S q tienen anotacion )
###---------------------------------------------------------------------###

###Se crea una tabla con los identificadores de los genes y se guarda como TSV y con el siguiente código se lee. 
genes_HPRvS = read.table("C:/Users/Alejandro/Desktop/Segundo artículo Transcriptómica Aedes/After_bam-R/2. Codigo R/RNAseq_definitivo_Ttm-S/5. High-R Permet v Susceptible/Genes_Anotaciones_P_Resistente75v_Susceptible.txt", sep="\t", header=T)

### Genes totales de Genes_HPRvS_anotados_log2FC_meanbase_PARAR: Esta corresponde a otra tabla con todos los genes de la comparación HPRvS y en una segunda columna se ponen los log2FC
ALL_genes_HPRvS = read.table("C:/Users/Alejandro/Desktop/Segundo artículo Transcriptómica Aedes/After_bam-R/2. Codigo R/RNAseq_definitivo_Ttm-S/5. High-R Permet v Susceptible/Genes_HPRvS_anotados_log2FC_meanbase_PARAR.txt", sep="\t", header=T)

###Combinar ambas listas
merged_genes_HPRvS = merge (genes_HPRvS, ALL_genes_HPRvS, by ="genes_HPRvS", all.x = T)
###Guardar la lista
#install.packages("openxlsx")
library(openxlsx)

write.xlsx(merged_genes_HPRvS, "C:/Users/Alejandro/Desktop/Segundo artículo Transcriptómica Aedes/After_bam-R/2. Codigo R/RNAseq_definitivo_Ttm-S/5. High-R Permet v Susceptible/merged_genes_HPRvS.xlsx", row.names = FALSE)

