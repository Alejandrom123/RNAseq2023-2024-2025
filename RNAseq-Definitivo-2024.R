#Definitive script 

#Packages
library("ggplot2")
library("gridExtra")
library("edgeR")
library("UpSetR")
library("dplyr")
library(gplots)
library(ggplot2)
library(RColorBrewer)
library(dendextend)

####1. Transcriptome Analysis with ONLY treated (Sample p3 eliminated)####
#The next analysis is performed with only treated samples. For improving resolution and to see differences, PCA is based ONLY on permethrin and ONLY on lambda-cyhalothrin DEG.

#Count table
countTable = read.table("C:/Users/alejandro/Desktop/Segundo artículo Transcriptómica Aedes/After_bam-R/1_Conteos/Counts_RNAseq_subsample_2024_V2.txt", sep="\t", header=T, row.names = 1)

#Remover columnas de counTable
countTable <-countTable[,9:26]
countTable$P3_filtered.bam <- NULL


colData = read.table("C:/Users/Alejandro/Desktop/Segundo artículo Transcriptómica Aedes/After_bam-R/Metadata/Factors_Metadata_especificos_subsample_V2.txt", sep="\t", header=T, row.names = 1)

#Mantener filas lambdacialotrina 4 a 12 
colData = colData[(4:21), , drop = F]
colData <- colData[-12, , drop = FALSE]

#Check if columns and rows names match
all(colnames(countTable) %in% rownames(colData))

all(colnames(countTable) == rownames(colData))

#Removing genes with read counts <32
countTable[, "max"]= apply(countTable[, 1:ncol(countTable)], 1, max)
countTable=countTable[countTable[,ncol(countTable)]>32,] 
countTable= countTable [,-ncol(countTable)]
#write.table(as.data.frame(countTable),file="countTableII.csv")

### guardar el universo de genes para GO
#gene_universe <-rownames(countTable)
#write.table(as.data.frame(gene_universe),file="gene_universe_countTable.csv")

#Perform the Differential Expression analysis
library(DESeq2)
dds<-DESeqDataSetFromMatrix(countData= countTable,colData= colData,design= ~ estado_susceptibilidad)
dds<-DESeq(dds)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
#plotDispEsts(dds)


library(dplyr)
#Lambda
##HLR75 vs LL25
res1<-results(dds,contrast=c("estado_susceptibilidad","L_resistente75","L_resistente25"))
Lambda_Res75vRes25= as.data.frame(res1)
#write.table(as.data.frame(res1), file = "res1_HLRvsLLR_1.csv")
Lambda_Res75vRes25= subset.data.frame(Lambda_Res75vRes25, padj<0.05 & abs(log2FoldChange)>1)
#write.table(as.data.frame(Lambda_Res75vRes25),file="Lambda_Resistente75vResistente25.csv")

#HLR75 vs HLS75
res2<-results(dds,contrast=c("estado_susceptibilidad","L_resistente75","L_susceptible75"))
Lambda_Res75vSusc75= as.data.frame(res2)
#write.table(as.data.frame(res2), file = "res2_HLRvsHLS.csv")
Lambda_Res75vSusc75= subset.data.frame(Lambda_Res75vSusc75, padj<0.05 & abs(log2FoldChange)>1)
write.table(as.data.frame(Lambda_Res75vSusc75),file="Lambda_Resistente75vSuseptible75.csv")

#Permet
#HPR75 vs LPR25
res3<-results(dds,contrast=c("estado_susceptibilidad","P_resistente75","P_resistente25"))
Perm_Res75vRes25= as.data.frame(res3)
#write.table(as.data.frame(res3),file="res3_HPRvsLPR.csv")
Perm_Res75vRes25= subset.data.frame(Perm_Res75vRes25, padj<0.05 & abs(log2FoldChange)>1)
write.table(as.data.frame(Perm_Res75vRes25),file="Perm_Resistente75vResistente25.csv")

#HPR75 vs HPS75
res4<-results(dds,contrast=c("estado_susceptibilidad","P_resistente75","P_susceptible75"))
Perm_Res75vSusc75= as.data.frame(res4)
#write.table(as.data.frame(res4),file="res4_HPRvsHPS.csv")
Perm_Res75vSusc75= subset.data.frame(Perm_Res75vSusc75, padj<0.05 & abs(log2FoldChange)>1)
write.table(as.data.frame(Perm_Res75vSusc75),file="Perm_Resistente75vSusceptible75.csv")


###comparations Susceptible 75 vs resistent 25###
#Permethrin
res5<-results(dds,contrast=c("estado_susceptibilidad","P_susceptible75","P_resistente25"))
Perm_Susc75vRes25= as.data.frame(res5)
#write.table(as.data.frame(res5),file="res5_HPSvsLPR.csv")
Perm_Susc75vRes25= subset.data.frame(Perm_Susc75vRes25, padj<0.05 & abs(log2FoldChange)>1)
write.table(as.data.frame(Perm_Susc75vRes25),file="Perm_Susceptible75vResistente25.csv")

#Lambda-cyhalothrin
res6<-results(dds,contrast=c("estado_susceptibilidad","L_susceptible75","L_resistente25"))
Lambda_Susc75vRes25= as.data.frame(res6)
#write.table(as.data.frame(res6),file="res6_HLSvsLPS.csv")
Lambda_Susc75vRes25= subset.data.frame(Lambda_Susc75vRes25, padj<0.05 & abs(log2FoldChange)>1)
write.table(as.data.frame(Lambda_Susc75vRes25),file="Lambda_Susceptible75vResistente25.csv")

###Comparison between insecticides
res7<-results(dds,contrast=c("estado_susceptibilidad","P_resistente75","L_resistente75"))
Lambda_Res75vPerm_Res75= as.data.frame(res7)
#write.table(as.data.frame(res7),file="Lambda_Res75vPerm_Res75.csv")
Lambda_Res75vPerm_Res75= subset.data.frame(Lambda_Res75vPerm_Res75, padj<0.05 & abs(log2FoldChange)>1)
#write.table(as.data.frame(Lambda_Res75vPerm_Res75),file="Lambda_Res75vPerm_Res75_0.05.csv")
###

res8<-results(dds,contrast=c("estado_susceptibilidad","P_resistente25","L_resistente25"))
Perm_Res25vLambda_Res25= as.data.frame(res7)
#write.table(as.data.frame(res7),file="Perm_Res25vLambda_Res25.csv")
Perm_Res25vLambda_Res25= subset.data.frame(Perm_Res25vLambda_Res25, padj<0.05 & abs(log2FoldChange)>1)
#write.table(as.data.frame(Perm_Res25vLambda_Res25),file="Perm_Res25vLambda_Res25_0.05.csv")

###Obtaining row names of each of the comparations
genes_nameres1 = rownames(Lambda_Res75vRes25)
genes_nameres2 = rownames(Lambda_Res75vSusc75)
genes_nameres3 = rownames(Perm_Res75vRes25)
genes_nameres4 = rownames(Perm_Res75vSusc75)
genes_nameres5 = rownames(Perm_Susc75vRes25)
genes_nameres6 = rownames(Lambda_Susc75vRes25)

#Merge of the rownames to obtain subsets of each genes and each insecticide
names_genes= c(genes_nameres1,genes_nameres2,genes_nameres3,genes_nameres4,genes_nameres5,genes_nameres6)
normalized_subset <- normalized_counts[rownames(normalized_counts) %in% names_genes, ]
normalized_subset_Permetrina <- normalized_subset[,10:17]
normalized_subset_Lambdacialotrina <- normalized_subset[,1:9]

#write.table(normalized_subset, file="normalized_counts_DEG.txt", sep="\t", quote=F)

##Correlación spearman y PCA con DEGs de comoparaciones con solo Permetrina##

#PCA of Normalized counts of Permethrine
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
heatmap.2(correlation, main= "Spearman DEG´s Permthrin", trace="none", col = colors , margin=c(10, 10), scale = "none")

###Distancia Euclidiana### Da horrible ni la hice

##Correlación spearman y PCA con DEGs de comoparaciones con solo Lambdacialotrina##

#Normalized counts lambda-cyhalothrin
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

###Correlation lambda-cyhalothrin 
correlation= cor(normalized_subset_Lambdacialotrina, method= "spearman")
colors = colorRampPalette(brewer.pal(9, "RdBu"))(100)
heatmap.2(correlation, main= "Spearman DEG´s Lambda", trace="none", col = colors , margin=c(10, 10), scale = "none")

###All treated
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
  ggtitle("PCA DEG´s solo Tratados.")

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



####2. Transcriptome Analysis with Treated and Susceptible (Sample p3 eliminated) FDR <0.01 and log2FC####
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

#write.table(as.data.frame(normalized_counts),file="normalized_counts.csv")
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
#write.table(as.data.frame(res12), file= "res12_HPRvLPR.csv")
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


####3. Union de Anotaciones con tabla con log2FC para MA plot en graphpad####

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

