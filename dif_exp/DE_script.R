### Taller de RNA-seq, Expresi?n diferencial con DESeq2

## Instalaci?n de paquetes desde Bioconductor

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("edgeR")

## Cargamos los paquetes a utilizar

library(DESeq2)
library(edgeR)

## Cargamos la matriz de conteos de RSEM.

data= read.table("ruta_del_archivo/RSEM_count_matrix.txt")

## Reordenamos las columnas para reflejar los tratamientos Contaminado vs No contaminado

col_ordering = c(4,5,6,1,2,3)
rnaseqMatrix = data[,col_ordering]
rnaseqMatrix = round(rnaseqMatrix)

## Filtramos aquellos genes que tengan una baja expresi?n en general.
rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 1) >= 2,]

## Formateamos las condiciones.
conditions = data.frame(conditions=factor(c(rep("Polluted", 3), rep("Unpolluted", 3))))
rownames(conditions) = colnames(rnaseqMatrix)

## Generamos un dataset del tipo DESeq2 a partir de nuestra matriz de conteo.
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = rnaseqMatrix,
  colData = conditions,
  design = ~ conditions)

## Generamos el modelo de expresion.

dds = DESeq(ddsFullCountTable)

## Realizamos el analisis de expresion diferencial
## Se entender?n como genes diferencialmente expresados aquellos con un FDR 
## menor a 0.05 y un valor absoluto de Log2FoldChange mayor a 1.

contrast=c("conditions","Polluted","Unpolluted")
res = results(dds, contrast, alpha=0.05, lfcThreshold=1)

## Veamos un resumen de los resultados
summary(res)

sum(res$padj < 0.05, na.rm=TRUE)


## Volcano plot para observar la distribucion de la expresi?n de genes.
DESeq2::plotMA(res)

## graficamos la expresi?n de un solo gen para cada individuo en los dos grupos.
plotCounts(dds, gene=which.min(res$padj), intgroup="conditions")

## Calculamos estadisticas para generar una tabla final con todos nuestros resultados.
baseMean_polluted <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "Polluted"])
baseMean_unpolluted <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "Unpolluted"])
res = cbind(baseMean_polluted, baseMean_unpolluted, as.data.frame(res))
res = cbind(sample_polluted="Polluted", sample_unpolluted="Unpolluted", as.data.frame(res))
res$padj[is.na(res$padj)]  <- 1
res = as.data.frame(res[order(res$pvalue),])


## Podemos exportar esta tabla si deseamos para analisis subsecuentes.

write.table(res, file='Differential_Exp_DESeq2_results', sep='	', quote=FALSE)
write.table(rnaseqMatrix, file='Normalized_count_matrix', sep='	', quote=FALSE)

