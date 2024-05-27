
library(Seurat)
library(monocle)
library(dplyr)
library(RColorBrewer)

DefaultAssay(scRNA) ='RNA'
data <- as(as.matrix(scRNA@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNA@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocds <- newCellDataSet(data,
								  phenoData = pd,
								  featureData = fd,
								  lowerDetectionLimit = 0.5,
								  expressionFamily = negbinomial.size())
print("format data done , filter select genes ")
#pData(monocds)$Cluster<-as.factor(pData(monocds)$celltype) 
pData(monocds)['Cluster']=scRNA@active.ident	
monocds <- estimateSizeFactors(monocds)
monocds <- estimateDispersions(monocds)
monocds <- detectGenes(monocds, min_expr = 0.1)
print(head(fData(monocds)))
expressed_genes <- row.names(subset(fData(monocds), num_cells_expressed >= 10)) # nolint
monocds <- monocds[expressed_genes, ]
disp_table <- dispersionTable(monocds)
express_genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
dir=paste0('monocle_output')
trajectory_cluster_dir=paste('./',dir,sep="")
monocds <- setOrderingFilter(monocds, express_genes)
monocds <- reduceDimension(
monocds,
max_components = 2,
method = "DDRTree")
monocds <- orderCells(monocds)

