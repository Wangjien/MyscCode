
library(Seurat)
library(monocle)
library(dplyr)
library(RColorBrewer)
library(qs)

# 读入数据
setwd('/root/wangje/Project/刘老师/NK_T/new_Result/Data/001_cytotrace2_result')
scRNA <- qread('CD4_celltype2聚类结果.qs')
Idents(scRNA) <- scRNA$celltype2
outfile_Prefix <- ''
# 设置Idents
DefaultAssay(scRNA) ='RNA'
message(paste0(Sys.time(),'构建Monocle2分析对象......'))	
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
# 过滤低质量的基因
message(paste0(Sys.time(),'过滤低质量的细胞......'))
expressed_genes <- row.names(subset(fData(monocds), num_cells_expressed >= 10)) # nolint
monocds <- monocds[expressed_genes, ]
monocds.raw <- monocds
# 高变基因筛选方式1 使用monocle2进行筛选*********************************************************************************
message(paste0(Sys.time(),'高变基因筛选方式1......'))
disp_table <- dispersionTable(monocds)
express_genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
# dir=paste0('monocle_output')
# trajectory_cluster_dir=paste('./',dir,sep="")
monocds <- setOrderingFilter(monocds, express_genes)
p1 <- plot_ordering_genes(monocds) # 查看选择的基因的分布
ggsave(filename = paste0(outfile_Prefix,"01直接根据阈值进行筛选.png"),plot = p1,height = 4,width = 4.5, bg = "white")
monocds <- reduceDimension(
	monocds,
	max_components = 2,
	method = "DDRTree")
monocds <- orderCells(monocds)
qsave(monocds,file = paste0(outfile_Prefix,"01直接根据阈值进行筛选.qs"))
p2 <- plot_cell_trajectory(monocds,color_by = "State")
p3 <- plot_cell_trajectory(monocds,color_by = "celltype2")
p4 <- plot_cell_trajectory(monocds,color_by = "Treat_assess")
p5 <- plot_cell_trajectory(monocds,color_by="Pseudotime")
ggsave(filename = paste0(outfile_Prefix,"01直接根据阈值进行筛选02.png"),plot = cowplot::plot_grid(plotlist = list(p1,p2,p3,p4,p5),align = 'hv',nrow = 1,ncol = 5),height = 4,width = 21, bg = "white")
# 高变基因筛选方式2 使用高变基因*************************************************************************************
message(paste0(Sys.time(),'高变基因筛选方式2......'))
monocds <- monocds.raw
express_genes <- Seurat::VariableFeatures(scRNA)
monocds <- setOrderingFilter(monocds, express_genes)
p1 <- plot_ordering_genes(monocds) # 查看选择的基因的分布
ggsave(filename = paste0(outfile_Prefix,"02高变基因筛选进行筛选.png"),plot = p1,height = 4,width = 4.5, bg = "white")
monocds <- reduceDimension(
	monocds,
	max_components = 2,
	method = "DDRTree")
monocds <- orderCells(monocds)
qsave(monocds,file = paste0(outfile_Prefix,"02高变基因筛选进行筛选.qs"))

p2 <- plot_cell_trajectory(monocds,color_by = "State")
p3 <- plot_cell_trajectory(monocds,color_by = "celltype2")
p4 <- plot_cell_trajectory(monocds,color_by = "Treat_assess")
p5 <- plot_cell_trajectory(monocds,color_by="Pseudotime")
ggsave(filename = paste0(outfile_Prefix,"cd8_02高变基因筛选进行筛选02.png"),plot = cowplot::plot_grid(plotlist = list(p1,p2,p3,p4,p5),align = 'hv',nrow = 1,ncol = 5),height = 4,width = 21, bg = "white")

# 高变基因筛选方式3 cluster差异基因***************************************************************************************
message(paste0(Sys.time(),'高变基因筛选方式3......'))
monocds <- monocds.raw
Idents(scRNA) <- "celltype2"
deg.cluster <- FindAllMarkers(scRNA)
express_genes <- subset(deg.cluster,p_val_adj <0.05)$gene
monocds <- setOrderingFilter(monocds, express_genes)
p1 <- plot_ordering_genes(monocds) # 查看选择的基因的分布
ggsave(filename = paste0(outfile_Prefix,"03cluster差异基因进行筛选.png"),plot = p1,height = 4,width = 4.5, bg = "white")
monocds <- reduceDimension(
	monocds,
	max_components = 2,
	method = "DDRTree")
monocds <- orderCells(monocds)
qsave(monocds,file = paste0("03cluster差异基因进行筛选.qs"))
p2 <- plot_cell_trajectory(monocds,color_by = "State")
p3 <- plot_cell_trajectory(monocds,color_by = "celltype2")
p4 <- plot_cell_trajectory(monocds,color_by = "Treat_assess")
p5 <- plot_cell_trajectory(monocds,color_by="Pseudotime")
ggsave(filename = paste0(outfile_Prefix,"cd8_03cluster差异基因进行筛选02.png"),plot = cowplot::plot_grid(plotlist = list(p1,p2,p3,p4,p5),align = 'hv',nrow = 1,ncol = 5),height = 4,width = 21, bg = "white")

# 高变基因筛选方式4 使用dpFeature进行筛选***********************************************************************************
# Description:

#      Tests each gene for differential expression as a function of
#      pseudotime or according to other covariates as specified.
#      ‘differentialGeneTest’ is Monocle's main differential analysis
#      routine.  It accepts a CellDataSet and two model formulae as
#      input, which specify generalized lineage models as implemented by
#      the ‘VGAM’ package.

# Usage:

#      differentialGeneTest(
#        cds,
#        fullModelFormulaStr = "~sm.ns(Pseudotime, df=3)",
#        reducedModelFormulaStr = "~1",
#        relative_expr = TRUE,
#        cores = 1,
#        verbose = FALSE
#      )
message(paste0(Sys.time(),'高变基因筛选方式4......'))
monocds <- monocds.raw
diff <- differentialGeneTest(monocds,
	fullModelFormulaStr = "~celltype2",
	core=50
	)
head(diff)
degs <- subset(diff,qval<0.01)
degs <- degs[order(degs$qval,decreasing = T)]
express_genes <- rownames(degs)
monocds <- setOrderingFilter(monocds, express_genes)
p1 <- plot_ordering_genes(monocds) # 查看选择的基因的分布
ggsave(filename = paste0(outfile_Prefix,"04dpFeature差异基因选进行筛选.png"),plot = p1,height = 4,width = 4.5, bg = "white")
monocds <- reduceDimension(
	monocds,
	max_components = 2,
	method = "DDRTree")
monocds <- orderCells(monocds)
qsave(monocds,file = paste0(outfile_Prefix,"04dpFeature差异基因选进行筛选.qs"))
p2 <- plot_cell_trajectory(monocds,color_by = "State")
p3 <- plot_cell_trajectory(monocds,color_by = "celltype2")
p4 <- plot_cell_trajectory(monocds,color_by = "Treat_assess")
p5 <- plot_cell_trajectory(monocds,color_by="Pseudotime")
ggsave(filename = paste0(outfile_Prefix,"04dpFeature差异基因选进行筛选02.png"),plot = cowplot::plot_grid(plotlist = list(p1,p2,p3,p4,p5),align = 'hv',nrow = 1,ncol = 5),height = 4,width = 21, bg = "white")



