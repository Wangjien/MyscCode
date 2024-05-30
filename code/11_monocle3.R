

# library(monocle3)
library(Seurat)
library(ggplot2)
library(dplyr)
library(qs)
library(Matrix)
# 1 构建CDS对象************************************************************************************************************
message("Prepare CDS object ......")
setwd('/root/wangje/Project/刘老师/NK_T/new_Result/Data/001_cytotrace2_result')
scRNA <- qread('CD4_celltype2聚类结果.qs')
Idents(scRNA) <- scRNA$celltype2
DefaultAssay(scRNA) ='RNA'
expression_matrix <- as(as.matrix(scRNA@assays$RNA@counts), 'sparseMatrix')
cell_metadata <- scRNA@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix))
rownames(gene_annotation) <- rownames(expression_matrix)

#构建cds对象
# library(monocle3,lib.loc = '/root/wangje/miniconda3/lib/R/library')
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 100)#method默认为PCA
p1 <- plot_pc_variance_explained(cds)#展示PC数，和seurat降维一摸一样