

# library(monocle3)
library(Seurat)
library(ggplot2)
library(dplyr)
library(qs)
library(Matrix)
library(monocle3)
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
cds <- reduce_dimension(cds,reduction_method='UMAP',
                        preprocess_method = 'PCA')
#将monocle对象与我们的seurat结合。转化为我们的聚类信息。
cds.embed <- cds@int_colData$reducedDims$UMAP#monocle3中UMAP信息
int.embed <- Embeddings(scRNA, reduction = "scVIUMAP2D")#seurat中UMAP降维信息
int.embed <- int.embed[rownames(cds.embed),]#替换
cds@int_colData$reducedDims$UMAP <- int.embed #替换

#画图看一下，分群就是我们seurat中的了
# color_cells_by <- "celltype2"
# p2 <- plot_cells(cds, color_cells_by=color_cells_by,
#            cell_size=0.5,group_label_size=4) 

mycds <- cds
mycds <- learn_graph(mycds,
                     verbose=T,
                     learn_graph_control=list(minimal_branch_len=20,#在图修剪过程中要保留的分支直径路径的最小长度。默认值是10。
                                              euclidean_distance_ratio=10#生成树中两个末端节点的欧氏距离与生成树上允许连接的任何连接点之间的最大距离之比。默认值为1。
                     ))

#推断轨迹树
p3 <- plot_cells(mycds, 
           color_cells_by = color_cells_by,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size=0.5,group_label_size=4)


p4 <- plot_cells(mycds, 
           color_cells_by = color_cells_by, 
           label_groups_by_cluster=FALSE,
           label_leaves=TRUE, 
           label_branch_points=TRUE,
           graph_label_size=4)
qsave(mycds, file = "./cd4_monocle3.qs")
ggsave(filename = "./cd4_monocle3.png",plot = p1|p2|p3|p4,width = 21,height = 4,bg = "white",limitsize = F)