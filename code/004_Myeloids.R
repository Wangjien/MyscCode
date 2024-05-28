
suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(Matrix)
    library(ggplot2)
    library(SCP)
    library(tidydr)
    library(cowplot)
    library(RColorBrewer)
    library(patchwork)
    library(ggpubr)
    library(qs)
})
# 读入数据
setwd('/root/wangje/Project/OveryArtical')
scRNA <- qread('04_大群数据_scp.qs')

# 1 提取髓系亚群进行聚类分析**************************************************************************************************
scRNA <- CreateSeuratObject(scRNA@assays$RNA, meta.data = scRNA@meta.data)
scRNA <- subset(scRNA,subset=celltype == "Myeloid cells")
scRNA <- CreateSeuratObject(scRNA@assays$RNA@counts,meta.data = scRNA@meta.data[,c('sample','group')])
scRNA <- SCP::Integration_SCP(scRNA, batch = 'sample', integration_method = 'Harmony',cluster_resolution = seq(0.1,1.5,0.1))
qsave(scRNA, file = './06_髓系数据_scp.qs')
# BBKNN聚类
scRNA <- SCP::Integration_SCP(scRNA, batch = 'sample', integration_method = 'BBKNN',cluster_resolution = seq(0.1,1.5,0.1))
qsave(scRNA, file = './06_髓系数据_scp.qs')
# scVI聚类
scRNA <- SCP::Integration_SCP(scRNA, batch = 'sample', integration_method = 'scVI',cluster_resolution = seq(0.1,1.5,0.1))
qsave(scRNA, file = './06_髓系数据_scp.qs')

# 2 使用singleR粗注释*********************************************************************************************************
source("SingleRAnno.R")
scRNA <- Run_singleR(
    scRNA,
    cluster = scRNA$BBKNN_res.1,
    resolution = "BBKNN_res.1"
)
table(scRNA$singleR)
# Adipocytes        Basophils               DC  Dendritic cells Epithelial cells      Macrophages        Monocytes      Neutrophils
#     1200              418               49              170               23             4090             4978             2341
# 分布保存数据
qsave(scRNA[,scRNA$singleR == "Adipocytes"],file = "06_髓系细胞中提取脂肪细胞.qs")
qsave(scRNA[,scRNA$singleR == "Epithelial cells"],file = "06_髓系细胞中提取上皮细胞.qs")
qsave(scRNA[,!scRNA$singleR %in%c('Adipocytes','Epithelial cells')],file = "06_髓系数据_scp.qs")
# 使用singleR进行粗注释，可以发现在提取的髓系亚群中是有其他细胞存在的，可以使用FeaturePlot查看注释的准确性
# 3 去除其他细胞后再次进行聚类*************************************************************************************************
scRNA <- qread('06_髓系数据_scp.qs')
scRNA <- CreateSeuratObject(scRNA@assays$RNA@counts)
scRNA$sample <- stringr::str_split_fixed(rownames(scRNA@meta.data), "_[A|T|G|C].*", n = 2)[, 1]
scRNA <- SCP::Integration_SCP(scRNA, batch = "sample", integration_method = "Harmony", cluster_resolution = seq(0.1, 1.5, 0.1))
qsave(scRNA, file = "./06_髓系数据_scp.qs")
# BBKNN聚类
scRNA <- SCP::Integration_SCP(scRNA, batch = "sample", integration_method = "BBKNN", cluster_resolution = seq(0.1, 1.5, 0.1))
qsave(scRNA, file = "./06_髓系数据_scp.qs")
# scVI聚类
scRNA <- SCP::Integration_SCP(scRNA, batch = "sample", integration_method = "scVI", cluster_resolution = seq(0.1, 1.5, 0.1))
qsave(scRNA, file = "./06_髓系数据_scp.qs")
# # Seurat聚类
# scRNA <- SCP::Integration_SCP(scRNA, batch = "sample", integration_method = "Seurat", cluster_resolution = seq(0.1, 1.5, 0.1))
# qsave(scRNA, file = "./06_髓系数据_scp.qs")
# ComBat聚类
scRNA <- SCP::Integration_SCP(scRNA, batch = "sample", integration_method = "ComBat", cluster_resolution = seq(0.1, 1.5, 0.1))
qsave(scRNA, file = "./06_髓系数据_scp.qs")
# # Conos聚类
# scRNA <- SCP::Integration_SCP(scRNA, batch = "sample", integration_method = "Conos", cluster_resolution = seq(0.1, 1.5, 0.1))
# qsave(scRNA, file = "./06_髓系数据_scp.qs")
if (dir.exists("02MyeloidsCluster")) {
  setwd("02MyeloidsCluster")
} else {
  dir.create("./02MyeloidsCluster")
  setwd("./02MyeloidsCluster")
}

# 3 区分出大分群celltype******************************************************************************************************
myeloids = list(
    Mac=c("C1QA","C1QB","C1QC","SELENOP","RNASE1","DAB2","LGMN","PLTP","MAF","SLCO2B1"),
    Mono=c("VCAN","FCN1","CD300E","S100A12","EREG","APOBEC3A","STXBP2","ASGR1","CCR2","NRG1"),
    Neutrophils =  c("FCGR3B","CXCR2","SLC25A37","G0S2","CXCR1","ADGRG3","PROK2","STEAP4","CMTM2" ),
    Basophils = c('CLC','GATA2','FCER1A'),
    pDC = c("GZMB","SCT","CLIC3","LRRC26","LILRA4","PACSIN1","CLEC4C","MAP1A","PTCRA","C12orf75"),
    DC1 = c("CLEC9A","XCR1","CLNK","CADM1","ENPP1","SNX22","NCALD","DBN1","HLA-DOB","PPY"),
    DC2=c( "CD1C","CD1E","AL138899.1","CD2","GPAT3","CCND2","ENHO","PKIB","CD1B"),
    DC3 =  c("HMSD","ANKRD33B","LAD1","CCR7","LAMP3","CCL19","CCL22","INSM1","TNNT2","TUBB2B"),
    Mast = c('CPA3','TPSAB1','TPSB2')
)
library(RColorBrewer)

p1 <- DotPlot(object = scRNA, features = myeloids ,scale = T,group.by = "BBKNN_res.1.5") + ##celltype_l3. ###seurat_clusters
  scale_colour_gradientn(colors=brewer.pal(9, "YlGnBu"))+
  annotate(geom = "segment", y = Inf, yend = Inf, color = "black", x = -Inf, xend = Inf, linewidth = 1) +
  annotate(geom = "segment", x = Inf, xend = Inf, color = "black", y = -Inf, yend = Inf, linewidth = 1) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 14),
    strip.text.x = element_text(size = 14, angle = 0), axis.text.y = element_text(size = 15)
        ) +
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 0.7, linetype = "solid")) +
    labs(x = "Gene Marker", y = "Cluster") +
  theme(
      legend.position = "bottom",
      panel.spacing.x = unit(0.1, "cm"),
      panel.spacing.y = unit(0.1, "cm"),
      # strip.text = element_blank(),
      # strip.text.x.top = element_blank(),
      panel.grid = element_line(color = "grey", linetype = "dashed", linewidth = unit(0.1, "cm")))
p2 <- DimPlot(scRNA,reduction = 'BBKNNUMAP2D',group.by = 'BBKNN_res.1.5',label = T, repel = T, label.size = 6,
    cols = c(RColorBrewer::brewer.pal(n=12,name = "Paired"),RColorBrewer::brewer.pal(n=8,'Pastel2'),RColorBrewer::brewer.pal(n=6,'Set2'))) + 
  labs(title = 'resolution:1.5', x= 'umap1',y = 'umap2')
ggsave(filename = "./002_髓系细胞大群分群.png", height = 7,width = 24 , plot = p2+p1+patchwork::plot_layout(widths = c(1,4)), bg = "white")

library(magrittr)
scRNA@meta.data %<>% dplyr::mutate(
    celltype2 = dplyr::case_when(
        BBKNN_res.1.5 %in% c(15) ~ 'pDC',
        BBKNN_res.1.5 %in% c(16) ~ 'DC1', 
        BBKNN_res.1.5 %in% c(7) ~ 'DC2',
        BBKNN_res.1.5 %in% c(18) ~ 'DC3',
        BBKNN_res.1.5 %in% c(9) ~ 'Mast',
        BBKNN_res.1.5 %in% c(0,13) ~ 'Neutrophils',
        BBKNN_res.1.5 %in% c(5,6,11) ~ 'Monocyte',
        BBKNN_res.1.5 %in% c(17) ~ "Basophil",
        BBKNN_res.1.5 %in% c(1,2,3,4,8,10,12,14) ~ 'Macrophage',
    )
)

p1 <- DimPlot(scRNA,
  group.by = "celltype2", label = F, repel = T, raster = F,reduction = 'BBKNNUMAP2D',cols = RColorBrewer::brewer.pal(n = 12, name = "Paired")) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme(plot.title = element_blank()) + labs(title = "") +
  theme_blank(xlab = "UMAP1", ylab = "UMAP2")+ggtitle(paste0('nCells:',dim(scRNA)[2]))
