
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

# 4 绘制DimPlot*******************************************************************************************************************************
p1 <- DimPlot(scRNA,
  group.by = "celltype2", label = T, repel = T, raster = F,reduction = 'BBKNNUMAP2D',cols = RColorBrewer::brewer.pal(n = 12, name = "Paired")) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme(plot.title = element_blank()) + labs(title = "") +
  theme_blank(xlab = "UMAP1", ylab = "UMAP2")+ggtitle(paste0('nCells:',dim(scRNA)[2]))
ggsave(filename = "./002_髓系大群DimPlot.png", height = 4, width = 6, plot = p1, bg = "white")
ggsave(filename = "./002_髓系大群DimPlot.pdf", height = 4, width = 6, plot = p1, bg = "white", family = "ArialMT")

p2 <- DimPlot(scRNA,
  group.by = "sample", label = F, repel = T, raster = F,reduction = 'BBKNNUMAP2D',
  cols = c(RColorBrewer::brewer.pal(n = 12, name = "Paired"), 
  RColorBrewer::brewer.pal(n = 8, name = "Set2"),RColorBrewer::brewer.pal(n = 6, name = "Dark2"))) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme(plot.title = element_blank()) + labs(title = "") +
  theme_blank(xlab = "UMAP1", ylab = "UMAP2")+ggtitle('Patient')
ggsave(filename = "./002_髓系大群DimPlotSample.png", height = 4, width = 6, plot = p2, bg = "white")
ggsave(filename = "./002_髓系大群DimPlotSample.pdf", height = 4, width = 6, plot = p2, bg = "white", family = "ArialMT")

p3 <- DimPlot(scRNA, group.by = "group1",reduction = 'BBKNNUMAP2D',label = F,repel = T,cols = c(RColorBrewer::brewer.pal(n = 12, name = "Paired"), 
  RColorBrewer::brewer.pal(n = 8, name = "Set2"),RColorBrewer::brewer.pal(n = 6, name = "Dark2"))) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme_blank(xlab = "UMAP1", ylab = "UMAP2")+ labs(title ="Group") +
  theme(plot.title = element_text(hjust = 0.5,size = 18,color="black"), legend.text = element_text(size = 10, color = "black"),
  legend.position = "right")
ggsave(filename = "./002_髓系大群DimPlotGroup.png", height = 4, width = 5, plot = p3, bg = "white")
ggsave(filename = "./002_髓系大群DimPlotGroup.pdf", height = 4, width = 5, plot = p3, bg = "white", family = "ArialMT")

p4 <- DimPlot(scRNA,group.by = "celltype2", label = F, repel = T, raster = F,split.by = 'group1',ncol = 6,reduction = 'BBKNNUMAP2D',
  cols = RColorBrewer::brewer.pal(n = 12, name = "Paired")) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme(plot.title = element_blank()) + labs(title = "") +
  theme_blank(xlab = "UMAP1", ylab = "UMAP2")
p4 <- p4 + theme(strip.text = element_text(size = 20, color = "black", face = "bold"),strip.background = element_rect(fill = '#fffbe9'))
ggsave(filename = "./002_髓系大群SplitDimPlot06.png", height = 4, width = 26, plot = p4, bg = "white",limitsize = FALSE)
ggsave(filename = "./002_髓系大群SplitDimPlot06.pdf", height = 4, width = 26, plot = p4, bg = "white", family = "ArialMT",limitsize = FALSE)

mark_list <- list(
    Macrophage = c("C1QA", "C1QB"),
    Monocyte = c("VCAN", "FCN1"),
    Neutrophils = c("FCGR3B", "CXCR2"),
    Basophils = c("CLC", "GATA2"),
    pDC = c("GZMB", "SCT"),
    DC1 = c("CLEC9A", "XCR1"),
    DC2 = c("CD1C", "CD1E"),
    DC3 = c("HMSD", "ANKRD33B"),
    Mast = c("CPA3", "TPSAB1")
)
plist <- list()
for (i in 1:length(mark_list)) {
    for (j in mark_list[[i]]) {
        plist[[paste0(i, "|", j)]] <-
            FeaturePlot(scRNA,
                features = j, reduction = "BBKNNUMAP2D", pt.size = 0.0001, max.cutoff = 1.5,
                cols = c("#FFEFD5", "#E6E6FA", "#87CEFA", "#6495ED", "#4169E1", "#0000CD", "#000080")
            ) +
            scale_x_continuous("") + scale_y_continuous("") +
            theme_bw() +
            theme(
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.ticks = element_blank(), axis.text = element_blank(),
                legend.position = "none",
                plot.title = element_text(hjust = 0.5, size = 18)
            ) + ggtitle(paste0(names(mark_list[i]), " (", j, ")"))
    }
}
ggsave(filename = "./002_髓系大群FeaturePlot.png", height = 9, width = 18,plot = patchwork::wrap_plots(plist, ncol = 6, nrow = 3))
ggsave(filename = "./002_髓系大群FeaturePlot.pdf", height = 9, width = 18,plot = patchwork::wrap_plots(plist, ncol = 6, nrow = 3))

message(paste0(Sys.time()," 绘制marker密度热图..."))
reduction <- 'BBKNNUMAP2D'
plist <- list()
for(i in names(mark_list)){
  for(j in mark_list[[i]]){
    result <- tryCatch({
       scCustomize::Plot_Density_Custom(scRNA, features = j,reduction = reduction, viridis_palette = "viridis")+
        labs(title = paste0(i,":",j))+
        Seurat::NoGrid()+
        Seurat::NoAxes()+
        theme_test(base_size = 15,
       base_line_size = 1.3,
       base_rect_size = 1.3)+
    theme(plot.title = element_text(hjust =0.5),axis.text =element_blank(),axis.ticks = element_blank(), axis.title = element_blank(),legend.title = element_text(size = 15), legend.text =element_text(size =12))
    },error = function(e){
      message(paste0(i,"|",j," error!"))
    },finally={
      message(paste0(Sys.time()," ",i,"|",j, " finish..."))
    })
    plist[[paste0(i,'|',j)]] <- result
  }
}
ggsave(filename = "./002_髓系大群基因表达密度图.png", height = 9, width = 21,plot = patchwork::wrap_plots(plist, ncol = 6, nrow = 3))
ggsave(filename = "./002_髓系大群基因表达密度图.pdf", height = 9, width = 21,plot = patchwork::wrap_plots(plist, ncol = 6, nrow = 3))
# 5 绘制比例变化箱线图***************************************************************************************************************************







# 7 提取出其中的Macrophage和Monocytes再次聚类分群*************************************************************************************************
scRNA <- scRNA[,scRNA$celltype2 %in% c('Macrophage','Monocyte')]
scRNA <- CreateSeuratObject(scRNA@assays$RNA@counts)
scRNA$sample <- stringr::str_split_fixed(rownames(scRNA@meta.data), "_[A|T|G|C].*", n = 2)[, 1]
scRNA <- SCP::Integration_SCP(scRNA, batch = "sample", integration_method = "Harmony", cluster_resolution = seq(0.1, 1.5, 0.1))
qsave(scRNA, file = "./06_髓系数据MacroMono_scp.qs")
# BBKNN聚类
scRNA <- SCP::Integration_SCP(scRNA, batch = "sample", integration_method = "BBKNN", cluster_resolution = seq(0.1, 1.5, 0.1))
qsave(scRNA, file = "./06_髓系数据MacroMono_scp.qs")
# scVI聚类
scRNA <- SCP::Integration_SCP(scRNA, batch = "sample", integration_method = "scVI", cluster_resolution = seq(0.1, 1.5, 0.1))
qsave(scRNA, file = "./06_髓系数据MacroMono_scp.qs")
# # Seurat聚类
# scRNA <- SCP::Integration_SCP(scRNA, batch = "sample", integration_method = "Seurat", cluster_resolution = seq(0.1, 1.5, 0.1))
# qsave(scRNA, file = "./06_髓系数据_scp.qs")
# ComBat聚类
scRNA <- SCP::Integration_SCP(scRNA, batch = "sample", integration_method = "ComBat", cluster_resolution = seq(0.1, 1.5, 0.1))
qsave(scRNA, file = "./06_髓系数据MacroMNO_scp.qs")

# 8 Macrophage和Monocytes进行分群**************************************************************************************************************
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
ggsave(filename = "./002_MacrophageMonocytes大群分群.png", height = 7,width = 24 , plot = p2+p1+patchwork::plot_layout(widths = c(1,4)), bg = "white")

## Top10差异基因绘制热图 
library(viridis)
DefaultAssay(scRNA) = "RNA"
Seurat::Idents(scRNA) = "BBKNN_res.1.5"
scRNA <- Seurat::ScaleData(scRNA)
scRNA.Findmarkers <- FindAllMarkers(
  scRNA, 
  only.pos = TRUE, 
  # min.pct = 0.25, 
  # logfc.threshold = 0.25, 
  BPPARAM = MulticoreParam(20))
Top10.genes <- scRNA.Findmarkers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
p3 <- DoHeatmap(scRNA,features=Top10.genes$gene)+scale_fill_viridis()
write.csv(Top10.genes, file ="MacroMono分辨率1.5Top10基因.csv", quote = F, row.names = T)