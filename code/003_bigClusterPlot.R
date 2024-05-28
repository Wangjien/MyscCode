
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



# 或使用tidydr中的theme_dr()

# 读入数据
setwd("/root/wangje/Project/OveryArtical")
scRNA <- qread("./04_大群数据_scp.qs")

# # 去除中六文章中的3个中间组样本
# scRNA_sub <- scRNA[,scRNA$sample %in% c("HRS421447","HRS421451","HRS421452","HRS421453",
#     "young1","young3","young4","middle2","middle3","middle4","na1129","old1","old2","old3","na-1010","na-426",
#     "na1111","na1122","na1117","na1128","NA_728","na-412")]

# 1. 绘制不同分组的DimPlot ********************************************************************************************************************
source("./plotFunc.R")
if (dir.exists("01AllCluster")) {
  setwd("01AllCluster")
} else {
  dir.create("./01AllCluster")
  setwd("./01AllCluster")
}
# 设置factor
scRNA$celltype <- factor(scRNA$celltype, levels = sort(names(sort(table(scRNA$celltype), decreasing = F))))
# Stromal & theca cells" "SMC" "Endothelial cells" "NK & T cells" "Myeloid cells" "Epithelial cells" "Oocyte|granule cells"
# "Adipocytes" "B cells"
p1 <- DimPlot(scRNA, group.by = "celltype", cols = c(RColorBrewer::brewer.pal(n = 7, name = "Dark2"), RColorBrewer::brewer.pal(n = 7, name = "Accent"))) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme_set
p1 <- p1 + labs(x = "UMAP1", y = "UMAP2")
ggsave(filename = "./001_AllClutserCelltypeDimPlot.png", height = 4, width = 6, plot = p1, bg = "white")
ggsave(filename = "./001_AllClutserCelltypeDimPlot.pdf", height = 4, width = 6, plot = p1, bg = "white", family = "ArialMT")
p1 <- DimPlot(scRNA, group.by = "celltype", cols = c(RColorBrewer::brewer.pal(n = 7, name = "Dark2"), RColorBrewer::brewer.pal(n = 7, name = "Accent"))) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme_test(base_rect_size = 1.1) + Seurat::NoAxes() + Seurat::NoGrid() + labs(title = paste0("nCells:", dim(scRNA)[2])) +
  theme(plot.title = element_text(hjust = 0.96, vjust = -8.5), legend.text = element_text(size = 14, color = "black"))
ggsave(filename = "./001_AllClutserCelltypeDimPlot02.png", height = 4, width = 6, plot = p1, bg = "white")
ggsave(filename = "./001_AllClutserCelltypeDimPlot02.pdf", height = 4, width = 6, plot = p1, bg = "white", family = "ArialMT")
p1 <- DimPlot(scRNA, group.by = "celltype",label = T,repel = T,cols = c(RColorBrewer::brewer.pal(n = 7, name = "Dark2"), RColorBrewer::brewer.pal(n = 7, name = "Accent"))) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme_test(base_rect_size = 1.1) + Seurat::NoAxes() + Seurat::NoGrid() + labs(title = paste0("nCells:", dim(scRNA)[2])) +
  theme(plot.title = element_text(hjust = 0.96, vjust = -8.5), legend.text = element_text(size = 14, color = "black"))
ggsave(filename = "./001_AllClutserCelltypeDimPlot03.png", height = 4, width = 6, plot = p1, bg = "white")
ggsave(filename = "./001_AllClutserCelltypeDimPlot03.pdf", height = 4, width = 6, plot = p1, bg = "white", family = "ArialMT")

p1 <- DimPlot(scRNA, group.by = "celltype",label = T,repel = T,cols = c(RColorBrewer::brewer.pal(n = 7, name = "Dark2"), RColorBrewer::brewer.pal(n = 7, name = "Accent"))) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme_test(base_rect_size = 1.1) + Seurat::NoAxes() + Seurat::NoGrid() + labs(title = paste0("nCells:", dim(scRNA)[2])) +
  theme(plot.title = element_text(hjust = 0.96, vjust = -8.5), 
  legend.text = element_text(size = 14, color = "black"),legend.position = "left")
ggsave(filename = "./001_AllClutserCelltypeDimPlot04.png", height = 4, width = 6, plot = p1, bg = "white")
ggsave(filename = "./001_AllClutserCelltypeDimPlot04.pdf", height = 4, width = 6, plot = p1, bg = "white", family = "ArialMT")

p2 <- DimPlot(scRNA, group.by = "sample",label = F,repel = T,cols = c(RColorBrewer::brewer.pal(n = 12, name = "Paired"), 
  RColorBrewer::brewer.pal(n = 8, name = "Set2"),RColorBrewer::brewer.pal(n = 6, name = "Dark2"))) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme_test(base_rect_size = 1.1) + Seurat::NoAxes() + Seurat::NoGrid() + labs(title ="Patient") +
  theme(plot.title = element_text(hjust = 0.5,size = 18,color="black"), legend.text = element_text(size = 10, color = "black"),
  legend.position = "right")
ggsave(filename = "./001_AllClutserCelltypeDimPlot05.png", height = 4,width = 5, plot = p2, bg = "white")
ggsave(filename = "./001_AllClutserCelltypeDimPlot05.pdf", height = 4, width = 5, plot = p2, bg = "white", family = "ArialMT")

p3 <- DimPlot(scRNA, group.by = "group1",label = F,repel = T,cols = c(RColorBrewer::brewer.pal(n = 12, name = "Paired"), 
  RColorBrewer::brewer.pal(n = 8, name = "Set2"),RColorBrewer::brewer.pal(n = 6, name = "Dark2"))) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme_test(base_rect_size = 1.1) + Seurat::NoAxes() + Seurat::NoGrid() + labs(title ="Group") +
  theme(plot.title = element_text(hjust = 0.5,size = 18,color="black"), legend.text = element_text(size = 10, color = "black"),
  legend.position = "right")
ggsave(filename = "./001_AllClutserCelltypeDimPlot06.png", height = 4, width = 6, plot = p3, bg = "white")
ggsave(filename = "./001_AllClutserCelltypeDimPlot06.pdf", height = 4, width = 6, plot = p3, bg = "white", family = "ArialMT")

# split dimplot
p4 <- DimPlot(scRNA,group.by = "celltype", label = F, repel = T, raster = F,split.by = 'group1',ncol = 6,
  cols= c(RColorBrewer::brewer.pal(n = 7, name = "Dark2"), RColorBrewer::brewer.pal(n = 7, name = "Accent"))
) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme(plot.title = element_blank()) + labs(title = "") +
  theme_blank(xlab = "UMAP1", ylab = "UMAP2")
p4 <- p4 + theme(strip.text = element_text(size = 20, color = "black", face = "bold"),strip.background = element_rect(fill = '#fffbe9'))

ggsave(filename = "./001_AllClutserCelltypeSplitDimPlot06.png", height = 4, width = 26, plot = p4, bg = "white",limitsize = FALSE)
ggsave(filename = "./001_AllClutserCelltypeSplitDimPlot06.pdf", height = 4, width = 26, plot = p4, bg = "white", family = "ArialMT",limitsize = FALSE)

# 2 绘制marker基因表达 ************************************************************************************************************************
mark_list <- list(
  NK = c("KLRD1"),
  "T cells" = c("CD3E"),
  "Theca cells" = c("DCN"),
  Myeloids = c("CD68"),
  Epithelial = c("EPCAM"),
  Endothelial = c("CLDN5"),
  Bcells = c("CD79A"),
  "Smooth muscle" = c("PLAC9"),
  Neutrophil = c("FCGR3B"),
  Mast = c("CPA3")
)
message(paste0(Sys.time()," 绘制marker基因FeaturePlot..."))
plist <- list()
for (i in 1:length(mark_list)) {
for (j in mark_list[[i]]) {
  plist[[paste0(i, "|", j)]] <-
    FeaturePlot(scRNA,
      features = j, reduction = "HarmonyUMAP2D", max.cutoff = 1.5,raster = TRUE,
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
ggsave(filename = "./001_AllClutserFeaturePlot.png", height = 6, width = 15,plot = patchwork::wrap_plots(plist, ncol = 5, nrow = 2))
ggsave(filename = "./001_AllClutserFeaturePlot.pdf", height = 6, width = 15,plot = patchwork::wrap_plots(plist, ncol = 5, nrow = 2))

# 3 绘制marker基因表达密度表达散点图 **********************************************************************************************************
message(paste0(Sys.time()," 绘制marker密度热图..."))
reduction <- 'HarmonyUMAP2D'
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
ggsave(filename = "./001_AllClutserDensityPlot.png", height = 6, width = 18,plot = patchwork::wrap_plots(plist, ncol = 5, nrow = 2))
ggsave(filename = "./001_AllClutserDensityPlot.pdf", height = 6, width = 18,plot = patchwork::wrap_plots(plist, ncol = 5, nrow = 2))

# 4 绘制变化比例箱线图 ***********************************************************************************************************************
# 绘制细胞比例柱形图
library(rstatix)
library(ggplot2)
library(ggpubr)

df1 <- scRNA@meta.data %>%
  dplyr::select(group1, celltype, sample) %>%
  dplyr::group_by(group1, celltype,sample) %>%
  dplyr::count() %>%
  dplyr::rename("n1" = "n")
df2 <- scRNA@meta.data %>%
  dplyr::group_by(group1) %>%
  dplyr::count()
df <- df1 %>%
  left_join(df2, by = "group1") %>%
  data.table::setDT()
# 计算比例
df[, Prop := n1 / n]
write.table(df,file = "001_AllClutser大群细胞比例.txt",sep = '\t',row.names = F,quote = F)
# 比较组合
compair <- list(c('18-29','37-39'),c('18-29','39-44'),c('18-29','32-35'),c('18-29','47-49'),c('18-29','55-60'),c('37-39','39-44'),
  c('37-39','32-35'),c('37-39','47-49'),c('37-39','55-60'),c('39-44','32-35'),c('39-44','47-49'),c('39-44','55-60'),c('32-35','47-49'),
  c('32-35','55-60'),c('47-49','55-60'))  
# df$group <- factor(df$group, levels = c("20_35", "38_44", "50_60"))
df$celltype <- varhandle::unfactor(df$celltype)

lapply(unique(df$celltype), FUN = function(x) {
  tryCatch(
    {
      df_sub <- df[df$celltype == x, ]
      # plot_boxplot
      p1 <- ggplot(df_sub) +
        geom_boxplot(
          aes(x = celltype, y = Prop, color = group1, fill = group1),
          position = position_dodge(width = 0.8),
          outlier.shape = NA, color = "black"
        ) +
        geom_jitter(aes(x = celltype, y = Prop, fill = group1),
          color = "black", position = position_dodge(width = 0.8), pch = 21
        )
      # 计算显著性
      stat.test <- df_sub %>%
        group_by(celltype) %>%
        wilcox_test(Prop ~ group1, comparisons = compair) %>%
        add_xy_position(x = "celltype", dodge = 0.8)
      # 添加新的显著性
      stat.test$new_signif <- case_when(
        0 <= stat.test$p & stat.test$p < 0.01 ~ "***",
        0.01 <= stat.test$p & stat.test$p < 0.05 ~ "**",
        0.05 <= stat.test$p & stat.test$p < 0.1 ~ "*",
        TRUE ~ "ns"
      )
      stat.test <- stat.test %>% filter(new_signif != "ns")
      p2 <- p1 + stat_pvalue_manual(
        stat.test,
        label = "new_signif", tip.length = 0.00, size = 6,
        hide.ns = FALSE
      ) + labs(y = "Fraction", x = "", title = x)
      p2 <- p2 +
        theme_classic(base_size = 20, base_line_size = 1) +
        scale_fill_manual(values = c(RColorBrewer::brewer.pal(n = 6, name = "Set3"))) +
        theme(
          axis.title.x = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 15, color = "black"),
          axis.text.y = element_text(size = 15, colour = "black"),
          axis.text.x = element_blank(),
          axis.title = element_text(size = 15, colour = "black"),
          axis.ticks.y = element_line(linewidth = 1.2),
          axis.ticks.x = element_blank(),
          axis.ticks.length = unit(0.4, "cm"),
          plot.title = element_text(hjust = 0.5, colour = "black", size = 18, angle = 0)
        ) +
        scale_x_discrete(expand = c(0.45, 0))
      return(p2)
    },
    error = function(e) {
      message(e$message)
      return(NULL)
    }
  )
}) -> plist

p5 <- ggpubr::ggarrange(plotlist = plist, ncol = 5, nrow = 2, common.legend = TRUE, align = "v", legend = "bottom")
ggsave(filename = "./001_AllClutser比例箱线图.png", height = 8, width = 15,plot = p5,bg = "white")
ggsave(filename = "./001_AllClutser比例箱线图.pdf", height = 8, width = 15,plot = p5,bg="white")

# 5 绘制细胞比例柱形图*************************************************************************************************************************
prop <- as.data.frame(prop.table(table(scRNA$celltype, scRNA$group1), margin = 2))
p1 <- ggplot(prop, aes(x = Var2, y = Freq, fill = Var1, stratum = Var1, alluvium = Var1)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", curve_type = "linear", alpha = 0.5, width = 0.7, color = "white") +
  geom_stratum(alpha = 1, color = "black", width = 0.7) +
  # geom_text(aes(label = Prop),position = position_stack(vjust = .5),size =3.5, color = 'black') +
  theme_classic(base_size = 20, base_line_size = 1) +
  guides(fill = guide_legend(override.aes = list(size = 6))) +
  # scale_x_discrete(expand = c(0.04,0))+
  labs(x = "", y = "% Fraction") +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 1, vjust = 0.5, size = 22, color = "black"),
    # panel.border = element_rect(fill = NA, color = "black", size = 1.3, linetype = "solid"),
    legend.text = element_text(size = 18, color = "black"),
    axis.title.y = element_text(size = 22, colour = "black"),
    axis.text.y = element_text(size = 22, colour = "black"),
    legend.position = 'top',
    # axis.text.x = element_text(size = 20, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
    # axis.text.x = element_text(size = 20, colour = "black", angle = 45, hjust = 1, vjust = 1),
    axis.ticks.y = element_line(linewidth = 1.3),
    axis.ticks.length.y = unit(0.4, "cm"),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
  ) +
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(n = 7, name = "Dark2"), RColorBrewer::brewer.pal(n = 7, name = "Accent")))