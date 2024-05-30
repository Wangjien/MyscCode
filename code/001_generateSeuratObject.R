suppressMessages(suppressWarnings({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(magrittr)
  library(ggpubr)
  library(stringr)
  library(qs)
}))

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#  合并两个文章的数据与自测数据，并保存raw Seurat对象
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#**********************************************************************************************
# Single-Cell Atlas of Human Ovaries Reveals The Role Of The Pyroptotic Macrophage in Ovarian Aging
## 从每个样本目录中读入表达矩阵文件并生成Seurt对象**********************************************
message(paste0(Sys.time(), "Single-Cell Atlas of Human Ovaries Reveals The Role Of The Pyroptotic Macrophage in Ovarian Aging"))
path <- "/root/wangje/Project/OveryArtical/Aritical2"
flist <- setNames(lapply(seq(421447, 421453), FUN = function(x) {
  message(paste0("HRS", x))
  RenameCells(CreateSeuratObject(Read10X(file.path(path, paste0("HRS", x), paste0("HRS", x), "outs", "filtered_feature_bc_matrix")),
    min.cell = 3, min.feature = 200
  ), add.cell.id = paste0("HRS", x))
}), paste0("HRS", seq(421447, 421453)))
## 合并Seurat对象列表生成一个Seurat对象*********************************************************
for (name in names(flist)) {
  flist[[name]][["sample"]] <- name
}
# scRNA.raw <- merge(flist[[1]],flist[-1])
# message(sprintf(c("基因为:%s","细胞为:%s"),dim(scRNA.raw)))

#**********************************************************************************************
# Spatiotemporal transcriptomic changes of human ovarian aging and the regulatory role of FOXP1
message(paste0(Sys.time(), "Spatiotemporal transcriptomic changes of human ovarian aging and the regulatory role of FOXP1"))
path <- "/root/wangje/Project/OveryArtical/GSE255690/scRNA/ExpressionData"
young <- setNames(lapply(file.path(path, "young", c("young1", "young3", "young4")),
  FUN = function(x) CreateSeuratObject(Read10X(x), min.cell = 3, min.feature = 200)
), c("young1", "young3", "young4"))

middle <- setNames(lapply(file.path(path, "middle", c("middle2", "middle3", "middle4")),
  FUN = function(x) CreateSeuratObject(Read10X(x), min.cell = 3, min.feature = 200)
), c("middle2", "middle3", "middle4"))

old <- setNames(lapply(file.path(path, "old", c("old1", "old2", "old3")),
  FUN = function(x) CreateSeuratObject(Read10X(x), min.cell = 3, min.feature = 200)
), c("old1", "old2", "old3"))
# 合并列表
flist1 <- c(young, middle, old)
# 添加sample列
for (name in names(flist1)) {
  flist1[[name]][["sample"]] <- name
  flist1[[name]] <- RenameCells(flist1[[name]], add.cell.id = name)
}

#**********************************************************************************************
# 自测数据
message(paste0(Sys.time(), "自测数据"))
inputDir <- "/root/wangje/Project/临床_薛老师_卵巢/Data/rawData"
setwd(inputDir)
dirs <- list.files("./")
flist2 <- list()
for (i in dirs) {
  for (j in list.files(i)) {
    tryCatch(
      {
        message(paste0(Sys.time(), " ", i, " ", j))
        print(file.path(inputDir, i, j, "filtered_feature_bc_matrix/"))
        tmp <- Seurat::Read10X(file.path(inputDir, i, j, "filtered_feature_bc_matrix/"))
        tmp <- Seurat::CreateSeuratObject(tmp, min.cell = 3, min.feature = 200)
        tmp <- Seurat::RenameCells(tmp, add.cell.id = j)
        flist2[[paste0(i, "|", j)]] <- tmp
      },
      error = function(e) {
        message(e$message)
        return(NULL)
      }
    )
  }
}
for (i in names(flist2)) {
  flist2[[i]]$sample <- stringr::str_split_fixed(names(flist2), "\\|", n = 2)[, 2]
}
#***********************************************************************************************
# 合并数据
AllList <- c(flist, flist1, flist2)
scRNA <- merge(AllList[[1]], AllList[-1])
message(sprintf(c("基因为:%s", "细胞为:%s"), dim(scRNA)))
scRNA <- PercentageFeatureSet(scRNA, pattern = '^MT', col.name = 'percent.mt')
#***********************************************************************************************
setwd("/root/wangje/Project/OveryArtical")
qsave(scRNA, file = "./01_mergeRawSeurat.qs")
# 绘制小提琴图
scRNA$group <- ifelse(scRNA$sample %in% c("HRS421447", "HRS421448", "HRS421449", "HRS421450", "HRS421451", "HRS421452", "HRS421453"), "Artical1", ifelse(
  scRNA$sample %in% c("young1", "young3", "young4", "middle2", "middle3", "middle4", "old1", "old2", "old3"), "Artical2", "Self"
))
message(paste0(Sys.time(), " 绘制质控小提琴图... ..."))
p1 <- Seurat::VlnPlot(scRNA, group.by = "sample", features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3, pt.size = 0.001, raster = TRUE) &
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.title.x = element_blank())

p2 <- Seurat::VlnPlot(scRNA, group.by = "sample", features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3, pt.size = 0) &
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.title.x = element_blank())
# 保存图片
ggsave(
  filename = file.path('./',"01_rawSeuratOnjectVlnplot.png"), height = 8, width = 16,
  plot = patchwork::wrap_plots(list(p1, p2), nrow = 2), limitsize = FALSE, bg = "white"
)

# 写出表格信息
message(paste0(Sys.time(), " 绘制统计结果... ..."))
df <- as.data.frame(scRNA@meta.data %>% dplyr::select(sample, group) %>% distinct())
rownames(df) <- 1:nrow(df)
colnames(df) <- c("Patient", "Group")
df$nCells <- sapply(unique(df$Patient), FUN = function(x) {
  ncol(scRNA[, scRNA$sample == x])
})
df$nGenes <- sapply(unique(df$Patient), FUN = function(x) {
  nrow(scRNA[, scRNA$sample == x])
})
df$`nCount_RNA(min->max)` <- sapply(unique(df$Patient), FUN = function(x) {
  paste0(min(scRNA[, scRNA$sample == x]$nCount_RNA), "->", max(scRNA[, scRNA$sample == x]$nCount_RNA))
})
df$nCount_RNA_median <- sapply(unique(df$Patient), FUN = function(x) {
  median(scRNA[, scRNA$sample == x]$nCount_RNA)
})
df$`nFeature_RNA(min->max)` <- sapply(unique(df$Patient), FUN = function(x) {
  paste0(min(scRNA[, scRNA$sample == x]$nFeature_RNA), "->", max(scRNA[, scRNA$sample == x]$nFeature_RNA))
})
df$nFeature_RNA_median <- sapply(unique(df$Patient), FUN = function(x) {
  median(scRNA[, scRNA$sample == x]$nFeature_RNA)
})
df$`percent.mt(min->max)` <- sapply(unique(df$Patient), FUN = function(x) {
  paste0(round(min(scRNA[, scRNA$sample == x]$percent.mt), digits = 4), "->", round(max(scRNA[, scRNA$sample == x]$percent.mt), digits = 4))
})
df$percent.mt_median <- sapply(unique(df$Patient), FUN = function(x) {
  median(scRNA[, scRNA$sample == x]$percent.mt)
})

# 绘制表格
tbody.style <- tbody_style(
  color = "black",
  fill = c("#e8f3de", "#d3e8bb"), hjust = 1, x = 0.9
)
ggtexttable(df,
  rows = NULL,
  theme = ttheme(
    colnames.style = colnames_style(color = "white", fill = "#8cc257"),
    tbody.style = tbody.style
  )
) -> t_table
# 保存统计表格
ggsave(
  filename = file.path('./', "01_rawSeuratOnjectCount.png"), height = 8, width = 16,
  plot = t_table, limitsize = FALSE, bg = "white"
)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#  数据粗过滤并使用DoubletFinder检测双细胞
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
filter_info <- scRNA$nCount_RNA <30000 & scRNA$nFeature_RNA > 200 & scRNA$nFeature_RNA < 7000 & scRNA$percent.mt < 20
scRNA.filter <- scRNA[,filter_info]
scRNA.filter$group <- ifelse(scRNA.filter$sample %in% c("HRS421447", "HRS421448", "HRS421449", "HRS421450", "HRS421451", "HRS421452", "HRS421453"), "Artical1", ifelse(
  scRNA$sample %in% c("young1", "young3", "young4", "middle2", "middle3", "middle4", "old1", "old2", "old3"), "Artical2", "Self"
))
message(paste0(Sys.time(), " 绘制质控小提琴图... ..."))
p1 <- Seurat::VlnPlot(scRNA.filter, group.by = "sample", features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3, pt.size = 0.001, raster = TRUE) &
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.title.x = element_blank())

p2 <- Seurat::VlnPlot(scRNA.filter, group.by = "sample", features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3, pt.size = 0) &
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.title.x = element_blank())
# 保存图片
ggsave(
  filename = file.path('./',"02_rawFilterSeuratOnjectVlnplot.png"), height = 8, width = 16,
  plot = patchwork::wrap_plots(list(p1, p2), nrow = 2), limitsize = FALSE, bg = "white"
)
# scRNA.filter 33047 features across 227793 samples within 1 assay
qsave(scRNA.filter, file = "./02_rawFilterSeurat.qs")
#进行双细胞检测*************************************************************************************
rm(scRNA);gc()
Merge_obj <- NormalizeData(scRNA.filter, normalization.method = "LogNormalize") %>% 
  FindVariableFeatures(., selection.method = "vst", nfeatures = 3000)
Merge_obj <- ScaleData(Merge_obj, vars.to.regress = c("percent.mt","nCount_RNA"), verbose = T)
Merge_obj <- RunPCA(Merge_obj, npcs = 50, verbose = FALSE)

#去批次
library(harmony)
Merge_obj <- RunHarmony(Merge_obj,"sample")
Merge_obj <- RunUMAP(Merge_obj,  dims = 1:20, reduction = "harmony")
Merge_obj <- FindNeighbors(Merge_obj, reduction = "harmony",dims = 1:20) 
Merge_obj  <- FindClusters(object = Merge_obj , resolution = seq(from = 0.1, to = 1.0,  by = 0.1))

#单个分开，用来做DoubletFinder
sce_list <- SplitObject(Merge_obj, split.by = "sample")

detectDoublet <- function(
  obj, #seurat obj
  dims, #Pc数目
  estDubRate, #期望双细胞率，该值最好通过细胞在10X/Drop-Seq装置上的负载密度来估计
  ncores=1,#线程数
  SCTransform,#True or False
  Homotypic=F,#是否需要优化同源双细胞
  annotation){ #听说最好是celltype
  #use DoubletFinder packages
  require(DoubletFinder)#2.0.4
  
  #select pK
  sweep.res.list <- paramSweep(obj, PCs=dims, sct=SCTransform, num.cores=ncores)
  sweep.stats <- summarizeSweep(sweep.res.list, GT=FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
  message(sprintf("Using pK = %s...", pK))
  
  
  #Doublet Proportion Estimate
  if(Homotypic==F){
    
    nExp_poi <- round(estDubRate * length(Cells(obj)))
    
  }else{
    
    homotypic.prop <- modelHomotypic(obj@meta.data[,annotation])
    nExp_poi <- round(estDubRate * length(Cells(obj)))
    nExp_poi  <- round(nExp_poi*(1-homotypic.prop))
  }
  
  # DoubletFinder:
  obj <- doubletFinder(obj, PCs = dims, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = SCTransform)
  
  # Rename results into more useful annotations
  pann <- grep(pattern="^pANN", x=names(obj@meta.data), value=TRUE)
  message(sprintf("Using pANN = %s...", pann))
  classify <- grep(pattern="^DF.classifications", x=names(obj@meta.data), value=TRUE)
  obj$pANN <- obj[[pann]]
  obj$DF.classify <- obj[[classify]]
  obj[[pann]] <- NULL
  obj[[classify]] <- NULL
  
  return(obj)
}

# 判断预估的双细胞比例（10x的表格）
estDubRateFun <- function(srt){
    if(!inherits(srt,"Seurat")){stop("input file is not Seurat Object...")}
    rate <- tryCatch({
        if(dim(srt)[2] > 500 && dim(srt)[2] <=1000){
            tmp <- 0.004 + ((dim(srt)[2]-500) * round(0.004/500,digits = 7))
        }else if (dim(srt)[2] > 1000 && dim(srt)[2] <= 2000) {
            tmp <- 0.008 + ((dim(srt)[2]-1000) * round(0.008/1000,digits = 7))
        }else if (dim(srt)[2] > 2000 && dim(srt)[2] <= 3000) {
            tmp <- 0.016 + ((dim(srt)[2]-1000) * round(0.007/1000,digits = 7))
        }else if (dim(srt)[2] > 3000 && dim(srt)[2] <= 4000) {
            tmp <- 0.023 + ((dim(srt)[2]-1000) * round(0.008/1000,digits = 7))
        }else if (dim(srt)[2] > 4000 && dim(srt)[2] <= 5000) {
           tmp <- 0.031 + ((dim(srt)[2]-1000) * round(0.008/1000,digits = 7))
        }else if (dim(srt)[2] > 5000 && dim(srt)[2] <= 6000) {
           tmp <- 0.039 + ((dim(srt)[2]-1000) * round(0.008/1000,digits = 7))
        }else if (dim(srt)[2] > 6000 && dim(srt)[2] <= 7000) {
           tmp <- 0.046 + ((dim(srt)[2]-1000) * round(0.008/1000,digits = 7))
        }else if (dim(srt)[2] > 7000 && dim(srt)[2] <= 8000) {
           tmp <- 0.054 + ((dim(srt)[2]-1000) * round(0.008/1000,digits = 7))
        }else if (dim(srt)[2] > 8000 && dim(srt)[2] <= 9000) {
           tmp <- 0.061 + ((dim(srt)[2]-1000) * round(0.008/1000,digits = 7))
        }else if (dim(srt)[2] > 9000 && dim(srt)[2] <= 10000) {
           tmp <- 0.069 + ((dim(srt)[2]-1000) * round(0.008/1000,digits = 7))
        }else if(dim(srt)[2] >= 10000){
           tmp <- 0.076 + ((dim(srt)[2]-1000) * round(0.008/1000,digits = 7))
        }else {
           return(NULL)
        }
        tmp
    },error = function(e){
        message(e$error)
        return(NULL)
    })
    return(rate)
}

# 对每个样本进行分析
message(paste0(Sys.time(),' 开始进行双细胞检测...'))
sce_list <- lapply(sce_list, FUN = function(x){
    detectDoublet(x,dims = 1:30,estDubRate=estDubRateFun(x),
                                       ncores = 2, SCTransform=F, Homotypic=F, annotation="seurat_clusters")
})
message(paste0(Sys.time(),' 保存DoubletFinder结果文件...'))
saveRDS(sce_list, file = file.path(outputDir,'03_doubletFinderResult.rds'))

# 进行绘图
plist <- lapply(names(sce_list),FUN = function(x){
    Seurat::DimPlot(sce_list[[x]], group.by = "DF.classify") + labs(title = x)
})

png(file.path('./','03_doubletFinderResult.png'),res = 300, units = 'in', height = 4, width = 36)
print(patchwork::wrap_plots(plist, ncol = 9))
dev.off()

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#  使用R包SCP进行去批次和聚类
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# 去除预测的双细胞
scRNA  <- merge(sce_list[[1]],sce_list[-1])
scRNA$new_group <- paste0(scRNA$sample,'_',scRNA$DF.classify)
message(paste0(Sys.time(), " 绘制质控小提琴图... ..."))
p1 <- Seurat::VlnPlot(scRNA.filter, group.by = "new_group", features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3, pt.size = 0.001, raster = TRUE) &
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.title.x = element_blank())

p2 <- Seurat::VlnPlot(scRNA.filter, group.by = "new_group", features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3, pt.size = 0) &
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), axis.title.x = element_blank())
# 保存图片
ggsave(
  filename = file.path('./',"03_DoubletFinderSeuratOnjectVlnplot.png"), height = 8, width = 26,
  plot = patchwork::wrap_plots(list(p1, p2), nrow = 2), limitsize = FALSE, bg = "white"
)
scRNA <- scRNA[,scRNA$DF.classify=="Singlet"] 
# harmony重新聚类
scRNA <- SCP::Integration_SCP(scRNA.filter, batch = 'sample', integration_method = 'Harmony',cluster_resolution = seq(0.1,1.5,0.1))
qsave(scRNA, file = './05_大群数据_scp.qs')
# BBKNN聚类
scRNA <- SCP::Integration_SCP(scRNA.filter, batch = 'sample', integration_method = 'BBKNN',cluster_resolution = seq(0.1,1.5,0.1))
qsave(scRNA, file = './05_大群数据_scp.qs')
# scVI聚类
scRNA <- SCP::Integration_SCP(scRNA.filter, batch = 'sample', integration_method = 'scVI',cluster_resolution = seq(0.1,1.5,0.1))
qsave(scRNA, file = './05_大群数据_scp.qs')
