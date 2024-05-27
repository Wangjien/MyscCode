
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

#坐标轴绘制，修改SCP中的theme_blank函数
theme_blank <- function(
    add_coord = TRUE, xlen_npc = 0.25, ylen_npc = 0.25,
    xlab = "", ylab = "", lab_size = 12, ...) {
  if (isTRUE(add_coord) && isTRUE(xlab != "")) {
    x_space <- lab_size + 2
  } else {
    x_space <- 0
  }
  if (isTRUE(add_coord) && isTRUE(ylab != "")) {
    y_space <- lab_size + 2
  } else {
    y_space <- 0
  }
  args1 <- list(
    panel.border = element_blank(), panel.grid = element_blank(),
    axis.title = element_blank(), axis.line = element_blank(),
    axis.ticks = element_blank(), axis.text = element_blank(),
    legend.background = element_blank(), legend.box.margin = margin(
      0,
      0, 0, 0
    ), legend.margin = margin(0, 0, 0, 0), plot.margin = margin(0,
      0, x_space, y_space,
      unit = "points"
    ), complete = FALSE
  )
  args2 <- as.list(match.call())[-1]
  call.envir <- parent.frame(1)
  args2 <- lapply(args2, function(arg) {
    if (is.symbol(arg)) {
      eval(arg, envir = call.envir)
    } else if (is.call(arg)) {
      eval(arg, envir = call.envir)
    } else {
      arg
    }
  })
  for (n in names(args2)) {
    args1[[n]] <- args2[[n]]
  }
  args <- args1[names(args1) %in% formalArgs(theme)]
  out <- do.call(what = theme, args = args)
  if (isTRUE(add_coord)) {
    g <- grobTree(gList(linesGrob(x = unit(
      c(0, xlen_npc),
      "npc"
    ), y = unit(c(0, 0), "npc"), arrow = arrow(angle = 25, length = unit(
      0.045,
      "npc"
    ), type = "open"), gp = gpar(lwd = 2,fill = "black")), textGrob(
      label = xlab,
      x = unit(0, "npc"), y = unit(0, "npc"), vjust = 4 / 3,
      hjust = 0, gp = gpar(fontsize = lab_size)
    ), linesGrob(x = unit(c(
      0,
      0
    ), "npc"), y = unit(c(0, ylen_npc), "npc"), arrow = arrow(angle = 25, length = unit(
      0.045,
      "npc"
    ), type = "open"), gp = gpar(lwd = 2,fill="black")), textGrob(
      label = ylab,
      x = unit(0, "npc"), y = unit(0, "npc"), vjust = -2 / 3,
      hjust = 0, rot = 90, gp = gpar(fontsize = lab_size)
    )))
    return(list(list(annotation_custom(g)), list(theme_scp() +
      out), list(coord_cartesian(clip = "off"))))
  } else {
    return(list(list(theme_scp() + out)))
  }
}

# 或使用tidydr中的theme_dr()

# 读入数据
setwd('/root/wangje/Project/OveryArtical')
scRNA <- qread('./04_大群数据_scp.qs')

# # 去除中六文章中的3个中间组样本
# scRNA_sub <- scRNA[,scRNA$sample %in% c("HRS421447","HRS421451","HRS421452","HRS421453", 
#     "young1","young3","young4","middle2","middle3","middle4","na1129","old1","old2","old3","na-1010","na-426",
#     "na1111","na1122","na1117","na1128","NA_728","na-412")]

# 1. 绘制不同分组的DimPlot
