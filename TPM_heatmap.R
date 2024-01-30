# 安装和加载必要的包
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
library(pheatmap)
library(readxl)
library(dplyr)

# 设置工作路径
setwd("E:/sigma因子的研究/组学数据/衍生处理数据")

# 1. 导入数据
data <- read_excel("E:/sigma因子的研究/组学数据/衍生处理数据/转录组热图聚类TPM.xlsx", sheet = 1)
# 删掉gene列
data <- data[-1]
# 对每个值取对数转换：log10(TPM + 1)
data_log_transformed <- log10(data + 1)

# 确保breaks是唯一的
breaks <- seq(1, 3, length.out = 100)

png(filename = "TPM_heatmap.png", width = 7, height = 12, units = 'in', res = 300, bg = "transparent")

pheatmap(data_log_transformed, 
         breaks = breaks,
         show_rownames = FALSE,
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(c("#3638F0", "white", "#F02E35"))(length(breaks) - 1),
         border_color = NA) 

dev.off()