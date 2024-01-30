library(pheatmap)
library(readxl)
library(dplyr)
 
# 设置工作路径
setwd("E:/sigma因子的研究/导出图片")

# 导入数据
data <- read_excel("E:/sigma因子的研究/组学数据/衍生处理数据/转录组通路图.xlsx", sheet = 2)

# 对每个值取对数转换：log10(TPM + 1)
data_log_transformed <- log10(data[-1] + 1)

rownames(data_log_transformed) <- data$Index


# 确保breaks是唯一的
breaks <- seq(1, 3, length.out = 100)

# 创建热图
png(filename = "figure3_heatmap_with_row_gap.png", width = 8, height = 20, units = 'in', res = 300, bg = "transparent")

pheatmap(data_log_transformed, 
         breaks = breaks,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         cellwidth = 20,
         cellheight = 20,
         gaps_row = seq(1, nrow(data_log_transformed)),  # 设置为大于0的值，以在每行之间添加间隔
         color = colorRampPalette(c("#3638F0", "white", "#F02E35"))(length(breaks) - 1),
         border_color = "black",
         fontsize = 20,
         fontfamily= "serif"
         #fontface = "italic"
         ) 

dev.off()
