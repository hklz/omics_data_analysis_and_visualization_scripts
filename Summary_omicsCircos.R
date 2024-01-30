# Load the OmicCircos library
library(OmicCircos)
library(colorspace)
library(Cairo)
library(RColorBrewer)
library(fields)
library(shiny)
library(OmicCircos)
setwd("E:/sigma因子的研究/数据处理作图")

# Read the data
all_TPM <- read.table("E:/sigma因子的研究/组学数据/衍生处理数据/all_TPM_for_circos_lg_size_4.txt", sep='\t', header = T)
Genes <- read.table("E:/sigma因子的研究/组学数据/衍生处理数据/Genes_for_circos.txt", sep='\t', header = T)
sigmafactor_labels <- read.table("E:/sigma因子的研究/组学数据/衍生处理数据/DE_origin/labels_for_circos.txt", sep='\t', header = T)
hot_labels <- read.table("E:/sigma因子的研究/组学数据/衍生处理数据/DE_origin/hot_labels_circos.txt", sep='\t', header = T)
MGI_up <- read.table("E:/sigma因子的研究/组学数据/衍生处理数据/DE_origin/MGI_up.txt", sep='\t', header = T)
MGH_up <- read.table("E:/sigma因子的研究/组学数据/衍生处理数据/DE_origin/MGH_up.txt", sep='\t', header = T)
MGE_up <- read.table("E:/sigma因子的研究/组学数据/衍生处理数据/DE_origin/MGE_up.txt", sep='\t', header = T)
MGF_up <- read.table("E:/sigma因子的研究/组学数据/衍生处理数据/DE_origin/MGF_up.txt", sep='\t', header = T)
MGD_up <- read.table("E:/sigma因子的研究/组学数据/衍生处理数据/DE_origin/MGD_up.txt", sep='\t', header = T)
MGS_up <- read.table("E:/sigma因子的研究/组学数据/衍生处理数据/DE_origin/MGS_up.txt", sep='\t', header = T)
MGN_up <- read.table("E:/sigma因子的研究/组学数据/衍生处理数据/DE_origin/MGH_up.txt", sep='\t', header = T)



db_Temp <- segAnglePo(Genes,seg = Genes$GeneId);
colors <- c('#F94144', '#F3722C', '#F8961E', 
            '#F9844A','#F9C74F','#90BE6D',
            '#43AA8B', '#4D908E','#577590', '#277DA1')
color_palette <- colorRampPalette(colors)(length(Genes$GeneId))
###########################

# Define the color scale
#pdffile <- "all_condition_direct.pdf";
pdffile <- "All_omics_up.pdf";
pdf(pdffile, 20, 20);
plot(c(-1000,1000), c(-1000,1000), type="n", axes=FALSE, xlab="", ylab="", main="");
#par(mar=c(5, 5, 5, 5)); #sets the margins of the plot 


#相关性基因
circos(R=475, cir=db_Temp, W=20, mapping=MGI_up, type="link", lwd=0.1, col='#CDDC39');
circos(R=475, cir=db_Temp, W=20, mapping=MGH_up, type="link", lwd=0.1, col='#9C27B0');
circos(R=475, cir=db_Temp, W=20, mapping=MGS_up, type="link", lwd=0.1, col='#E91E63');
circos(R=475, cir=db_Temp, W=20, mapping=MGE_up, type="link", lwd=0.1, col='#2196F3');
circos(R=475, cir=db_Temp, W=20, mapping=MGD_up, type="link", lwd=0.1, col='#FF9800');
circos(R=475, cir=db_Temp, W=20, mapping=MGF_up, type="link", lwd=0.1, col='#F44336');
circos(R=475, cir=db_Temp, W=20, mapping=MGN_up, type="link", lwd=0.1, col='#4CAF50');

#外层基因组
circos(R=578, cir=db_Temp, type="chr", col=color_palette, print.chr.lab=FALSE, W=100, scale=FALSE);

#标签
circos(R = 580, cir =db_Temp, W = 40, mapping = sigmafactor_labels, type = "label", side = "out", col = "black");
circos(R = 573, cir =db_Temp, W = 200, mapping = hot_labels, type = "label", side = "out",lwd=0.5, col = "red");
#表达量
circos(R=480,cir=db_Temp,W=100,mapping=all_TPM,type="heatmap",lwd=0.5, col.bar = FALSE);

#区域标明

dev.off()