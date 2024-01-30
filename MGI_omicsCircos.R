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
#Temp_exp <- read.table("Temp_cond_all.txt", sep='\t', header = T)
all_TPM <- read.table("E:/sigma因子的研究/组学数据/衍生处理数据/all_TPM_for_circos_lg.txt", sep='\t', header = T)
Genes <- read.table("E:/sigma因子的研究/组学数据/衍生处理数据/Genes_for_circos.txt", sep='\t', header = T)
labels <- read.table("E:/sigma因子的研究/组学数据/衍生处理数据/DE_origin/Label_MGI.txt", sep='\t', header = T)
MGI_TPM <- read.table("E:/sigma因子的研究/组学数据/衍生处理数据/DE_origin/MGI_TPM.txt", sep='\t', header = T, fileEncoding="UTF-16LE")
MGI_up <- read.table("E:/sigma因子的研究/组学数据/衍生处理数据/DE_origin/MGI_down.txt", sep='\t', header = T)
MGI_down <- read.table("E:/sigma因子的研究/组学数据/衍生处理数据/DE_origin/MGI_up.txt", sep='\t', header = T)
db_Temp <- segAnglePo(Genes,seg = Genes$GeneId);
colors <- c('#F94144', '#F3722C', '#F8961E', 
            '#F9844A','#F9C74F','#90BE6D',
            '#43AA8B', '#4D908E','#577590', '#277DA1')
color_palette <- colorRampPalette(colors)(length(Genes$GeneId))
###########################

# Define the color scale
#pdffile <- "all_condition_direct.pdf";
pdffile <- "MGI.pdf";
pdf(pdffile, 40, 40);
plot(c(-2000,2000), c(-2000,2000), type="n", axes=FALSE, xlab="", ylab="", main="");
par(mar=c(5, 5, 5, 5)); #sets the margins of the plot 

#表达量
circos(R=400,cir=db_Temp,W=400,mapping=MGI_TPM,col.v=4,type="b",B=FALSE,col='#CDDC39',lwd=0.5, scale = TRUE);
circos(R=0,cir=db_Temp,W=600,mapping=MGI_TPM,col.v=4, type="b3", B=FALSE,col='white',lwd=2, scale = FALSE);
#相关性基因
circos(R=475, cir=db_Temp, W=20, mapping=MGI_up, type="link", lwd=0.1, col='#F5655F');
circos(R=475, cir=db_Temp, W=20, mapping=MGI_down, type="link", lwd=0.1, col='#3AA0F2');
#外层基因组
circos(R=475, cir=db_Temp, type="chr", col=color_palette, print.chr.lab=FALSE, W=100, scale=FALSE);

#标签
circos(R = 680, cir =db_Temp, W = 40, mapping =labels, type = "label", side = "out", col = "black");

dev.off()