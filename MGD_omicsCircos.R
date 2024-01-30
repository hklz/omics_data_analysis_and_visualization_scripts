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
Genes <- read.table("E:/sigma因子的研究/组学数据/衍生处理数据/Genes_for_circos.txt", sep='\t', header = T)
labels <- read.table("E:/sigma因子的研究/组学数据/衍生处理数据/DE_origin/Label_MGD.txt", sep='\t', header = T)
MGD_TPM <- read.table("E:/sigma因子的研究/组学数据/衍生处理数据/DE_origin/MGD_TPM.txt", sep='\t', header = T)
MGD_up <- read.table("E:/sigma因子的研究/组学数据/衍生处理数据/DE_origin/MGD_down.txt", sep='\t', header = T)
MGD_down <- read.table("E:/sigma因子的研究/组学数据/衍生处理数据/DE_origin/MGD_up.txt", sep='\t', header = T)
db_Temp <- segAnglePo(Genes,seg = Genes$GeneId);
colors <- c('#F94144', '#F3722C', '#F8961E', 
            '#F9844A','#F9C74F','#90BE6D',
            '#43AA8B', '#4D908E','#577590', '#277DA1')
color_palette <- colorRampPalette(colors)(length(Genes$GeneId))
###########################

# Define the color scale
#pdffile <- "all_condition_direct.pdf";
pdffile <- "MGD.pdf";
pdf(pdffile, 40, 40);
plot(c(-2000,2000), c(-2000,2000), type="n", axes=FALSE, xlab="", ylab="", main="");
par(mar=c(5, 5, 5, 5)); #sets the margins of the plot 

#表达量
circos(R=400,cir=db_Temp,W=400,mapping=MGD_TPM,col.v=4,type="b",B=FALSE,col='#FF9800',lwd=0.5, scale = TRUE);
circos(R=0,cir=db_Temp,W=600,mapping=MGD_TPM,col.v=4, type="b3", B=FALSE,col='white',lwd=2, scale = FALSE);
#相关性基因
circos(R=475, cir=db_Temp, W=20, mapping=MGD_up, type="link", lwd=0.1, col='#F5655F');
circos(R=475, cir=db_Temp, W=20, mapping=MGD_down, type="link", lwd=0.1, col='#3AA0F2');
#外层基因组
circos(R=475, cir=db_Temp, type="chr", col=color_palette, print.chr.lab=FALSE, W=100, scale=FALSE);

#标签
circos(R = 680, cir =db_Temp, W = 40, mapping =labels, type = "label", side = "out", col = "black");

dev.off()