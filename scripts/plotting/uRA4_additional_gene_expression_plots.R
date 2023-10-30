library(Seurat)
library(xlsx)
library(ggplot2)

############################################
###---Additional Gene Expression Plots---###
############################################

uRA4.1 = readRDS("data/uRA4_w_subs")

#---PCs and Other Genes---#
DoHeatmap(uRA4.1, features = c("PRDM1", "IRF4", "XBP1", "IL6R", "CD79A", "SYK", 
	"TNFRSF1B", "ITGB2", "CD79B", "IL21R", "IFNGR1", "IFNGR2"), raster = FALSE, 
	group.colors = c("#F8766D", "#B79F00", "#00BA38", "#00C08B", "#00B4F0", "#C77CFF", 
	"#F564E3", "#619CFF"))

#---MHC Genes---#

DoHeatmap(uRA4.1, features = c("HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DPA1", "HLA-DPB1", 
	"HLA-DQA1", "HLA-DQB1", "CD74", "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-A", 
	"HLA-B", "HLA-C", "HLA-E", "HLA-F", "B2M", CD80", "CD86"), raster = FALSE, 
	group.colors = c("#F8766D", "#B79F00", "#00BA38", "#00C08B", "#00B4F0", "#C77CFF", 
	"#F564E3", "#619CFF"))

#---Cluster 3.1 Vs. 3.0 Markers---#

genes = read.xlsx("reports/uRA4_3.1_vs_3.0_DEGs.xlsx", 
	sheetName = "DEGs", col.names = TRUE, row.names = FALSE)	

DoHeatmap(uRA4.1, features = genes$gene, raster = FALSE, 
	group.colors = c("#F8766D", "#B79F00", "#00BA38", "#00C08B", "#00B4F0", 
	"#C77CFF", "#F564E3", "#619CFF")) + 
	ggtitle("up DEGs in 3.1 vs. 3.0 (n=95)") + 
	theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_blank())
 
