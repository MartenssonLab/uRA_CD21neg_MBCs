library(Seurat)
library(ggplot2)
library(Nebulosa)
library(gridExtra)

####################################
###---Gene Expression Analysis---###
####################################

ueRA4 = readRDS("data/ueRA4_base")

genes = c("CR2","CD27","TBX21","ITGAX","FCRL5","FCRL3","FAS",
	"CD19","MS4A1","SELL","CD24","CD38","CXCR5","CXCR4") # predefined markers

#---Heatmap---#
DoHeatmap(ueRA4, features = genes, raster = FALSE, 
	group.colors = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#F564E3"))

#---Density Plots---#
for(i in 1:14){
  p = plot_density(ra2.1, features = genes[i], reduction = "umap") 
  assign(paste0("p", i), p)
} # Make density plots

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, 
             p10, p11, p12, p13, p14, ncol = 7)

#---Feature Plots---#
cd21 = WhichCells(ueRA4, expression = CR2 > 0)) # Extract cells expressing CD21
cd27 = WhichCells(ueRA4, expression = CD27 > 0)) # Extract cells expressing CD27

# Plot CD21+ cells
DimPlot(ueRA4, cells.highlight = cd21, cols.highlight = "red", label = TRUE, pt.size = 2, 
	label.size = 7, sizes.highlight = 2, order = TRUE) +
	ggtitle("CD21") + theme(plot.title = element_text(hjust = 0.5)) + NoLegend()

# Plot CD27+ cells
DimPlot(ueRA4, cells.highlight = cd27, cols.highlight = "palecioletred3", label = TRUE, pt.size = 2, 
	label.size = 7, sizes.highlight = 2, order = TRUE) + 
	ggtitle("CD27") + theme(plot.title = element_text(hjust = 0.5)) + NoLegend()
 
