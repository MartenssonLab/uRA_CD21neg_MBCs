library(Seurat)
library(xlsx)
library(ggplot2)
library(ggtree)
library(clustree)
library(gridExtra)
library(plotrix)

#######################
###---Basic Plots---###
#######################
ueRA4 = readRDS("data/ueRA4_base")

DimPlot(ueRA4, label = TRUE, pt.size = 2, label.size = 7, reduction = "umap",
	cols = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#F564E3")) # Plot clusters in UMAP space

#---Pie Chart of Cell Numbers---#
table(Idents(ueRA4)) # Number of cells per cluster

slices = c(889, 827, 713, 443, 37)
pie(slices, labels = NA, col = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#F564E3"), 
	clockwise = TRUE, main = "number of cells per cluster")

draw.circle(0, 0, 0.35, col = "white")
text(x = 0, y = 0, "2 909", cex = 2)
text(x = 0.45, y = 0.35, "889", cex = 2)
text(x = 0.2, y = -0.55, "827", cex = 2)
text(x = -0.53, y = -0.1, "713", cex = 2)
text(x = -0.3, y = 0.45, "443", cex = 2)
text(x = -0.05, y = 0.9, "37", cex = 2)

mylab = c("0", "1", "2", "3", "4")
myorder = matrix(c(1, 2, 3, 4, 5), nrow = 1, ncol = 5, byrow = F)
legend(-0.9, 1.15, mylab[myorder], cex = 1.5, 
	fill = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#F564E3"), bty = "n", ncol = 5)


###################################
###---Hierarchical Clustering---###
###################################

tree = BuildClusterTree(ueRA4, assay = "SCT")
myphytree = Tool(tree, slot = "BuildClusterTree")
ggtree(myphytree) + 
	geom_tiplab() +
	theme_tree() +
	xlim(NA, 200) + 
	ggtitle("Hierarchical Clustering") # Plot dendrogram

##########################################
###---Differentially Expressed Genes---###
##########################################

degs = read.xlsx("data/ueRA4_DEGs.xlsx", sheetName = "Upregulated", col.names = TRUE, row.names = FALSE) # Read in DEGs

DoHeatmap(ueRA4, features = degs$gene, raster = FALSE, 
	group.colors = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#F564E3")) + 
	theme(axis.text.y = element_blank())

#---Pie Chart---#
# Identify number of significantly upregulated DEGs by cluster
for (i in levels(ueRA4)) {
  d = nrow(subset(degs, subset = cluster == i))
  print(paste("Cluster", i, sep = ""))
  print(d)
}

slices = c(73, 68, 73, 245, 68)
pie(slices, labels = NA, col = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#F564E3"), 
	clockwise = TRUE, main = "upregulated DEGs by cluster", cex = 3)
text(x = 0.25, y = 0.5, "73", cex = 2.5, col = "black")
text(x = 0.5, y = 0.13, "68", cex = 2.5, col = "black")
text(x = 0.5, y = -0.3, "73", cex = 2.5, col = "black")
text(x = -0.35, y = -0.3, "245", cex = 2.5, col = "black")
text(x = -0.22, y = 0.5, "68", cex = 2.5, col = "black")

mylab = c("0", "1", "2", "3", "4")
myorder = matrix(c(1, 2, 3, 4, 5), nrow = 1, ncol = 5, byrow = F)
legend(-0.9, 1.1, mylab[myorder], cex = 1.5, fill = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#F564E3"), bty = "n", ncol = 5)