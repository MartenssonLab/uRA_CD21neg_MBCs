library(Seurat)
library(xlsx)
library(ggplot2)
library(ggtree)
library(clustree)

###################################
###---Hierarchical Clustering---###
###################################

ueRA4 = readRDS("data/ueRA4_base")

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

degs = read.xlsx("data/ueRA4_DEGs.xlsx", sheetName = "Upregulated", col.names = TRUE, row.names = FALSE)
DoHeatmap(ueRA4, features = degs$gene, raster = FALSE, 
	group.colors = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#F564E3")) # Plot significant upregulated DEGs

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
