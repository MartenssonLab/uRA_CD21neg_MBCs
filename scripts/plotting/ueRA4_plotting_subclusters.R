library(Seurat)
library((ggtree)
library(clustree)
library(ggplot2)

#######################
###---Basic Plots---###
#######################

ueRA4.1 = readRDS("data/ueRA4_w_subs")

DimPlot(ueRA4.1, label = TRUE, pt.size = 2, label.size = 7, 
	cols = c("#F8766D", "#B79F00", "#00BA38", "#00C08B", "#00B4F0", "#C77CFF", "#F564E3", "#619CFF")) # Plot clusters

###################################
###---Hierarchical Clustering---###
###################################

## C3 Subclusters
ueRA4.3 = readRDS("ueRA4_C3_subs")

tree = BuildClusterTree(ueRA4.3, assay = "SCT")
myphytree = Tool(tree, slot = "BuildClusterTree")
ggtree(myphytree) + 
	geom_tiplab() + 
	theme_tree() + 
	xlim(NA, 200) + 
	ggtitle("Hierarchical Clustering") # Plot dendrogram

## All clusters
tree = BuildClusterTree(ueRA4.1, assay = "SCT")
myphytree = Tool(tree, slot = "BuildClusterTree")
ggtree(myphytree) + 
	geom_tiplab() + 
	theme_tree() + 
	xlim(NA, 200) + 
	ggtitle("Hierarchical Clustering") # Plot dendrogram


####################################
###---Gene Expression Analysis---###
####################################

#---Top10 Differentially Expressed Genes---#
top10_up_degs = read.xlsx("reports/ueRA4_C3_subs.xlsx", 
	sheetName = "Top10_Up", col.names = TRUE, row.names = FALSE)
top10_down_degs = read.xlsx("reports/ueRA4_C3_subs.xlsx", 
	sheetName = "Top10_Down", col.names = TRUE, row.names = FALSE)

DoHeatmap(ueRA4.3, features = top10_up_degs$gene, raster = FALSE, 
	group.colors = c("#00C08B", "#00B4F0", "#C77CFF")) # Plot top10 up DEGs

DoHeatmap(ueRA4.3, features = top10_down_degs$gene, raster = FALSE, 
	group.colors = c("#00C08B", "#00B4F0", "#C77CFF")) # Plot top10 down DEGs

#---Predefined Markers---#

genes = c("CR2","CD27","TBX21","ITGAX","FCRL5","FCRL3","FAS",
	"CD19","MS4A1","SELL","CD24","CD38","CXCR5","CXCR4") # predefined markers

DoHeatmap(ueRA4.3, features = genes, raster = FALSE, 
	group.colors = c("#00C08B", "#00B4F0", "#C77CFF"))

###########################
###---% of Total MBCs---###
###########################

## Calculate the % of total MBCs contained in cluster 3.1
p = round(table(Idents(ueRA4.1)) / sum(table(Idents(ueRA4.1))

slices = c(p[5], 100 - p[5])
pie(slices, labels = NA, col = c("lightblue3", "lightgray"), clockwise = TRUE, 
	main = "Cluster 3.1", cex = 2)
text(x = 0.37, y = 0.35, paste(p[5], "%", sep = ""), cex = 2)
text(x = -0.05, y = -0.25, "RNA Seq", cex = 2)