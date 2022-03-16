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

###################################################
###---Distribution of Immunoglobulin Isotypes---###
###################################################

#---Determine Proportions of Main Ig Isotypes by Cluster---#
## Exclude IgE and IgD (poorly detected)
tab = table(Idents(ueRA4.1), 
	ueRA4.1$Isotype)[,c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHM")]

## Calculate proportions
for (i in seq(rownames(tab))) {
	# Calculate percent of each isotype for each cluster
	props = c(
		IgM = tab[i, "IGHM"] / sum(tab[i,]) * 100,
		IgG = sum(tab[i, c("IGHG1", "IGHG2", "IGHG3", "IGHG4")]) / sum(tab[i,]) * 100, 
		IgA = sum(tab[i, c("IGHA1", "IGHA2")]) / sum(tab[i,]) * 100)
	print(paste("Cluster", rownames(tab)[i], sep = " "))
	print(round(props, 1))
}

# Plot for cluster 3.1
slices = c(5, 52, 43.9)
pie(slices, labels = NA, col = c("cornflowerblue", "brown4", "tan1"), 
	clockwise = TRUE, main = "Cluster 3.1 (VDJ-Seq)")
text(x = 0.4, y = -0.12, "52.0", cex = 2, col = "white")
text(x = -0.4, y = 0.15, "42", cex = 2, col = "white")
text(x = 0.09, y = 0.55, "5", cex = 2, col = "white")
legend(0.95, 0.35, c("IgM", "IgG", "IgA1"), cex = 2, 
	fill = c("cornflowerblue", "brown4", "tan1"), bty = "n")

#---Determine Proportions of IgG Subtypes---#
for(i in seq(rownames(tab1))) {
	props = c(
		IgG1 = tab[i, "IGHG1"] / sum(tab1[i, c["IGHG1", "IGHG2", "IGHG3", "IGHG4")] * 100,
		IgG2 = tab[i, "IGHG2"] / sum(tab1[i, c["IGHG1", "IGHG2", "IGHG3", "IGHG4")] * 100,
		IgG3 = tab[i, "IGHG3"] / sum(tab1[i, c["IGHG1", "IGHG2", "IGHG3", "IGHG4")] * 100,
		IgG4 = tab[i, "IGHG4"] / sum(tab1[i, c["IGHG1", "IGHG2", "IGHG3", "IGHG4")] * 100)
	print(paste("Cluster", rownames(tab1)[i], sep = " "))
	print(round(props, 1))
}

# Plot for Cluster 3.1
slices = c(73, 12, 15) # No cells in 3.1 express IgG4
pie(slices, labels = NA, col = c("brown4", "indianred3", "indianred1"), 
	clockwise = TRUE, main = "Cluster 3.1 (VDJ-Seq)")
text(x = 0.27, y = -0.18, "73", cex = 2, col = "white")
text(x = -0.5, y = 0.15, "12", cex = 2, col = "white")
text(x = -0.25, y = 0.52, "15", cex = 2, col = "white")
legend(0.95, 0.35, c("IgG1", "IgG2", "IgG3"), cex = 2, 
	fill = c("brown4", "indianred3", "indianred1"), bty = "n")
 
