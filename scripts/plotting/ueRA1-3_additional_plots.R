library(Seurat)
library(ggtree)
library(clustree)
library(ggplot2)
library(Nebulosa)
library(xlsx)

#######################
###---Basic Plots---###
#######################

ueRA1 = readRDS("data/ueRA1")
ueRA2 = readRDS("data/ueRA2")
ueRA3 = readRDS("data/ueRA3")

DimPlot(ueRA1.1, label = TRUE, pt.size = 2, label.size = 7, 
	cols = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#F564E3", "#619CFF")) # Plot clusters

DimPlot(ueRA2.1, label = TRUE, pt.size = 2, label.size = 7, 
	cols = c("#F8766D", "#B79F00", "#00BA38", "#00C08B", "#00B4F0", "#C77CFF", "#F564E3", "#619CFF")) # Plot clusters

DimPlot(ueRA3.1, label = TRUE, pt.size = 2, label.size = 7, 
	cols = c("#F8766D", "#B79F00", "#00BA38", "#00C08B", "#00B4F0", "#C77CFF", "#F564E3")) # Plot clusters

#############################################
###--Examine Expression of Flow Markers---###
#############################################

#---ueRA1---#
DoHeatmap(ueRA1, features = c("TBX21", "ITGAX", "FCRL5", "FCRL3", "FAS", "CD19", "MS4A1", "SELL", "CD24", "CD38", "CXCR5", "CXCR4"), 
	raster = FALSE)

#---ueRA2---#
DoHeatmap(ueRA2, features = c("TBX21", "ITGAX", "FCRL5", "FCRL3", "FAS", "CD19", "MS4A1", "SELL", "CD24", "CD38", "CXCR5", "CXCR4"), 
	raster = FALSE)

#---ueRA3---#
DoHeatmap(ueRA3, features = c("TBX21", "ITGAX", "FCRL5", "FCRL3", "FAS", "CD19", "MS4A1", "SELL", "CD24", "CD38", "CXCR5", "CXCR4"), 
	raster = FALSE)

###################################
###---Hierarchical Clustering---###
###################################

#---ueRA1---#
ueRA1_tree = BuildClusterTree(ueRA1, assay = "SCT")
ueRA1_phytree = Tool(object = ueRA1_tree, slot = "BuildClusterTree")
ggtree(ueRA1_phytree) + geom_tiplab() +
	theme_tree() + 
	xlim(NA, 200) +
	ggtitle("Hierarchical Clustering ueRA1")

#---ueRA2---#
ueRA2_tree = BuildClusterTree(ueRA2, assay = "SCT")
ueRA2_phytree = Tool(object = ueRA2_tree, slot = "BuildClusterTree")
ggtree(ueRA2_phytree) + geom_tiplab() +
	theme_tree() + 
	xlim(NA, 200) +
	ggtitle("Hierarchical Clustering ueRA2")

#---ueRA3---#
ueRA3_tree = BuildClusterTree(ueRA3, assay = "SCT")
ueRA3_phytree = Tool(object = ueRA3_tree, slot = "BuildClusterTree")
ggtree(ueRA3_phytree) + geom_tiplab() +
	theme_tree() + 
	xlim(NA, 200) +
	ggtitle("Hierarchical Clustering ueRA3")

##################################################
###---Examine Differentially Expressed Genes---###
##################################################

#---ueRA1---#
ueRA1$order = factor(Idents(ueRA1))
ueRA1$order = factor(ueRA1$order, levels = c("3", "4", "2", "5", "0", "1")) # set order based on HC
Idents(ueRA1) = ueRA1$order

ueRA1_marks = FindAllMarkers(ueRA1, only.pos = TRUE, min.pct = 0.25) 
ueRA1_marks = subset(ueRA1_marks, subset = p_val_adj < 0.05)
DoHeatmap(ueRA1, features = ueRA1_marks$gene, raster = FALSE)

#---ueRA2---#
ueRA2$order = factor(Idents(ueRA2))
ueRA2$order = factor(ueRA2$order, levels = c("3", "4", "5", "2", "0", "1")) # set order based on HC
Idents(ueRA2) = ueRA2$order

ueRA2_marks = FindAllMarkers(ueRA2, only.pos = TRUE, min.pct = 0.25) 
ueRA2_marks = subset(ueRA2_marks, subset = p_val_adj < 0.05)
DoHeatmap(ueRA2, features = ueRA2_marks$gene, raster = FALSE)

#---ueRA3---#
ueRA3$order = factor(Idents(ueRA3))
ueRA3$order = factor(ueRA3$order, levels = c("3", "4", "2", "0", "1")) # set order based on HC
Idents(ueRA3) = ueRA3$order

ueRA3_marks = FindAllMarkers(ueRA3, only.pos = TRUE, min.pct = 0.25) 
ueRA3_marks = subset(ueRA3_marks, subset = p_val_adj < 0.05)
DoHeatmap(ueRA3, features = ueRA3_marks$gene, raster = FALSE)

###################################################################
###---Examine Combined Expression of Upregulated Flow Markers---###
###################################################################

#---ueRA1---#
plot_density(ueRA1, features = c("TBX21", "ITGAX", "FCRL5", "FCRL3", "FAS", "CD19", "MS4A1"), 
	reduction = "umap", joint = TRUE, combine = FALSE, size = 2)[8] # plot combined density

#---ueRA2---#
plot_density(ueRA2, features = c("TBX21", "ITGAX", "FCRL5", "FCRL3", "FAS", "CD19", "MS4A1"), 
	reduction = "umap", joint = TRUE, combine = FALSE, size = 2)[8] # plot combined density

#---ueRA3---#
plot_density(ueRA3, features = c("TBX21", "ITGAX", "FCRL5", "FCRL3", "FAS", "CD19", "MS4A1"), 
	reduction = "umap", joint = TRUE, combine = FALSE, size = 2)[8] # plot combined density

###########################################################
###---Examine Expression of ueRA4 Cluster 3.1 Markers---###
###########################################################

ueRA1.1 = readRDS("data/ueRA1_w_subs")
ueRA2.1 = readRDS("data/ueRA2_w_subs")
ueRA3.1 = readRDS("data/ueRA3_w_subs")
ueRA4.1 = readRDS("data/ueRA4_w_subs")

ueRA4_3.1_marks = read.xlsx("reports/ueRA4_3.1_vs_3.0_DEGs.xlsx", sheetName = "DEGs", 
	col.names = TRUE, row.names = FALSE) # read in list of genes

#---ueRA1---#
ueRA1.1$order = factor(Idents(ueRA1.1))
ueRA1.1$order = factor(ueRA1.1$order, levels = c("3.1", "3.0", "0", "1", "2", "4", "5")) # set order
Idents(ueRA1.1) = ueRA1.1$order

DoHeatmap(ueRA1.1, features = ueRA4_3.1_marks$gene, 
	group.colors = c("#00B4F0", "#00C08B", "#F8766D", "#B79F00", "#00BA38", "#F564E3", "#619CFF"), 
	raster = FALSE)

#---ueRA2---#
ueRA2.1$order = factor(Idents(ueRA2.1))
ueRA2.1$order = factor(ueRA2.1$order, levels = c("3.0", "3.1", "3.2", "0", "1", "2", "4", "5")) # set order based on HC

DoHeatmap(ueRA2.1, features = ueRA4_3.1_marks$gene, 
	group.colors = c("#00C08B", "#00B4F0", "#C77CFF", "#F8766D", "#B79F00", "#00BA38", "#619CFF", "#C77CFF"), 
	raster = FALSE, group.by = ueRA2.1$order)

#---ueRA3---#
ueRA3.1$order = factor(Idents(ueRA3.1))
ueRA3.1$order = factor(ueRA3.1$order, levels = c("3.2", "3.1", "3.0", "0", "1", "2", "4")) # set order based on HC

DoHeatmap(ueRA3.1, features = ueRA4_3.1_marks$gene, 
	group.colors = c("#C77CFF", "#00B4F0", "#00C08B", "#F8766D", "#B79F00", "#00BA38", "#F564E3"), 
	raster = FALSE, group.by = ueRA3.1$order)