library(Seurat)
library(ggtree)
library(clustree)
library(ggplot2)
library(Nebulosa)
library(xlsx)

#######################
###---Basic Plots---###
#######################

uRA1 = readRDS("data/uRA1_base")
uRA2 = readRDS("data/uRA2_base")
uRA3 = readRDS("data/uRA3_base")

DimPlot(uRA1.1, label = TRUE, pt.size = 2, label.size = 7, 
	cols = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#F564E3", "#619CFF")) # Plot clusters

DimPlot(uRA2.1, label = TRUE, pt.size = 2, label.size = 7, 
	cols = c("#F8766D", "#B79F00", "#00BA38", "#00C08B", "#00B4F0", "#C77CFF", "#F564E3", "#619CFF")) # Plot clusters

DimPlot(uRA3.1, label = TRUE, pt.size = 2, label.size = 7, 
	cols = c("#F8766D", "#B79F00", "#00BA38", "#00C08B", "#00B4F0", "#C77CFF", "#F564E3")) # Plot clusters

#############################################
###--Examine Expression of Flow Markers---###
#############################################

#---uRA1---#
DoHeatmap(uRA1, features = c("TBX21", "ITGAX", "FCRL5", "FCRL3", "FAS", "CD19", "MS4A1", "SELL", "CD24", "CD38", "CXCR5", "CXCR4"), 
	raster = FALSE)

#---uRA2---#
DoHeatmap(uRA2, features = c("TBX21", "ITGAX", "FCRL5", "FCRL3", "FAS", "CD19", "MS4A1", "SELL", "CD24", "CD38", "CXCR5", "CXCR4"), 
	raster = FALSE)

#---uRA3---#
DoHeatmap(uRA3, features = c("TBX21", "ITGAX", "FCRL5", "FCRL3", "FAS", "CD19", "MS4A1", "SELL", "CD24", "CD38", "CXCR5", "CXCR4"), 
	raster = FALSE)

###################################
###---Hierarchical Clustering---###
###################################

#---uRA1---#
uRA1_tree = BuildClusterTree(uRA1, assay = "SCT")
uRA1_phytree = Tool(object = uRA1_tree, slot = "BuildClusterTree")
ggtree(uRA1_phytree) + geom_tiplab() +
	theme_tree() + 
	xlim(NA, 200) +
	ggtitle("Hierarchical Clustering uRA1")

#---uRA2---#
uRA2_tree = BuildClusterTree(uRA2, assay = "SCT")
uRA2_phytree = Tool(object = uRA2_tree, slot = "BuildClusterTree")
ggtree(uRA2_phytree) + geom_tiplab() +
	theme_tree() + 
	xlim(NA, 200) +
	ggtitle("Hierarchical Clustering uRA2")

#---uRA3---#
uRA3_tree = BuildClusterTree(uRA3, assay = "SCT")
uRA3_phytree = Tool(object = uRA3_tree, slot = "BuildClusterTree")
ggtree(uRA3_phytree) + geom_tiplab() +
	theme_tree() + 
	xlim(NA, 200) +
	ggtitle("Hierarchical Clustering uRA3")

##################################################
###---Examine Differentially Expressed Genes---###
##################################################

#---uRA1---#
uRA1$order = factor(Idents(uRA1))
uRA1$order = factor(uRA1$order, levels = c("3", "4", "2", "5", "0", "1")) # set order based on HC
Idents(uRA1) = uRA1$order

uRA1_marks = FindAllMarkers(uRA1, only.pos = TRUE, min.pct = 0.25) 
uRA1_marks = subset(uRA1_marks, subset = p_val_adj < 0.05)
DoHeatmap(uRA1, features = uRA1_marks$gene, raster = FALSE)

#---uRA2---#
uRA2$order = factor(Idents(uRA2))
uRA2$order = factor(uRA2$order, levels = c("3", "4", "5", "2", "0", "1")) # set order based on HC
Idents(uRA2) = uRA2$order

uRA2_marks = FindAllMarkers(uRA2, only.pos = TRUE, min.pct = 0.25) 
uRA2_marks = subset(uRA2_marks, subset = p_val_adj < 0.05)
DoHeatmap(uRA2, features = uRA2_marks$gene, raster = FALSE)

#---uRA3---#
uRA3$order = factor(Idents(uRA3))
uRA3$order = factor(uRA3$order, levels = c("3", "4", "2", "0", "1")) # set order based on HC
Idents(uRA3) = uRA3$order

uRA3_marks = FindAllMarkers(uRA3, only.pos = TRUE, min.pct = 0.25) 
uRA3_marks = subset(uRA3_marks, subset = p_val_adj < 0.05)
DoHeatmap(uRA3, features = uRA3_marks$gene, raster = FALSE)

###################################################################
###---Examine Combined Expression of Upregulated Flow Markers---###
###################################################################

#---uRA1---#
plot_density(uRA1, features = c("TBX21", "ITGAX", "FCRL5", "FCRL3", "FAS", "CD19", "MS4A1"), 
	reduction = "umap", joint = TRUE, combine = FALSE, size = 2)[8] # plot combined density

#---uRA2---#
plot_density(uRA2, features = c("TBX21", "ITGAX", "FCRL5", "FCRL3", "FAS", "CD19", "MS4A1"), 
	reduction = "umap", joint = TRUE, combine = FALSE, size = 2)[8] # plot combined density

#---uRA3---#
plot_density(uRA3, features = c("TBX21", "ITGAX", "FCRL5", "FCRL3", "FAS", "CD19", "MS4A1"), 
	reduction = "umap", joint = TRUE, combine = FALSE, size = 2)[8] # plot combined density

###########################################################
###---Examine Expression of uRA4 Cluster 3.1 Markers---###
###########################################################

uRA1.1 = readRDS("data/uRA1_w_subs")
uRA2.1 = readRDS("data/uRA2_w_subs")
uRA3.1 = readRDS("data/uRA3_w_subs")
uRA4.1 = readRDS("data/uRA4_w_subs")

uRA4_3.1_marks = read.xlsx("reports/uRA4_3.1_vs_3.0_DEGs.xlsx", sheetName = "DEGs", 
	col.names = TRUE, row.names = FALSE) # read in list of genes

#---uRA1---#
uRA1.1$order = factor(Idents(uRA1.1))
uRA1.1$order = factor(uRA1.1$order, levels = c("3.1", "3.0", "0", "1", "2", "4", "5")) # set order
Idents(uRA1.1) = uRA1.1$order

DoHeatmap(uRA1.1, features = uRA4_3.1_marks$gene, 
	group.colors = c("#00B4F0", "#00C08B", "#F8766D", "#B79F00", "#00BA38", "#F564E3", "#619CFF"), 
	raster = FALSE)

#---uRA2---#
uRA2.1$order = factor(Idents(uRA2.1))
uRA2.1$order = factor(uRA2.1$order, levels = c("3.0", "3.1", "3.2", "0", "1", "2", "4", "5")) # set order based on HC

DoHeatmap(uRA2.1, features = uRA4_3.1_marks$gene, 
	group.colors = c("#00C08B", "#00B4F0", "#C77CFF", "#F8766D", "#B79F00", "#00BA38", "#619CFF", "#C77CFF"), 
	raster = FALSE, group.by = uRA2.1$order)

#---uRA3---#
uRA3.1$order = factor(Idents(uRA3.1))
uRA3.1$order = factor(uRA3.1$order, levels = c("3.2", "3.1", "3.0", "0", "1", "2", "4")) # set order based on HC

DoHeatmap(uRA3.1, features = uRA4_3.1_marks$gene, 
	group.colors = c("#C77CFF", "#00B4F0", "#00C08B", "#F8766D", "#B79F00", "#00BA38", "#F564E3"), 
	raster = FALSE, group.by = uRA3.1$order)
