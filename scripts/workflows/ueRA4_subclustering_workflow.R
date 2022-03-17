library((Seurat)
library(dplyr)
library(xlsx)

##########################
###---Sub-Clustering---###
##########################

ueRA4.3 = subset(ueRA4, subset = seurat_clusters == "3") # Isolate cells in cluster 3

#---Single-Cell Transform---#
ueRA4.3 = SCTransform(ueRA4.3, variable.features.n = dim(ueRA4.3)[1], assay = "RNA")

#---Principal Component Analysis---#
ueRA4.3 = RunPCA(ueRA4.3, features = VariableFeatures(ueRA4.3))
# ElbowPlot(ueRA4.3, reduction = "pca", ndims = 50)

#---Unsupervised Re-Clustering---#
ueRA4.3 = FindNeighbors(ueRA4.3, dims = 1:40)
ueRA4.3 = FindClusters(ueRA4.3, resolution = 0.9)

#---Running UMAP---#
ueRA4.3 = RunUMAP(ueRA4.3, dims = 1:40)
# UMAPPlot(ueRA4.3, label = TRUE)

#---Rename Clusters Accordingly---#
c3.0 = WhichCells(ueRA4.3, idents = "0")
c3.1 = WhichCells(ueRA4.3, idents = "1")
c3.2 = WhichCells(ueRA4.3, idents = "2")

Idents(ueRA4.3, cells = c3.0) = "3.0"
Idents(ueRA4.3, cells = c3.1) = "3.1"
Idents(ueRA4.3, cells = c3.2) = "3.2"

########################################
###---Differential Gene Expression---###
########################################

ueRA4.3.marks = subset(FindAllMarkers(ueRA4.3, only.pos = FALSE, min.pct = 0.25), 
	subset = p_val_adj < 0.05)
pos.marks = subset(ueRA4.3.marks, subset = avg_log2FC > 0.0)
neg.marks = subset(ueRA4.3.marks, subset = avg_log2FC < 0.0)

write.xlsx(pos.marks %>% group_by(cluster) %>% top_n(-10, p_val_adj), 
	"reports/ueRA4_C3_subs.xlsx", sheetName = "Top10_Up", col.names = TRUE, row.names = FALSE)
write.xlsx(neg.marks %>% group_by(cluster) %>% top_n(-10, p_val_adj), )
	"reports/ueRA4_C3_subs.xlsx", sheetName = "Top10_Down", col.names = TRUE, row.names = FALSE, append = TRUE)

saveRDS(ueRA4.3, file = "data/ueRA4_C3_subs") # Save sub-cluster object

#########################################
###---Map Sub-Clusters to Base Data---###
#########################################

ueRA4.1 = ueRA4
Idents(ueRA4.1, cells = c3.0) = "3.0"
Idents(ueRA4.1, cells = c3.1) = "3.1"
Idents(ueRA4.1, cells = c3.2) = "3.2"
Idents(ueRA4.1) = factor(Idents(ueRA4.1), 
	levels = c("0", "1", "2", "3.0", "3.1", "3.2", "4")

###############################
###---Cluster 3.1 Vs. 3.0---###
###############################

marks = subset(FindMarkers(ueRA4.1, ident.1 = "3.1", ident.2 = "3.0", min.pct = 0.25, only.pos = TRUE), 
	subset = p_val_adj <0.05)
write.xlsx(marks, file = "reports/ueRA4_3.1_vs_3.0_DEGs.xlsx", sheetName = "DEGs", col.names = TRUE, row.names = FALSE)

saveRDS(ueRA4.1, file = "data/ueRA4_w_subs") # Save Seurat object
