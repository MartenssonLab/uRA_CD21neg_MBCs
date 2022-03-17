library((Seurat)
library(dplyr)
library(xlsx)

##########################
###---Sub-Clustering---###
##########################

ueRA2.3 = subset(ueRA2, subset = seurat_clusters == "3") # Isolate cells in cluster 3

#---Single-Cell Transform---#
ueRA2.3 = SCTransform(ueRA2.3, variable.features.n = dim(ueRA2.3)[1], assay = "RNA")

#---Principal Component Analysis---#
ueRA2.3 = RunPCA(ueRA2.3, features = VariableFeatures(ueRA2.3))
# ElbowPlot(ueRA2.3, reduction = "pca", ndims = 50)

#---Unsupervised Re-Clustering---#
ueRA2.3 = FindNeighbors(ueRA2.3, dims = 1:40)
ueRA2.3 = FindClusters(ueRA2.3, resolution = 0.95)

#---Running UMAP---#
ueRA2.3 = RunUMAP(ueRA2.3, dims = 1:40)
# UMAPPlot(ueRA2.3, label = TRUE)

#---Rename Clusters Accordingly---#
c3.0 = WhichCells(ueRA2.3, idents = "0")
c3.1 = WhichCells(ueRA2.3, idents = "1")
c3.2 = WhichCells(ueRA2.3, idents = "2")

Idents(ueRA2.3, cells = c3.0) = "3.0"
Idents(ueRA2.3, cells = c3.1) = "3.1"
Idents(ueRA2.3, cells = c3.2) = "3.2"

########################################
###---Differential Gene Expression---###
########################################

ueRA2.3.marks = subset(FindAllMarkers(ueRA2.3, only.pos = FALSE, min.pct = 0.25), 
	subset = p_val_adj < 0.05)
pos.marks = subset(ueRA2.3.marks, subset = avg_log2FC > 0.0)
neg.marks = subset(ueRA2.3.marks, subset = avg_log2FC < 0.0)

write.xlsx(pos.marks %>% group_by(cluster) %>% top_n(-10, p_val_adj), 
	"reports/ueRA2_C3_subs.xlsx", sheetName = "Top10_Up", col.names = TRUE, row.names = FALSE)
write.xlsx(neg.marks %>% group_by(cluster) %>% top_n(-10, p_val_adj), )
	"reports/ueRA2_C3_subs.xlsx", sheetName = "Top10_Down", col.names = TRUE, row.names = FALSE, append = TRUE)

saveRDS(ueRA2.3, file = "data/ueRA2_C3_subs") # Save sub-cluster object

#########################################
###---Map Sub-Clusters to Base Data---###
#########################################

ueRA2.1 = ueRA2
Idents(ueRA2.1, cells = c3.0) = "3.0"
Idents(ueRA2.1, cells = c3.1) = "3.1"
Idents(ueRA2.1, cells = c3.2) = "3.2"
Idents(ueRA2.1) = factor(Idents(ueRA2.1), 
	levels = c("0", "1", "2", "3.0", "3.1", "4", "5")
