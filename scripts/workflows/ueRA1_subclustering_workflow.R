library(Seurat)
library(dplyr)
library(xlsx)

##########################
###---Sub-Clustering---###
##########################

ueRA1.3 = subset(ueRA1, subset = seurat_clusters == "3") # Isolate cells in cluster 3

#---Single-Cell Transform---#
ueRA1.3 = SCTransform(ueRA1.3, variable.features.n = dim(ueRA1.3)[1], assay = "RNA")

#---Principal Component Analysis---#
ueRA1.3 = RunPCA(ueRA1.3, features = VariableFeatures(ueRA1.3))
# ElbowPlot(ueRA1.3, reduction = "pca", ndims = 50)

#---Unsupervised Re-Clustering---#
ueRA1.3 = FindNeighbors(ueRA1.3, dims = 1:40)
ueRA1.3 = FindClusters(ueRA1.3, resolution = 0.95)

#---Running UMAP---#
ueRA1.3 = RunUMAP(ueRA1.3, dims = 1:40)
# UMAPPlot(ueRA1.3, label = TRUE)

#---Rename Clusters Accordingly---#
c3.0 = WhichCells(ueRA1.3, idents = "0")
c3.1 = WhichCells(ueRA1.3, idents = "1")

Idents(ueRA1.3, cells = c3.0) = "3.0"
Idents(ueRA1.3, cells = c3.1) = "3.1"

########################################
###---Differential Gene Expression---###
########################################

ueRA1.3.marks = subset(FindAllMarkers(ueRA1.3, only.pos = FALSE, min.pct = 0.25), 
	subset = p_val_adj < 0.05)
pos.marks = subset(ueRA1.3.marks, subset = avg_log2FC > 0.0)
neg.marks = subset(ueRA1.3.marks, subset = avg_log2FC < 0.0)

write.xlsx(pos.marks %>% group_by(cluster) %>% top_n(-10, p_val_adj), 
	"reports/ueRA1_C3_subs.xlsx", sheetName = "Top10_Up", col.names = TRUE, row.names = FALSE)
write.xlsx(neg.marks %>% group_by(cluster) %>% top_n(-10, p_val_adj), )
	"reports/ueRA1_C3_subs.xlsx", sheetName = "Top10_Down", col.names = TRUE, row.names = FALSE, append = TRUE)

saveRDS(ueRA1.3, file = "data/ueRA1_C3_subs") # Save sub-cluster object

#########################################
###---Map Sub-Clusters to Base Data---###
#########################################

ueRA1.1 = ueRA1
Idents(ueRA1.1, cells = c3.0) = "3.0"
Idents(ueRA1.1, cells = c3.1) = "3.1"
Idents(ueRA1.1) = factor(Idents(ueRA1.1), 
	levels = c("0", "1", "2", "3.0", "3.1", "4", "5")

saveRDS(ueRA1.1, file = "data/ueRA1_w_subs") # Save Seurat object
