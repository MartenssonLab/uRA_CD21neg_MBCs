library((Seurat)
library(dplyr)
library(xlsx)

##########################
###---Sub-Clustering---###
##########################

uRA2.3 = subset(uRA2, subset = seurat_clusters == "3") # Isolate cells in cluster 3

#---Single-Cell Transform---#
uRA2.3 = SCTransform(uRA2.3, variable.features.n = dim(uRA2.3)[1], assay = "RNA")

#---Principal Component Analysis---#
uRA2.3 = RunPCA(uRA2.3, features = VariableFeatures(uRA2.3))
# ElbowPlot(uRA2.3, reduction = "pca", ndims = 50)

#---Unsupervised Re-Clustering---#
uRA2.3 = FindNeighbors(uRA2.3, dims = 1:40)
uRA2.3 = FindClusters(uRA2.3, resolution = 0.95)

#---Running UMAP---#
uRA2.3 = RunUMAP(uRA2.3, dims = 1:40)
# UMAPPlot(uRA2.3, label = TRUE)

#---Rename Clusters Accordingly---#
c3.0 = WhichCells(uRA2.3, idents = "0")
c3.1 = WhichCells(uRA2.3, idents = "1")
c3.2 = WhichCells(uRA2.3, idents = "2")

Idents(uRA2.3, cells = c3.0) = "3.0"
Idents(uRA2.3, cells = c3.1) = "3.1"
Idents(uRA2.3, cells = c3.2) = "3.2"

########################################
###---Differential Gene Expression---###
########################################

uRA2.3.marks = subset(FindAllMarkers(uRA2.3, only.pos = FALSE, min.pct = 0.25), 
	subset = p_val_adj < 0.05)
pos.marks = subset(uRA2.3.marks, subset = avg_log2FC > 0.0)
neg.marks = subset(uRA2.3.marks, subset = avg_log2FC < 0.0)

write.xlsx(pos.marks %>% group_by(cluster) %>% top_n(-10, p_val_adj), 
	"reports/uRA2_C3_subs.xlsx", sheetName = "Top10_Up", col.names = TRUE, row.names = FALSE)
write.xlsx(neg.marks %>% group_by(cluster) %>% top_n(-10, p_val_adj), )
	"reports/uRA2_C3_subs.xlsx", sheetName = "Top10_Down", col.names = TRUE, row.names = FALSE, append = TRUE)

saveRDS(uRA2.3, file = "data/uRA2_C3_subs") # Save sub-cluster object

#########################################
###---Map Sub-Clusters to Base Data---###
#########################################

uRA2.1 = uRA2
Idents(uRA2.1, cells = c3.0) = "3.0"
Idents(uRA2.1, cells = c3.1) = "3.1"
Idents(uRA2.1, cells = c3.2) = "3.2"
Idents(uRA2.1) = factor(Idents(uRA2.1), 
	levels = c("0", "1", "2", "3.0", "3.1", "3.2","4", "5")
