library((Seurat)
library(dplyr)
library(xlsx)

##########################
###---Sub-Clustering---###
##########################

uRA4.3 = subset(uRA4, subset = seurat_clusters == "3") # Isolate cells in cluster 3

#---Single-Cell Transform---#
uRA4.3 = SCTransform(uRA4.3, variable.features.n = dim(uRA4.3)[1], assay = "RNA")

#---Principal Component Analysis---#
uRA4.3 = RunPCA(uRA4.3, features = VariableFeatures(uRA4.3))
# ElbowPlot(uRA4.3, reduction = "pca", ndims = 50)

#---Unsupervised Re-Clustering---#
uRA4.3 = FindNeighbors(uRA4.3, dims = 1:40)
uRA4.3 = FindClusters(uRA4.3, resolution = 0.9)

#---Running UMAP---#
uRA4.3 = RunUMAP(uRA4.3, dims = 1:40)
# UMAPPlot(uRA4.3, label = TRUE)

#---Rename Clusters Accordingly---#
c3.0 = WhichCells(uRA4.3, idents = "0")
c3.1 = WhichCells(uRA4.3, idents = "1")
c3.2 = WhichCells(uRA4.3, idents = "2")

Idents(uRA4.3, cells = c3.0) = "3.0"
Idents(uRA4.3, cells = c3.1) = "3.1"
Idents(uRA4.3, cells = c3.2) = "3.2"

########################################
###---Differential Gene Expression---###
########################################

uRA4.3.marks = subset(FindAllMarkers(uRA4.3, only.pos = FALSE, min.pct = 0.25), 
	subset = p_val_adj < 0.05)
pos.marks = subset(uRA4.3.marks, subset = avg_log2FC > 0.0)
neg.marks = subset(uRA4.3.marks, subset = avg_log2FC < 0.0)

write.xlsx(pos.marks %>% group_by(cluster) %>% top_n(-10, p_val_adj), 
	"reports/uRA4_C3_subs.xlsx", sheetName = "Top10_Up", col.names = TRUE, row.names = FALSE)
write.xlsx(neg.marks %>% group_by(cluster) %>% top_n(-10, p_val_adj), )
	"reports/uRA4_C3_subs.xlsx", sheetName = "Top10_Down", col.names = TRUE, row.names = FALSE, append = TRUE)

saveRDS(uRA4.3, file = "data/uRA4_C3_subs") # Save sub-cluster object

#########################################
###---Map Sub-Clusters to Base Data---###
#########################################

uRA4.1 = uRA4
Idents(uRA4.1, cells = c3.0) = "3.0"
Idents(uRA4.1, cells = c3.1) = "3.1"
Idents(uRA4.1, cells = c3.2) = "3.2"
Idents(uRA4.1) = factor(Idents(uRA4.1), 
	levels = c("0", "1", "2", "3.0", "3.1", "3.2", "4")

###############################
###---Cluster 3.1 Vs. 3.0---###
###############################

marks = subset(FindMarkers(uRA4.1, ident.1 = "3.1", ident.2 = "3.0", min.pct = 0.25, only.pos = TRUE), 
	subset = p_val_adj <0.05)
write.xlsx(marks, file = "reports/uRA4_3.1_vs_3.0_DEGs.xlsx", sheetName = "DEGs", col.names = TRUE, row.names = FALSE)

saveRDS(uRA4.1, file = "data/uRA4_w_subs") # Save Seurat object
