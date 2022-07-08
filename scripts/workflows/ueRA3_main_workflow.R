library(Seurat)
library(dplyr)
library(DoubletFinder)
library(ggplot2)
library(scPred)
library(xlsx)

########################
###---Read in Data---###
########################
set.seed(1001)
ueRA3.data = Read10X("data/gene_expression/ueRA3") # Read in CellRanger output files
ueRA3 = CreateSeuratObject(counts = ueRA3.data, project = "ueRA3", min.cells = 3, min.features = 200) # Create Seurat object

###########################
###---Quality Control---###
###########################

ueRA3[["percent.mt"]] = PercentageFeatureSet(ueRA2, pattern = "^MT-") # Add % mitochondrial genes expressed as metadata
ueRA3[["percent.ribo"]] = PercentageFeatureSet(ueRA2, pattern = "^RP[SL]") # Add % ribosomal genes expressed as metadata

# VlnPlot(ueRA3, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), ncol = 2) # Plot relevant quality control metrics

ueRA3 = subset(ueRA3, subset = nFeature_RNA > 750 & nFeature_RNA < 3000 & percent.mt < 11 & percent.ribo >5) # Filter cells

##################################################
###---Normalization, Scaling, Gene Filtering---###
##################################################

#-----Normalize Data-----#
ueRA3 = NormalizeData(ueRA3, normalization.method = "LogNormalize", scale.factor = 10000) # Log-normalize the data
ueRA3 = FindVariableFeatures(ueRA3, selection.method = "vst") # Identify highly variable genes

#-----Remove Immunoglobulin Genes from Matrix-----#
ueRA3 = ueRA3[!grepl("^IG[HKL][VDJ]|^IGH[A][1-2]|^IGH[G][1-4]|^IGHD$|^IGHM$|^IGHE$|^IGKC$|^IGLC[1-7]|^AC233755.1$", rownames(ueRA3)),] # See Sundell et al. (2022) for details. 

#---Scaling Data---#
ueRA3 = ScaleData(ueRA3, features = rownames(ueRA3))

###################################################
###---Dimensionality Reduction and Clustering---###
###################################################

#---Principal Component Analysis---#
ueRA3 = RunPCA(object = ueRA3, features = VariableFeatures(object = ueRA3))

#---Unsupervised Clustering---#
ueRA3 = FindNeighbors(object = ueRA3, dims = 1:25) # Construct SNN graph, using 25 PCs
ueRA3 = FindClusters(object = ueRA3, resolution = 0.8) # Identify clusters

#-----Running UMAP-----#
ueRA3 = RunUMAP(object = ueRA3, dims = 1:25) # Run UMAP
UMAPPlot(ueRA3, label = TRUE) # Plot clusters in UMAP space

############################################
###---Doublet Prediction and Filtering---###
############################################

annotations = ueRA3@meta.data$seurat_clusters # Define clusters
homotypic.prop = modelHomotypic(annotations) # Model homotypic proportions
nExp = round(ncol(ueRA3)* 0.031) # Expect 3.1% doublets for 4 000 cells
nExp_adj = round(nExp*(1-homotypic.prop)) # Adjust expectations
ueRA3 = doubletFinder_v3(ueRA3, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:25) # Perform doublet prediction
DF.name = colnames(ueRA3@meta.data)[grepl("DF.classification", colnames(ueRA3@meta.data))]

# DimPlot(ueRA3, group.by = DF.name) + ggtitle("Predicted Doublets") + theme(plot.title = element_text(hjust = 0.5)) # Plot predicted doublets
# VlnPlot(ueRA3, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1) +
	ggtitle("Predicted Doublets") + theme(plot.title = element_text(hjust = 0.5)) # Plot number of detected genes in cells

dim(ueRA3[,ueRA3@meta.data[,DF.name] == "Doublet"])[2] # Number of doublets detected
ueRA3 = ueRA3[,ueRA3@meta.data[,DF.name] == "Singlet"] # Exclude predicted doublets

################################
###---Cell Type Prediction---###
################################

reference = scPred::pbmc_1 # Load reference
reference = reference %>% 
	NormalizeData() %>%
	FindVaraibleFeatures() %>%
	ScaleData() %>%
	RunPCA(verbose = FALSE) %>%
	RunUMAP(dims = 1:25) # Perform all processing steps on the reference data
# DimPlot(reference, group.by = "cell_type", label = TRUE, repel = TRUE) + 
#	NoAxes() + 
#	ggtitle("Reference") + 
#	theme(plot.title = element_text(hjust=0.5)) # Plot predicted cell types in reference in UMAP space

transfer.anchors = FindTransferAnchors(reference = reference, query = ueRA3, dims = 1:25) # Identify transfer anchors
predictions = TransferData(anchorset = transfer.anchors, refdata = reference$cell_type, dims = 1:25) # Transfer predictions
ueRA3 = AddMetaData(ueRA3, metadata = predictions) # Add predictions as metadata

# DimPlot(ueRA3, group.by = "predicted.id", label = TRUE, repel = TRUE) +
#	ggtitle("Predicted Cell Types") + 
#	theme(plot.title = element_text(hjust = 0.5)) # Plot predicted 

###################################
###---Remove Irrelevant Cells---###
###################################

Idents(ueRA3) = ueRA3$predicted.id # Set predicted cell types as main identities
vdj_doublets = read.table("data/VDJ_data/ueRA3_VDJ-seq_doublets.txt", header = FALSE, row.names = NULL) # Read in list of VDJ-seq doublets
Idents(object = ueRA3, cells = vdj_doublets) = "VDJDubs" # Indicate VDJ-seq data doublets
ueRA3 = subset(ueRA3, idents = "B cell") # Only retain predicted B cells, and discard VDJ-seq doublets

#################################
###---Single-Cell Transform---###
#################################

ueRA3 = SCTransform(object = ueRA3, variable.features.n = dim(ueRA3)[1]) # Due to high homegenity in the data (all memory B cells), all genes are retained as variable features

#---Principal Component Analysis---#
ueRA3 = RunPCA(object = ueRA3, features = VariableFeatures(object = ueRA3), assay = "SCT")
# ElbowPlot(ueRA3, reduction = "pca", ndims = 50)

#---Unsupervised Clustering---#
ueRA3 = FindNeighbors(object = ueRA3, dims = 1:35)
ueRA3 = FindClusters(object = ueRA3, resolution = 0.35)

#---Running UMAP---#
ueRA3 = RunUMAP(object = ueRA3, dims = 1:35)
DimPlot(ueRA3, label = T, pt.size = 2, label.size = 7, 
	cols = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#F564E3", "#619CFF"))

##########################################
###---Differentially Expressed Genes---###
##########################################

ueRA3.marks = subset(FindAllMarkers(ueRA3, only.pos = FALSE, min.pct = 0.25),
	subset = p_val_adj < 0.05) # Identify DEGs and retain only significant DEGs 
pos.marks = subset(ueRA3.marks, subset = avg_log2FC > 0.0) # Filter out upregulated DEGs
neg.marks = subset(ueRA3.marks, subset = avg_log2FC < 0.0) # Filter out downregulated DEGs

write.xlsx(pos.marks, "data/ueRA3_DEGs.xlsx", sheetName = "Upregulated", col.names = TRUE, row.names = FALSE)
write.xlsx(neg.marks, "data/ueRA3_DEGs.xlsx", sheetName = "Downregulated", col.names = TRUE, row.names = FALSE, append = TRUE)

#---Save Seurat Object---#
saveRDS(ueRA3, file = "data/ueRA3_base")
