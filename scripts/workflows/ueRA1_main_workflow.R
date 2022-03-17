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

#---Gene Expression Data---#
ueRA1.data = Read10X("data/gene_expression/ueRA1") # Read in CellRanger output files
ueRA1 = CreateSeuratObject(counts = ueRA1.data, project = "ueRA1", min.cells = 3, min.features = 200) # Create Seurat object

#---VDJ-Seq-Data---#
vdj1 = read.csv("data/VDJ_data/ueRA1_filtered_contig_annotations.csv, sep = ",") %>%
	filter(productive == "true", chain == "IGH", is_cell == "true") # Only retain productive heavy chains in confirmed cells

vdjdubs1 = vdj1[(duplicated(vdj1$barcode) | duplicated(vdj1$barcode, fromLast = TRUE)),] # Save duplicated barcodes
vdj1 = vdj1[!(duplicated(vdj1$barcode) | duplicated(vdj1$barcode, fromLast = TRUE)),] # Remove cells with reoccurring barcodes (duplicates)
vdj1 = vdj1[!(vdj1$c_gene == ""),] # Remove cells with no annotated C gene

vdj1 = vdj1[,c("barcode", "c_gene", "reads", "umis")] # Trim the contents for relevance
names(vdj1)[names(vdj1) == "c_gene"] = "Isotype" # Rename for simplicity
rownames(vdj1) = vdj1[,1] # Set barcodes as rownames
ueRA1 = AddMetaData(ueRA1, metadata = vdj1) # Add IgH isotype data to Seurat object

###########################
###---Quality Control---###
###########################

ueRA1[["percent.mt"]] = PercentageFeatureSet(ueRA1, pattern = "^MT-") # Add % mitochondrial genes expressed as metadata
ueRA1[["percent.ribo"]] = PercentageFeatureSet(ueRA1, pattern = "^RP[SL]") # Add % ribosomal genes expressed as metadata

# VlnPlot(ueRA1, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), ncol = 2) # Plot relevant quality control metrics

ueRA1 = subset(ueRA1, subset = nFeature_RNA > 700 & nFeature_RNA < 3000 & percent.mt < 11 & percent.ribo >5) # Filter cells

##################################################
###---Normalization, Scaling, Gene Filtering---###
##################################################

#-----Normalize Data-----#
ueRA1 = NormalizeData(ueRA1, normalization.method = "LogNormalize", scale.factor = 10000) # Log-normalize the data
ueRA1 = FindVariableFeatures(ueRA1, selection.method = "vst") # Identify highly variable genes

#-----Remove Immunoglobulin Genes from Matrix-----#
ueRA1 = ueRA1[!grepl("^IG[HKL][VDJ]|^IGH[A][1-2]|^IGH[G][1-4]|^IGHD$|^IGHM$|^IGHE$|^IGKC$|^IGLC[1-7]|^AC233755.1$", rownames(ueRA1)),] # See Sundell et al. (2022) for details. 

#---Scaling Data---#
ueRA1 = ScaleData(ueRA1, features = rownames(ueRA1))

###################################################
###---Dimensionality Reduction and Clustering---###
###################################################

#---Principal Component Analysis---#
ueRA1 = RunPCA(object = ueRA1, features = VariableFeatures(object = ueRA1))

#---Unsupervised Clustering---#
ueRA1 = FindNeighbors(object = ueRA1, dims = 1:25) # Construct SNN graph, using 25 PCs
ueRA1 = FindClusters(object = ueRA1, resolution = 0.8) # Identify clusters

#-----Running UMAP-----#
ueRA1 = RunUMAP(object = ueRA1, dims = 1:25) # Run UMAP
UMAPPlot(ueRA1, label = TRUE) # Plot clusters in UMAP space

############################################
###---Doublet Prediction and Filtering---###
############################################

annotations = ueRA1@meta.data$seurat_clusters # Define clusters
homotypic.prop = modelHomotypic(annotations) # Model homotypic proportions
nExp = round(ncol(ueRA1)* 0.031) # Expect 3.1% doublets for 4 000 cells
nExp_adj = round(nExp*(1-homotypic.prop)) # Adjust expectations
ueRA1 = doubletFinder_v3(ueRA1, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:25) # Perform doublet prediction
DF.name = colnames(ueRA1@meta.data)[grepl("DF.classification", colnames(ueRA1@meta.data))]

# DimPlot(ueRA1, group.by = DF.name) + ggtitle("Predicted Doublets") + theme(plot.title = element_text(hjust = 0.5)) # Plot predicted doublets
# VlnPlot(ueRA1, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1) +
	ggtitle("Predicted Doublets") + theme(plot.title = element_text(hjust = 0.5)) # Plot number of detected genes in cells

dim(ueRA1[,ueRA1@meta.data[,DF.name] == "Doublet"])[2] # Number of doublets detected
ueRA1 = ueRA1[,ueRA1@meta.data[,DF.name] == "Singlet"] # Exclude predicted doublets

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

transfer.anchors = FindTransferAnchors(reference = reference, query = ueRA1, dims = 1:25) # Identify transfer anchors
predictions = TransferData(anchorset = transfer.anchors, refdata = reference$cell_type, dims = 1:25) # Transfer predictions
ueRA1 = AddMetaData(ueRA1, metadata = predictions) # Add predictions as metadata

# DimPlot(ueRA1, group.by = "predicted.id", label = TRUE, repel = TRUE) +
#	ggtitle("Predicted Cell Types") + 
#	theme(plot.title = element_text(hjust = 0.5)) # Plot predicted 

###################################
###---Remove Irrelevant Cells---###
###################################

Idents(ueRA1) = ueRA1$predicted.id # Set predicted cell types as main identities
Idents(object = ueRA1, cells = vdjdubs1$barcode) = "VDJDubs" # Indicate VDJ-seq data doublets
ueRA1 = subset(ueRA1, idents = "B cell") # Only retain predicted B cells, and discard VDJ-seq doublets

#################################
###---Single-Cell Transform---###
#################################

ueRA1 = SCTransform(object = ueRA1, variable.features.n = dim(ueRA1)[1]) # Due to high homegenity in the data (all memory B cells), all genes are retained as variable features

#---Principal Component Analysis---#
ueRA1 = RunPCA(object = ueRA1, features = VariableFeatures(object = ueRA1), assay = "SCT")
# ElbowPlot(ueRA1, reduction = "pca", ndims = 50)

#---Unsupervised Clustering---#
ueRA1 = FindNeighbors(object = ueRA1, dims = 1:35)
ueRA1 = FindClusters(object = ueRA1, resolution = 0.23)

#---Running UMAP---#
ueRA1 = RunUMAP(object = ueRA1, dims = 1:35)
DimPlot(ueRA1, label = T, pt.size = 2, label.size = 7, 
	cols = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#F564E3", "#619CFF"))

##########################################
###---Differentially Expressed Genes---###
##########################################

ueRA1.marks = subset(FindAllMarkers(ueRA1, only.pos = FALSE, min.pct = 0.25),
	subset = p_val_adj < 0.05) # Identify DEGs and retain only significant DEGs 
pos.marks = subset(ueRA1.marks, subset = avg_log2FC > 0.0) # Filter out upregulated DEGs
neg.marks = subset(ueRA1.marks, subset = avg_log2FC < 0.0) # Filter out downregulated DEGs

write.xlsx(pos.marks, "reports/ueRA1_DEGs.xlsx", sheetName = "Upregulated", col.names = TRUE, row.names = FALSE)
write.xlsx(neg.marks, "reports/ueRA1_DEGs.xlsx", sheetName = "Downregulated", col.names = TRUE, row.names = FALSE, append = TRUE)

#---Save Seurat Object---#
saveRDS(ueRA1, file = "data/ueRA1_base")library(Seurat)
library(dplyr)
library(DoubletFinder)
library(ggplot2)
library(scPred)
library(xlsx)

########################
###---Read in Data---###
########################
set.seed(1001)

#---Gene Expression Data---#
ueRA1.data = Read10X("Rdata/gene_expression/ueRA1") # Read in CellRanger output files
ueRA1 = CreateSeuratObject(counts = ueRA1.data, project = "ueRA1", min.cells = 3, min.features = 200) # Create Seurat object

#---VDJ-Seq-Data---#
vdj1 = read.csv("data/VDJ_data/ueRA1_filtered_contig_annotations.csv, sep = ",") %>%
	filter(productive == "true", chain == "IGH", is_cell == "true") # Only retain productive heavy chains in confirmed cells

vdjdubs1 = vdj1[(duplicated(vdj1$barcode) | duplicated(vdj1$barcode, fromLast = TRUE)),] # Save duplicated barcodes
vdj1 = vdj1[!(duplicated(vdj1$barcode) | duplicated(vdj1$barcode, fromLast = TRUE)),] # Remove cells with reoccurring barcodes (duplicates)
vdj1 = vdj1[!(vdj4$c_gene == ""),] # Remove cells with no annotated C gene

vdj1 = vdj1[,c("barcode", "c_gene", "reads", "umis")] # Trim the contents for relevance
names(vdj1)[names(vdj1) == "c_gene"] = "Isotype" # Rename for simplicity
rownames(vdj1) = vdj1[,1] # Set barcodes as rownames
ueRA1 = AddMetaData(ueRA1, metadata = vdj1) # Add IgH isotype data to Seurat object

###########################
###---Quality Control---###
###########################

ueRA1[["percent.mt"]] = PercentageFeatureSet(ueRA1, pattern = "^MT-") # Add % mitochondrial genes expressed as metadata
ueRA1[["percent.ribo"]] = PercentageFeatureSet(ueRA1, pattern = "^RP[SL]") # Add % ribosomal genes expressed as metadata

# VlnPlot(ueRA1, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), ncol = 2) # Plot relevant quality control metrics

ueRA1 = subset(ueRA1, subset = nFeature_RNA > 700 & nFeature_RNA < 3000 & percent.mt < 11 & percent.ribo >5) # Filter cells

##################################################
###---Normalization, Scaling, Gene Filtering---###
##################################################

#-----Normalize Data-----#
ueRA1 = NormalizeData(ueRA1, normalization.method = "LogNormalize", scale.factor = 10000) # Log-normalize the data
ueRA1 = FindVariableFeatures(ueRA1, selection.method = "vst") # Identify highly variable genes

#-----Remove Immunoglobulin Genes from Matrix-----#
ueRA1 = ueRA1[!grepl("^IG[HKL][VDJ]|^IGH[A][1-2]|^IGH[G][1-4]|^IGHD$|^IGHM$|^IGHE$|^IGKC$|^IGLC[1-7]|^AC233755.1$", rownames(ueRA1)),] # See Sundell et al. (2022) for details. 

#---Scaling Data---#
ueRA1 = ScaleData(ueRA1, features = rownames(ueRA4))

###################################################
###---Dimensionality Reduction and Clustering---###
###################################################

#---Principal Component Analysis---#
ueRA1 = RunPCA(object = ueRA1, features = VariableFeatures(object = ueRA1))

#---Unsupervised Clustering---#
ueRA1 = FindNeighbors(object = ueRA1, dims = 1:25) # Construct SNN graph, using 25 PCs
ueRA1 = FindClusters(object = ueRA1, resolution = 0.8) # Identify clusters

#-----Running UMAP-----#
ueRA1 = RunUMAP(object = ueRA1, dims = 1:25) # Run UMAP
UMAPPlot(ueRA1, label = TRUE) # Plot clusters in UMAP space

############################################
###---Doublet Prediction and Filtering---###
############################################

annotations = ueRA1@meta.data$seurat_clusters # Define clusters
homotypic.prop = modelHomotypic(annotations) # Model homotypic proportions
nExp = round(ncol(ueRA1)* 0.031) # Expect 3.1% doublets for 4 000 cells
nExp_adj = round(nExp*(1-homotypic.prop)) # Adjust expectations
ra1 = doubletFinder_v3(ueRA1, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:25) # Perform doublet prediction
DF.name = colnames(ueRA1@meta.data)[grepl("DF.classification", colnames(ueRA1@meta.data))]

# DimPlot(ueRA1, group.by = DF.name) + ggtitle("Predicted Doublets") + theme(plot.title = element_text(hjust = 0.5)) # Plot predicted doublets
# VlnPlot(ueRA1, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1) +
	ggtitle("Predicted Doublets") + theme(plot.title = element_text(hjust = 0.5)) # Plot number of detected genes in cells

dim(ueRA1[,ueRA1@meta.data[,DF.name] == "Doublet"])[2] # Number of doublets detected
ueRA1 = ueRA1[,ueRA1@meta.data[,DF.name] == "Singlet"] # Exclude predicted doublets

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

transfer.anchors = FindTransferAnchors(reference = reference, query = ueRA1, dims = 1:25) # Identify transfer anchors
predictions = TransferData(anchorset = transfer.anchors, refdata = reference$cell_type, dims = 1:25) # Transfer predictions
ueRA1 = AddMetaData(ueRA1, metadata = predictions) # Add predictions as metadata

# DimPlot(ueRA1, group.by = "predicted.id", label = TRUE, repel = TRUE) +
#	ggtitle("Predicted Cell Types") + 
#	theme(plot.title = element_text(hjust = 0.5)) # Plot predicted 

###################################
###---Remove Irrelevant Cells---###
###################################

Idents(ueRA1) = ueRA1$predicted.id # Set predicted cell types as main identities
Idents(object = ueRA1, cells = vdjdubs1$barcode) = "VDJDubs" # Indicate VDJ-seq data doublets
ueRA1 = subset(ueRA1, idents = "B cell") # Only retain predicted B cells, and discard VDJ-seq doublets

#################################
###---Single-Cell Transform---###
#################################

ueRA1 = SCTransform(object = ueRA1, variable.features.n = dim(ueRA1)[1]) # Due to high homegenity in the data (all memory B cells), all genes are retained as variable features

#---Principal Component Analysis---#
ueRA1 = RunPCA(object = ueRA1, features = VariableFeatures(object = ueRA1), assay = "SCT")
# ElbowPlotueRA1, reduction = "pca", ndims = 50)

#---Unsupervised Clustering---#
ueRA1 = FindNeighbors(object = ueRA1, dims = 1:35)
ueRA1 = FindClusters(object = ueRA1, resolution = 0.23)

#---Running UMAP---#
ueRA1 = RunUMAP(object = ueRA1, dims = 1:35)
DimPlot(ueRA1, label = T, pt.size = 2, label.size = 7, 
	cols = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#F564E3", "#619CFF"))

##########################################
###---Differentially Expressed Genes---###
##########################################

ueRA1.marks = subset(FindAllMarkers(ueRA1, only.pos = FALSE, min.pct = 0.25),
	subset = p_val_adj < 0.05) # Identify DEGs and retain only significant DEGs 
pos.marks = subset(ueRA1.marks, subset = avg_log2FC > 0.0) # Filter out upregulated DEGs
neg.marks = subset(ueRA1.marks, subset = avg_log2FC < 0.0) # Filter out downregulated DEGs

write.xlsx(pos.marks, "reports/ueRA1_DEGs.xlsx", sheetName = "Upregulated", col.names = TRUE, row.names = FALSE)
write.xlsx(neg.marks, "reports/ueRA1_DEGs.xlsx", sheetName = "Downregulated", col.names = TRUE, row.names = FALSE, append = TRUE)

#---Save Seurat Object---#
saveRDS(ueRA1, file = "data/ueRA1_base")
