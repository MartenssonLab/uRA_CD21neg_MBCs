library(Seurat)
library(dplyr)
library(DoubletFinder)
library(ggplot2)
library(scPred)

########################
###---Read in Data---###
########################
set.seed(1001)

#---Gene Expression Data---#
ueRA4.data = Read10X("data/gene_expression/ueRA4")
ueRA4 = CreateSeuratObject(counts = ueRA4.data, project = "ueRA4", min.cells = 3, min.features = 200)

#---VDJ-Seq Data---#
vdj4 = read.csv("data/VDJ_data/ueRA4/ueRA4_filtered_contig_annotations.csv, sep = ",") %>%
	filter(productive == "true", chain == "IGH", is_cell == "true") # Only retain productive heavy chains in confirmed cells


vdjdubs4 = vdj4[(duplicated(vdj4$barcode) | duplicated(vdj4$barcode, fromLast = TRUE)),] # Save duplicated barcodes
vdj4 = vdj4[!(duplicated(vdj4$barcode) | duplicated(vdj4$barcode, fromLast = TRUE)),] # Remove cells with reoccurring barcodes (duplicates)
vdj4 = vdj[!(vdj4$c_gene == ""),] # Remove cells with no annotated C gene

vdj4 = vdj4[,c("barcode", "c_gene", "reads", "umis")] # Trim the contents for relevance
names(vdj4)[names(vdj4) == "c_gene"] = "Isotype" # Rename for simplicity
rownames(vdj4) = vdj4[,1] # Set barcodes as rownames
ueRA4 = AddMetaData(ueRA4, metadata = vdj4) # Add IgH isotype data to SEurat object

###########################
###---Quality Control---###
###########################

ueRA4[["percent.mt"]] = PercentageFeatureSet(ueRA4, pattern = "^MT-") # Add percent mitochondrial genes expressed as metadata
ueRA4[["percent.ribo"]] = PercentageFeatureSet(ueRA4, pattern = "^RP[SL]") # Add percent ribosomal genes expressed as metadata

# VlnPlot(ueRA4, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), ncol = 2) # Plot relevant quality control metrics

ueRA4 = subset(ueRA4, subset = nFeature_RNA > 500 & nFeature_RNA < 3000 & percent.mt < 12 & percent.ribo > 5) # Filter cells

##################################################
###---Normalization, Scaling, Gene Filtering---###
##################################################

#---Normalize Data---#
ueRA4 = NormalizeData(ueRA4, normalization.method = "LogNormalize", scale.factor = 10000) # Log-normalize the data
ueRA4 = FindVariableFeatures(ueRA4, selection.method = "vst") # Identify highly variable genes

#---Remove Immunoglobulin Genes from Gene Expression Matrix---#
ueRA4 = ueRA4[!grepl("^IG[HKL][VDJ]|^IGH[A][1-2]|^IGH[G][1-4]|^IGHD$|^IGHM$|^IGHE$|^IGKC$|^IGLC[1-7]|^AC233755.1$", rownames(ueRA4)),] # See Sundell et al. (2022) for details.

#---Scaling Data---#
ueRA4 = ScaleData(ueRA4, features = rownames(ueRA4))

###################################################
###---Dimensionality Reduction and Clustering---###
###################################################

#---Principal Component Analysis---#
ueRA4 = RunPCA(ueRA4, features = VariableFeatures(ueRA4)) # Perform PCA

#---Unsupervised Clustering---#
ueRA4 = FindNeighbors(ueRA4, dims = 1:25) # Construct SNN graph, using 25 PCs
ueRA = FindClusters(ueRA4, resolution = 0.8) # Identify clusters

#---Running UMAP---#
ueRA4 = RunUMAP(ueRA4, dims = 1:25)
UMAPPlot(ueRA4, label = TRUE)

############################################
###---Doublet Prediction and Filtering---###
############################################

annotations = ueRA4@meta.data$seurat_clusters # Define clusters
homotypic.prop = modelHomotypic(annotations) # Model homotypic proportions
nExp = round(ncol(ueRA4)* 0.031) # Expect 3.1% doublets based on 4 000 cells
nExp_adj = round(nExp*(1-homotypic.prop)) # Adjust expectations
ueRA4 = doubletFinder_v3(ueRA4, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:25) # Perform doublet prediction
DF.name = colnames(ueRA4@meta.data)[grepl("DF.classification", colnames(ueRA4@meta.data))]

# DimPlot(ueRA4, group.by = DF.name) + ggtitle("Predicted Doublets") + theme(plot.title = element_text(hjust = 0.5)) # Plot predicted doublets
# VlnPlot(ueRA4, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1) +
	ggtitle("Predicted Doublets") + theme(plot.title = element_text(hjust = 0.5)) # Plot number of detected genes in cells

dim(ueRA4[,ueRA4@meta.data[,DF.name] == "Doublet"])[2] # Number of doublets detected
ueRA4 = ueRA4[,ueRA4@meta.data[,DF.name] == "Singlet"] # Exclude predicted doublets

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

transfer.anchors = FindTransferAnchors(reference = reference, query = ueRA4, dims = 1:25) # Identify transfer anchors
predictions = TransferData(anchorset = transfer.anchors, refdata = reference$cell_type, dims = 1:25) # Transfer predictions
ueRA4 = AddMetaData(ueRA4, metadata = predictions) # Add predictions as metadata
# DimPlot(ueRA4, group.by = "predicted.id", label = TRUE, repel = TRUE) +
#	ggtitle("Predicted Cell Types") + 
#	theme(plot.title = element_text(hjust = 0.5)) # Plot predicted 

###################################
###---Remove Irrelevant Cells---###
###################################

Idents(ueRA4) = ueRA4$predicted.id
Idents(ueRA4, cells = vdjdubs4$barcode) = "VDJdubs" # Indicate VDJ-seq data doublets
ueRA4 = subset(ueRA4, idents = "B cell") # Only retain predicted B cells, and discard VDJ-seq doublets

#################################
###---Single-Cell Transform---###
#################################
ueRA4 = SCTransform(ueRA4, variable.features.n = dim(ueRA4)[1]) # Due to high homegenity in the data (all memory B cells), all genes are retained as variable features

#---Principal Component Analysis---#
ueRA4 = RunPCA(ueRA4, features = VariableFeatures(ueRA4), assay = "SCT") # Perform PCA
# ElbowPlot(ueRA4.1, reduction = "pca", ndims = 50)

#---Unsupervised Clustering---#
ueRA4 = FindNeighbors(ueRA4, dims = 1:35)
ueRA4 = FindClusters(ueRA4, resolution = 0.35)

#---Running UMAP---#
ueRA4 = RunUMAP(ueRA4, dims = 1:35)
DimPlot(ueRA4, label = TRUE, pt.size = 2, label.size = 7)

##########################################
###---Differentially Expressed Genes---###
##########################################

ueRA4.marks = subset(FindAllMarkers(ueRA4, only.pos = FALSE, min.pct = 0.25),
	subset = p_val_adj < 0.05) # Identify DEGs and retain only significant DEGs 
pos.marks = subset(ueRA4.marks, subset = avg_log2FC > 0.0) # Filter out upregulated DEGs
neg.marks = subset(ueRA4.marks, subset = avg_log2FC < 0.0) # Filter out downregulated DEGs

write.xlsx(pos.marks, "reports/ueRA4_DEGs.xlsx", sheetName = "Upregulated", col.names = TRUE, row.names = FALSE)
write.xlsx(neg.marks, "reports/ueRA4_DEGs.xlsx", sheetName = "Downregulated", col.names = TRUE, row.names = FALSE, append = TRUE)

#---Save Seurat Object---#
saveRDS(ueRA4, file = "data/ueRA4_base")
