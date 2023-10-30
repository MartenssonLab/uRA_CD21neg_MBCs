library(Seurat)
library(ggplot2)
library(monocle3)
library(SeuratWrappers)
library(gridExtra)

# Read in and process data
uRA4.1 = readRDS("data/uRA4_w_subs.rds")

## Extract cluster information
uRA4_traj = as.cell_data_set(uRA4.1) # Convert to cell dataset
uRA4_traj = cluster_cells(uRA4_traj) # Create initial clustering
plot_cells(uRA4_traj, show_trajectory_graph = F, cell_size = 0.75, 
           color_cells_by = "cluster")

colData(uRA4_full)
fData(uRA4_traj) # Access row data table in monocle object
rownames(fData(uRA4_traj))[1:10] # Genes
fData(uRA4_traj)$gene_short_name = rownames(fData(uRA4_traj))

# Extract cluster information from Seurat object
recreate.partition = c(rep(1, length(uRA4_traj@colData@rownames)))
names(recreate.partition) = uRA4_traj@colData$rownames
recreate.partition = as.factor(recreate.partition)
uRA4_traj@clusters$UMAP$UMAP$clusters = list_cluster

## Assign cluster information
list_cluster = uRA4.1@active.ident
uRA4_traj@reductions$umap@cell.embeddings

## Assign UMAP coordinate cell embeddings
uRA4_traj@int_colData@listData$reducedDims$UMAP = 
  uRA4.1@reductions$umap@cell.embeddings

# Learn trajectory graph
uRA4_traj = learn_graph(uRA4_traj, use_partition = T)

plot_cells(uRA4_traj, 
           color_cells_by = "cluster", 
           label_groups_by_cluster = F, 
           group_label_size = 7, 
           cell_size = 2, 
           trajectory_graph_segment_size = 1.5, 
           trajectory_graph_color = "black") # Plot trajectory

# Pseudo-temporal ordering of cells
uRA4_traj = order_cells(uRA4_traj) # Select node in cluster 0 (unswitched MBCs)

plot_cells(uRA4_traj, 
           color_cells_by = "pseudotime", 
           label_groups_by_cluster = T, 
           group_label_size = 7, 
           cell_size = 2, 
           trajectory_graph_segment_size = 1.5, 
           trajectory_graph_color = "black") # Plot trajectory with cells colored by pseudotime
