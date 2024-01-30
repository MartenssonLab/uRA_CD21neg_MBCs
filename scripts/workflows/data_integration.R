
# Packages used
```{r}
library(Seurat)
library(ggplot2)
library(Nebulosa)
library(gridExtra)
library(xlsx)
library(ComplexHeatmap)
library(scales)
```

# Integrate uRA datasets 1 - 4
```{r}
# Read in processed data
uRA1.1 = readRDS("/data/uRA1_w_subs.rds")
uRA2.1 = readRDS("/data/uRA2_w_subs.rds")
uRA3.1 = readRDS("/data/uRA3_w_subs.rds")
uRA4.1 = readRDS("/data/uRA4_w_subs.rds")

uRA_files = list(uRA1.1, uRA2.1, uRA3.1, uRA4.1)

# Data integration
## Identify 10 000 integration features
int_features = SelectIntegrationFeatures(uRA_files, 
                                         nfeatures = 10000)

## Prep SCT-processed data for integration
uRA_int = PrepSCTIntegration(uRA_files, 
                             anchor.features = int_features)

## Identify integration anchors
anchors = FindIntegrationAnchors(uRA_files, 
                                 dims = 1:35, 
                                 reduction = "pca", 
                                 anchor.features = int_features)

## Integrate the datasets
uRA_int = IntegrateData(anchors, 
                        dims = 1:35, 
                        new.assay.name = "CCA")

## Scale integrated data
uRA_int = ScaleData(uRA_int)

# Cluster integrated data
uRA_int = uRA_int %>%
  RunUMAP(dims = 1:35) %>%
  FindNeighbors(dims = 1:35) %>%
  FindClusters(resolution = 0.6)
```

# Examine Expression of Marker Genes from Flow Cytometry
```{r}
# Main markers
plot_density(uRA_int, 
             features = c("CR2", "TBX21", "ITGAX", "CD27"))

VlnPlot(uRA_int, 
        features = "CD27", 
        assay = "RNA")

# Upregulated flow markers
plot_density(uRA_int, 
             features = c("FCRL5", "FCRL3", "FAS", "CXCR3", "CD19", "MS4A1"))

# Downregulated flow markers
plot_density(uRA_int, 
             features = c("SELL", "CD24", "CD38", "CXCR5", "CXCR4"))
```

# Examine Ig Isotype Distribution
```{r}
# IgM
DimPlot(uRA_int, 
        cells.highlight = WhichCells(uRA_int, 
                                     expression = isotype == "IgM"), 
        cols.highlight = "cornflowerblue", 
        pt.size = 1.5, 
        sizes.highlight = 1.5, 
        order = T) + 
  ggtitle("IgM") + 
  theme(plot.title = element_text(hjust = 0.5, 
                                  face = 1, 
                                  size = 18)) + 
  NoLegend()

# IgG
DimPlot(uRA_int, 
        cells.highlight = WhichCells(uRA_int, 
                                     expression = isotype == "IgG"), 
        cols.highlight = "brown4", 
        pt.size = 1.5, 
        sizes.highlight = 1.5, 
        order = T) + 
  ggtitle("IgG") + 
  theme(plot.title = element_text(hjust = 0.5, 
                                  face = 1, 
                                  size = 18)) + 
  NoLegend()

# IgA
DimPlot(uRA_int, 
        cells.highlight = WhichCells(uRA_int, 
                                     expression = isotype == "IgA"), 
        cols.highlight = "tan1", 
        pt.size = 1.5, 
        sizes.highlight = 1.5, 
        order = T) + 
  ggtitle("IgA") + 
  theme(plot.title = element_text(hjust = 0.5, 
                                  face = 1, 
                                  size = 18)) + 
  NoLegend()
```

# Subcluster Cells In the Cluster Equivalent to uRA4 Cluster 3
```{r}
int_c6 = subset(uRA_int, idents = "6")

int_c6 = int_c6 %>% 
  ScaleData() %>%
  FindVariableFeatures(nfeatures = 10000) %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:40) %>%
  FindClusters(resolution = 0.5) %>%
  RunUMAP(dims = 1:40)


int_c6.1 = int_c6 %>% 
  ScaleData() %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 10000) %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:40) %>%
  FindClusters(resolution = 0.5) %>%
  RunUMAP(dims = 1:40)

## Map subclusters back to integrated main data
Idents(uRA_int, cells = WhichCells(int_c6, 
                                   idents = "0")) = "6.0"
Idents(uRA_int, cells = WhichCells(int_c6, 
                                   idents = "1")) = "6.1"
Idents(uRA_int, cells = WhichCells(int_c6, 
                                   idents = "2")) = "6.2"
Idents(uRA_int, cells = WhichCells(int_c6, 
                                   idents = "3")) = "6.3"

uRA_int$seurat_subclusters = factor(Idents(uRA_int), levels = c("0", "1", "2", "3", "4", "5", 
                                                                "6.0", "6.1", "6.2", "6.3", 
                                                                "7", "8", "9", "10", "11", "12"))

Idents(uRA_int) = uRA_int$seurat_subclusters

## Rename subclusters to reflect uRA4 subclusters
#Idents(uRA_int, cells = WhichCells(uRA_int, idents = "6.0")) = "3.0a"
#Idents(uRA_int, cells = WhichCells(uRA_int, idents = "6.1")) = "3.0b"
#Idents(uRA_int, cells = WhichCells(uRA_int, idents = "6.2")) = "3.1"
#Idents(uRA_int, cells = WhichCells(uRA_int, idents = "6.3")) = "3.2"

# Save Seurat object
saveRDS(uRA_int, "data/uRA_integrated.rds")
```

# Examine Isotype Proportions in uRA4 3.1- and 3.2-Equivalent Subclusters
```{r}
# Calculate
## NOTE: Due to poor annotation of IgD, and negligible expression of IgE, only IgM, IgG, and IgA expression is considered.
tab = table(Idents(uRA_int), uRA_int$isotype)[, c(1, 3, 4)]
for(i in seq_len(nrow(tab))){
  percentages = c(
    IgM = sum(tab[i, 1]) / sum(tab[i,]) * 100,
    IgG = sum(tab[i, 2]) / sum(tab[i,]) * 100,
    IgA = sum(tab[i, 3]) / sum(tab[i,]) * 100
  )
  results_df = rbind(results_df, cbind(Cluster = rownames(tab)[i], round(t(as.data.frame(percentages)), 1)))
}

# Plot isotype proportions in uRA4 3.1- and 3.2-equivalent subclusters
cluster_df = percent_df[percent_df$Cluster %in% c("6.2", "6.3"),]
for(i in seq_len(nrow(cluster_df))){
  df_for_plot = data.frame(
    Isotype = c("IgM", "IgG", "IgA"), 
    levels = c("IgM", "IgG", "IgA"),
    Value = as.integer(c(cluster_df$IgM[i], 
                         cluster_df$IgG[i], 
                         cluster_df$IgA[i])), 
    Labels = paste(c(cluster_df$IgM[i], 
                     cluster_df$IgG[i], 
                     cluster_df$IgA[i]), 
                   "%", 
                   sep = "")
  )
  
  # Calculate the cumulative percentages for label placement
  df_for_plot$Position = cumsum(as.integer(df_for_plot$Value)) - 0.5 * as.integer(df_for_plot$Value)
  
  p = ggplot(df_for_plot, aes(x = "", 
                              y = Value, 
                              fill = Isotype)) +
    geom_bar(width = 1, stat = "identity", colour = "black") +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = c("IgM" = "cornflowerblue", 
                                 "IgG" = "brown4", 
                                 "IgA" = "tan1"), 
                      breaks = c("IgM", "IgG", "IgA")) +
    geom_text(aes(y = Position, 
                  label = Labels), 
              color = "white", 
              size = 5) +
    theme_void() +
    labs(title = paste("Cluster", cluster_df$Cluster[i])) +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5))
  
  plot_list[[i]] = p
}

gridExtra::grid.arrange(grobs = plot_list, 
                        ncol = 2)
```

# Examine Expresison Patterns for DN Tbet+ CD11c+ and CD27+ Tbet+ CD11c+ Signatures
```{r}
# Read in gene lists
DN_sign = read.xlsx("/reports/uRA4_3.1_vs_3.0_DEGs.xlsx", sheetIndex = 1)
CD27_sign = read.xlsx("/reports/uRA4_3.2_vs_3.0_DEGs.xlsx", sheetIndex = 1)

# Calculate and extract scaled average expression of the markers genes for each cluster
avg_DN = t(scale(t(as.matrix(AverageExpression(uRA_int, features = DN_sign$gene)$RNA))))
avg_CD27 = t(scale(t(as.matrix(AverageExpression(uRA_int, features = CD27_sign$gene)$RNA))))

# Define annotation identities and colors
annotations = data.frame(levels(Idents(uRA_int)))
names(annotations) = "Cluster"
rownames(annotations) = levels(Idents(uRA_int))

anno_cols = scales::hue_pal()(length(levels(Idents(uRA_int))))

# Plot average expression heatmaps
## DN Tbet+ CD11c+ signature
ComplexHeatmap::pheatmap(avg_DN, 
                         color = colorRampPalette(c("#FF00FF", "black", "#FFFF00"))(100), 
                         cluster_rows = F, 
                         cluster_cols = F, 
                         fontsize_row = 15, 
                         fontsize_col = 15,
                         border_color = "black",
                         show_rownames = F, 
                         show_colnames = T, 
                         annotation_col = annotations, 
                         annotation_colors = anno_cols, 
                         annotation_legend = F, 
                         main = "DN signature", 
                         column_names_side = "top", 
                         angle_col = "45")

## CD27+ Tbet+ CD11c+ signature
ComplexHeatmap::pheatmap(avg_CD27, 
                         color = colorRampPalette(c("#FF00FF", "black", "#FFFF00"))(100), 
                         cluster_rows = F, 
                         cluster_cols = F, 
                         fontsize_row = 15, 
                         fontsize_col = 15,
                         border_color = "black",
                         show_rownames = F, 
                         show_colnames = T, 
                         annotation_col = annotations, 
                         annotation_colors = anno_cols, 
                         annotation_legend = F, 
                         main = "CD27 signature", 
                         column_names_side = "top", 
                         angle_col = "45")
```
