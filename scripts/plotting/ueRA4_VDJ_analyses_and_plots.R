library(Seurat)
library(dplyr)
library(ggplot2)
library(xlsx)

###################################################
###---Distribution of Immunoglobulin Isotypes---###
###################################################

ueRA4.1 = readRDS("data/ueRA4_w_subs") # read in subclustered data

igh = read.xlsx("data/VDJ_data/ueRA4_VDJ-seq_data.xlsx") # read in relevant VDJ-seq/IgH data for cluster 3.1
igh$cluster = forcats::as_factor(igh$cluster) # convert to factor

rownames(igh) = igh(1) # set barcodes as rownames
isotypes = igh[,"c_gene"] # Extract only isotypes

ueRA4.1 = AddMetaData(ueRA4.1, metadata = igh) # Add isotypes as metadata

#---Examine Distribution of Isotypes---#

igm = WhichCells(object = ueRA4.1, expression = c_gene == "IGHM")

igg = c(WhichCells(object = ueRA4.1, expression = c_gene == "IGHG1"), 
	WhichCells(object = ueRA4.1, expression = c_gene == "IGHG2"), 
	WhichCells(object = ueRA4.1, expression = c_gene == "IGHG3")) 

## NOTE: no cells in cluster 3.1 express IgG4

iga = c(WhichCells(object = ueRA4.1, expression = c_gene == "IGHA1"), 
	WhichCells(object = ueRA4.1, expression = c_gene == "IGHA2"))

p1 = DimPlot(ueRa4.1, cells.highlight = igm, cols.highlight = "cornflowerblue", label = TRUE, 
	pt.size = 0.85, sizes.highlight = 0.85, label.size = 7) + 
	ggtitle("IgM") +
	theme(plot.title = element_text(hjust = 0.5) +
	NoLegend()

p2 = DimPlot(ueRa4.1, cells.highlight = igg, cols.highlight = "brown4", label = TRUE, 
	pt.size = 0.85, sizes.highlight = 0.85, label.size = 7) + 
	ggtitle("IgA") +
	theme(plot.title = element_text(hjust = 0.5) +
	NoLegend()

p3 = DimPlot(ueRa4.1, cells.highlight = iga, cols.highlight = "tan1", label = TRUE, 
	pt.size = 0.85, sizes.highlight = 0.85, label.size = 7) + 
	ggtitle("IgA") +
	theme(plot.title = element_text(hjust = 0.5) +
	NoLegend()

grid.arrange(p1, p2, p3, nrow = 1)

#---Determine Proportions of Main Ig Isotypes by Cluster---#
## Exclude IgE and IgD (poorly detected)
tab = table(Idents(ueRA4.1), 
	ueRA4.1$Isotype)[,c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHM")]

## Calculate proportions
for (i in seq(rownames(tab))) {
	# Calculate percent of each isotype for each cluster
	props = c(
		IgM = tab[i, "IGHM"] / sum(tab[i,]) * 100,
		IgG = sum(tab[i, c("IGHG1", "IGHG2", "IGHG3", "IGHG4")]) / sum(tab[i,]) * 100, 
		IgA = sum(tab[i, c("IGHA1", "IGHA2")]) / sum(tab[i,]) * 100)
	print(paste("Cluster", rownames(tab)[i], sep = " "))
	print(round(props, 1))
}

# Plot for cluster 3.1
slices = c(5, 52, 43.9)
pie(slices, labels = NA, col = c("cornflowerblue", "brown4", "tan1"), 
	clockwise = TRUE, main = "Cluster 3.1 (VDJ-Seq)")
text(x = 0.4, y = -0.12, "52.0", cex = 2, col = "white")
text(x = -0.4, y = 0.15, "42", cex = 2, col = "white")
text(x = 0.09, y = 0.55, "5", cex = 2, col = "white")
legend(0.95, 0.35, c("IgM", "IgG", "IgA1"), cex = 2, 
	fill = c("cornflowerblue", "brown4", "tan1"), bty = "n")

#---Determine Proportions of IgG Subtypes---#
for(i in seq(rownames(tab1))) {
	props = c(
		IgG1 = tab[i, "IGHG1"] / sum(tab1[i, c["IGHG1", "IGHG2", "IGHG3", "IGHG4")] * 100,
		IgG2 = tab[i, "IGHG2"] / sum(tab1[i, c["IGHG1", "IGHG2", "IGHG3", "IGHG4")] * 100,
		IgG3 = tab[i, "IGHG3"] / sum(tab1[i, c["IGHG1", "IGHG2", "IGHG3", "IGHG4")] * 100,
		IgG4 = tab[i, "IGHG4"] / sum(tab1[i, c["IGHG1", "IGHG2", "IGHG3", "IGHG4")] * 100)
	print(paste("Cluster", rownames(tab1)[i], sep = " "))
	print(round(props, 1))
}

# Plot for Cluster 3.1
slices = c(73, 12, 15) # No cells in 3.1 express IgG4
pie(slices, labels = NA, col = c("brown4", "indianred3", "indianred1"), 
	clockwise = TRUE, main = "Cluster 3.1 (VDJ-Seq)")
text(x = 0.27, y = -0.18, "73", cex = 2, col = "white")
text(x = -0.5, y = 0.15, "12", cex = 2, col = "white")
text(x = -0.25, y = 0.52, "15", cex = 2, col = "white")
legend(0.95, 0.35, c("IgG1", "IgG2", "IgG3"), cex = 2, 
	fill = c("brown4", "indianred3", "indianred1"), bty = "n")

#######################################
###---Examine V- and J-Gene Usage---###
#######################################

# igh = read.xlsx("data/VDJ_data/ueRA4_VDJ-seq_data.xlsx") # read in relevant VDJ-seq/IgH data for cluster 3.1
# igh$cluster = forcats::as_factor(igh$cluster) # convert to factor

#---Examine V-Gene Usage---#
ighv = igh %>% select(v_gene, cluster) # extract V genes
ighv = ighv %>% group_by(cluster) %>% count(v_gene) # calculate frequency of each identified V gene

## Calculate percentage usage
ighv_perc = ighv %>%
	group_by(cluster, v_gene) %>%
	summarise(n) %>%
	mutate(percentage = (n/sum(n))*100)

## Plot V-gene usage
ggplot(data = ighv_perc, aes(x = v_gene, y = percentage, fill = "#00B4F0")) +
	geom_col(position = position_dodge()) +
	theme_bw() +
	scale_fill_manual(values = "#00B4F0") +
	facet_grid(rows = vars(cluster)) +
	xlab("V gene") +
	ylab("Percentage (%)") +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
	NoLegend()

#---Examine J-Gene Usage---#
ighj = igh %>% select(j_gene, cluster) # extract J genes
ighj = ighj %>% group_by(cluster) %>% count(j_gene) # calculate frequency of each identified J gene

## Calculate percentage usage
ighj_perc = ighj %>%
	group_by(cluster, j_gene) %>%
	summarise(n) %>%
	mutate(percentage = (n/sum(n))*100)

## Plot J-gene usage
ggplot(data = ighj_perc, aes(x = j_gene, y = percentage, fill = "#00B4F0")) + 
	geom_col(position = position_dodge()) +
	theme_bw() +
	scale_fill_manual(values = "#00B4F0") +
	facet_grid(rows = vars(cluster)) +
	xlab("J gene") +
	ylab("Percentage (%)") +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
	NoLegend()

#############################################
###---Examine HCDR3 Length Distribution---###
#############################################

# Calculate HCDR3 lengths
igh$cdr3_length = nchar(igh$cdr3)

# Plot length distribution
ggplot(data = igh, aes(cdr3_length)) +
  geom_bar(fill = "#00B4F0") +
  theme_bw() +
  ylab("Frequency") + 
  xlab("H_CDR3 length") + 
  theme(panel.border = element_blank(), 
        axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 
