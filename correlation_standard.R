library(dplyr)
library(pheatmap)
library(ggraph)
library(tidygraph) 

install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap)
library(circlize)


# Expression Correlation Analysis in pan-cancer data - Figure 5---------------------

data <- read.csv('Batch_corrected_Expression_Public_24Q4_subsetted.csv', stringsAsFactors = FALSE)
data <- data[, -c(1, 4:7)]

# Our potential metastasis markers
mechgenes <- c("ACTN4","ANXA2", "ANXA6", "CFL1", "PFN1", "EZR", "ECM1", "FSCN1", "TAGLN",  "MYH9", "GSN")

# Standard markers devided by functional groups 
invasion <- c("MMP1", "MMP2","MMP9", "CDC42", "RAC1","RHOA","BCAR3")
adhesion <- c("EPCAM", "CLDN1", "CLDN7")
EMT <- c("S100A4","CDH1","VIM", "MTDH", "NOTCH3","LGALS1")
stemness <- c("CD44", "CD24", "ALDH2", "PROM1", "LGR5", "POU5F1",  "NANOG", "SOX2", "BMI1")
signaling <-  c("STAT3", "YAP1", "VEGFC", "VEGFD", "FLT4")

comparison_genes <- unique(c(invasion, adhesion, EMT, stemness, signaling))
comparison_genes <- setdiff(comparison_genes, mechgenes)  

results <- data.frame(
  mech_gene = character(),
  comparison_gene = character(),
  pearson_cor = numeric(),
  pearson_p = numeric(),
  stringsAsFactors = FALSE)

for (g1 in mechgenes) {
  for (g2 in comparison_genes) {
    if (g1 %in% names(data) && g2 %in% names(data)) {
      pearson_test <- suppressWarnings(cor.test(data[[g1]], data[[g2]], method = "pearson"))
      
      results <- rbind(results, data.frame(
        mech_gene = g1,
        comparison_gene = g2,
        pearson_cor = pearson_test$estimate,
        pearson_p = pearson_test$p.value))}}}

results <- results[order(-abs(results$pearson_cor)), ]
head(results, 10)
# write_xlsx(results, "Correlation_mech_standard.xlsx")

significant_p <- subset(results, (pearson_p < 0.05 & abs(pearson_cor) > 0.5))

available_genes <- intersect(c(mechgenes, comparison_genes), names(data))
data_sub <- data[, available_genes]

cor_matrix_p <- cor(data_sub[, mechgenes], data_sub[, comparison_genes], method = "pearson", use = "pairwise.complete.obs")

strong_corr_mask <- abs(cor_matrix_p) >= 0.5
labels_matrix <- ifelse(strong_corr_mask, round(cor_matrix_p, 2), "")

gene_types <- sapply(comparison_genes, function(gene) {
  if (gene %in% invasion) return("Invasion & Motility")
  if (gene %in% adhesion) return("Adhesion")
  if (gene %in% EMT) return("EMT")
  if (gene %in% stemness) return("Stemness - associated")
  if (gene %in% signaling) return("Signaling")})

annotation_col <- data.frame(Marker = factor(gene_types))
rownames(annotation_col) <- comparison_genes

ann_colors <- c(
  `Invasion & Motility` = "#4C72B0",
  Adhesion = "#55A868",
  EMT = "#C44E52",
  `Stemness - associated` = "#8172B2",
  Signaling = "#CCB974")

col_annot <- HeatmapAnnotation(
  Marker = annotation_col$Marker,
  col = list(Marker = ann_colors),
  annotation_legend_param = list(
    title = "Marker",
    at = names(ann_colors),
    labels = names(ann_colors),
    labels_gp = gpar(fontsize = 10),
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    grid_height = unit(0.5, "cm"),  
    grid_width = unit(0.5, "cm")))
 
breaks <- unique(c(
  seq(-0.65, -0.5001, length.out = 60),  
  seq(-0.5, 0.5, length.out = 50),       
  seq(0.5001, 0.65, length.out = 60)))

col_fun <- colorRamp2(
  breaks = breaks,
  colors = colorRampPalette(c("blue", "white", "red"))(length(breaks)))

Heatmap(
  cor_matrix_p,
  name = "Pearson\nCorrelation",
  top_annotation = col_annot,
  col = col_fun, 
  cluster_rows = TRUE,  
  cluster_columns = TRUE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  column_title = "Correlation between expression of standard and potential markers in pan-cancer dataset",
  column_title_gp = gpar(fontsize = 16, fontface = "bold", col = "black"),
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    r <- cor_matrix_p[i, j]
    if (!is.na(r) && abs(r) > 0.5) {
      grid.rect(
        x = x, y = y, width = width, height = height,
        gp = gpar(col = "black", fill = NA, lwd = 1.5))}},
  
  heatmap_legend_param = list(
    title = "Correlation", 
    title_position = "leftcenter-rot",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 10),
    legend_height = unit(4, "cm"), 
    legend_width = unit(1.2, "cm"),
    direction = "vertical"))


# Copy Number correlation analysis  --------------------------------------------------

cn <- read.csv('Copy_Number_Public_24Q4_subsetted.csv', stringsAsFactors = FALSE)
cn <- cn[, -c(1, 4:7)]

results_cn <- data.frame(
  mech_gene = character(),
  comparison_gene = character(),
  pearson_cor = numeric(),
  pearson_p = numeric(),
  stringsAsFactors = FALSE)

for (g1 in mechgenes) {
  for (g2 in comparison_genes) {
    if (g1 %in% names(cn) && g2 %in% names(cn)) {
      pearson_test <- suppressWarnings(cor.test(cn[[g1]], cn[[g2]], method = "pearson"))
      results_cn <- rbind(results_cn, data.frame(
        mech_gene = g1,
        comparison_gene = g2,
        pearson_cor = pearson_test$estimate,
        pearson_p = pearson_test$p.value))}}}

results_cn <- results_cn[order(-abs(results_cn$pearson_cor)), ]
significant_p_cn <- subset(results_cn, (pearson_p < 0.05 & abs(pearson_cor) > 0.5))

available_genes_cn <- intersect(c(mechgenes, comparison_genes), names(cn))
data_sub <- cn[, available_genes_cn]

comparison_genes <- comparison_genes[-which(comparison_genes == "NANOG")]
cor_matrix_p_cn <- cor(data_sub[, mechgenes], data_sub[, comparison_genes], method = "pearson", use = "pairwise.complete.obs")

strong_corr_mask_cn <- abs(cor_matrix_p_cn) >= 0.5
labels_matrix_cn <- ifelse(strong_corr_mask_cn, round(cor_matrix_p_cn, 2), "")
an_col_cn <-  annotation_col[rownames(annotation_col) != "NANOG", , drop = FALSE]

col_annot <- HeatmapAnnotation(
  Marker = an_col_cn$Marker,
  col = list(Marker = ann_colors),
  annotation_legend_param = list(
    title = "Marker",
    at = names(ann_colors),
    labels = names(ann_colors),
    labels_gp = gpar(fontsize = 10),
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    grid_height = unit(0.5, "cm"),  
    grid_width = unit(0.5, "cm")))

Heatmap(
  cor_matrix_p_cn,
  name = "Pearson\nCorrelation",
  top_annotation = col_annot,
  col = col_fun, 
  cluster_rows = TRUE,  
  cluster_columns = TRUE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  column_title = "Correlation between Copy Number of standard and potential markers in pan-cancer dataset",
  column_title_gp = gpar(fontsize = 16, fontface = "bold", col = "black"),
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    r <- cor_matrix_p_cn[i, j]
    if (!is.na(r) && abs(r) > 0.5) {
      grid.rect(
        x = x, y = y, width = width, height = height,
        gp = gpar(col = "black", fill = NA, lwd = 1.5))}},
  
  heatmap_legend_param = list(
    title = "Correlation", 
    title_position = "leftcenter-rot",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 10),
    legend_height = unit(4, "cm"), 
    legend_width = unit(1.2, "cm"),
    direction = "vertical"))

# Graphs  - Figure 6 --------------------------------------------------------------
edges_pearson <- significant_p %>%
  transmute(
    from = mech_gene,
    to = comparison_gene,
    cor = pearson_cor)

graph_pearson <- tbl_graph(edges = edges_pearson, directed = FALSE)
graph_pearson <- graph_pearson %>%
  mutate(Marker = ifelse(name %in% mechgenes, "Mechbio", "Standard"))

set.seed(3)
ggraph(graph_pearson, layout = "fr") +
  geom_edge_link(
    aes(
      edge_width = abs(cor),
      edge_alpha = abs(cor),
      color = cor),
    show.legend = c(edge_width = TRUE, edge_color = FALSE, edge_alpha = FALSE)) +
  geom_node_point(
    aes(color = Marker),
    size = 5,
    show.legend = FALSE) +
  geom_node_text(aes(label = name), repel = TRUE, size = 5) +
  scale_edge_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    guide = "none") +
  scale_edge_alpha(range = c(0.2, 1)) +
  scale_color_manual(
    values = c(Mechbio = "dodgerblue", Standard = "forestgreen")) +
  guides(
    color = "none",
    edge_alpha = "none") +
  theme_void() +
  labs(edge_width = "Pearson |r|") +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12))