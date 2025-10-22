library(dplyr)
library(readxl)
library(openxlsx)
library(tidyr)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(patchwork)
library(grid)
library(forcats)
library(purrr)

cell_lines <- read.csv('/data/Model.csv', stringsAsFactors = FALSE)
cell_lines <- cell_lines %>% filter(GrowthPattern != 'Suspension') 

# cell lines only cancer
cell_lines_canc <- cell_lines %>% filter(OncotreeLineage %in% c('Breast', 'Pancreas', 'CNS/Brain') & OncotreePrimaryDisease != 'Non-Cancerous')
cell_lines_canc %>% group_by(OncotreeLineage) %>% summarise(n=n())

# cell lines normal and non-cancerous across all organs
lines_noncanc <- cell_lines %>% filter(OncotreeLineage == 'Normal' | OncotreePrimaryDisease == 'Non-Cancerous') 

metmap <- read_xlsx('/data/Supplementary Table 04 MetMap 500 met potential.xlsx', sheet=6)  

metmap$...1 <- sub("_.*", "", metmap$...1)

cell_lines_canc <- cell_lines_canc %>%
  mutate(MetMap = ifelse(StrippedCellLineName %in% metmap$"...1", metmap$mean[match(StrippedCellLineName, metmap$"...1")], NA))

lines_noncanc <- lines_noncanc %>%
  mutate(MetMap = ifelse(StrippedCellLineName %in% metmap$"...1", metmap$mean[match(StrippedCellLineName, metmap$"...1")], NA))

cell_lines_canc <- cell_lines_canc %>% mutate(Metastatic_MetMap = ifelse(cell_lines_canc$MetMap >= -2, 'High confidence',ifelse(cell_lines_canc$MetMap <= -4, 'Non-metastatic', 'Low confidence')))

lines_noncanc <- lines_noncanc %>% mutate(lines_noncanc = ifelse(lines_noncanc$MetMap >= -2, 'High confidence',ifelse(lines_noncanc$MetMap <= -4, 'Non-metastatic', 'Low confidence')))

# write.xlsx(cell_lines_canc, file = "/data/Cell_lines_canc_metmap.xlsx", overwrite = TRUE)

# Get Data  ------------------------------------------------------------

# --- Brain ---
Brain_mechbio <- read.csv('/data/CNS_Brain.csv', stringsAsFactors = FALSE)

Brain_mechbio <- Brain_mechbio %>%
  left_join(cell_lines_canc %>% select(ModelID, PrimaryOrMetastasis, Metastatic_MetMap, OncotreeSubtype), by = c("X" = "ModelID"))

Brain_mechbio <- Brain_mechbio[!is.na(Brain_mechbio$OncotreeSubtype), ]

Brain_mechbio$aggressivness <- ifelse(Brain_mechbio$OncotreeSubtype %in% c('Glioblastoma', 'Medulloblastoma', 'Atypical Teratoid/Rhabdoid Tumor', 'Primitive Neuroectodermal Tumor', 'Diffuse Intrinsic Pontine Glioma', 'Gliosarcoma'), 'High-Aggressive', 'Low-Aggressive')

vec1 <- cell_lines_canc %>% filter(OncotreeLineage == "CNS/Brain") %>% pull(ModelID)
vec2 <- Brain_mechbio$X
setequal(vec1, vec2)
setdiff(vec1, vec2)  
setdiff(vec2, vec1)
Brain_mechbio <- Brain_mechbio %>% filter(X %in% vec1)
Brain_mechbio$Organ <- 'CNS/Brain'

# --- Pancreas ---
Pancreas_mechbio <- read.csv('/data/Pancreas.csv', stringsAsFactors = FALSE)
Pancreas_mechbio <- Pancreas_mechbio %>%
  left_join(cell_lines_canc %>% select(ModelID, PrimaryOrMetastasis, Metastatic_MetMap), by = c("X" = "ModelID"))

vec1 <- cell_lines_canc %>% filter(OncotreeLineage == "Pancreas") %>% pull(ModelID)
vec2 <- Pancreas_mechbio$X
setequal(vec1, vec2)
setdiff(vec1, vec2)  
setdiff(vec2, vec1)
Pancreas_mechbio <- Pancreas_mechbio %>% filter(X %in% vec1)
Pancreas_mechbio$Organ <- 'Pancreas'

#  --- Breast ---
Breast_mechbio <- read.csv('/data/Breast.csv', stringsAsFactors = FALSE)

Breast_mechbio <- Breast_mechbio %>%
  left_join(cell_lines_canc %>% select(ModelID, OncotreePrimaryDisease, OncotreeSubtype,  PrimaryOrMetastasis, Metastatic_MetMap, GrowthPattern), by = c("X" = "ModelID"))
vec1 <- cell_lines_canc %>% filter(OncotreeLineage == "Breast") %>% pull(ModelID)
vec2 <- Breast_mechbio$X
setequal(vec1, vec2)
setdiff(vec1, vec2)  
setdiff(vec2, vec1)
Breast_mechbio <- Breast_mechbio %>% filter(X %in% vec1)
Breast_mechbio$PrimaryOrMetastasis[Breast_mechbio$PrimaryOrMetastasis == ""] <- NA
Breast_mechbio$Organ <- 'Breast'

# --- Normal ---
Non_canc_mechbio <- read.csv('/data/Non-canc.csv', stringsAsFactors = FALSE)

Non_canc_mechbio <- Non_canc_mechbio %>%
  left_join(lines_noncanc %>% select(ModelID, OncotreeLineage), by = c("X" = "ModelID"))
vec1 <- lines_noncanc %>% pull(ModelID)
vec2 <- Non_canc_mechbio$X
setequal(vec1, vec2)
setdiff(vec1, vec2)  
setdiff(vec2, vec1)
Non_canc_mechbio <- Non_canc_mechbio %>% filter(X %in% vec1)
str(Non_canc_mechbio)

summary_stats <- function(x) {
  x <- x[!is.na(x)]
  mean_val <- round(mean(x, na.rm = TRUE), 2)
  sd_val <- round(sd(x, na.rm = TRUE), 2)
  count_val <- length(x)
  return(c(mean = mean_val, sd = sd_val, count_val = count_val))
}
results <- as.data.frame(t(apply(Non_canc_mechbio[, -c(1, 46:49)], 2, summary_stats)))
results$mean_sd <- paste0(results$mean, " ± ", results$sd)

results2 <- as.data.frame(t(apply(Brain_mechbio %>% filter(aggressivness == 'Low-Aggressive') %>% select(2:45), 2, summary_stats)))
results2$mean_sd <- paste0(results2$mean, " ± ", results2$sd)

# count statistics --------------------------------------------------------

num_cols <- sapply(Brain_mechbio, is.numeric)
gene_data <- Brain_mechbio[, num_cols]

safe_wilcox <- function(x, group) {
  df <- data.frame(x = x, group = group)
  df <- df[!is.na(df$x) & !is.na(df$group), ]
  
  if (length(unique(df$group)) < 2 || any(table(df$group) < 2)) {
    return(NA)
  }
  
  tryCatch(wilcox.test(x ~ group, data = df)$p.value, error = function(e) NA)
}

p_aggr <- sapply(gene_data, function(col) safe_wilcox(col, Brain_mechbio$aggressivness))
p_prim <- sapply(gene_data, function(col) safe_wilcox(col, Breast_mechbio$PrimaryOrMetastasis))
p_meta <- sapply(gene_data, function(col) safe_wilcox(col, Brain_mechbio$Metastatic_MetMap))

pvals <- data.frame(
  Feature = names(gene_data),
  p_aggressivness = p_aggr,
  # p_PrimMet = p_prim,
  p_metmap = p_meta)

# FDR correction
pvals$p_aggressivness_adj <- p.adjust(pvals$p_aggressivness, method = "fdr")
pvals$p_PrimMet_adj <- p.adjust(pvals$p_PrimMet, method = "fdr")
pvals$p_metmap_adj <- p.adjust(pvals$p_metmap, method = "fdr")

pvals <- pvals[order(pvals$p_metmap_adj), ]

# Kruskal- Wallis test -------------------------------------------------------------------
all_data_pan <- Breast_mechbio[, c(2:45, 48, 49, 51)] 

Brain_mechbio_anova <- Brain_mechbio[, c(2:47, 49, 50)] 
Brain_mechbio_anova$PrimaryOrMetastasis <- Brain_mechbio_anova$aggressivness
Brain_mechbio_anova <- Brain_mechbio_anova[, -47]
Non_canc_mechbio$PrimaryOrMetastasis <- 'normal'
Non_canc_mechbio$Metastatic_MetMap <- 'normal'
Non_canc_mechbio$Organ <- 'normal'
all_data_pan <- rbind(all_data_pan, Brain_mechbio_anova, Pancreas_mechbio[, c(2:48)], Non_canc_mechbio[,-c(1, 46)])

str(all_data_pan)
all_data_pan$PrimaryOrMetastasis <- as.factor(all_data_pan$PrimaryOrMetastasis)
all_data_pan$Metastatic_MetMap <- as.factor(all_data_pan$Metastatic_MetMap)
all_data_pan$Organ <- as.factor(all_data_pan$Organ)

run_kruskal_all <- function(df, grouping_column, target_columns, alpha = 0.05, adjust_method = "BH") {
  results <- purrr::map_dfr(target_columns, function(param) {
    temp_df <- df %>%
      select(all_of(grouping_column), all_of(param)) %>%
      drop_na() %>%
      rename(group = all_of(grouping_column), value = all_of(param))
    if (length(unique(temp_df$group)) < 2) return(NULL)
    test <- kruskal.test(value ~ group, data = temp_df)
    
    tibble(
      Parameter = param,
      p = test$p.value)})
  
  results <- results %>%
    mutate(p.adj = p.adjust(p, method = adjust_method)) %>%
    filter(p.adj < alpha)
  return(results)}

test_columns <- colnames(all_data_pan)[grepl("Copy.Number|RNAi|CRISPR|Expression", colnames(all_data_pan))]

results_BN_Met <- run_kruskal_all(
  df = all_data_pan %>% filter(Organ %in% c('Breast', 'normal')),
  grouping_column = "Metastatic_MetMap",
  target_columns = test_columns)

results_BN_Prim <- run_kruskal_all(
  df = all_data_pan %>% filter(Organ %in% c('Breast', 'normal')),
  grouping_column = "PrimaryOrMetastasis",
  target_columns = test_columns)

results_PN_Met <- run_kruskal_all(
  df = all_data_pan %>% filter(Organ %in% c('Pancreas', 'normal')),
  grouping_column = "Metastatic_MetMap",
  target_columns = test_columns)

results_PN_Prim <- run_kruskal_all(
  df = all_data_pan %>% filter(Organ %in% c('Pancreas', 'normal')),
  grouping_column = "PrimaryOrMetastasis",
  target_columns = test_columns)

results_BCN_Met <- run_kruskal_all(
  df = all_data_pan %>% filter(Organ %in% c('CNS/Brain', 'normal')),
  grouping_column = "Metastatic_MetMap",
  target_columns = test_columns)

results_BCN_Prim <- run_kruskal_all(
  df = all_data_pan %>% filter(Organ %in% c('CNS/Brain', 'normal')),
  grouping_column = "PrimaryOrMetastasis",
  target_columns = test_columns)

# Heatmaps with p-value - Figure 1 -------------------------------------------------------

names(all_data_pan) <- gsub("^Batch\\.corrected\\.Expression\\.Public\\.24Q4\\.", "Expression ", names(all_data_pan))
names(all_data_pan) <- gsub("^Copy\\.Number\\.Public\\.24Q4\\.", "Copy Number ", names(all_data_pan))
names(all_data_pan) <- gsub("^CRISPR\\.\\.DepMap\\.Public\\.24Q4\\.Score\\.\\.Chronos\\.\\.", "Gene Effect (CRISPR) ", names(all_data_pan))
names(all_data_pan) <- gsub("^RNAi\\.\\.Achilles\\.DRIVE\\.Marcotte\\.\\.DEMETER2\\.\\.", "RNAi ", names(all_data_pan))

all_data_pan <- all_data_pan %>%
  mutate(
    PrimaryOrMetastasis = case_when(
      PrimaryOrMetastasis == "High-Aggressive" ~ "Highly aggressive",
      PrimaryOrMetastasis == "Low-Aggressive"  ~ "Lowly aggressive",
      TRUE ~ as.character(PrimaryOrMetastasis)
    ),
    PrimaryOrMetastasis = factor(PrimaryOrMetastasis))

run_kruskal_heatmap <- function(df, group_var = "Metastatic_MetMap", organ) {
  df[[group_var]] <- as.factor(df[[group_var]])
  df_filtered <- df %>%
    filter(Organ %in% c(organ, 'normal') & !is.na(.data[[group_var]]))
  
  data_cols <- colnames(df_filtered)[str_detect(colnames(df_filtered),
                                                "Copy Number|RNAi|CRISPR|Expression")]
  
  results <- map_dfr(data_cols, function(col) {
    temp_df <- df_filtered %>%
      select(all_of(group_var), !!sym(col)) %>%
      rename(Value = !!sym(col)) %>%
      drop_na()
    
    if (n_distinct(temp_df[[group_var]]) < 2) return(NULL)
    
    test <- kruskal.test(as.formula(paste("Value ~", group_var)), data = temp_df)
    
    parts <- str_match(col, "^(.*?)\\s+(\\S+)$")
    param <- parts[, 2]
    gene  <- parts[, 3]
    
    tibble(
      Gene = gene,
      Parameter = param,
      p_value = test$p.value
    )
  })
  
  results <- results %>%
    mutate(Significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    ))
  
  results$Significance <- factor(results$Significance, levels = c("***", "**", "*", "ns"))
  
  sig_colors <- c("***" = "#990000", "**" = "#d7301f", "*" = "#fc9272", "ns" = "grey90")
  
  ggplot(results, aes(x = Gene, y = fct_rev(Parameter), fill = Significance)) +
    geom_tile(color = "white") +
    scale_fill_manual(values = sig_colors, name = "p-value") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    labs(x = NULL, y = NULL)
}

p1 <- run_kruskal_heatmap(df = all_data_pan, group_var = "Metastatic_MetMap", organ = "Breast")
p2 <- run_kruskal_heatmap(df = all_data_pan, group_var = "Metastatic_MetMap", organ = "Pancreas")
p3 <- run_kruskal_heatmap(df = all_data_pan, group_var = "Metastatic_MetMap", organ = "CNS/Brain")
p1_prim <- run_kruskal_heatmap(df = all_data_pan, group_var = "PrimaryOrMetastasis", organ = "Breast")
p2_prim <- run_kruskal_heatmap(df = all_data_pan, group_var = "PrimaryOrMetastasis", organ = "Pancreas")
p3_prim <- run_kruskal_heatmap(df = all_data_pan, group_var = "PrimaryOrMetastasis", organ = "CNS/Brain")

strip_y <- theme(axis.text.y = element_blank(),
                 axis.ticks.y = element_blank(),
                 axis.title.y = element_blank())

strip_x <- theme(axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.title.x = element_blank())

p1 <- p1 + theme(axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.title.x = element_blank()) + theme(legend.position = "none")
p2 <- p2 + theme(axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.title.x = element_blank()) + theme(legend.position = "none")
p3 <- p3 + theme(axis.title.x = element_blank()) + theme(legend.position = "none")

p1_prim <- p1_prim + theme(axis.text.y = element_blank(),
                           axis.ticks.y = element_blank(),
                           axis.title.y = element_blank(),
                           axis.text.x = element_blank(),
                           axis.ticks.x = element_blank(),
                           axis.title.x = element_blank()) + theme(legend.position = "none")
p2_prim <- p2_prim + theme(axis.text.y = element_blank(),
                           axis.ticks.y = element_blank(),
                           axis.title.y = element_blank(),
                           axis.text.x = element_blank(),
                           axis.ticks.x = element_blank(),
                           axis.title.x = element_blank())
p3_prim <- p3_prim + theme(axis.text.y = element_blank(),
                           axis.ticks.y = element_blank(),
                           axis.title.y = element_blank()) + theme(legend.position = "none")

col_title_metmap <- wrap_elements(grid::textGrob("MetMap level", gp = gpar(fontsize = 16, fontface = "bold")))
col_title_tumor <- wrap_elements(grid::textGrob("Tumor Type", gp = gpar(fontsize = 16, fontface = "bold")))

organ_title_b <- wrap_elements(grid::textGrob("Breast", gp = gpar(fontsize = 14, fontface = "bold")))
organ_title_p <- wrap_elements(grid::textGrob("Pancreas", gp = gpar(fontsize = 14, fontface = "bold")))
organ_title_c <- wrap_elements(grid::textGrob("CNS/Brain", gp = gpar(fontsize = 14, fontface = "bold")))

x_axis_label <- wrap_elements(grid::textGrob("Gene", gp = gpar(fontsize = 14, fontface = "bold")))
y_axis_label <- wrap_elements(grid::textGrob("Parameter", rot = 90, gp = gpar(fontsize = 14, fontface = "bold")))

empty <- ggplot() + theme_void()

final_plot <- 
  (  col_title_metmap | col_title_tumor ) /
  (empty | organ_title_b | empty) /
  (  p1 | p1_prim) /
  (empty | organ_title_p | empty) /
  (  p2 | p2_prim) /
  (empty | organ_title_c | empty) /
  ( p3 | p3_prim)  +
  plot_layout(
    heights = c(0.7, 0.55, 5, 0.55, 5, 0.55, 5),
    widths = c(0.7, 5, 5, 0.7)  
  ) &
  theme(plot.margin = margin(10, 10, 10, 10))

final_plot

# Difference plots only significant - Figure 3 ---------------------------------------

long_data <- all_data_pan %>%
  filter(Organ != 'normal') %>%
  pivot_longer(
    cols = matches("^(Copy Number|Gene Effect \\(CRISPR\\)|Expression|RNAi)"),
    names_to = "FullName",
    values_to = "Value"
  ) %>%
  mutate(
    DataType = case_when(
      str_detect(FullName, "^Copy Number") ~ "Copy Number",
      str_detect(FullName, "^Gene Effect \\(CRISPR\\)") ~ "Gene Effect (CRISPR)",
      str_detect(FullName, "^Expression") ~ "Expression",
      str_detect(FullName, "^RNAi") ~ "RNAi",
      TRUE ~ NA_character_
    ),
    Gene = case_when(
      DataType == "Copy Number" ~ str_remove(FullName, "Copy Number "),
      DataType == "Gene Effect (CRISPR)" ~ str_remove(FullName, "Gene Effect \\(CRISPR\\) "),
      DataType == "Expression" ~ str_remove(FullName, "Expression "),
      DataType == "RNAi" ~ str_remove(FullName, "RNAi "),
      TRUE ~ FullName
    ),
    DataType = factor(DataType, levels = c("Gene Effect (CRISPR)", "Expression", "RNAi", "Copy Number"))
  ) %>%
  drop_na(Value)

get_significant_wilcox <- function(data, group_var, threshold = 0.05) {
  data %>%
    filter(!is.na(.data[[group_var]])) %>%
    group_by(Organ, Gene, DataType) %>%
    filter(n_distinct(.data[[group_var]]) == 2) %>%  
    summarise(
      p_value = tryCatch(
        wilcox.test(Value ~ .data[[group_var]])$p.value,
        error = function(e) NA
      ),
      .groups = "drop"
    ) %>%
    filter(!is.na(p_value) & p_value < threshold) %>%
    mutate(ComparisonType = group_var)}

sig_metmap <- get_significant_wilcox(long_data, "Metastatic_MetMap")
sig_pom    <- get_significant_wilcox(long_data, "PrimaryOrMetastasis")
sig_all <- bind_rows(sig_metmap, sig_pom)

group_levels_extended <- function(comparison, organ) {
  if (comparison == "PrimaryOrMetastasis" && organ == "CNS/Brain") {
    return(c("Low-Aggressive", "High-Aggressive"))
  } else if (comparison == "PrimaryOrMetastasis") {
    return(c("Primary", "Metastatic"))
  } else if (comparison == "Metastatic_MetMap") {
    return(c("Low confidence", "High confidence"))
  } else {
    return(NULL) }}

plot_data <- long_data %>%
  pivot_longer(
    cols = c("PrimaryOrMetastasis", "Metastatic_MetMap"),
    names_to = "ComparisonType",
    values_to = "Group") %>%
  filter(!is.na(Group))

summary_data <- plot_data %>%
  group_by(Organ, Gene, DataType, ComparisonType, Group) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    N = n(),
    .groups = "drop") %>%
  pivot_wider(
    names_from = Group,
    values_from = c(Mean, SD, N),
    names_sep = "_")

bar_data <- summary_data %>%
  rowwise() %>%
  mutate(
    levels = list(group_levels_extended(ComparisonType, Organ)),
    Group1 = levels[1],
    Group2 = levels[2],
    Facet = paste(Gene, as.character(DataType), sep = ", "),
    col_mean_1 = paste0("Mean_", Group1),
    col_mean_2 = paste0("Mean_", Group2),
    col_sd_1 = paste0("SD_", Group1),
    col_sd_2 = paste0("SD_", Group2),
    col_n_1 = paste0("N_", Group1),
    col_n_2 = paste0("N_", Group2),
    Diff = if (col_mean_1 %in% names(cur_data()) && col_mean_2 %in% names(cur_data())) {
      cur_data()[[col_mean_2]] - cur_data()[[col_mean_1]]
    } else {
      NA_real_},
    SE = if (
      col_sd_1 %in% names(cur_data()) &&
      col_sd_2 %in% names(cur_data()) &&
      col_n_1 %in% names(cur_data()) &&
      col_n_2 %in% names(cur_data())) {
      sqrt(
        (cur_data()[[col_sd_1]]^2 / cur_data()[[col_n_1]]) +
          (cur_data()[[col_sd_2]]^2 / cur_data()[[col_n_2]]))} else {NA_real_},
    signif = NA_character_) %>%
  ungroup()

sig_combos <- sig_all %>%
  mutate(
    Facet = paste(Gene, as.character(DataType), sep = ", "),
    signif = case_when(
      p_value < 0.001 ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ NA_character_)) %>%
  select(Organ, Facet, ComparisonType, signif)

bar_data <- bar_data %>%
  left_join(sig_combos, by = c("Organ", "Facet", "ComparisonType")) %>%
  mutate(signif = ifelse(is.na(signif.y), NA, signif.y)) %>%
  select(-signif.x, -signif.y)

panels_with_significance <- sig_combos %>%
  distinct(Facet, ComparisonType)

bar_data_filtered <- bar_data %>%
  semi_join(panels_with_significance, by = c("Facet", "ComparisonType")) %>%
  mutate(
    DataType = factor(DataType, levels = c("Gene Effect (CRISPR)", "Expression", "RNAi", "Copy Number")),
    ComparisonType = recode(ComparisonType,
                            "PrimaryOrMetastasis" = "Tumor Type",
                            "Metastatic_MetMap" = "MetMap level"),
    Facet = Gene) %>% 
  arrange(DataType) %>% 
  mutate(Facet = factor(Facet, levels = unique(Facet)))  

facet_structure <- list(
  "Gene Effect (CRISPR)" = list(
    c("ANXA2", "Tumor Type"),
    c("ANXA6", "Tumor Type"),
    c("ECM1",  "Tumor Type")
  ),
  "Expression" = list(
    c("FSCN1", "Tumor Type"),
    c("FSCN1", "MetMap level"),
    c("EZR",   "MetMap level")
  ),
  "RNAi" = list(
    c("ANXA2", "MetMap level"),
    c("PFN1",  "MetMap level"),
    c("EZR",   "MetMap level")
  ),
  "Copy Number" = list(
    c("PFN1",  "Tumor Type"),
    c("ANXA6", "MetMap level"),
    c("EZR",   "MetMap level"),
    c("MYH9",  "MetMap level")))

data_type_order <- c("Gene Effect (CRISPR)", "Expression", "RNAi", "Copy Number")

make_bar_plot <- function(df) {
  ggplot(df, aes(x = Organ, y = Diff, fill = Organ)) +
    geom_bar(stat = "identity", color = "black") +
    geom_errorbar(aes(ymin = Diff - SE, ymax = Diff + SE), width = 0.2) +
    geom_text(
      data = df %>% filter(!is.na(signif)),
      aes(label = signif, y = Diff + SE + 0.2 * SE),
      size = 5, color = "red"
    ) +
    labs(title = paste(df$Gene[1], df$ComparisonType[1])) +
    theme_bw() +
    scale_fill_brewer(palette = "Set2") +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(size = 12),
      axis.text.x = element_text(size = 12))}

empty_plot <- ggplot() + 
  theme_void()
plot_rows <- list()

for (dtype in names(facet_structure)) {
  items <- facet_structure[[dtype]]

  panel_list <- lapply(items, function(x) {
    gene <- x[1]
    comp <- x[2]
    
    df <- bar_data_filtered %>%
      filter(Gene == gene, DataType == dtype, ComparisonType == comp)
    if (nrow(df) == 0) return(empty_plot)
    make_bar_plot(df) })

  while (length(panel_list) < 4) {
    panel_list <- append(panel_list, list(empty_plot))}

  label <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = dtype, angle = 90, size = 3.7) +
    theme_void()
  row_plot <- wrap_plots(c(list(label), panel_list), nrow = 1, widths = c(0.1, rep(1, 4)))
  plot_rows <- append(plot_rows, list(row_plot))}

final_plot <- wrap_plots(plot_rows, ncol = 1)
final_plot

# Heatmaps for means - Figure 2------------------------------------------------------------

analyze_and_plot_from_all_data_pan <- function(all_data_pan_df, mode = c("PrimaryOrMetastasis", "Metastatic_MetMap")) {
  mode <- match.arg(mode)
  get_group_col <- function(cancer_type_str) { 
    if (mode == "PrimaryOrMetastasis") {
      return("PrimaryOrMetastasis") 
    } else { # mode == "Metastatic_MetMap"
      return("Metastatic_MetMap")}}
  
  get_valid_groups <- function(cancer_type_str) {
    if (cancer_type_str == "CNS/Brain" && mode == "PrimaryOrMetastasis") {
      return(c("Highly aggressive", "Lowly aggressive"))
    } else if (mode == "PrimaryOrMetastasis") {
      return(c("Primary", "Metastatic"))
    } else {
      return(c("High confidence", "Low confidence")) }}

  safe_wilcox <- function(x_values, group_values) {
    group_as_factor <- factor(group_values)
    df_for_test <- data.frame(x = x_values, group = group_as_factor)
    df_for_test <- df_for_test[!is.na(df_for_test$x) & !is.na(df_for_test$group), ]
    df_for_test$group <- droplevels(df_for_test$group)
    group_counts <- table(df_for_test$group)
    if (length(group_counts) < 2 || any(group_counts < 2)) {
      return(NA_real_)}
    
    tryCatch(
      wilcox.test(x ~ group, data = df_for_test, exact = FALSE)$p.value, 
      error = function(e) NA_real_)}
  
  extract_and_compare <- function(df_subset_for_cancer, current_cancer_type_name) {
    group_col_name <- get_group_col(current_cancer_type_name)
    valid_group_levels <- get_valid_groups(current_cancer_type_name)
    
    df_subset_for_cancer[[group_col_name]] <- factor(df_subset_for_cancer[[group_col_name]])
    processed_data_for_cancer <- df_subset_for_cancer %>%
      filter(.data[[group_col_name]] %in% valid_group_levels) %>%
      
      mutate(!!sym(group_col_name) := droplevels(.data[[group_col_name]]))
    
    processed_data_for_cancer %>%
      pivot_longer( cols = matches("^(Copy Number|Expression|RNAi|Gene Effect \\(CRISPR\\))"),
                   names_to = "RawParameter", values_to = "Value", values_drop_na = TRUE) %>% 
      mutate(
        Marker = str_extract(RawParameter, "[^ ]+$"),
        Parameter = str_extract(RawParameter, "^(Copy Number|Expression|RNAi|Gene Effect \\(CRISPR\\))")
      ) %>%
      filter(Parameter != "Other", !is.na(Marker)) %>% 
      group_by(CancerType = current_cancer_type_name, Parameter, Marker) %>% 
      group_modify(~ {
        p_val <- safe_wilcox(.x$Value, .x[[group_col_name]])
        
        .x %>%
          group_by(Status = .data[[group_col_name]]) %>% 
          summarise(
            Mean = mean(Value, na.rm = TRUE),
            SD = sd(Value, na.rm = TRUE),
            N = sum(!is.na(Value)), 
            .groups = "drop"       
          ) %>%
          mutate(p_value = p_val)}) %>%
      ungroup()}

  all_data_pan_prepared <- all_data_pan_df %>%
    filter(Organ != "normal", !is.na(Organ)) %>%
    mutate(CancerType_standardised = case_when(
      Organ == "CNS/Brain" ~ "CNS/Brain",
      TRUE ~ as.character(Organ)))
  
  list_of_df_by_cancer <- split(all_data_pan_prepared, all_data_pan_prepared$CancerType_standardised)
  data_summary <- imap_dfr(list_of_df_by_cancer, function(df_subset_for_cancer, current_cancer_type_name) {
    extract_and_compare(df_subset_for_cancer = df_subset_for_cancer, 
                        current_cancer_type_name = current_cancer_type_name)})
  
  data_summary <- data_summary %>%
    mutate(
      Group = paste(CancerType, Status),
      ValueLabel = ifelse(N > 0, sprintf("n=%d",N), NA),
      SigLabel = case_when(
        is.na(p_value) ~ "",
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "" ))
  
  plot_title_mode_specific <- if (mode == "PrimaryOrMetastasis") {
    "Tumor Type"
  } else { 
    "MetMap level"}
  
  plot_list <- data_summary %>%
    split(.$Parameter) %>%
    map(~ {
      plot_data <- .x
      if(nrow(plot_data) == 0) return(NULL) 
      
      current_midpoint <- median(plot_data$Mean, na.rm = TRUE)
      if(is.na(current_midpoint) || !is.finite(current_midpoint)) current_midpoint <- 0 
      plot_data$Marker <- factor(plot_data$Marker, levels = sort(unique(plot_data$Marker)))
      plot_data$Group <- factor(plot_data$Group, levels = unique(plot_data$Group[order(plot_data$CancerType, plot_data$Status)]))
      
  pretty_parameter <- case_when(
        unique(plot_data$Parameter) == "Gene Effect (CRISPR)" ~ "CRISPR",
        TRUE ~ unique(plot_data$Parameter))    
      
      ggplot(plot_data, aes(x = Marker, y = fct_rev(Group), fill = Mean)) +
        geom_tile(color = "white", linewidth = 0.5) + 
        geom_text(aes(label = ValueLabel), size = 4, na.rm = TRUE) + 
        geom_text(aes(label = SigLabel), vjust = 4, size = 6, color = "black", nudge_y = 0.05, na.rm = TRUE) + 
        scale_fill_gradientn(
          colours = scales::col_numeric(palette = c("red", "white", "green"), domain = NULL)(seq(0, 1, length.out = 20)),
          name = paste(unique(pretty_parameter), "Mean"),
          na.value = "grey80" 
        ) +
        theme_minimal(base_size = 12) + 
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 14),
          strip.text = element_text(face = "bold", size = 14), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 14, hjust = 0.5) 
        ) +
        labs(
          # title = paste0(unique(plot_data$Parameter), " for ", plot_title_mode_specific),
          x = "Gene",
          y = "")}) %>%
    compact() 
  
    walk(plot_list, print)
  return(invisible(plot_list)) }

analyze_and_plot_from_all_data_pan(all_data_pan, mode = "PrimaryOrMetastasis")
analyze_and_plot_from_all_data_pan(all_data_pan, mode = "Metastatic_MetMap")


# Z-score all genes - Figure 4 -------------------------------------------------------

calculate_zscore <- function(df, type = c("Expression", "CopyNumber", "CRISPR", "RNAi")) {
  type <- match.arg(type)
  pattern_list <- list(
    Expression = "^Expression\\s",
    CopyNumber = "^Copy Number\\s",
    CRISPR     = "^Gene Effect \\(CRISPR\\)\\s",
    RNAi       = "^RNAi\\s")
  pattern <- pattern_list[[type]]
  cols <- grep(pattern, colnames(df), value = TRUE)
  z_matrix <- scale(df[, cols]) 
  cumulative_z <- rowMeans(z_matrix, na.rm = TRUE)
  df[[paste0(type, "_zscore")]] <- cumulative_z
  return(df)}

prepare_zscore_data <- function(df,
                                selected_genes = NULL,
                                group_var = "Metastatic_MetMap") {
  df <- df |>
    select(-any_of(c("Expression_zscore", "CopyNumber_zscore", "CRISPR_zscore", "RNAi_zscore")))
  
  if (is.null(selected_genes)) {
    df <- df |>
      calculate_zscore("Expression") |>
      calculate_zscore("CopyNumber") |>
      calculate_zscore("CRISPR") |>
      calculate_zscore("RNAi")
  } else {
    select_gene_columns <- function(prefix) {
      pattern <- paste0("^", prefix, "\\s(", paste0(selected_genes, collapse = "|"), ")$")
      grep(pattern, colnames(df), value = TRUE)}
    
    expr_cols   <- select_gene_columns("Expression")
    cnv_cols    <- select_gene_columns("Copy Number")
    crispr_cols <- select_gene_columns("Gene Effect \\(CRISPR\\)")
    rnai_cols   <- select_gene_columns("RNAi")
    
    df$Expression_zscore <- if (length(expr_cols) > 0) {
      rowMeans(scale(df[, expr_cols, drop = FALSE]), na.rm = TRUE)
    } else { NA }
    
    df$CopyNumber_zscore <- if (length(cnv_cols) > 0) {
      rowMeans(scale(df[, cnv_cols, drop = FALSE]), na.rm = TRUE)
    } else { NA }
    
    df$CRISPR_zscore <- if (length(crispr_cols) > 0) {
      rowMeans(scale(df[, crispr_cols, drop = FALSE]), na.rm = TRUE)
    } else { NA }
    
    df$RNAi_zscore <- if (length(rnai_cols) > 0) {
      rowMeans(scale(df[, rnai_cols, drop = FALSE]), na.rm = TRUE)
    } else { NA }}

  long_data <- df |>
    select(all_of(group_var),
           Expression_zscore, CopyNumber_zscore, CRISPR_zscore, RNAi_zscore) |>
    pivot_longer(
      cols = c(Expression_zscore, CopyNumber_zscore, CRISPR_zscore, RNAi_zscore),
      names_to = "Parameter",
      values_to = "Value"
    ) |>
    drop_na(Parameter, !!sym(group_var), Value) |>
    mutate(
      Parameter = recode(Parameter,
                         "Expression_zscore" = "Expression",
                         "CopyNumber_zscore" = "Copy Number",
                         "CRISPR_zscore"     = "Gene Effect (CRISPR)",
                         "RNAi_zscore"       = "RNAi"
      ),
      group_var = .data[[group_var]]
    )
  
  return(long_data)}

plot_zscore_summary <- function(long_data,
                                group_label = "group_var",
                                show_p = TRUE) {
  
  ggplot(long_data, aes(x = group_var, y = Value, fill = group_var)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3) +
    facet_wrap(~Parameter, scales = "free_y", ncol = 2) +
    labs(
      title = "Cumulative Z-score Comparison",
      x = "Group",
      y = "Z-score",
      fill = "Group"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 12),
      legend.position = "bottom",
      plot.title = element_text(size = 14, hjust = 0.5)
    ) +
    if (show_p) {
      ggpubr::stat_compare_means(
        method = "wilcox.test",
        label = "p.signif"
      )
    } else {
      NULL}}

data_z <- all_data_pan %>% filter(Organ != 'normal')
long_data_all <- prepare_zscore_data(data_z, group_var = "Metastatic_MetMap")
plot_zscore_summary(long_data_all)

# finding significant combination in pan-cancer dataset

all_genes <- c("ANXA2", "MYH9", "ECM1", "EZR", "FSCN1",
               "ACTN4", "PFN1", "TAGLN", "CFL1", "GSN", "ANXA6")

# only breast data

find_successful_combinations <- function(data_z,
                                         all_genes,
                                         k_min,
                                         k_max,
                                         group_var) {
  successful_combinations <- list()
  
  for (k in k_min:k_max) {
    gene_combos <- combn(all_genes, k, simplify = FALSE)
    
    for (genes in gene_combos) {
      long_data <- prepare_zscore_data(
        df = data_z,
        selected_genes = genes,
        group_var = group_var)
      
      if (nrow(long_data) == 0 || length(unique(long_data[[group_var]])) < 2) next
      test_formula <- reformulate(group_var, response = "Value")
      result <- long_data %>% 
        group_by(Parameter) %>% 
        wilcox_test(formula = test_formula) %>% 
        ungroup()
      
      if (any(result$p < 0.05, na.rm = TRUE)) {
        combo_name <- paste(genes, collapse = ", ")
        successful_combinations[[combo_name]] <- result}}}
  
  param_hits_table <- successful_combinations %>% 
    imap_dfr(function(result_tbl, combo_name) {
      result_tbl %>% 
        filter(p < 0.05) %>% 
        distinct(Parameter) %>% 
        mutate(combo = combo_name)}) %>% 
    count(Parameter, name = "n_combos") %>% 
    arrange(desc(n_combos))
  list(
    combinations = successful_combinations,
    parameter_summary = param_hits_table)}

combos_breast <- find_successful_combinations(
  data_z = data_z %>% filter(Organ == 'Breast'),
  all_genes = all_genes,
  k_min = 2,
  k_max = 11,
  group_var = "PrimaryOrMetastasis")

combos_breast_metmap <- find_successful_combinations(
  data_z = data_z %>% filter(Organ == 'Breast'),
  all_genes = all_genes,
  k_min = 2,
  k_max = 11,
  group_var = "Metastatic_MetMap")

combos_breast$parameter_summary
combos_breast$combinations %>% length()
combos_breast$combinations %>% View()

# Pancreas
combos_pancreas <- find_successful_combinations(
  data_z = data_z %>% filter(Organ == 'Pancreas'),
  all_genes = all_genes,
  k_min = 2,
  k_max = 11,
  group_var = "PrimaryOrMetastasis")

combos_pancreas_metmap <- find_successful_combinations(
  data_z = data_z %>% filter(Organ == 'Pancreas'),
  all_genes = all_genes,
  k_min = 2,
  k_max = 11,
  group_var = "Metastatic_MetMap")

combos_pancreas$parameter_summary
combos_pancreas$combinations %>% length()
combos_pancreas$combinations %>% View()


# cns/brain
combos_cns <- find_successful_combinations(
  data_z = data_z %>% filter(Organ == 'CNS/Brain'),
  all_genes = all_genes,
  k_min = 2,
  k_max = 11,
  group_var = "PrimaryOrMetastasis")

combos_cns_metmap <- find_successful_combinations(
  data_z = data_z %>% filter(Organ == 'CNS/Brain'),
  all_genes = all_genes,
  k_min = 2,
  k_max = 11,
  group_var = "Metastatic_MetMap")

combos_cns$parameter_summary
combos_cns$combinations %>% length()
combos_cns$combinations %>% View()

# for Pan-cancer Tumor Type classification we'll rename CNS/Brain levels. Although it is not exactly the same comparison, for this analysis it is appropriate by sense 

data_z_ch <- data_z %>% 
  mutate(
    PrimaryOrMetastasis = case_when(
      PrimaryOrMetastasis == "Highly aggressive" ~ "Metastatic",
      PrimaryOrMetastasis == "Lowly aggressive"  ~ "Primary",
      TRUE ~ PrimaryOrMetastasis))

combos_all_changed_tumtype <- find_successful_combinations(
  data_z = data_z_ch,
  all_genes = all_genes,
  k_min = 2,
  k_max = 11,
  group_var = "PrimaryOrMetastasis")

combos_all_changed_metmap <- find_successful_combinations(
  data_z = data_z_ch,
  all_genes = all_genes,
  k_min = 2,
  k_max = 11,
  group_var = "Metastatic_MetMap")

combos_all_changed_tumtype$parameter_summary
combos_all_changed_tumtype$combinations %>% length()
combos_all_changed_tumtype$combinations %>% View()

# finding intersection in lists of combinations and parameters across organs and pan-cancer dataset

TumType_brall <- intersect(BrTumType_all_param, AllTumType_all_param)
TumType_panall <- intersect(PaTumType_all_param, AllTumType_all_param)
TumType_cnsall <- intersect(CNTumType_all_param, AllTumType_all_param)
Tum_Br_pan_all <- intersect(TumType_brall, TumType_panall)
intersect(Tum_Br_pan_all, TumType_cnsall)
TumType_br_pan <-  intersect(BrTumType_all_param, PaTumType_all_param)
intersect(TumType_br_pan, CNTumType_all_param)
intersect(TumType_brall, CNTumType_all_param)
intersect(TumType_panall, CNTumType_all_param)
intersect(Tum_Br_pan_all, CNTumType_all_param)

MetMap_brall <- intersect(BrMetMap_all_param, AllMetMap_all_param)
MetMap_panall <- intersect(PaMetMap_all_param, AllMetMap_all_param)
MetMap_cnsall <- intersect(CNMetMap_all_param, AllMetMap_all_param)
Met_Br_pan_all <- intersect(MetMap_brall, MetMap_panall)
intersect(Met_Br_pan_all, MetMap_cnsall)
MetMap_br_pan <-  intersect(BrMetMap_all_param, PaMetMap_all_param)
intersect(MetMap_br_pan, CNMetMap_all_param)
intersect(MetMap_brall, CNMetMap_all_param)
intersect(MetMap_panall, CNMetMap_all_param)

# checking significant combos for significant combination 
combos_breast$combinations$`MYH9, ECM1, EZR, FSCN1, ACTN4, PFN1, TAGLN, GSN`
combos_pancreas$combinations$`MYH9, ECM1, EZR, FSCN1, PFN1, TAGLN, GSN`
combos_cns$combinations$`MYH9, ECM1, EZR, FSCN1, ACTN4, PFN1, TAGLN, GSN`
combos_all_changed_tumtype$combinations$`MYH9, ECM1, EZR, FSCN1, ACTN4, PFN1, TAGLN, GSN`

combos_breast_metmap$combinations$`EZR, FSCN1, ACTN4, PFN1`
combos_pancreas_metmap$combinations$`MYH9, EZR, FSCN1, ACTN4, PFN1, TAGLN, CFL1`
combos_cns_metmap$combinations$`MYH9, EZR, FSCN1, ACTN4, PFN1, TAGLN, CFL1`
combos_all_changed_metmap$combinations$`EZR, ACTN4, PFN1, TAGLN, CFL1, ANXA6`


# MetMap level - expression was the best
## for pan-cancer dataset 
gene_sets <- list(
  c("EZR", "FSCN1"),
  c("ANXA2", "EZR", "FSCN1"),
  c("MYH9", "EZR", "FSCN1"),
  c("ECM1", "EZR", "FSCN1"),
  c("EZR", "FSCN1", "PFN1"),
  c("EZR", "FSCN1", "CFL1"))

prepare_cumulative_zscore <- function(df, gene_set, group_name = NULL) {
  expr_cols <- paste0("Expression ", gene_set)
  expr_cols <- intersect(expr_cols, colnames(df))
  if (length(expr_cols) == 0) return(NULL)
  scaled_expr <- df[, expr_cols, drop = FALSE] %>% 
    scale(center = TRUE, scale = TRUE)
  zscore <- rowMeans(scaled_expr, na.rm = TRUE)
  tibble(
    zscore = zscore,
    group = df$Metastatic_MetMap,
    gene_set = paste(gene_set, collapse = ", "))}

zscore_non_normal <- map_dfr(gene_sets, \(genes) {
  prepare_cumulative_zscore(data_z, genes)})

zscore_data <- zscore_non_normal %>% 
  filter(group %in% c("Low confidence", "High confidence")) %>% 
  mutate(group = factor(group, levels = c("Low confidence", "High confidence")))

p <- ggplot(zscore_data, aes(x = group, y = zscore, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  facet_wrap(~gene_set, scales = "free_y") +
  ggpubr::stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("High confidence", "Low confidence")),
    label = "p.signif",
    vjust = 0.5
  ) +
  labs(
    title = "Expression cumulative z-score comparison for MetMap level classification",
    x = "",
    y = "Expression cumulative z-score"
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.position = "none") +
   scale_fill_manual(
    values = c(
      "High confidence" = "#F8766D",  
      "Low confidence"  = "#00BFC4" ))

ggdraw(p) +
  draw_plot_label(label = "A", x = 0.01, y = 0.99, hjust = 0, vjust = 1, size = 20, fontface = "bold")


# Supplementary Figure 2 - only Breast/ Pancreas / CNS/Brain data for MetMap Expression 

zscore_non_normal <- map_dfr(gene_sets, \(genes) {
  prepare_cumulative_zscore(data_z %>% filter(Organ == 'CNS/Brain'), genes)})

zscore_data <- zscore_non_normal %>% 
  filter(group %in% c("Low confidence", "High confidence")) %>% 
  mutate(group = factor(group, levels = c("Low confidence", "High confidence")))

ggplot(zscore_data, aes(x = group, y = zscore, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  facet_wrap(~gene_set, scales = "free_y") +
  ggpubr::stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("High confidence", "Low confidence")),
    label = "p.signif",
    vjust = 0.5
  ) +
  labs(
    title = "Expression cumulative z-score comparison for MetMap level classification for CNS/Brain data",
    x = "",
    y = "Expression cumulative z-score"
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.position = "none") +
  scale_fill_manual(
    values = c(
      "High confidence" = "#F8766D",  
      "Low confidence"  = "#00BFC4" ))

# Tumor Type - Copy Number was the best
## for pan-cancer dataset 
gene_sets <- list(
  c("MYH9", "PFN1", "TAGLN"),
  c("MYH9", "PFN1", "CFL1"),
  c("MYH9", "FSCN1", "PFN1", "TAGLN"))

prepare_cumulative_zscore <- function(df, gene_set, group_name = NULL) {
  copy_cols <- paste0("Copy Number ", gene_set)
  copy_cols <- intersect(copy_cols, colnames(df))
  if (length(copy_cols) == 0) return(NULL)
  scaled_copy <- df[, copy_cols, drop = FALSE] %>% 
    scale(center = TRUE, scale = TRUE)
  zscore <- rowMeans(scaled_copy, na.rm = TRUE)
  tibble(
    zscore = zscore,
    group = df$PrimaryOrMetastasis,
    gene_set = paste(gene_set, collapse = ", "))}

zscore_non_normal <- map_dfr(gene_sets, \(genes) {
  prepare_cumulative_zscore(data_z %>% mutate(
    PrimaryOrMetastasis = case_when(
      PrimaryOrMetastasis == "Highly aggressive" ~ "Metastatic",
      PrimaryOrMetastasis == "Lowly aggressive"  ~ "Primary",
      TRUE ~ PrimaryOrMetastasis)), genes)})

zscore_data <-  zscore_non_normal %>% 
  filter(group %in% c("Primary", "Metastatic")) %>% 
  mutate(group = factor(group, levels = c("Primary", "Metastatic")))

p1 <- ggplot(zscore_data, aes(x = group, y = zscore, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  facet_wrap(~gene_set, scales = "free_y") +
  ggpubr::stat_compare_means(
    method = "wilcox.test",
    comparisons = list(
      c("Metastatic", "Primary")),
    label = "p.signif") +
  labs(
    title = "Copy Number cumulative z-score comparison for Tumor Type classification",
    x = "",
    y = "Copy Number cumulative z-score") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size=12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.position = "none") +
  scale_fill_manual(
    values = c(
      "Metastatic" = "#F8766D",  
      "Primary"  = "#00BFC4" ))

ggdraw(p1) +
  draw_plot_label(label = "B", x = 0.01, y = 0.99, hjust = 0, vjust = 1, size = 20, fontface = "bold")

# Supplementary Figure 2 - only Breast/ Pancreas / CNS/Brain data for Tumor Type Copy Number 

zscore_non_normal <- map_dfr(gene_sets, \(genes) {
  prepare_cumulative_zscore(data_z %>%  filter(Organ == 'Breast'), genes)}) # or Pancreas

# or CNS/Brain: 
zscore_non_normal <- map_dfr(gene_sets, \(genes) {
  prepare_cumulative_zscore(data_z  %>%  filter(Organ == 'CNS/Brain') %>% mutate(
    PrimaryOrMetastasis = case_when(
      PrimaryOrMetastasis == "Highly aggressive" ~ "Metastatic",
      PrimaryOrMetastasis == "Lowly aggressive"  ~ "Primary")), genes)})


zscore_data <-  zscore_non_normal %>% 
  filter(group %in% c("Primary", "Metastatic")) %>% 
  mutate(group = factor(group, levels = c("Primary", "Metastatic")))

ggplot(zscore_data, aes(x = group, y = zscore, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  facet_wrap(~gene_set, scales = "free_y") +
  ggpubr::stat_compare_means(
    method = "wilcox.test",
    comparisons = list(
      c("Metastatic", "Primary")),
    label = "p.signif") +
  labs(
    title = "Copy Number cumulative z-score comparison for Tumor Type classification for CNS/Brain data",
    x = "",
    y = "Copy Number cumulative z-score") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size=12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.position = "none") +
  scale_fill_manual(
    values = c(
      "Metastatic" = "#F8766D",  
      "Primary"  = "#00BFC4" ))
