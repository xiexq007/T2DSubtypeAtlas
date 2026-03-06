
# Load libraries

library(Seurat) 
library(dplyr)
library(ggplot2)
library(mascarade)
library(ggalluvial)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(mdp)
library(stringr)

# Figure 2A ----------

# HPAP Fluidigm C1
fluidigm <- readRDS("./fluidigm.rds")
filtered_fluidigm <- subset(fluidigm, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 25) 
filtered_fluidigm <- NormalizeData(filtered_fluidigm, normalization.method = "LogNormalize", scale.factor = 10000)
hormones <- c("INS","GCG")
RidgePlot(filtered_fluidigm,features = hormones,group.by = "method", ncol = 1)

# HPAP 10X-Chromium
chromium <- readRDS("./10xchromium.rds")
filtered_chromium <- subset(chromium, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 25) 
filtered_chromium <- NormalizeData(filtered_chromium, normalization.method = "LogNormalize", scale.factor = 10000)
RidgePlot(filtered_chromium,features = hormones,group.by = "method", ncol = 1)

# Figure 2B ----------

panc_integrated <- readRDS("./integrated_data.rds")

plot_df <- cbind(panc_integrated@reductions$umap.harmony@cell.embeddings, panc_integrated@meta.data)

# Calculate the distance from the center
maskTable <- generateMask(
  dims = plot_df[, c("umapharmony_1", "umapharmony_2")], 
  cluster = plot_df$cell_type, 
  minDensity = 1.0, 
  smoothSigma = 0.05
)

plot_df$cell_type <- factor(plot_df$cell_type, levels = c("alpha", "beta", "delta", "PP/gamma", "acinar", "stellate", "ductal", "endothelial", "macrophage", "mast"))
color_mapping <-  c("#D58583", "#EBAEA9", "#BCD5B4","#DBC9B3","#AED0DF","#F8BE88", "#EED0E0", "#CBE5DE", "#cdb6da","#b3b0ae")
ggplot(plot_df, aes(x = umapharmony_1, y = umapharmony_2)) +
  geom_point(aes(color = cell_type), size = 0.05) +  
  geom_path(data = maskTable, aes(group = group), 
            linewidth = 0.3, linetype = 2, color = "black") +
  coord_fixed(ratio = 1) +
  theme_classic() +
  labs(x = "UMAP1", y = "UMAP2")+
  scale_color_manual(values = color_mapping, guide = guide_legend(override.aes = list(size = 4)) ) +  
  theme_dr() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",  
    legend.title = element_text(face = "bold")  
  )

# Figure 2C ----------
# Example for INS & GCG

INS <- FeaturePlot(panc_integrated, features = "INS",reduction =  "umap.harmony") + 
  scale_color_gradient(low = "#c6c8c7",high = "#e2161c")+theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(color = "black", size = 1)) 

GCG <- FeaturePlot(panc_integrated, features = "GCG",reduction =  "umap.harmony") + 
  scale_color_gradient(low = "#c6c8c7",high = "#e2161c")+theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(color = "black", size = 1)) 

# Figure 2D-E ----------

donor_data <- data.frame(panc_integrated@meta.data)
cell_order <- c("alpha", "beta", "delta", "PP/gamma", "acinar", "stellate", "ductal", "endothelial", "macrophage", "mast")
donor_data$cell_type <- factor(donor_data$cell_type, levels = cell_order)

# Calculate the frequency of each cell type per donor
cell_ratio <- donor_data %>% group_by(donor, cell_type, subgroup) %>%
  summarise(count = n())%>%
  group_by(donor) %>%
  mutate(total = sum(count)) %>% 
  mutate(Freq = count/total)

comparisons <- list( c("ND","MARD"),c("ND","MOD"),c("ND","SIDD"),c("MARD","MOD"), c("MARD","SIDD"),c("MOD","SIDD"))
cols <- c('#49ae48', '#2965ae', '#805188','#fda402','#d22f42')
cell <- c("beta","alpha","delta","PP/gamma")
subdata <- cell_ratio[which(cell_ratio $cell_type %in% cell),]

ggplot(subdata, aes(x = subgroup, y = Freq, fill = subgroup)) + 
  geom_bar(stat = "summary", fun = mean, width = 0.8) + 
  geom_jitter(width = 0.2,size = 2,pch = 20,color="black")+
  stat_summary(geom = "errorbar", fun.data = 'mean_se', width = 0.2, size = 0.5) +  
  scale_fill_manual(values = cols) +
  stat_compare_means(comparisons = comparisons, method = 'wilcox.test', label = "p.signif") +
  scale_y_continuous(expand = c(0, 0), limits = c(NA, 1)) + 
  labs(x = "", y = "Percentage") +
  theme_classic() +
  facet_grid(~cell_type, switch = "x") +  
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank(),  
    axis.text.y = element_text(color = "black", size = 12),
    axis.title.y = element_text(size = 12),
    strip.placement = "outside"  
  )

# Figure 2F ----------

ratio <- donor_data %>% 
  group_by(cell_type, subgroup) %>%
  count() %>%
  group_by(cell_type) %>% 
  mutate(Freq = n/sum(n)) %>%
  mutate(subgroup = factor(subgroup, levels = c("ND","MARD", "MOD", "SIDD", "SIRD")))

# Stacked bar with flow alluvium
ggplot(ratio, aes(x = cell_type, y = Freq, fill = subgroup, stratum = subgroup, alluvium = subgroup)) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(expand = c(0,0),limits = c(0, 1.08)) + 
  geom_col(width = 0.7, color = NA, size = 0.5) + 
  geom_flow(width = 0.7, alpha = 0.22, knot.pos = 0.35, color = 'white', linewidth = 0.5) +
  geom_alluvium(width = 0.7, alpha = 1, knot.pos = 0.35,fill = NA, color = 'white', linewidth = 0.5) +
  theme_classic() + 
  labs(x = "", y = "Percentage") +
  theme(axis.text.x = element_text(color = "black", size = 12,angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(color = "black", size = 12))
  
# Figure 2G ----------

run_cor_heatmap <- function(seurat_obj, cell_type_target, genelist, method = "pearson") {
  
  subtype_colors = c("MARD" = "#66C2A5", "MOD" = "#8DA0CB", "SIDD" = "#FC8D62", "SIRD" = "#E78AC3")
  target_order = c("ND", "MARD", "MOD", "SIDD", "SIRD")
  
  if ("ALL" %in% cell_type_target) {
    message("Running correlation for ALL cells...")
    obj_to_proc = seurat_obj
  } else {
    message(paste0("Running correlation for cell type: ", cell_type_target))
    obj_to_proc = subset(seurat_obj, subset = cell_type %in% cell_type_target)
  }
  
  pb = AggregateExpression(obj_to_proc, 
                           group.by = c("subgroup", "donor"), 
                           assays = "RNA", 
                           slot = "data", 
                           return.seurat = FALSE)
  mat = as.matrix(pb$RNA)
  mat = mat[rownames(mat) %in% genelist, ]
  mat = mat[rowSums(mat) > 0, ]
  
  gene_sds = apply(mat, 1, sd)
  mat = mat[gene_sds > 0, ]
  
  mat_scaled = t(scale(t(mat)))
  mat_scaled = mat_scaled[complete.cases(mat_scaled), ]
  cor_res = cor(mat_scaled, method = method)
  
  current_subtypes = sub("_.*", "", colnames(cor_res))
  sample_order = order(factor(current_subtypes, levels = target_order))
  cor_res = cor_res[sample_order, sample_order]
  
  subtype_labels = sub("_.*", "", colnames(cor_res))
  col_anno = HeatmapAnnotation(Subtype = subtype_labels, col = list(Subtype = subtype_colors), show_annotation_name = FALSE)
  row_anno = rowAnnotation(Subtype = subtype_labels, col = list(Subtype = subtype_colors), show_annotation_name = FALSE)
  
  min_cor = min(cor_res, na.rm = TRUE)
  col_fun = colorRamp2(c(min_cor, (min_cor + 1)/2, 1), c("#5a99d2", "white", "red"))
  
  hp = Heatmap(cor_res, 
               name = method,
               top_annotation = col_anno,left_annotation = row_anno,
               cluster_rows = FALSE, cluster_columns = FALSE,
               show_row_names = TRUE, show_column_names = FALSE,
               col = col_fun,
               column_split = factor(subtype_labels, levels = target_order),
               row_split = factor(subtype_labels, levels = target_order))
  
  return(hp)
}

# Select the top 100 highly expressed genes for similarity calculation
gene_means = rowMeans(as.matrix(panc_integrated[["RNA"]]$data))
genes = names(sort(gene_means, decreasing = TRUE)[1:100])
genes = genes[!grepl("^RP[SL]", genes)]
panc_t2d = subset(panc_integrated, subset = subgroup != "ND")
run_cor_heatmap(panc_t2d, "ALL",genes,'pearson')

# Figure 2H ----------

# Calculate the molecular degree of perturbation
# Example for alpha

panc_integrated$disease_group <- paste(panc_integrated$disease_status,panc_integrated$subgroup,sep = "_")
panc_alpha <- subset(panc_integrated, cell_type == 'alpha')
Idents(panc_alpha) <- 'disease_group'

# Extract log-normalized data
alpha_sc <- as.data.frame(panc_alpha@assays$RNA$data)
colnames(alpha_sc) <- paste(panc_alpha$disease_group,c(1:dim(alpha_sc)[2]),sep = "")

# Create phenotype dataFrame directly from Seurat metadata
pheno_df <- data.frame(Sample = colnames(alpha_sc))
pheno_df$Class <- ifelse(grepl("MOD", pheno_df$Sample), "MOD", 
                        ifelse(grepl("ND", pheno_df$Sample), "ND",
                               ifelse(grepl("MARD", pheno_df$Sample), "MARD",
                                      ifelse(grepl("SIDD", pheno_df$Sample), "SIDD",
                                             ifelse(grepl("SIRD", pheno_df$Sample), "SIRD", NA)))))

# Run MDP analysis
mdp_sc_five <- mdp(
  data = alpha_sc, 
  pdata = pheno_df, 
  control_lab = "ND",
  directory = "./mdp", 
  file_name = "02_alpha_five_"
)


# Visualization
plot_mdp_five_groups <- function(file_path, y_limit) {
  
  # Load MDP results
  mdp_data <- read.table(file_path, header = TRUE)
  
  # Ensure factor levels for plotting
  mdp_data$perturbedgenes.Class <- factor(
    mdp_data$perturbedgenes.Class, 
    levels = c("ND", "MARD", "MOD", "SIDD", "SIRD")
  )
  
  # Custom color palette
  cols <- c('#49ae48', '#2965ae', '#805188', '#fda402', '#d22f42')
  
  # Generate Boxplot
  p <- ggplot(mdp_data, aes(x = perturbedgenes.Class, y = perturbedgenes.Score, fill = perturbedgenes.Class)) +
    geom_boxplot(na.rm = T,outlier.shape = NA) +
    scale_fill_manual(values = cols) +
    ylim(0, y_limit) +
    theme_classic() +
    labs(title = "", x = "", y = "Score") + 
    theme(legend.position = "none", 
          axis.text = element_text(size = 12, color = "black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  
  # Statistical Analysis (One-way ANOVA)
  # Image of one-way ANOVA assumptions and interpretation
  anova_test <- aov(perturbedgenes.Score ~ perturbedgenes.Class, data = mdp_data)
  print(summary(anova_test))
  
  return(p)
}

mdp_file <- "./02_alpha_five_sample_scores.tsv"
p_mdp <- plot_mdp_five_groups(mdp_file, y_limit = 0.5)

# Figure 2I ----------
# Summarize donor-level perturbation scores across pancreatic cell types

data_path <- "./mdp/"

alpha <- read.delim(paste0(data_path, "02_alpha_five_sample_scores.tsv"), sep = "\t", header = TRUE)
beta  <- read.delim(paste0(data_path, "01_beta_five_sample_scores.tsv"), sep = "\t", header = TRUE)
delta <- read.delim(paste0(data_path, "03_delta_five_sample_scores.tsv"), sep = "\t", header = TRUE)
gamma <- read.delim(paste0(data_path, "04_gamma_five_sample_scores.tsv"), sep = "\t", header = TRUE)

summarize_donor_mdp <- function(df, score_col = "perturbedgenes.Score", class_col = "perturbedgenes.Class") {
  df_processed = df %>%
    mutate(donor = str_extract(perturbedgenes.Sample, "HPAP\\d+"))
  donor_summary = df_processed %>%
    group_by(donor, !!sym(class_col)) %>%
    summarise(
      mean_perturbation_score = mean(!!sym(score_col), na.rm = TRUE),
      cell_count = n(), 
      .groups = 'drop') %>%
    rename(Subtype = !!sym(class_col))
  return(donor_summary)
}

# Data transformation
alpha_donor <- summarize_donor_mdp(alpha)
beta_donor  <- summarize_donor_mdp(beta)
delta_donor <- summarize_donor_mdp(delta)
gamma_donor <- summarize_donor_mdp(gamma)

all_donor <- bind_rows(
  mutate(beta_donor,  CellType = "Beta"),
  mutate(alpha_donor, CellType = "Alpha"),
  mutate(delta_donor, CellType = "Delta"),
  mutate(gamma_donor, CellType = "Gamma")
)

target_order <- c("ND", "MARD", "MOD", "SIDD", "SIRD")
all_donor$Subtype <- factor(all_donor$Subtype, levels = target_order)
donor_info <- all_donor %>%
  select(donor, Subtype) %>%
  distinct() %>%
  arrange(Subtype)

# Heatmap generation
subtype_colors <- c("ND" = "#66C2A5", "MARD" = "#FC8D62", "MOD" = "#8DA0CB", "SIDD" = "#E78AC3", "SIRD" = "#A6D854")
cell_types <- c("Alpha", "Beta", "Delta", "Gamma")
ht_list <- NULL 

for (ct in cell_types) {
  ct_data = all_donor %>%
    filter(CellType == ct) %>%
    select(donor, mean_perturbation_score) %>%
    right_join(donor_info, by = "donor") %>%
    arrange(Subtype)
  mat = as.matrix(ct_data$mean_perturbation_score)
  rownames(mat) = ct_data$donor
  colnames(mat) = ct

  min_val = min(mat, na.rm = TRUE)
  max_val = max(mat, na.rm = TRUE)
  mid_val = median(mat, na.rm = TRUE) 
  col_fun = colorRamp2(c(min_val, mid_val, max_val), c("#5a99d2", "white", "red"))
  
  if (is.null(ht_list)) {
    left_anno = rowAnnotation(
      Subtype = donor_info$Subtype,
      col = list(Subtype = subtype_colors),
      show_annotation_name = FALSE)
  } else {left_anno <- NULL}
  
  current_ht = Heatmap(mat,
                       name = paste0(ct, "_Score"),
                       column_title = ct,
                       left_annotation = left_anno,
                       col = col_fun,
                       cluster_rows = FALSE,
                       cluster_columns = FALSE,
                       show_row_names = (ct == "Gamma"), 
                       row_split = donor_info$Subtype,
                       row_gap = unit(1, "mm"),
                       border = TRUE,
                       width = unit(1.5, "cm")) 
  if (is.null(ht_list)) {
    ht_list <- current_ht} else {ht_list <- ht_list + current_ht}
}

draw(ht_list, ht_gap = unit(3, "mm"))
