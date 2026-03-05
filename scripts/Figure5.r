
# Load libraries

library(Seurat) 
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggalluvial)
library(scales) 
library(ggpubr)
library(circlize)
library(ComplexHeatmap)
library(clusterProfiler)
library(org.Hs.eg.db) 
library(stringr)
library(ggnewscale)
library(CytoTRACE2)
library(slingshot)
library(SingleCellExperiment)
library(monocle)
library(cowplot)
library(tidyr)


beta <- readRDS("./sub_beta.rds")


# Figure 5A (The left panel) ----------

mycol <- c('#a2d792','#a7cee3','#949494','#f29900')
DimPlot(beta, reduction = "umap.harmony", label = T)+
  theme_classic()+
  labs(title = "",x = 'UMAP1', y = 'UMAP2')+
  scale_color_manual(values = mycol)+
  theme(panel.grid = element_blank(),axis.text = element_text(size = 12, color = "black"))


# Figure 5A (The right panel) ----------

DimPlot(beta, reduction = "umap.harmony", label = T,split.by = 'subgroup')+
  theme_classic()+
  labs(title = "",x = 'UMAP1', y = 'UMAP2')+
  scale_color_manual(values = mycol)+
  theme(panel.grid = element_blank(),axis.text = element_text(size = 12, color = "black"))


# Figure 5B ----------
# markers <- FindAllMarkers(beta, only.pos = TRUE, min.pct = 0.20)

sig_markers <- markers %>% 
  filter(p_val_adj < 0.001, avg_log2FC > 0.5)    

aver_exp <- AverageExpression(beta,  layer = 'data' )
aver_exp <- as.data.frame(aver_exp$RNA)
aver_exp <- t(scale(t(aver_exp)))

# Filter for markers with high specificity (pct.1 > 0.6) and take the top 5 per cluster
top_markers <- sig_markers %>%
  group_by(cluster) %>%
  filter(pct.1 > 0.6) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 5)

aver_exp2 <- aver_exp[top_markers$gene,]
mycol <- colorRamp2(c(-1, 0, 1), c("#fffdfe", "#eb9499", '#9d0506'))
Heatmap(aver_exp2,
        cluster_columns = F,
        cluster_rows = F,
        show_row_names = T,
        col = mycol,
        row_names_gp = gpar(fontsize = 10),
        row_names_side = "left",
        column_names_rot = 0,
        column_names_gp = gpar(fontsize = 10) ,
        column_names_side = c('top'),
        name = 'expression',
        cell_fun = function(j, i, x, y, width, height, fill) {grid.rect(x, y, width, height, gp = gpar(fill = fill, col = "white", lwd = 0.7))})


# Figure 5C ----------

# Batch GO enrichment calculation
go_list <- list()
for (i in 0:3) {
  gene <- sig_markers[which(sig_markers$cluster == i),'gene']
  gene <- AnnotationDbi::select(org.Hs.eg.db, keys = gene,column = "ENTREZID", keytype="SYMBOL")
  go_list[[i+1]] <- enrichGO(gene = gene$ENTREZID,OrgDb = org.Hs.eg.db,ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = T)
  go_list[[i+1]] <- as.data.frame(go_list[[i+1]])
  go_list[[i+1]] <- go_list[[i+1]][go_list[[i+1]]$pvalue<0.05,]
}
names(go_list) <- c("cluster 0","cluster 1","cluster 2", "cluster 3")

plot_cluster_go <- function(cluster_name, category_map, go_list) {
  res_df <- go_list[[cluster_name]]
  plot_df <- data.frame()
  for (cat_name in names(category_map)) {
    terms <- category_map[[cat_name]]
    tmp <- res_df %>% 
      filter(tolower(Description) %in% tolower(terms)) %>% 
      mutate(Category = cat_name)
    plot_df <- rbind(plot_df, tmp)
  }
  plot_df <- plot_df %>%
    mutate(
      GeneRatio_val = as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio)),
      logP = -log10(pvalue)
    ) %>%
    arrange(desc(Category), pvalue) %>% 
    mutate(Description = factor(Description, levels = unique(Description)))
  
  p <- ggplot(plot_df, aes(x = GeneRatio_val, y = Description)) +
    geom_rect(aes(ymin = as.numeric(Description) - 0.45, 
                  ymax = as.numeric(Description) + 0.45, 
                  xmin = -Inf, xmax = 0, fill = Category), alpha = 0.3) +
    geom_segment(aes(x = 0, xend = GeneRatio_val, y = Description, yend = Description, color = Category), size = 1) +
    scale_fill_manual(values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99")) +
    scale_color_manual(values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99")) +
    new_scale_color() +
    geom_point(aes(size = Count, color = logP)) +
    scale_color_gradientn(colors = colorRampPalette(c("#3daeb7", "#eeeeee",'#ff615d', "#ff5743"))(100), name = "-Log10(pvalue)") +
    theme_bw() +
    labs(x = "GeneRatio", y = "", title = paste("GO Enrichment:", cluster_name)) +
    theme(
      axis.text.y = element_text(size = 9, color = "black"),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 8)
    )
  return(p)
}

# Manually curated map for meaningful biological interpretation
category_map_c1 <- list(
  "Energy metabolism" = c("oxidative phosphorylation", "aerobic respiration", 
                          "ATP biosynthetic process", "electron transport chain"),
  "Protein processing" = c("glycosylation","positive regulation of protein processing", 
                           "glycoprotein biosynthetic process", "glycoprotein metabolic process",
                           "protein processing"),
  "Membrane transportation" = c("protein transmembrane transport", "protein targeting to membrane", 
                                "membrane raft organization"),
  "Metal ion homeostasis" = c("detoxification","stress response to copper ion","detoxification of copper ion"),
  "Glucose metabolism" = c("glucose metabolic process", "hexose metabolic process", "gluconeogenesis", 
                           "glucose 6-phosphate metabolic process")
)

# Plotting for Cluster 1
plot_cluster_go("cluster 1", category_map_c1, go_list)


# Figure 5D ----------

features <- list(
  INS_response = c('SRSF4','AKT3','FOXO1'),
  Maturity = c('INS','ERO1B','PCSK1N','FXYD2','G6PC2'),
  ATP = c('NDUFA3','DNAJC15','NDUFA11'),
  Ribosomal = c('RPL12','RPL21','RPL3'),
  Immaturity = c('CHGA','CHGB','CLU','RBP4'))

DotPlot(beta, features = features, dot.scale = 5.6, col.min = -1, col.max = 1.0)+ 
  scale_color_gradient2(low = "#3c1c49", high = "#fde427", mid = "#1da686")+
  labs(y=" ", x="")+
  RotatedAxis()+
  theme(text = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.position = "top")

# Figure 5E ----------

mature_genes <- c('CD63','GAPDH','UBB','CALM1','CFAP126','GJD2','SLC2A2','PDX1','NKX6-1',
                 'UCN3','PCSK1','MAFA','INS','PFKFB2','G6PC2','MDH1','NEUROD1','RFX6',
                 'CHL1','PFKM','GPD2','GCK','SIX2','SIX3','CREB1','DNMT3A',
                 'SYTL4','ERO1B','GLRB','RXRA','PTF1A','ESRRG', 'BCL6', 'HNF4G','ATF5','PPARG',
                 'FAM159B','MNX1','HOPX','FXYD2','MT1X','SAT1','BMP5','NEFM')
immature_genes <- c('CD81', 'TMEM27','CLTRN','PYY','RPL3','EEF1A1' ,
                   'MAFB','CHGB','RBP4','ALDH1A3','MKI67','CHGA',
                   'PCSK2','PAM','IAPP','CLU',
                   'POU3F4', 'NEUROG3','GATA4', 'SOX4', 'NR4A1', 'EGR1', 'KLF11', 'NR1D1', 'CEBPB', 'HES1',
                   'DBP', 'RELB', 'ELF3', 'KLF10', 'PLAGL1','ATF3', 'PLAG1','MAFF','MYC','SOX9','SMAD3','ID1',
                   'PRDM16','FEV')
ER_genes <- c('PDIA3','DNAJB9','ERN1','EIF2AK3','ATF6','XBP1','EIF2A','ATF4','DDIT3',
             'HSPA5','SEL1L','SYVN1','DERL3','DNAJB1','DNAJB2','DNAJB11','HERPUD1','SDF2L1',
             'JUND','CREB3L2')

# Calculate module scores
# AddModuleScore calculates the average expression of each gene set relative to control genes
beta <- AddModuleScore(beta, features = list(mature_genes), name = "Mature_score")
beta <- AddModuleScore(beta, features = list(immature_genes), name = "Immature_score")
beta <- AddModuleScore(beta, features = list(ER_genes), name = "ER_score")

# Extract UMAP coordinates and merge with scores
umap_data <- as.data.frame(Embeddings(beta, "umap.harmony")) 

# AddModuleScore appends "1" to the specified name
umap_data$maturity_score <- beta@meta.data$Mature_score1

# Standardize scores (Z-score normalization) for cross-sample comparison
umap_data$maturity_score_standardized <- (umap_data$maturity_score - mean(umap_data$maturity_score)) / sd(umap_data$maturity_score)

# Visualization
ggplot(umap_data, aes(x = umapharmony_1, y = umapharmony_2, color = maturity_score_standardized)) +
  geom_point(size = 0.4, alpha = 0.4) +
  scale_color_gradientn(colours = colorRampPalette(c('#4e82bb','#d7f0f7','#fefdbd','#f97b4b','#e53038'))(100), 
                        limits = c(0, 1), oob = scales::squish) +
# Add 2D density contours to highlight cell clusters
  geom_density_2d(color = 'black', size = 0.5) +
  theme_minimal()+
  labs( x = "", y = "", color = "Score") +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),        
    axis.ticks = element_blank(),        
    axis.title = element_blank()
  )


# Figure 5F ----------

plot_data <- beta@meta.data %>%
  dplyr::select(seurat_clusters, Mature_score1, Immature_score1, ER_score1) %>%
  pivot_longer(cols = starts_with(c("Mature", "Immature", "ER")), 
               names_to = "Module", values_to = "Score") %>%
  mutate(Module = gsub("_score1", "", Module)) %>%
  mutate(Module = factor(Module, levels = c("Mature", "Immature", "ER")))

ggplot(plot_data, aes(x = seurat_clusters, y = Score, fill = Module)) +
  geom_violin(trim = FALSE, width = 1,alpha = 0.90, color = "white", show.legend = FALSE) +
  geom_boxplot(width = 0.03, fill = "white", color = "white", outlier.shape = NA, fatten = NULL, coef = 0) + 
  stat_summary(fun = median, geom = "crossbar", width = 0.1, aes(color = Module), show.legend = FALSE) +
  facet_wrap(~Module, scales = "free_y", ncol = 3) + 
  scale_fill_manual(values = c("Mature" = "#f39800", "ER" = "#b6141b", "Immature" = "#4972a6")) +
  scale_color_manual(values = c("Mature" = "#f39800", "ER" = "#b6141b","Immature" = "#4972a6")) +
  theme_classic() +
  theme(
    strip.background = element_blank(), 
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    panel.grid = element_blank()) +
  labs(x = "Beta cell clusters", y = "Module score")


# Figure 5G ----------

plot_data <- beta@meta.data %>%
  dplyr::select(seurat_clusters, subgroup, CytoTRACE2_Score)

# Potency score comparison in beta clusters
vln_cols <- c("#49ae48", "#93d6e7", '#949494',"#f19800") 
ggviolin(plot_data, x="seurat_clusters", y="CytoTRACE2_Score", 
         color="black",add ='mean_sd',fill='seurat_clusters',alpha = 0.7, 
         add.params = list(color="black"),trim = TRUE) +
  scale_fill_manual(values = vln_cols) +
  NoLegend() + 
  labs(x = 'Clusters',y = 'CytoTRACE2 potency score') + 
  ylim(0, 0.3) 


# Figure 5H ----------

metadata <- as.data.frame(beta@meta.data)

# ND subgroup comparison
# Calculate the frequency of each Seurat cluster within clinical subgroups
ratio <- metadata %>% group_by(subgroup, seurat_clusters) %>%
  dplyr::tally() %>%
  dplyr::group_by(subgroup) %>%
  dplyr::mutate(Freq = n / sum(n)) %>%
  dplyr::ungroup()

# Stacked bar with alluvial flows
cols <- c('#a2d792', '#a7cee3', '#949494','#f29900')
ggplot(ratio, aes(x = subgroup, y = Freq, fill = seurat_clusters, stratum = seurat_clusters, alluvium = seurat_clusters)) +
  scale_fill_manual(values = cols)+
  scale_y_continuous(expand = c(0,0), labels = percent_format(), limits = c(0, 1.05)) + 
  geom_col(width = 0.7, color = NA) +
  geom_flow(width = 0.7, alpha = 0.22, knot.pos = 0.35, color = 'white', linewidth = 0.5) +
  geom_alluvium(width = 0.7, alpha = 1, knot.pos = 0.35, fill = NA, color = 'white', linewidth = 0.5) +
  theme_classic() + 
  theme(axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        legend.title = element_blank(),
        panel.grid = element_blank()) +
  labs(x = "", y = "Beta cells (%)")


# Figure 5I ----------

# Calculate the percentage of each cluster for every individual donor
ratio <- metadata %>% group_by(donor,subgroup, seurat_clusters) %>%
  dplyr::tally() %>%
  dplyr::group_by(donor) %>%
  dplyr::mutate(Freq = n / sum(n)) %>%
  dplyr::ungroup()

# Pairwise comparisons for group-level differences
comparisons <- list( c("ND","MARD"),c("ND","MOD"),c("ND","SIDD"),c("MARD","MOD"), c("MARD","SIDD"),c("MOD","SIDD"))

# Faceted barplots with error bars
cols  <- c("#66c2a4", "#fb8d61", '#8d9fca',"#e789c3","#a6d753") 
ggplot(ratio, aes(x = subgroup, y = Freq, fill = subgroup)) + 
  geom_bar(stat = "summary", fun = mean, width = 0.7,color = 'black',alpha = 0.9) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, size = 0.5) + 
  geom_jitter(width = 0.2, size = 2, pch = 20, color="black", alpha = 0.8) +
  stat_compare_means(comparisons = comparisons, method = "t.test", label = "p.format",method.args = list(exact = FALSE)) +
  facet_wrap(~seurat_clusters, scales = "free_y", ncol = 4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = cols) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text.x = element_text(color = "black", size = 10, angle = 45, hjust = 1),
    strip.placement = "outside",
    legend.position = "none") +
  labs(x = "", y = "Percentage")


# Figure 5J ----------

# Subgroup mature score comparison

mycol <- c("#939599", "#49ae48", "#9680bd", "#e06083", "#09adc4")
ggplot(metadata, aes(x = subgroup, y = Mature_score1, fill = subgroup)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.9, size = 0.4,width = 0.85, aes(color = subgroup)) +
  geom_boxplot(width = 0.1, fill = "black", color = "black", outlier.shape = NA) +      
  stat_summary(fun = "median", geom = "point", shape = 21, fill = "white", color = "white", size = 2, stroke = 0.8) +
  scale_fill_manual(values = mycol) +
  scale_color_manual(values = mycol) +
  theme_pubr() +
  theme(legend.position = "none") +
  labs(x = '', y = 'Mature score') +
  stat_compare_means(comparisons = list(c('ND', 'MARD'), c('ND', 'MOD'),c('ND', 'SIDD'),c('ND', 'SIRD')), 
                     label = "p.signif", 
                     method = "t.test",
                     step.increase = 0.1)


# Figure 5K ----------
# Pseudotime density distribution across subgroups
# Using Slingshot trajectory results to visualize subtype shifts

density_df <- data.frame(
  Pseudotime = beta$Pseudotime,
  Subgroup = beta$subgroup,
  Donor = beta$donor) %>% 
  filter(!is.na(Pseudotime)) 

density_df$Subgroup <- as.factor(density_df$Subgroup, 
                                 levels = c("ND", "MARD", "MOD", "SIDD", "SIRD"))

# Calculate median pseudotime per subgroup for vertical indicators
mu <- density_df %>%
  dplyr::group_by(Subgroup) %>%
  dplyr::summarise(grp.median = median(Pseudotime, na.rm = TRUE))

# Density plot
ggplot(density_df, aes(x = Pseudotime, fill = Subgroup)) +
  geom_density(alpha = 0.8, color = "white") + 
  geom_vline(data = mu, aes(xintercept = grp.median, color = Subgroup),linetype = "dashed", linewidth = 1) +
  scale_fill_manual(values = c("ND" = "#d5d5d5", "MARD" = '#9c94c6','MOD'='#fae9a2','SIDD' = '#e5a283','SIRD' = '#cfdee6')) +
  scale_color_manual(values = c("ND" = "#d5d5d5", "MARD" = '#9c94c6','MOD'='#fae9a2','SIDD' = '#e5a283','SIRD' = '#cfdee6')) +
  theme_classic() +
  labs(x = "Pseudotime",  y = "Cell density") +
  theme(legend.position = "right",
        axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold"))


# Figure 5L ----------
# Donor-level pseudotime statistical comparison
# Aggregate pseudotime scores per donor for statistical robustness

target_order <- c("ND", "MARD", "MOD", "SIDD", "SIRD")
donor_summary <- density_df %>%
  dplyr::mutate(Subgroup = factor(Subgroup, levels = target_order)) %>% 
  dplyr::group_by(Donor, Subgroup) %>%
  dplyr::summarise(
    Median_Pseudotime = median(Pseudotime, na.rm = TRUE),
    Cell_Count = n(), .groups = 'drop') %>%
  dplyr::ungroup()

# Calculate the global median of the ND group as a reference baseline
nd_median = donor_summary %>%
  filter(Subgroup == "ND") %>%
  summarise(m = median(Median_Pseudotime, na.rm = TRUE)) %>%
  pull(m)

# Visualization
colors <- c("black", "#fb8d61","#8d9fca", "red", "#a6d753")
ggplot(donor_summary, aes(x = Subgroup, y = Median_Pseudotime, fill = Subgroup)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black") +
  geom_jitter(width = 0.2, size = 3, shape = 21, color = "black") +
  scale_fill_manual(values = colors) +
  stat_compare_means(comparisons = list(c("ND", "MARD"), c("ND", "MOD"), c("ND", "SIDD")),
                     label = "p.format", method = "wilcox.test") + 
  geom_hline(yintercept = nd_median, linetype = "dashed", color = "red", size = 0.8) +
  theme_classic() +
  labs(x = '', y = 'Pseudotime') +
  theme(legend.position = "none")


# Figure 5M ----------
# Monocle2 trajectory analysis

# Define Monocle pipeline function

run_monocle_trajectory <- function(seurat_obj, group_name, colors) {
  message(">>> Running Monocle 2 for Subgroup: ", group_name)
  
  # 1. Create CellDataSet
  # Extract raw counts and metadata
  counts_matrix <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
  
  pd <- new('AnnotatedDataFrame', data = seurat_obj@meta.data)
  fData <- data.frame(gene_short_name = row.names(seurat_obj), row.names = row.names(seurat_obj))
  fd <- new('AnnotatedDataFrame', data = fData)
  
  cds <- newCellDataSet(as(counts_matrix, 'sparseMatrix'),
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())
  
  # 2. Preprocessing: Normalization and Dispersion estimation
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  
  # 3. Feature Selection (Ordering Genes)
  # Selecting genes with high dispersion relative to the mean expression
  cds <- detectGenes(cds, min_expr = 0.1)
  disp_table <- dispersionTable(cds)
  ordering_genes <- subset(disp_table, mean_expression >= 0.1 & 
                             dispersion_empirical >= 1 * dispersion_fit)$gene_id
  
  cds <- setOrderingFilter(cds, ordering_genes)
  
  # 4. Dimensionality Reduction & Cell Ordering
  cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
  cds <- orderCells(cds)
  
  # 5. Visualization
  p <- plot_cell_trajectory(cds, show_branch_points = TRUE, 
                            color_by = "seurat_clusters", cell_size = 1) + 
       scale_color_manual(values = colors) +
       theme_classic() +
       labs(title = paste("Trajectory:", group_name))
  
  return(list(cds = cds, plot = p))
}


cluster_cols <- c('#a2d792', '#a7cee3', '#949494', '#f29900')

# Run Analysis for ND (Normal/Control)
beta_nd <- subset(beta, subgroup == "ND")
res_nd <- run_monocle_trajectory(beta_nd, "ND", cluster_cols)

# Run Analysis for SIDD (Severe Insulin-Deficient Diabetes)
beta_sidd <- subset(beta, subgroup == "SIDD")
res_sidd <- run_monocle_trajectory(beta_sidd, "SIDD", cluster_cols)