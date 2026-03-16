
# Load libraries

library(Seurat) 
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggrepel)
library(patchwork)
library(circlize)
library(ComplexHeatmap)
library(pheatmap)
library(DESeq2)
library(ggpubr)
library(tibble)
library(GSVA)
library(fmsb)


# Figure 4A ----------

deg_subgroup_ovo <- readRDS("./deg_subgroup_ovo.rds")

# Define criteria for Significant DEGs: adj. P-value < 0.001 and |log2FC| > 1
filter_deg <- function(df) {
  result <- df[df$p_val_adj < 0.001 & abs(df$avg_log2FC) > 1, ]
  result$gene <- rownames(result)
  return(result)
}

sig_subgroup_ovo <- lapply(deg_subgroup_ovo, function(group) {
    filter_deg(group)
})

# Number statistics
cell_types <- c('alpha','beta','delta','PP/gamma','acinar','stellate','ductal','endothelial','macrophage','mast')
comparisons <- c("SIRD_SIDD", "SIRD_MOD","SIRD_MARD","SIRD_ND",  "SIDD_MOD","SIDD_MARD", "SIDD_ND", "MOD_MARD","MOD_ND","MARD_ND")

# Initialize an empty matrix for DEG counts
deg_num <- data.frame(matrix(0, nrow = length(cell_types), ncol = length(comparisons)))
rownames(deg_num) <- comparisons
colnames(deg_num) <- cell_types

for (name in names(sig_subgroup_ovo)) {
  parts <- strsplit(name, "_")[[1]]
  cell_type <- parts[1]  
  comparison <- paste(parts[c(2,4)], collapse = "_")  
  if (cell_type %in% cell_types && comparison %in% comparisons) {
    count <- nrow(sig_subgroup_ovo[[name]])
    deg_num[comparison,cell_type] <- count
  }}

# Reorder columns and rows for prioritized biological display
target_cell_order <- c('acinar','ductal','delta','beta','alpha','endothelial','stellate','macrophage','mast','PP/gamma')
target_comp_order <- c('SIRD_SIDD','SIRD_MARD','SIRD_MOD','SIRD_ND','SIDD_MOD','MOD_MARD','SIDD_ND','MOD_ND','SIDD_MARD','MARD_ND')

deg_num <- deg_num[target_comp_order, target_cell_order]

# Visualization: Heatmap
color_palette <- colorRampPalette(c("white", "#fb5e49", "#b82531"))(100)
pheatmap(as.matrix(deg_num),
         color = color_palette,
         display_numbers = T,
         number_format = "%.0f",
         cluster_rows = F,
         cluster_cols = F,
         number_color = "black",
         name = 'Number of DEGs')

# Visualization: Bar charts for marginal totals

# Comparison-wise totals
comparison_df <- data.frame(
  Comparison = rownames(deg_num), 
  count = rowSums(deg_num)
)

ggplot(comparison_df, aes(x = reorder(Comparison, count), y = count)) +
  geom_bar(stat = "identity", fill = "#ffb150", width = 0.8) +
  coord_flip() + # Flip for better readability of comparison names
  theme_bw() +
  labs(x = NULL, y = "Number of DEGs", title = "DEGs by Comparison") +
  theme(panel.grid = element_blank(), axis.text = element_text(color = "black"))

# Cell-type-wise totals
cell_type_df <- data.frame(
  CellType = colnames(deg_num), 
  count = colSums(deg_num)
)

ggplot(cell_type_df, aes(x = reorder(CellType, count), y = count)) +
  geom_bar(stat = "identity", fill = "#9ed8a6", width = 0.8) +
  coord_flip() + # Consistent flip orientation
  theme_bw() +
  labs(x = NULL, y = "Number of DEGs", title = "DEGs by Cell Type") +
  theme(panel.grid = element_blank(), axis.text = element_text(color = "black"))


# Figure 4B ----------
# Scatter plot

plot_scatter <- function(df,selected_gene,title,left_x,right_x,left_y,right_y){
  df <- df %>%
    mutate(Difference = pct.1 - pct.2) %>% 
    rownames_to_column("gene") %>%
    mutate(change = ifelse(p_val_adj < 0.001 & avg_log2FC > 1, "up",
                      ifelse(p_val_adj < 0.001 & avg_log2FC < -1, "down", "ns"))) 
  sub_df <- df[match(selected_gene,df$gene),]
  
  ggplot(df, aes(x = Difference, y = avg_log2FC, color = change)) + 
    geom_point(size = 2) +
    geom_point(data = df[df$change == "up", ], aes(x = Difference, y = avg_log2FC), color = "#dc7720", size = 2) +
    geom_point(data = df[df$change == "down", ], aes(x = Difference, y = avg_log2FC), color = "#0d1b46", size = 2) +
    scale_color_manual('change',labels=c(paste0("down(",table(df$change)[[1]],')'),
                                         'ns',
                                         paste0("up(",table(df$change)[[3]],')' )),
                       values=c("#0d1b46", "grey","#dc7720" ))+
    geom_text_repel(data = sub_df[which(sub_df$avg_log2FC > 0),],aes(label = gene),color = "black", size = 3,
                    segment.size = 0.3, segment.color = "black", max.overlaps = 200,direction = 'y',
                    hjust = "left",nudge_x = right_x, nudge_y = right_y) +
    geom_text_repel(data = sub_df[which(sub_df$avg_log2FC < 0),],aes(label = gene),color = "black", size = 3,
                    segment.size = 0.3, segment.color = "black", max.overlaps = 200,direction = 'y',
                    hjust = "left",nudge_x = left_x, nudge_y = left_y) +
    geom_vline(xintercept = 0,linetype = 2) +
    geom_hline(yintercept = 0,linetype = 2) +
    labs(x="Percentage difference",y = "Log-Fold Change",title = title)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 12),
          axis.text = element_text(color = "black", size = 12))
}

# SIDD vs ND: alpha
plot_scatter(deg_subgroup_ovo$alpha_SIDD_vs_ND,
            c("XIST","SST","FXYD2","TMEM196","INS","UQCR10","USP9Y","RPS4Y1","TTTY14","CLPS","PLA2G1B","CELA3A","CTRC","CELA3B","CPA1","REG1A","HSPA5","DNAJA4"),
             'SIDD-ND: alpha',-0.4,0.4,-0.6,0.65)

# SIDD vs ND: beta
plot_scatter(deg_subgroup_ovo$beta_SIDD_vs_ND,
             c("XIST","PPP1R1A","TMED6","IAPP","IFI27L2","DNAJC15","COA3","SST","IER3","FXYD2","FKBP11","CHGA","IGFBP7","PDK4","NDUFAF3","ALDH1A1","MPG",
               "DDX3Y","USP9Y","UPY","CEL","CELA3A","CELA3B","CPA1","REG1A","TSPAN5","TMEM259","KAT6A","SYT11","SNHG25","PCSK2","CLPS"),
             'SIDD-ND: beta',-0.75,0.5,-0.65,0.5)

# SIRD vs ND: alpha
plot_scatter(deg_subgroup_ovo$alpha_SIRD_vs_ND,
             c("DDIT3","AGPAT5","SRSF4","SLC9A8","ATXN1","SLC38A2","CDK8","HERPUD1","PRKAA2","FOXO3","KLF12","ATG7","MAFG","SESTD1","AKT3","ETV6",
               "CRYBA2","NDUFB1","FKBP11","MAFB","NEUROD1","FKBP2","FOS","FXYD3","SSTR2","NDUFA1"),
             'SIRD-ND: alpha',-0.6,0.4,-0.7,0.65)

# SIRD vs ND: beta
plot_scatter(deg_subgroup_ovo$beta_SIRD_vs_ND,
            c("SNHG5","IP6K1","SETD7","SLC22A5","CUL1","SRSF4","FOXO3","MAFG","HERPUD1","SOGA3","PPP2R2A","COA1","MT1X","MT2A","AQP3","CDKN1C","G6PC2","EGR1","CLPS",
               "CKB","NDUFB1",),
             'SIRD-ND: beta',-0.6,0.4,-0.7,0.65)


# Figure 4C (The left panel: boxplot) ----------

plot_donor_concordance <- function(seurat_obj, genes, celltype, target_group, mode, ref_group,num_columns = 4) {
  panc_sub <- subset(seurat_obj, subset = cell_type == celltype)
  pb_res <- AggregateExpression(panc_sub,features = genes,group.by = c("donor", "subgroup"), 
                               assays   = "RNA",slot = "data",return.seurat = FALSE)$RNA
  pb_df <- as.data.frame(pb_res) %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "combined", values_to = "expression") %>%
    separate(combined, into = c("donor", "subgroup"), sep = "_")

  # Baseline calculation
  if (mode == "ovr") {
    ref_baseline <- pb_df %>%
      filter(subgroup != target_group) %>%
      group_by(gene) %>%
      summarise(baseline_mean = median(expression), .groups = "drop")
    plot_data = pb_df %>%
      mutate(plot_group = ifelse(subgroup == target_group, target_group, "The Rest"))
    title_suffix <- paste(target_group, "vs The Rest (OVR)")
    group_levels <- c(target_group, "The Rest")
    
  } else {
    
    ref_baseline <- pb_df %>%
      filter(subgroup == ref_group) %>%
      group_by(gene) %>%
      summarise(baseline_mean = median(expression), .groups = "drop")
    plot_data <- pb_df %>% 
      filter(subgroup %in% c(target_group, ref_group)) %>%
      mutate(plot_group = subgroup)
    title_suffix = paste(target_group, "vs", ref_group, "(OVO)")
    group_levels = c(target_group, ref_group)
  }

  # Calculate Log2 fold change for each donor
  donor_logFC_df <- plot_data %>%
    left_join(ref_baseline, by = "gene") %>%
    mutate(logFC = log2((expression + 0.1) / (baseline_mean + 0.1)))
  
  donor_logFC_df$plot_group <- factor(donor_logFC_df$plot_group, levels = group_levels)
  donor_logFC_df$gene <- factor(donor_logFC_df$gene, levels = genes)
  subgroup_colors <- c("ND"="#66C2A5", "MARD"="#FC8D62", "MOD"="#8DA0CB", "SIDD"="#E78AC3", "SIRD"="#A6D854")

  # Visualization
  p <- ggplot(donor_logFC_df, aes(x = plot_group, y = logFC)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5, color = "black", fill = "gray") +
    geom_jitter(aes(fill = subgroup), width = 0.2, size = 3, shape = 21, color = "black", alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.8) +
    facet_wrap(~gene, scales = "free_y",ncol = num_columns) +
    theme_bw() +
    scale_fill_manual(values = subgroup_colors) +
    labs(
      title = paste("Donor-level Consistency:", celltype),
      subtitle = paste("Comparison:", title_suffix),
      y = "log2(Fold Change vs Baseline)",
      x = "Comparison Groups",
      fill = "Original Subtype"
    )
 
  print(p)
  return(donor_logFC_df)
}

# ovo -----------
# target group vs ND
panc_integrated <- readRDS("./integrated_data.rds")
panc_integrated$cell_group <- paste(panc_integrated$cell_type, panc_integrated$subgroup, sep = "_")

# SIDD-beta
sidd_ovo_beta  plot_donor_concordance(
  seurat_obj   = panc_integrated,
  genes        =  'CTRB1',
  celltype     = "beta",
  target_group = "SIDD",
  ref_group = "ND",
  num_columns = 1,
  mode = "ovo" )

# SIDD-alpha
sidd_ovo_alpha <- plot_donor_concordance(
  seurat_obj   = panc_integrated,
  genes        = "PNLIP",
  celltype     = "alpha",
  target_group = "SIDD",
  ref_group = "ND",
  num_columns = 1,
  mode = "ovo" )


# Figure 4C (The right panel: hetamap) ----------

# Getting the subset includes alpha, beta, delta, and gamma cells
panc_sub <- subset(panc_integrated, cell_type %in% c("alpha",'beta','delta','PP/gamma'))
aver_exp <- AverageExpression(panc_sub, group.by = 'cell_group', layer = 'data')
aver_exp <- as.data.frame(aver_exp$RNA)
aver_exp <- t(scale(t(aver_exp)))

plot_heatmap_ovo <- function(aver_exp,cell_group,genes){
  aver_exp2 <- aver_exp[genes,cell_group]
  mycol <- colorRamp2(c(-2, 0, 2), c("#0e7b7d", "white", "#d2565e"))
  ht <- Heatmap(aver_exp2,
                cluster_columns = F,
                cluster_rows = F,
                show_row_names = T,
                col = mycol,
                row_names_gp = gpar(fontsize = 9) ,
                column_names_rot = 60,
                column_names_gp = gpar(fontsize = 9),
                column_names_side = c('top'),
                name = 'expression',
                cell_fun = function(j, i, x, y, width, height, fill) {grid.rect(x, y, width, height, gp = gpar(fill = fill, col = "white", lwd = 0.5))}
               )
  draw(ht)
}

# SIDD-ND
selected_ovo <- sig_subgroup_ovo[c("alpha_SIDD_vs_ND","beta_SIDD_vs_ND","delta_SIDD_vs_ND","PP/gamma_SIDD_vs_ND")]  
plot_heatmap_ovo(aver_exp,
                 c("alpha-SIDD","beta-SIDD","delta-SIDD","PP/gamma-SIDD","alpha-ND","beta-ND","delta-ND","PP/gamma-ND"),
                 c("CELA3A","CELA3B","CPA1","PRSS1","CTRB1", # Pancreatic secretion
                   "MTRNR2L1","MTRNR2L8","MTRNR2L12", # Apoptosis related
                   "PNLIP","PLA2G1B","PLCG2","TRIB3", # lipid metabolic
                   "REG1A","HDAC9","CHD7","RFX3","NR4A1","SPINK1" # hormone metabolic realted
                 )
                )


# Figure 4D ----------
# Dot plot

genes <- c("NDUFB1","NDUFAF3","UQCC2","NDUFA1","NDUFC1","ATP5ME", # respiration
          "MT1E","MT1X","MT2A","MT1F","MTLN", # Ion homeostasis
          "SLC30A8","NEUROD1","RFX6","CAMK2G","G6PC2","ISL1", # Insulin secretion
          "FKBP11","HSPE1","FKBP2","ERO1B","HSPA8", # protein folding
          "CTRB2","CTRB1","PRSS1","CLPS","CPA1", # Pancreatic secretion
          "TRIB3","HSPA5","CLGN","TMEM259","PPP1R15A","CREBRF", # ER stree response
          "MTRNR2L1","MTRNR2L8","MTRNR2L12","BNIP3L","IFI27L2", # Apoptosis related
          "RPL38","RPS29","RPLP2","RPL31","RPS10" # cytoplasmic translation
)

group_order <- c('alpha_ND','beta_ND','delta_ND','PP/gamma_ND','alpha_MARD','beta_MARD','delta_MARD','PP/gamma_MARD',
                 'alpha_MOD','beta_MOD','delta_MOD','PP/gamma_MOD','alpha_SIDD','beta_SIDD','delta_SIDD','PP/gamma_SIDD',
                 'alpha_SIRD','beta_SIRD','delta_SIRD','PP/gamma_SIRD')
Idents(panc_sub) <- factor(Idents(panc_sub), levels = group_order)

DotPlot(panc_sub, 
        features = genes,
        dot.scale = 5.6,
        col.min = -1, 
        col.max = 1.0) + 
  scale_color_gradient2(low = "#3c1c49", high = "#fde427", mid = "#1da686") +
  labs(y=" ", x="") + 
  RotatedAxis()+
  theme(text = element_text(size = 8), axis.text = element_text(size = 8))


# Figure 4E ----------
# GSEA analysis

run_gsea <- function(data, type = "ovr", group = NULL, cell = NULL, contrast = NULL) {

  if (type == "ovr") {
    df <- data[[group]][[cell]]
  } else if (type == "ovo") {
    df <- data[[contrast]] } 
  
  gene_list <- df$avg_log2FC
  names(gene_list) = rownames(df)
  gene_list <- sort(gene_list, decreasing = TRUE)
 
  gse_res <- gseGO(geneList = gene_list, 
                  OrgDb = org.Hs.eg.db, 
                  keyType = "SYMBOL",
                  ont = "BP", 
                  pAdjustMethod = "BH", 
                  pvalueCutoff = 1, 
                  eps = 0,
                  seed = 123)
  return(gse_res)
}

# SIDD
ovo_sidd_beta_gsea <- run_gsea(deg_subgroup_ovo, type="ovo",  contrast="beta_SIDD_vs_ND")
ovo_sidd_alpha_gsea <- run_gsea(deg_subgroup_ovo, type="ovo", contrast="alpha_SIDD_vs_ND")
ovo_sidd_delta_gsea <- run_gsea(deg_subgroup_ovo, type="ovo", contrast="delta_SIDD_vs_ND")

# SIRD
ovo_sird_beta_gsea <- run_gsea(deg_subgroup_ovo, type="ovo",  contrast="beta_SIRD_vs_ND")
ovo_sird_alpha_gsea <- run_gsea(deg_subgroup_ovo, type="ovo", contrast="alpha_SIRD_vs_ND")
ovo_sird_delta_gsea <- run_gsea(deg_subgroup_ovo, type="ovo", contrast="delta_SIRD_vs_ND")

# MOD
ovo_mod_beta_gsea <- run_gsea(deg_subgroup_ovo, type="ovo",  contrast="beta_MOD_vs_ND")
ovo_mod_alpha_gsea <- run_gsea(deg_subgroup_ovo, type="ovo", contrast="alpha_MOD_vs_ND")
ovo_mod_delta_gsea <- run_gsea(deg_subgroup_ovo, type="ovo", contrast="delta_MOD_vs_ND")

target_terms <-  c( 'oxidative phosphorylation',
                    'proton motive force-driven ATP synthesis',
                    'ATP biosynthetic process',
                    'response to copper ion',
                    'response to insulin',
                    'response to endoplasmic reticulum stress',
                    'cellular response to unfolded protein',
                    'response to topologically incorrect protein',
                    'transcription by RNA polymerase II',
                    'positive regulation of RNA splicing',
                    'cellular response to chemical stress',
                    'DNA damage response',
                    'response to mechanical stimulus',
                    'response to hormone',
                    'cellular response to external stimulus'
                   )

gsea_list <- list(
  "SIDD_alpha" = ovo_sidd_alpha_gsea,
  "SIDD_beta"  = ovo_sidd_beta_gsea,
  "SIDD_delta" = ovo_sidd_delta_gsea,
  "SIRD_alpha" = ovo_sird_alpha_gsea,
  "SIRD_beta"  = ovo_sird_beta_gsea,
  "SIRD_delta" = ovo_sird_delta_gsea,
  "MOD_alpha" = ovo_mod_alpha_gsea,
  "MOD_beta"  = ovo_mod_beta_gsea,
  "MOD_delta"  = ovo_mod_delta_gsea
)

# Extract NES, P-value, and FDR matrices
nes_mat <- matrix(0, nrow = length(target_terms), ncol = length(gsea_list))
p_mat   <- matrix(1, nrow = length(target_terms), ncol = length(gsea_list)) 
fdr_mat <- matrix(1, nrow = length(target_terms), ncol = length(gsea_list)) 

rownames(nes_mat) <- target_terms; colnames(nes_mat) = names(gsea_list)
rownames(p_mat)   <- target_terms; colnames(p_mat)   = names(gsea_list)
rownames(fdr_mat) <- target_terms; colnames(fdr_mat) = names(gsea_list)

for (name in names(gsea_list)) {
  
  res_df <- as.data.frame(gsea_list[[name]])
  matched_idx <- match(target_terms, res_df$Description)
  
  nes_val <- res_df$NES[matched_idx]
  nes_mat[!is.na(nes_val), name] <- nes_val[!is.na(nes_val)]
  
  p_val <- res_df$pvalue[matched_idx]
  p_mat[!is.na(p_val), name] <- p_val[!is.na(p_val)]
  
  fdr_val <- res_df$p.adjust[matched_idx]
  fdr_mat[!is.na(fdr_val), name] <- fdr_val[!is.na(fdr_val)]
}

# Heatmap visualization
col_fun <- colorRamp2(c(-4, 0, 4), c("#5a99d2", "white", "red"))
Heatmap(nes_mat, 
        name = "NES", 
        col = col_fun,
        cluster_rows = FALSE,       
        cluster_columns = FALSE,   
        rect_gp = gpar(col = "black", lwd = 0.5), 
        cell_fun = function(j, i, x, y, width, height, fill) {
          if (fdr_mat[i, j] < 0.05) {
            # Asterisk for FDR < 0.05
            grid.text("*", x, y, gp = gpar(fontsize = 15, fontface = "bold"))} 
          else if (p_mat[i, j] < 0.05) {
            # Plus sign for nominal p-value < 0.05
            grid.text("+", x, y, gp = gpar(fontsize = 12))}},
        column_split = factor(rep(c("SIDD", "SIRD","MOD"), each = 3), levels = c("SIDD", "SIRD","MOD")),
        row_names_side = "left",
        column_names_rot = 45,
        row_names_gp = gpar(fontsize = 10)
       )


# Figure 4F ----------

panc_beta <- subset(panc_integrated, cell_type == 'beta')
panc_bulk <- AverageExpression(panc_beta, group.by = 'cell_group', layer = 'data')
panc_bulk <- as.matrix(panc_bulk$RNA)

# GSVA analysis
gsva_results <- gsva(expr = panc_bulk, 
                     gset.idx.list = geneset_list, # genelist for enrichment
                     method = "gsva",
                     kcdf = "Gaussian" )
# Radar plot
radar_data <- as.data.frame(t(gsva_results))

# Reorder rows to match clinical progression: ND -> MARD -> MOD -> SIDD -> SIRD
target_rows <- c('ND', 'MARD', 'MOD', 'SIDD', 'SIRD')
radar_data <- radar_data[target_rows, ]

# Define axes limits (Max/Min) for the radar chart
max_min <- data.frame(matrix(rep(c(0.45, -0.7), ncol(radar_data)), 
                             nrow = 2, byrow = FALSE))
colnames(max_min) <- colnames(radar_data)
rownames(max_min) <- c("Max", "Min")
radar_data <- rbind(max_min, radar_data)

# Visualization
par(mar = c(1, 1, 1, 1))
mycolor <- c('#5470c6','#91cc75','#fac858','#ee6666','#73c0de')
radarchart(radar_data,
           plwd = 2, 
           plty = 1,
           pcol = mycolor,
           pfcol = scales::alpha(my9color,0.25),
           cglcol = "grey", cglty = 1, cglwd = 0.8,vlcex = 0)

# Add legend
legend(
  x = "left", legend = rownames(radar_data[-c(1,2),]), horiz = F,
  bty = "n", pch = 20 , col = my9color,
  text.col = "black", cex = 1, pt.cex = 1.5
)


# Figure 4G ----------

run_pseudobulk_deseq2 <- function(seurat_obj, cell_type_name, mode = "ovo", target_group = "SIDD",
                                 ref_group = "ND", min_counts = 1, min_donors = 2) {
 
  message("--- Processing: ", cell_type_name, " | Mode: ", mode, " ---")
  
  sub_obj <- subset(seurat_obj, subset = cell_type == cell_type_name)
  if (mode == "ovo") {
    sub_obj <- subset(sub_obj, subset = subgroup %in% c(target_group, ref_group))
    title_suffix <- paste0(target_group, "_vs_", ref_group)
  } else {
    sub_obj$condition <- ifelse(sub_obj$subgroup == target_group, target_group, "TheRest")
    ref_group <- "TheRest"
    title_suffix <- paste0(target_group, "_vs_TheRest")
  }
  
  n_target <- length(unique(sub_obj$donor[sub_obj$subgroup == target_group]))
  if (n_target < 2) {
    warning("Target group '", target_group, "' has only ", n_target, " donor. DESeq2 requires N >= 2.")
    return(NULL) 
  }
  
  # Aggregate counts per donor per subgroup
  pb_counts <- AggregateExpression(
    sub_obj, group.by = c("donor", "subgroup"), 
    assays = "RNA", slot = "counts", return.seurat = FALSE)$RNA

  sample_info <- data.frame(col_name = colnames(pb_counts)) %>%
    tidyr::separate(col_name, into = c("donor", "subgroup_orig"), sep = "_", remove = FALSE) %>%
    column_to_rownames("col_name")
  
  if (mode == "ovr") {
    sample_info$condition = factor(ifelse(sample_info$subgroup_orig == target_group, target_group, "TheRest"),
                                   levels = c("TheRest", target_group))
  } else {
    sample_info$condition = factor(sample_info$subgroup_orig, levels = c(ref_group, target_group))
  }
  
  # DESeq2 analysis
  dds <- DESeqDataSetFromMatrix(countData = round(pb_counts), colData = sample_info, design = ~ condition)
  keep <- rowSums(counts(dds) >= min_counts) >= min_donors
  dds <- dds[keep,]
  
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("condition", target_group, ref_group), alpha = 0.1)
  
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    mutate(cell_type = cell_type_name, comparison = title_suffix) %>%
    arrange(pvalue)
  
  return(res_df)
}

# SIDD-ovo-alpha ----
alpha_sidd_ovo <- run_pseudobulk_deseq2(panc_integrated, "alpha", mode = "ovo", target_group = "SIDD", ref_group = "ND")

cell_degs <- sig_subgroup_ovo$alpha_SIDD_vs_ND$gene
donor_degs <- alpha_sidd_ovo %>% filter(pvalue < 0.05) %>% pull(gene)
overlap_genes <- intersect(cell_degs, donor_degs)
length(overlap_genes)
length(overlap_genes) / dim(sig_subgroup_ovo$alpha_SIDD_vs_ND)[1]

# Cell-level (Wilcoxon)
cell_res <- sig_subgroup_ovo$alpha_SIDD_vs_ND %>%
  filter(gene %in% overlap_genes) %>%
  select(gene, avg_log2FC)

# Donor-level (DESeq2)
donor_res <- alpha_sidd_ovo %>%
  filter(gene %in% overlap_genes) %>%
  select(gene, log2FoldChange)

plot_df <- inner_join(cell_res, donor_res, by = "gene")
target_labels <- c("CELA3A", "CELA3B", "CLPS", "CPA1","CTRB1","CTRB2", # pancreatic secretion
                   "MTRNR2L12", "MTRNR2L8",  # apoptosis
                   "PLA2G1B", "PLCG2","PNLIP", # lipid metabolic
                   "TRIB3", # ER stress
                   "DACH1","SEMA3E"
                 )
ggplot(plot_df, aes(x = avg_log2FC, y = log2FoldChange)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(color = "#b9e3c2", alpha = 0.8, size = 3) +
  geom_smooth(method = "lm", color = "red", fill = "pink", alpha = 0.2) +
  stat_cor(method = "pearson") +
  geom_text_repel(
    aes(label = gene),
    data = filter(plot_df, gene %in% target_labels),
    size = 3.5,
    box.padding = 0.5,     
    point.padding = 0.3,   
    segment.color = 'black', 
    fontface = "italic",
    max.overlaps = Inf     
  ) +
  theme_bw() +
  labs(
    title = "Log2FC Consistency: SIDD vs ND (Alpha Cells)",
    subtitle = paste("Comparison of", length(overlap_genes), "Overlap Genes"),
    x = "Cell-level avg_log2FC (Wilcoxon)",
    y = "Donor-level log2FoldChange (DESeq2)"
  )

direction_check <- all(sign(plot_df$avg_log2FC) == sign(plot_df$log2FoldChange))
message("Do all overlap genes have the same direction? ", direction_check)
mismatched_genes = plot_df %>% filter(sign(avg_log2FC) != sign(log2FoldChange))

# SIDD-ovo-beta ----
beta_sidd_ovo <- run_pseudobulk_deseq2(panc_integrated, "beta", mode = "ovo", target_group = "SIDD", ref_group = "ND")

cell_degs <- sig_subgroup_ovo$beta_SIDD_vs_ND$gene
donor_degs <- beta_sidd_ovo %>% filter(pvalue < 0.05) %>% pull(gene)
overlap_genes <- intersect(cell_degs, donor_degs)
length(overlap_genes)
length(overlap_genes) / dim(sig_subgroup_ovo$beta_SIDD_vs_ND)[1]

# Cell-level (Wilcoxon)
cell_res <- sig_subgroup_ovo$beta_SIDD_vs_ND %>%
  filter(gene %in% overlap_genes) %>%
  select(gene, avg_log2FC)

# Donor-level (DESeq2)
donor_res <- beta_sidd_ovo %>%
  filter(gene %in% overlap_genes) %>%
  select(gene, log2FoldChange)

plot_df <- inner_join(cell_res, donor_res, by = "gene")
target_labels <- c("CEL", "CLPS", "CTRB1", "CTRC", # pancreatic secretion
                   "MTRNR2L12", "MTRNR2L8", "XKR6", # apoptosis
                   "PLA2G1B", "PLCG2","PNLIP", # lipid metabolic
                   "PPP1R1A","SST","CHGA","ALDH1A1","FKBP11","PTPRM","TUNAR", # beta-cell identity
                   "NDUFAF3","PDK4", "DNAJC15", # energy metabolism & mitochondrial dysfunction
                   "FXYD2","IGFBP7","KCNB2" # ion transport & cell signaling
)
ggplot(plot_df, aes(x = avg_log2FC, y = log2FoldChange)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(color = "#b9e3c2", alpha = 0.8, size = 3) +
  geom_smooth(method = "lm", color = "red", fill = "pink", alpha = 0.2) +
  stat_cor(method = "pearson", label.x = -8, label.y = 4.5) +
  geom_text_repel(
    aes(label = gene),
    data = filter(plot_df, gene %in% target_labels),
    size = 3.5,
    box.padding = 0.5,     
    point.padding = 0.3,   
    segment.color = 'black', 
    fontface = "italic",
    max.overlaps = Inf     
  ) +
  theme_bw() +
  labs(
    title = "Log2FC Consistency: SIDD vs ND (Beta Cells)",
    subtitle = paste("Comparison of", length(overlap_genes), "Overlap Genes"),
    x = "Cell-level avg_log2FC (Wilcoxon)",
    y = "Donor-level log2FoldChange (DESeq2)"
  )

direction_check <- all(sign(plot_df$avg_log2FC) == sign(plot_df$log2FoldChange))
message("Do all overlap genes have the same direction? ", direction_check)
mismatched_genes = plot_df %>% filter(sign(avg_log2FC) != sign(log2FoldChange))


# Figure 4H ----------

plot_scatter(deg_subgroup_ovo$beta_SIRD_vs_SIDD,
             c("RPS4Y1","EIF1AY","UTY","DDX3Y","AGPAT5","CPA1","PPP2R2A","FOXO3","DNAJB2","SYT11","MAPK6","KAT6A","CUL1","OPA1","REV1","SESTD1","AKT3",
               "XIST","AQP3","MT1E","MT1X","G6PC2","SST","PPP1R1A","UQCC2","FKBP11","CDKN1C","FXYD2","SERF2","IAPP"),
             'SIRD-SIDD: beta',-0.5,0.5,-0.9,0.9)


# Figure 4I ----------

do_enrich <- function(group){
  df <- sig_subgroup_ovo[[group]]
  genes <- list(up = df[which(df$avg_log2FC>0),6],down = df[which(df$avg_log2FC<0),6],all = df$gene)
  go_list <- list()
  kegg_list <- list()
  for (i in 1:length(genes)) {
    index <- names(genes)[i]
    print(index)
    gene_id <- AnnotationDbi::select(org.Hs.eg.db, keys = genes[[i]],column = "ENTREZID", keytype="SYMBOL")
    go_list[[index]] <- enrichGO(gene = gene_id$ENTREZID,OrgDb = org.Hs.eg.db,ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = T)
    go_list[[index]] <- as.data.frame(go_list[[index]])
    go_list[[index]] <- go_list[[index]][which(go_list[[index]]$pvalue<0.05),]
    
    kegg_list[[index]] <- enrichKEGG(gene = gene_id$ENTREZID,organism = "hsa",pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1)
    kegg_list[[index]] <- setReadable(kegg_list[[index]], OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    kegg_list[[index]] <- as.data.frame(kegg_list[[index]])  
    kegg_list[[index]] <- kegg_list[[index]][which(kegg_list[[index]]$pvalue<0.05),]
  }
  return(list(go = go_list,kegg = kegg_list))
}

rdDd_enrich <- do_enrich('beta_SIRD_vs_SIDD')


# Figure 4J ----------
# Figure 4J (Protein-protein interactions) was generated using STRING (https://string-db.org/) and Cytoscape (https://cytoscape.org/).
