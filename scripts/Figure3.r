
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


# Figure 3A ----------

deg_subgroup_ovr <- readRDS("./deg_subgroup_ovr.rds")

# Define criteria for Significant DEGs: adj. P-value < 0.001 and |log2FC| > 1
filter_deg <- function(df) {
  result <- df[df$p_val_adj < 0.001 & abs(df$avg_log2FC) > 1, ]
  result$gene <- rownames(result)
  return(result)
}

sig_subgroup_ovr <- lapply(deg_subgroup_ovr, function(group) {
  lapply(group, function(cell_type) {
    filter_deg(cell_type)
  })
})

# Number statistics
num_df <- data.frame()
for (group in names(sig_subgroup_ovr)) {
  cell_types <- sig_subgroup_ovr[[group]]
  for (cell_type in names(cell_types)) {
    df <- cell_types[[cell_type]]
    
    up_count <- sum(df$avg_log2FC > 0)
    down_count <- sum(df$avg_log2FC < 0)
    
    # Store results in a tidy format
    num_df <- rbind(num_df, data.frame(
      group = group,
      cell_type = cell_type,
      regulation = "Up-regulated",
      count = up_count
    ))
    num_df <- rbind(num_df, data.frame(
      group = group,
      cell_type = cell_type,
      regulation = "Down-regulated",
      count = down_count
    ))
  }
}

num_df$group <- factor(num_df$group, levels = c("ND", "MARD", "MOD", "SIDD", "SIRD"))
num_df$cell_type <- factor(num_df$cell_type, levels = c("alpha", "beta", "delta", "PP/gamma", "acinar",
                                                       "stellate","ductal","endothelial","macrophage","mast"))
mirror_df <- num_df %>%
  mutate(display_count = ifelse(regulation == "Down-regulated", -count, count))

ggplot(mirror_df, aes(x = cell_type, y = display_count, fill = regulation)) +
  geom_bar(stat = "identity", width = 0.7, color = "white", linewidth = 0.1) +
  facet_wrap(~ group, scales = "free_y", nrow = 1) + 
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_text(aes(label = ifelse(count > 0, count, ""), 
                vjust = ifelse(regulation == "Up-regulated", -0.3, 1.2)), 
            size = 3) +
  scale_fill_manual(values = c("Up" = "#d3e0f2", "Down" = "#cccccc")) +
  theme_classic() +
  labs(y = "Number of DEGs (Down < 0 > Up)", x = "") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    strip.background = element_rect(fill = "gray95", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "top")


# Figure 3B ----------
# Figure 3B was generated using a publicly available online bioinformatics visualization platform.


# Figure 3C (The left panel: heatmap) ----------

plot_cluster_heatmap <- function(seurat_obj, cell_type, deg_list, highlight_genes, col_fun) {
  
 """
  Generates a publication-quality heatmap for normalized gene expression.
  
  Args:
    seurat_obj: Integrated Seurat object containing 'cell_type' and 'subgroup' metadata.
    cell_type: Target cell population (e.g., 'beta').
    deg_list: Vector of gene names to be included in the heatmap.
    highlight_genes: Vector of specific genes to label on the right axis.
    col_fun: Color mapping function from circlize.
  """
  
  panc_sub <- subset(seurat_obj, cell_type == cell_type)
  aver_exp <- AverageExpression(panc_sub, group.by = 'subgroup', layer = 'data')
  aver_exp <- as.data.frame(aver_exp$RNA)
  aver_exp_scaled <- t(scale(t(aver_exp)))
  target_groups <- c('ND', 'MARD', 'MOD', 'SIDD', 'SIRD')
  plot_mat <- aver_exp_scaled[intersect(deg_list, rownames(aver_exp_scaled)), target_groups]

  label_indices <- which(rownames(plot_mat) %in% highlight_genes)
  labels <- rownames(plot_mat)[label_indices]
  
  right_anno <- rowAnnotation(
    link = anno_mark(
      at = label_indices, 
      labels = labels, 
      labels_gp = gpar(fontsize = 10)
    )
  )


  ht = Heatmap(aver_exp2,cluster_columns = F,cluster_rows = T,show_row_names = F,col = mycol,right_annotation = a,
               column_names_rot = 60,column_names_gp = gpar(fontsize = 10) ,column_names_side = c('top'),
               name = 'expression')
    
  ht <- Heatmap(
    plot_mat,
    name = "expression",
    cluster_columns = FALSE,    
    cluster_rows = TRUE,        
    show_row_names = FALSE,     
    col = col_fun,
    right_annotation = right_anno,
    column_names_rot = 60,
    column_names_side = "top",
    column_names_gp = gpar(fontsize = 10)
  )
  
  draw(ht)
}


mycol <- colorRamp2(c(-2, 0, 2), c("#f8f8d6", "#57bac2", "#202c5c"))
genes_to_label <- c(
  "CD44","UTY","RFX3", "FRMPD4", "TSPAN5", "PDE3A","RFX2","ARMC2","CTRB1",
  "CELA3B", "CELA3A", "PPP1R1A", "IAPP","TIMM13", "MPG", "CHGA", "FKBP11", "IGFBP7"
)

panc_integrated <- readRDS("./integrated_data.rds")
plot_cluster_heatmap(panc_integrated, 'beta', sidd_beta_degs, genes_to_label, mycol)


# Figure 3C (The right panel: boxplot) ----------

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

# ovr -----------
# SIDD-beta
# Beta cell homeostasis: PPP1R1A, IAPP, 
# Hormone-related: GCG, SST
# Energy metabolism-related: NDUFAF3, NDUFB2
# Pancreatic secretion: CELA3A, CELA3B, CTRB1, CLPS, CTRC, CEL,PLA2G1B,PNLIP 
# Apoptosis-related: MTRNR2L12, MTRNR2L8
# Cell adhesion & ECM remodeling: CD44, CTNND2
# Sex-linked: DDX3Y, UTY, USP9Y, RPS4Y1

sidd_ovr_beta = plot_donor_concordance(
  seurat_obj   = panc_integrated,
  genes        = c('CELA3A','PPP1R1A'),
  celltype     = "beta",
  target_group = "SIDD",
  mode         = "ovr" )


# Figure 3D ----------

mycol <- c("#500e57","#1d3263","#226c87","#0d6e29","#4f9759","#eab518","#e48e11","#c6790d","#e03215","#ab2321")

plot_multi_volcano <- function(group,dfbar1,dfbar2,title){
  
  cell_order <- c("alpha", "beta", "delta","PP/gamma","acinar","stellate","ductal","endothelial","macrophage","mast")
  
  # Combine DEG lists from sig_subgroup_ovr for the specific clinical group
  results_df <-  bind_rows(sig_subgroup_ovr[[group]], .id = "cell") %>%
    mutate(cell = factor(cell,levels = cell_order))  %>%
    mutate(label = ifelse(avg_log2FC>0,"up","down"))
  
  dfcol  <- data.frame(x = cell_order,y = 0,label=c(0:9))
  dfcol$x <- factor(dfcol$x, levels = cell_order)

  # Top genes selection  
  top_list <- list()
  for (cell_type in cell_order) {
    cell_data <- results_df %>%
      filter(cell == cell_type)
    
    num_up <- sum(cell_data$avg_log2FC > 0)
    num_down <- sum(cell_data$avg_log2FC < 0)
    
    top_up <- cell_data %>%
      filter(avg_log2FC > 0) %>%
      arrange(desc(avg_log2FC)) %>%
      slice_head(n = min(3, num_up)) 
        top_down = cell_data %>%
      filter(avg_log2FC < 0) %>%
      arrange(avg_log2FC) %>%
      slice_head(n = min(3, num_down))
    top_list[[cell_type]] = bind_rows(top_up, top_down)
  }
  
  top_df <- bind_rows(top_list)
  
  ggplot() +
    geom_col(data = dfbar1,mapping = aes(x = x,y = y),fill = "#dcdcdc",alpha = 0.6)+
    geom_col(data = dfbar2,mapping = aes(x = x,y = y),fill = "#dcdcdc",alpha = 0.6)+
    geom_jitter(data = results_df,aes(x = cell, y = avg_log2FC, color = label),size = 0.8,width =0.4)+
    geom_tile(data = dfcol,aes(x = x,y = y),height=1.6,color = "black",fill = mycol, size = 1, show.legend = F) + 
    geom_text_repel(data = top_df,aes(x = cell,y = avg_log2FC,label=gene),force = 0.5,max.overlaps = 20,arrow = arrow(length = unit(0.008, "npc"),type = "open", ends = "last")) + 
    scale_color_manual(name=NULL,values = c('#66ced5','#fb4e06')) +
    geom_text(data = dfcol,aes(x = x,y = y,label = x),size = 3,color ="white",fontface = "bold") +
    theme_minimal() + 
    labs(x = "", y = "average log2FC",title = title) +
    theme(
      axis.title = element_text(size = 12,color = "black"),
      axis.line.y = element_line(color = "black",size = 0.8),
      axis.line.x = element_blank(),
      axis.text.x = element_blank(),
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 12),
      legend.position = "top",
      legend.direction = "vertical",
      legend.justification = c(1,0),
      legend.text = element_text(size = 10) ) 
}

dfbar1 <- data.frame(x = cell_order,y = c(3.5,4,4,6,4,3.5,3.5,4,6,3))
dfbar2 <- data.frame(x = cell_order,y = c(-3.5,-3.5,-3.5,-6,-3,-5,-7,-3.5,-5,-8))
plot_multi_volcano('SIRD',dfbar1,dfbar2,'SIRD vs the rest as reference')


# Figure 3E ----------
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
ovr_sidd_beta_gsea <- run_gsea(deg_subgroup_ovr, type="ovr", group="SIDD", cell="beta")
ovr_sidd_alpha_gsea <- run_gsea(deg_subgroup_ovr, type="ovr", group="SIDD", cell="alpha")
ovr_sidd_delta_gsea <- run_gsea(deg_subgroup_ovr, type="ovr", group="SIDD", cell="delta")

# SIRD
ovr_sird_beta_gsea <- run_gsea(deg_subgroup_ovr, type="ovr", group ="SIRD", cell="beta")
ovr_sird_alpha_gsea <- run_gsea(deg_subgroup_ovr, type="ovr", group="SIRD", cell="alpha")
ovr_sird_delta_gsea <- run_gsea(deg_subgroup_ovr, type="ovr", group="SIRD", cell="delta")

target_terms <-  c('chromatin remodeling',
                  'DNA damage response',
                  'cellular response to stress',
                  'response to endoplasmic reticulum stress',
                  'lipid homeostasis',
                  'digestion',
                  'response to peptide hormone',
                  'response to insulin',
                  'cellular response to insulin stimulus',
                  'oxidative phosphorylation',
                  'ATP biosynthetic process',
                  'insulin secretion',
                  'hormone metabolic process',
                  'regulation of insulin secretion',
                  'regulation of peptide secretion',
                  'regulation of protein secretion' )

gsea_list <- list(
  "SIDD_alpha" = ovr_sidd_alpha_gsea,
  "SIDD_beta"  = ovr_sidd_beta_gsea,
  "SIDD_delta" = ovr_sidd_delta_gsea,
  "SIRD_alpha" = ovr_sird_alpha_gsea,
  "SIRD_beta"  = ovr_sird_beta_gsea,
  "SIRD_delta" = ovr_sird_delta_gsea
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
        cluster_rows = TRUE,       
        cluster_columns = FALSE,   
        rect_gp = gpar(col = "black", lwd = 0.5), 
        cell_fun = function(j, i, x, y, width, height, fill) {
          if (fdr_mat[i, j] < 0.05) {
            # Asterisk for FDR < 0.05
            grid.text("*", x, y, gp = gpar(fontsize = 15, fontface = "bold"))} 
          else if (p_mat[i, j] < 0.05) {
            # Plus sign for nominal p-value < 0.05
            grid.text("+", x, y, gp = gpar(fontsize = 12))}},
          column_split = factor(rep(c("SIDD", "SIRD"), each = 3), levels = c("SIDD", "SIRD")),
          row_names_side = "left",
          column_names_rot = 45,
          row_names_gp = gpar(fontsize = 10)
          )

# Figure 3F ----------
# GSEA Plot for oxidative phosphorylation

df <- as.data.frame(ovr_sidd_beta_gsea)
metabolic_id <- df %>% 
  filter(grepl("oxidative phosphorylation", Description, ignore.case = T)) %>% 
  slice(1) %>% pull(ID)

p1 <- gseaplot2(ovr_sidd_alpha_gsea, 
               geneSetID = metabolic_id, 
               title = "SIDD alpha Cells",
               color = "#86b35b", 
               pvalue_table = TRUE, pvalue_table_rownames = NULL,
               ES_geom = "line")
p2 <- gseaplot2(ovr_sidd_beta_gsea, 
               geneSetID = metabolic_id, 
               title = "SIDD beta Cells",
               color = "#86b35b", 
               pvalue_table = TRUE, pvalue_table_rownames = NULL,
               ES_geom = "line")
p3 <- gseaplot2(ovr_sidd_delta_gsea, 
               geneSetID = metabolic_id, 
               title = "SIDD delta Cells",
               color = "#86b35b", 
               pvalue_table = TRUE, pvalue_table_rownames = NULL,
               ES_geom = "line")
p4 <- gseaplot2(ovr_sird_alpha_gsea, 
               geneSetID = metabolic_id, 
               title = "SIRD alpha Cells",
               color = "firebrick", 
               pvalue_table = TRUE, pvalue_table_rownames = NULL,
               ES_geom = "line")
p5 <- gseaplot2(ovr_sird_beta_gsea, 
               geneSetID = metabolic_id, 
               title = "SIRD beta Cells",
               color = "firebrick", 
               pvalue_table = TRUE, pvalue_table_rownames = NULL,
               ES_geom = "line")
p6 <- gseaplot2(ovr_sird_delta_gsea, 
               geneSetID = metabolic_id, 
               title = "SIRD delta Cells",
               color = "firebrick", 
               pvalue_table = TRUE, pvalue_table_rownames = NULL,
               ES_geom = "line")

p1 <- as.ggplot(p1)
p2 <- as.ggplot(p2)
p3 <- as.ggplot(p3)
p4 <- as.ggplot(p4)
p5 <- as.ggplot(p5)
p6 <- as.ggplot(p6)

cowplot::plot_grid(p1, p2, p3, p4, p5, p6, ncol = 3)

# Figure 3G ----------

get_gene_fc <- function(gsea_obj) {
  genelist = gsea_obj@geneList 
  return(genelist)}

fold_change_list <- get_gene_fc(ovr_sidd_alpha_gsea)

selected_pathways <-  c("insulin secretion", 
                        "ATP biosynthetic process",
                        "cellular response to insulin stimulus")

cnetplot(ovr_sidd_alpha_gsea, 
         showCategory = selected_pathways, 
         foldChange = fold_change_list, 
         circular = TRUE,               
         colorEdge = TRUE,              
         node_label = "all") +           
  scale_color_gradient2(low = "#5a99d2", mid = "white", high = "red", midpoint = 0) 