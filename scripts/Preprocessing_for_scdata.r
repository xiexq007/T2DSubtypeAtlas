
# ------------------------------------------------------------------------------

# scRNA-seq meta-analysis of Human Pancreas Analysis Program (HPAP) data
# 1. HPAP Fluidigm C1
# 2. HPAP 10X-Chromium

# ------------------------------------------------------------------------------

# Load libraries
library(Seurat) 
library(dplyr)
library(tidyverse)
library(ggplot2)
library(stringr)
library(clustree)
library(harmony)
library(patchwork)
library(circlize)
library(tidydr)
library(cluster)
library(SoupX)


# 1. Ambient RNA correction and validation -------------------------------------------

# Load raw 10X data
# SoupX automatically finds 'raw' and 'filtered' matrices in the directory
data_dir <- "./raw_fastq/"
sc <- load10X(data_dir)

# Create temporary Seurat object for clustering
temp_obj <- CreateSeuratObject(counts = sc$toc)
temp_obj <- NormalizeData(temp_obj)
temp_obj <- FindVariableFeatures(temp_obj)
temp_obj <- ScaleData(temp_obj)
temp_obj <- RunPCA(temp_obj)
temp_obj <- FindNeighbors(temp_obj, dims = 1:20)
temp_obj <- FindClusters(temp_obj, resolution = 0.5)

# Pass cluster labels to SoupX
sc <- setClusters(sc, Idents(temp_obj))

# Run SoupX
sc <- autoEstCont(sc)

rho_val <- sc$fit$rhoEst
cat(paste0("Estimated contamination fraction: ", round(rho_val, 6), "\n"))

clean_counts <- adjustCounts(sc)

# Create Seurat objects for the validation comparison
orig_obj <- CreateSeuratObject(counts = sc$toc, project = "Original")
soup_obj <- CreateSeuratObject(counts = out, project = "SoupX_Cleaned")
orig_obj <- NormalizeData(orig_obj)
soup_obj <- NormalizeData(soup_obj)

# Visualization comparison
orig_obj <- RunPCA(orig_obj)
orig_obj <- RunUMAP(orig_obj, dims = 1:20)

soup_obj <- RunPCA(soup_obj)
soup_obj <- RunUMAP(soup_obj, dims = 1:20)

gene_to_plot <- "INS"
p1 <- FeaturePlot(orig_obj, features = gene_to_plot) +
      ggtitle("Before SoupX")
p2 <- FeaturePlot(soup_obj, features = gene_to_plot) +
      ggtitle("After SoupX")
p1 | p2

# Global expression stability
common_genes <- intersect(rownames(orig_obj), rownames(soup_obj))

avg_orig <- rowMeans(orig_obj@assays$RNA$data[common_genes, ])
avg_soup <- rowMeans(soup_obj@assays$RNA$data[common_genes, ])
correlation <- cor(avg_orig, avg_soup, method = "spearman")
df_cor <- data.frame(Original = avg_orig,
                     Cleaned  = avg_soup)

p_cor <- ggplot(df_cor, aes(x=Original, y=Cleaned)) +
  geom_point(alpha=0.3, size=0.5, color="darkblue") +
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed") +
  theme_bw() +
  labs(title = "Global Expression Stability",
       subtitle = paste("Spearman Correlation =", round(correlation, 6)),
       x = "Original Expression (Log-Normalized)",
       y = "SoupX Corrected Expression (Log-Normalized)")+
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


# 2. Load data -----------------------------------------------------------------------

# HPAP Fluidigm C1
fluidigm = readRDS("./fluidigm.rds")

# HPAP 10X-Chromium
chromium = readRDS("./10xchromium.rds")

dim(fluidigm) # 30595 genes 2176 cells
dim(chromium) # 33959 genes 146410 cells

# 3. QC, pre-processing  --------------------------------------------------------------

for (object in objects) {
    
    object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
    object[["percent.ribo"]] <- PercentageFeatureSet(object, pattern = "^RPS|^RPL")
    
    # Exploratory data analysis(gene, UMI, mito count distribution) was visualized
}

# Based on EDA , filter cells with genes < 200 and > 8000, mt percent > 28%
filtered_fluidigm <- subset(objects[[1]], subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 28) 
filtered_chromium <- subset(objects[[2]], subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 28) 

# 4. Normalize ------------------------------------------------------------------------

filtered_fluidigm = NormalizeData(filtered_fluidigm, normalization.method = "LogNormalize", scale.factor = 10000)
filtered_chromium = NormalizeData(filtered_chromium, normalization.method = "LogNormalize", scale.factor = 10000)

# 5. Supervised filtering to remove double positive hormonal cells in each dataset ----

# HPAP Fluidigm C1
hormones <- c("INS","GCG","PPY","SST","GHRL")
n <- length(hormones)
fluidigm_ridges <- list()

for (i in 1:5) {
  fluidigm_ridges[[i]] <- RidgePlot(filtered_fluidigm,features=hormones[i], group.by = "method")
}

fluidigm_scatters <- list()
for (i in 1:(n-1)) {
  for (j in (i+1):n){
    fluidigm_scatter = FeatureScatter(filtered_fluidigm, feature1 = hormones[i],feature2 = hormones[j], group.by = "method")+
      labs(title = NULL, caption = NULL)+
      theme(legend.position = "none",axis.text = element_blank(),axis.ticks = element_blank())
    fluidigm_scatters[[length(fluidigm_scatters)+1]] = fluidigm_scatter
  }
}

filtered_fluidigm <- subset(filtered_fluidigm, subset = ((INS>6 & GCG>5.5)|(INS>6 & PPY>4.5)|(INS>6 & SST>4.5)|(INS>6 & GHRL>4)|(GCG>5.5 & PPY>4.5)|(GCG>5.5 & SST>4.5)|(GCG>5.5 & GHRL>4)|(PPY>4.5 & SST>4.5)|(PPY>4.5 & GHRL>4)|(SST>4.5 & GHRL>4)), invert = T)

dim(filtered_fluidigm) # 30595 genes 1695 cells

# HPAP 10X-Chromium
chromium_ridges <- list()
for (i in 1:5) {
  chromium_ridges[[i]] <- RidgePlot(filtered_chromium,features=hormones[i],group.by = "method")
}

chromium_scatters <- list()
for (i in 1:(n-1)) {
  for (j in (i+1):n){
    chromium_scatter = FeatureScatter(filtered_chromium, feature1 = hormones[i],feature2 = hormones[j], group.by = "method")+
      labs(title = NULL, caption = NULL)+
      theme(legend.position = "none",axis.text = element_blank(),axis.ticks = element_blank())
    chromium_scatters[[length(chromium_scatters)+1]] = chromium_scatter
  }
}

filtered_chromium <- subset(filtered_chromium, subset = ((INS>6 & GCG>5.5)|(INS>6 & PPY>5)|(INS>6 & SST>4.5)|(INS>6 & GHRL>4)|(GCG>5.5 & PPY>5)|(GCG>5.5 & SST>4.5)|(GCG>5.5 & GHRL>4)|(PPY>5 & SST>4.5)|(PPY>5 & GHRL>4)|(SST>4.5 & GHRL>4)), invert = T)

dim(filtered_chromium) # 33959 genes 129388 cells

# 5. Combining datasets --------------------------------------------------------------

DefaultAssay(object <- filtered_fluidigm) = "RNA"
panc_combined <- merge(filtered_fluidigm, y = filtered_chromium, add.cell.ids = c("fluidigm", "chromium"), project = "AllPanc")
panc_combined # 46753 genes 131083 cells

# 6. Integrating datasets with harmony -----------------------------------------------

panc_combined <- NormalizeData(panc_combined)
panc_combined <- FindVariableFeatures(panc_combined, selection.method = "vst", nfeatures = 3000)
panc_combined <- ScaleData(panc_combined)
panc_combined <- RunPCA(panc_combined)
panc_combined <- RunUMAP(panc_combined, dims = 1:15, reduction = "pca", reduction.name = "umap.unintegrated")

# UMAP before harmony
DimPlot(panc_combined, reduction = "umap.unintegrated", group.by = "method")

# Harmony integration
panc_integrated <- IntegrateLayers(object = panc_combined, method = HarmonyIntegration,
                                  orig.reduction = "pca", new.reduction = "harmony",verbose = F)
panc_integrated[["RNA"]] <- JoinLayers(panc_integrated[["RNA"]])

# UMAP after harmony
DimPlot(panc_integrated, reduction = "umap.harmony", group.by = "method")

# Silhouette
calculate_asw <- function(obj, reduction, batch_col, n_sample = 5000, seed = 123) {
  set.seed(seed)
  embeddings = Embeddings(obj, reduction)[, 1:30]
  cells_to_sample = sample(rownames(embeddings), n_sample)
  sub_embeddings = embeddings[cells_to_sample, ]
  sub_labels = obj@meta.data[cells_to_sample, batch_col]
  
  dist_matrix = dist(sub_embeddings)
  sil = silhouette(as.numeric(as.factor(sub_labels)), dist_matrix)
  
  return(mean(sil[, 3]))
}

n_iterations <- 100
asw_pre_list <- c()
asw_post_list <- c()

for (i in 1:n_iterations) {
  message(paste("Running iteration", i, "..."))
  current_seed = i * 100
  
  # Before harmony (PCA embeddings)
  asw_pre = calculate_asw(panc_integrated, reduction = "pca", batch_col = "method", seed = current_seed)
  asw_pre_list = c(asw_pre_list, asw_pre)
  
  # After harmony (harmony embeddings)
  asw_post = calculate_asw(panc_integrated, reduction = "harmony", batch_col = "method", seed = current_seed)
  asw_post_list = c(asw_post_list, asw_post)
}

results_df <- data.frame(
  Value = c(asw_pre_list, asw_post_list),
  Group = rep(c("Before Harmony", "After Harmony"), each = n_iterations)
)

summary_stats <- results_df %>%
  group_by(Group) %>%
  summarise(Mean_ASW = mean(Value), SD_ASW = sd(Value))
print(summary_stats)

ggplot(results_df, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(width = 0.5) +
  theme_classic() +
  labs(y = "Average Silhouette Width") +
  scale_fill_manual(values = c("Before Harmony" = "#E64B35FF", "After Harmony" = "#4DBBD5FF"))+
  theme(legend.position = "none")+
  stat_compare_means(method = "wilcox.test",label = "p.signif")

# 7. Clustering ----------------------------------------------------------------------

# Resolution decision
seq <- seq(0.1, 1.5, by = 0.1)
for (res in seq) {
  panc_integrated <- FindClusters(panc_integrated, resolution = res)
}
clustree(panc_integrated, prefix = 'RNA_snn_res.')

panc_integrated <- FindNeighbors(panc_integrated, dims = 1:15, reduction = "harmony")
panc_integrated <- FindClusters (panc_integrated, resolution = 0.9)
panc_integrated <- RunUMAP(panc_integrated, reduction = "harmony", dims = 1:15, reduction.name = "umap.harmony")

# Clustering visuals 
DimPlot(panc_integrated, reduction = "umap.harmony", group.by = "seurat_clusters",label = T)

# 8. CellType Annotation--------------------------------------------------------------

# Plotting known cell type marker genes
panc_markers <- c("GCG","TTR","INS","IAPP","SST","HHEX","LEPR",
                 "PPY","GHRL","KRT19","CFTR","KRT7","PRSS1","CPA1",
                 "CPA2","REG1A","COL1A1","PDGFRB","COL1A2","TIMP1","FLT1",
                 "VWF","PECAM1","CD68","CD74","CD163","TPSAB1","CPA3")

DotPlot(panc_integrated, features = panc_markers, group.by = "seurat_clusters")+
  scale_x_discrete("")+scale_y_discrete("")+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3497d8","white","#e23232"))+
  coord_flip()

# Assigning cell type identity to clusters 
new.cluster.ids <- c("alpha","beta","alpha","beta","acinar","alpha","alpha","ductal","acinar",
                    "stellate","endothelial","beta","beta","delta","acinar","acinar","alpha",
                    "stellate","beta","ductal","acinar","acinar","beta","stellate","PP/gamma",
                    "ductal","macrophage","mast","acinar","alpha","stellate","endothelial","ductal",
                    "ductal","endothelial","endothelial","acinar","ductal","acinar","alpha"
                    )

names(new.cluster.ids) <- levels(panc_integrated)
panc_integrated <- RenameIdents(panc_integrated, new.cluster.ids)
panc_integrated@meta.data$cell_type <- Idents(object = panc_integrated)