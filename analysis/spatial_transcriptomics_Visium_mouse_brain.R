############################################################
# Here, will will carry out spatial transcriptomics analysis of
# a Visium HD mouse brain to gain insight into the 
# spatially restricted areas of gene expression in the brain. 
# We will follow the Seurat Visium HD analysis workflow. 
# The dataset can be found here: https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-mouse-brain-he
# We will be using only the 8um binned dataset. 

############################################################
# First, we load all needed libraries. 
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SeuratWrappers)
library(Banksy)  
library(ape)
library(hdf5r)
library(arrow)
library(tibble)
library(clusterProfiler)
library(org.Mm.eg.db)
library(msigdbr)

# Helper: saving all generated plots
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

save_plot_png <- function(path, plot, width = 2100, height = 1500, res = 300) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  png(filename = path, width = width, height = height, res = res)
  print(plot)
  dev.off()
  invisible(path)
}

# Now, we load in the 8um binned dataset. 
# We're using the filtered_feature_bc_matrix.h5 file and spatial/ directory.
localdir <- "/home/vg8/Visium_mouse_brain"
mouse_brain <- Load10X_Spatial(data.dir = localdir, bin.size = 8)

# We want to ensure we're using the 8um binned assay.
DefaultAssay(mouse_brain) <- "Spatial.008um"

# Next, we ensure our data has loaded in and we can plot our RNA counts onto the mouse brain image. 
# VlnPlot enables us to see spread of RNA counts.
vln.plot <- VlnPlot(mouse_brain, features = "nCount_Spatial.008um", pt.size = 0) + 
  theme(axis.text = element_text(size = 6)) + NoLegend()
# SpatialFeaturePlot allows us to see the RNA counts in relation to their localization (in each region of the brain).
# Note that low cellular density in specific tissue regions can be the reason for low RNA count
count.plot <- SpatialFeaturePlot(mouse_brain, features = "nCount_Spatial.008um") + 
  theme(legend.position = "right")
vln.plot | count.plot

p_qc <- vln.plot | count.plot
print(p_qc)
save_plot_png(
  "results/figures/qc_counts_violin_and_spatial.png",
  p_qc, width = 4200, height = 2100, res = 300)

############################################################
# Next, we use standard log-normalization to normalize.
DefaultAssay(mouse_brain) <- "Spatial.008um"
mouse_brain <- NormalizeData(mouse_brain)

# With normalization complete, we can begin to visualize gene expression and 
# ensure gene expression patterns are as expected.
p1 <- SpatialFeaturePlot(mouse_brain, features = "Snap25") + ggtitle("Snap25 expression")
p2 <- SpatialFeaturePlot(mouse_brain, features = "Hpca") + ggtitle("Hpca expression")
p3 <- SpatialFeaturePlot(mouse_brain, features = "Gfap") + ggtitle("Gfap expression")
p4 <- SpatialFeaturePlot(mouse_brain, features = "Tbr1") + ggtitle("Tbr1 expression")
p1 | p2 | p3 | p4

p_markers <- p1 | p2 | p3 | p4
print(p_markers)
save_plot_png(
  "results/figures/marker_expression_snap25_hpca_gfap_tbr1.png",
  p_markers, width = 4800, height = 2100, res = 300)
# Plotting Snap25 (marks all neurons), Mbp (marks white matter), 
# Hpca (marks hippocampus), and Tbr1 (marks cortex) demonstrates
# that expected gene expression patterns are indeed present. 

############################################################
# Next, we carry out unsupervised clustering. 
# The Seurat team has noted that using a sketch clustering workflow
# has an improved performance for Visium HD datasets, enabling
# better identification of rare and spatially restricted groups. 
# So here, we subsample 50,000 cells, perform clustering on them, and
# then project the cluster labels back onto the full dataset. 

# Find the variable genes
mouse_brain <- FindVariableFeatures(mouse_brain)
# Scale the data for PCA and neighbors
mouse_brain <- ScaleData(mouse_brain)
# Subsample 50,000 cells and create a new sketch assay 
mouse_brain <- SketchData(object = mouse_brain, ncells = 50000, method = "LeverageScore", sketched.assay = "sketch")
# Now use the sketch assay to carry out clustering workflow. 
DefaultAssay(mouse_brain) <- "sketch"
mouse_brain <- FindVariableFeatures(mouse_brain)
mouse_brain <- ScaleData(mouse_brain)
mouse_brain <- RunPCA(mouse_brain, assay = "sketch", reduction.name = "pca.sketch")
mouse_brain <- FindNeighbors(mouse_brain, assay = "sketch", reduction = "pca.sketch", dims = 1:50)
mouse_brain <- FindClusters(mouse_brain, cluster.name = "seurat_cluster.sketched", resolution = 3)
mouse_brain <- RunUMAP(mouse_brain, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = TRUE, dims = 1:50)

# Now we project the sketch-derived clusters and embeddings back to full dataset. 
mouse_brain <- ProjectData(
  object = mouse_brain,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched"))

# If desired, we can compare the clustering results for the sketched and full datasets.
# Sketched clustering 
DefaultAssay(mouse_brain) <- "sketch"
Idents(mouse_brain) <- "seurat_cluster.sketched"
p1 <- DimPlot(mouse_brain, reduction = "umap.sketch", label = FALSE) + 
  ggtitle("Sketched clustering (50,000 bins)") + theme(legend.position = "bottom")
# Full clustering 
DefaultAssay(mouse_brain) <- "Spatial.008um"
Idents(mouse_brain) <- "seurat_cluster.projected"
p2 <- DimPlot(mouse_brain, reduction = "full.umap.sketch", label = FALSE) + 
  ggtitle("Projected clustering (full dataset)") + theme(legend.position = "bottom")
p1 | p2

p_umap_compare <- p1 | p2
print(p_umap_compare)
save_plot_png(
  "results/figures/umap_sketched_vs_projected.png",
  p_umap_compare, width = 4200, height = 2100, res = 300)

# If desired, we can also visualize the unsupervised clusters based on their spatial location. 
# Note that running SpatialDimPlot(object, interactive = TRUE), also enables interactive visualization and exploration.
SpatialDimPlot(mouse_brain, group.by = "seurat_cluster.projected", label = TRUE, repel = TRUE, label.size = 4)
SpatialDimPlot(mouse_brain, group.by = "seurat_cluster.projected", interactive = TRUE, label = TRUE, repel = TRUE, label.size = 4)

p_spatial_clusters <- SpatialDimPlot(mouse_brain, group.by = "seurat_cluster.projected", label = TRUE, repel = TRUE, label.size = 4)
print(p_spatial_clusters)
save_plot_png(
  "results/figures/spatial_clusters_projected.png",
  p_spatial_clusters, width = 3000, height = 2400, res = 300)

# As we can see, while some clusters are clear and spatially restricted, some
# are mixed throughout the image. Therefore, it can be helpful to 
# plot only specific clusters to more easily see their spatial localization. 
Idents(mouse_brain) <- "seurat_cluster.projected"
cells <- CellsByIdentities(mouse_brain, idents = c(0, 9, 10, 29, 34, 35))
p1 <- SpatialDimPlot(mouse_brain,
                     cells.highlight = cells[setdiff(names(cells), "NA")],
                     cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T) + NoLegend()
p1

print(p1)
save_plot_png(
  "results/figures/spatial_highlight_selected_clusters_projected.png",
  p1, width = 4200, height = 3000, res = 300)

############################################################
# We can next utilize BANKSY to define and segment spatial tissue domains. 
# When performing clustering, BANKSY augments a spot’s expression pattern with 
# both the mean and the gradient of gene expression levels in a spot’s broader neighborhood.
# You must carefully consider 2 parameters: k_geom and lambda. 
# k_geom controls local neighborhood size. Larger values will yield larger domains.
# lambda controls influence of the neighborhood. Larger values yield more spatially coherent domains. 
# Running BANKSY and creating new BANKSY assay.
mouse_brain_BANKSY <- RunBanksy(mouse_brain, lambda = 0.8, verbose = TRUE,
                                assay = "Spatial.008um", slot = "data", features = "variable", k_geom = 50)
# Now use the BANKSY assay to carry out clustering workflow.
DefaultAssay(mouse_brain_BANKSY) <- "BANKSY"
mouse_brain_BANKSY <- RunPCA(mouse_brain_BANKSY, assay = "BANKSY", reduction.name = "pca.banksy", features = rownames(mouse_brain_BANKSY), npcs = 30)
mouse_brain_BANKSY <- FindNeighbors(mouse_brain_BANKSY, reduction = "pca.banksy", dims = 1:30)
mouse_brain_BANKSY <- FindClusters(mouse_brain_BANKSY, cluster.name = "banksy_cluster", resolution = 0.5)
Idents(mouse_brain_BANKSY) <- "banksy_cluster"
p_segmented_domains <- SpatialDimPlot(mouse_brain_BANKSY, group.by = "banksy_cluster", label = TRUE, repel = TRUE, label.size = 4)
p_segmented_domains

print(p_segmented_domains)
save_plot_png(
  "results/figures/banksy_segmented_domains.png",
  p_segmented_domains, width = 3000, height = 2400, res = 300)

# We can once again choose to plot specific clusters to better resolve their spatial localization.
# This time, domains should be clearer thanks to the BANKSY clustering. 
banksy_cells <- CellsByIdentities(mouse_brain_BANKSY)
p1 <- SpatialDimPlot(mouse_brain_BANKSY, cells.highlight = banksy_cells[setdiff(names(banksy_cells), "NA")], cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T) + NoLegend()
p1

print(p1)
save_plot_png(
  "results/figures/banksy_spatial_highlight_all_domains.png",
  p1, width = 4200, height = 3000, res = 300)

############################################################
# Using our BANKSY clustering, we can now compare gene expression in 
# distinct regions of the brain. Let's compare the white matter and 
# hippocampus and identify top differentially expression (DE) genes in each. 

# Before we begin DE, we need to identify where the white matter and hippocampus are using canonical markers. 
whitematter_markers <- c("Mbp", "Plp1", "Mobp", "Mag", "Mog", "Cnp", "Olig2")
hippocampus_markers  <- c("Prox1", "Wfs1", "Pcp4", "Grik1", "Zbtb20", "Calb1", "Neurod1")
markers_all <- unique(c(whitematter_markers, hippocampus_markers))
# We visualize expression of these canonical markers.
p_sp_wm <- SpatialFeaturePlot(
  mouse_brain_BANKSY,
  features = whitematter_markers,
  pt.size.factor = 1.2,
  ncol = 3
) + plot_annotation(title = "White Matter markers")
p_sp_hippo <- SpatialFeaturePlot(
  mouse_brain_BANKSY,
  features = hippocampus_markers,
  pt.size.factor = 1.2,
  ncol = 3
) + plot_annotation(title = "Hippocampus markers")
print(p_sp_wm)
print(p_sp_hippo)
# Zbtb20 clearly shows hippocampus region, overlaps with clusters 24, 25, 29.
SpatialFeaturePlot(mouse_brain_BANKSY, features = "Zbtb20") +
  SpatialDimPlot(mouse_brain_BANKSY, group.by = "banksy_cluster", label = TRUE)
# Plp1 clearly shows white matter region, overlaps with clusters 0, 26.
SpatialFeaturePlot(mouse_brain_BANKSY, features = "Plp1") +
  SpatialDimPlot(mouse_brain_BANKSY, group.by = "banksy_cluster", label = TRUE)

save_plot_png(
  "results/figures/white_matter_canonical_markers.png",
  p_sp_wm, width = 3600, height = 2700, res = 300)
save_plot_png(
  "results/figures/hippocampus_canonical_markers.png",
  p_sp_hippo, width = 3600, height = 2700, res = 300)
p_zbtb20_overlay <- SpatialFeaturePlot(mouse_brain_BANKSY, features = "Zbtb20") +
  SpatialDimPlot(mouse_brain_BANKSY, group.by = "banksy_cluster", label = TRUE)
print(p_zbtb20_overlay)
save_plot_png(
  "results/figures/zbtb20_featureplot_plus_banksy_clusters.png",
  p_zbtb20_overlay, width = 4200, height = 2100, res = 300)
p_plp1_overlay <- SpatialFeaturePlot(mouse_brain_BANKSY, features = "Plp1") +
  SpatialDimPlot(mouse_brain_BANKSY, group.by = "banksy_cluster", label = TRUE)
print(p_plp1_overlay)
save_plot_png(
  "results/figures/plp1_featureplot_plus_banksy_clusters.png",
  p_plp1_overlay, width = 4200, height = 2100, res = 300)

# First, define which clusters comprise the white matter and hippocampus. 
region_field <- "banksy_cluster"
table(mouse_brain_BANKSY[[region_field]][,1])[1:20]
clusters_vec <- mouse_brain_BANKSY[[region_field]][,1]
print(head(sort(unique(as.character(clusters_vec))), 50))
mouse_brain_BANKSY$anatomy <- "Other"
whitematter_clusters <- c("0", "26")
hippoc_clusters <- c("24", "25", "29")
mouse_brain_BANKSY$anatomy[as.character(clusters_vec) %in% whitematter_clusters] <- "White Matter"
mouse_brain_BANKSY$anatomy[as.character(clusters_vec) %in% hippoc_clusters] <- "Hippocampus"
table(mouse_brain_BANKSY$anatomy)

# Next, we want to keep only the white matter and hippocampus bins for DE
de_obj <- subset(mouse_brain_BANKSY, subset = anatomy %in% c("White Matter", "Hippocampus"))
Idents(de_obj) <- de_obj$anatomy
DefaultAssay(de_obj) <- "Spatial.008um"

# Finally, run DE and print the top DE genes 
markers_wm_vs_hip <- FindMarkers(
  object = de_obj,
  ident.1 = "White Matter",
  ident.2 = "Hippocampus",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.1)
markers_wm_vs_hip <- markers_wm_vs_hip %>%
  tibble::rownames_to_column("gene") %>%
  arrange(p_val_adj)
print(head(markers_wm_vs_hip, 20))
# This enables us to identify key expression differences between distinct regions of the brain. 
# As we can see, Plp1 has a 6.5 log2fold enrichment in the white matter cells, as expected!

############################################################
# Finally, we can compare the hippocampus and white matter regions of the brain
# in a different way: utilizing KEGG/GO terms to determine pathway enrichment. 
# We are using the objects generated in the previous section, so make sure to run it. 

# First, we determine which genes are more enriched in each tissue. 
wm_up <- markers_wm_vs_hip %>%
  filter(!is.na(p_val_adj), p_val_adj < 0.2, avg_log2FC > 0.1) %>%
  pull(gene) %>%
  unique()

hip_up <- markers_wm_vs_hip %>%
  filter(!is.na(p_val_adj), p_val_adj < 0.2, avg_log2FC < -0.1) %>%
  pull(gene) %>%
  unique()

# Next, we define the KEGG pathways we want to use from MSigDB (mouse, Legacy)
kegg_term2gene <- msigdbr(species = "mouse", category = "C2", subcategory = "CP:KEGG_LEGACY") %>%
  dplyr::select(term = gs_name, gene = gene_symbol) %>%
  dplyr::distinct()

# Finally, we run the pathway enrichment analysis 
kegg_wm_offline <- enricher(
  gene = wm_up,
  TERM2GENE = kegg_term2gene,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.2,
  qvalueCutoff = 0.25,
  minGSSize = 5,
  maxGSSize = 500)

kegg_hip_offline <- enricher(
  gene = hip_up,
  TERM2GENE = kegg_term2gene,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.2,
  qvalueCutoff = 0.25,
  minGSSize = 5,
  maxGSSize = 500)

# View top 15 KEGG results
wm_kegg_df <- as.data.frame(kegg_wm_offline)
hip_kegg_df <- as.data.frame(kegg_hip_offline)
print(head(wm_kegg_df, 15))
print(head(hip_kegg_df, 15))

# Plot the KEGG results
wm_kegg_df <- as.data.frame(kegg_wm_offline)
hip_kegg_df <- as.data.frame(kegg_hip_offline)
print(head(wm_kegg_df, 15))
print(head(hip_kegg_df, 15))

# Plot the KEGG results
p_kegg_wm <- dotplot(kegg_wm_offline, showCategory = 15) + ggtitle("KEGG Pathway Enrichment: White Matter enriched genes")
print(p_kegg_wm)

save_plot_png(
  "results/figures/kegg_dotplot_white_matter.png",
  p_kegg_wm, width = 3600, height = 2400, res = 300)

p_kegg_hip <- dotplot(kegg_hip_offline, showCategory = 15) + ggtitle("KEGG Pathway Enrichment: Hippocampus enriched genes")
print(p_kegg_hip)

save_plot_png(
  "results/figures/kegg_dotplot_hippocampus.png",
  p_kegg_hip, width = 3600, height = 2400, res = 300)

message("Done! All plots have been saved as PNGs in: results/figures/")
