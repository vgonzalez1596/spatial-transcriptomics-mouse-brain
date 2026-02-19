# Install required packages for this repo

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_pkgs <- c(
  "clusterProfiler",
  "org.Mm.eg.db"
)

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

# CRAN packages
cran_pkgs <- c(
  "Seurat",
  "SeuratObject",
  "SeuratWrappers",
  "Banksy",
  "ggplot2",
  "patchwork",
  "dplyr",
  "tibble",
  "hdf5r",
  "arrow",
  "ape",
  "msigdbr"
)

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

message("Packages installed. Next run: source('analysis/'mouse_brain_spatial_Visium.R')")
