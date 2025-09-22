
createSTobj <- function (gene_expr, spatial_loc, st_name, results_save_dir,
                                    min.cells = 400,
                                    min.features = 200,
                                    scale.factor = 10000,
                                    nfeatures = 2000,
                                    select_frac = 0.1) {
  # create results_save_dir if not exist
  if (!file.exists(results_save_dir)) {
    dir.create(results_save_dir)
  }
  results_save_dir <- paste0(results_save_dir, "/ST_raw_counts")
  if (!file.exists(results_save_dir)) {
    dir.create(results_save_dir)
  }
  results_save_dir <- paste0(results_save_dir, "/", st_name)
  if (!file.exists(results_save_dir)) {
    dir.create(results_save_dir)
  }

  # create seurat object and append location
  st_dataset <- CreateSeuratObject(counts = as.matrix(gene_expr), project = st_name, min.cells = min.cells, min.features = min.features)
  barcodes <- GetAssayData(st_dataset, "counts", "RNA") %>% colnames()
  spatial_loc <- spatial_loc[barcodes, ]
  #  https://divingintogeneticsandgenomics.com/post/how-to-construct-a-spatial-object-in-seurat/
  cents <- CreateCentroids(spatial_loc)
  segmentations.data <- list(
    "centroids" = cents,
    "segmentation" = NULL
  )
  coords <- CreateFOV(
    coords = segmentations.data,
    type = c("segmentation", "centroids"),
    molecules = NULL,
    assay = "RNA"
  )

  st_dataset[["coords"]] <- coords

  # Quality Control
  st_dataset[["percent.mt"]] <- PercentageFeatureSet(st_dataset, pattern = "^MT.")
  violin <- VlnPlot(st_dataset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  png(paste0(results_save_dir, "/QC_violin.png"))
  print(violin)
  dev.off()
  plot1 <- FeatureScatter(st_dataset, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(st_dataset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  combplot <- CombinePlots(plots = list(plot1, plot2))
  png(paste0(results_save_dir, "/QC_scatter.png"))
  print(combplot)
  dev.off()
  st_dataset <- subset(st_dataset, subset = nFeature_RNA >= quantile(st_dataset$nFeature_RNA, 0.05) &
                         percent.mt <= quantile(st_dataset$percent.mt, 0.95))
  # Normalizing the data
  st_dataset <- NormalizeData(st_dataset, normalization.method = "LogNormalize", scale.factor = scale.factor)

  # Scaling the data and then run PCA
  st_dataset <- FindVariableFeatures(st_dataset, selection.method = "vst", nfeatures = nfeatures)
  st_dataset <- ScaleData(st_dataset, features = rownames(st_dataset))
  st_dataset <- RunPCA(st_dataset, features = VariableFeatures(object = st_dataset))

  # Cluster and Run Umap
  st_dataset <- FindNeighbors(st_dataset, dims = 1:10)
  st_dataset <- FindClusters(st_dataset, resolution = 0.5)
  st_dataset <- RunUMAP(st_dataset, dims = 1:10)
  umap <- DimPlot(st_dataset, reduction = "umap")
  png(paste0(results_save_dir, "/umap.png"))
  print(umap)
  dev.off()
  # find marker genes for each cluster
  sc.markers <- FindAllMarkers(st_dataset, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25) %>% group_by(cluster) %>% top_frac(select_frac, wt = avg_log2FC) %>% top_frac(-select_frac, wt = p_val)
  write.csv(sc.markers, paste0(results_save_dir,  "/st_diff_expr.csv"))

  feature_map <- FeaturePlot(st_dataset, features = sc.markers$gene[1:9])
  png(paste0(results_save_dir, "/feature_map.png"))
  print(feature_map)
  dev.off()

  saveRDS(st_dataset, paste0(results_save_dir, "/", st_name, ".rds"))
}



createSCobj <- function(gene_expr, sc_name, sc_meta = c(), results_save_dir = "my_data",
                                   min.cells = 400,
                                   min.features = 200,
                                   scale.factor = 10000,
                                   nfeatures = 2000,
                                   select_frac = 0.1) {
  # note that sc_property function must has the same row name as barcodes.
  if (!file.exists(results_save_dir)) {
    dir.create(results_save_dir)
  }
  results_save_dir <- paste0(results_save_dir, "/SC_raw_counts")
  if (!file.exists(results_save_dir)) {
    dir.create(results_save_dir)
  }
  # create seurat object
  sc_dataset <- CreateSeuratObject(counts = as.matrix(gene_expr), project = sc_name, min.cells = min.cells, min.features = min.features)
  if (length(sc_meta) > 0) {
    meta_data <- sc_dataset@meta.data
    sc_meta <- sc_meta[rownames(meta_data), ]
    meta_data <- cbind(meta_data, sc_meta)
    # add other meta data here
    sc_dataset@meta.data <- meta_data
  }

  # Quality Control
  sc_dataset[["percent.mt"]] <- PercentageFeatureSet(sc_dataset, pattern = "^MT.")
  violin <- VlnPlot(sc_dataset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  png(paste0(results_save_dir, "/QC_violin.png"))
  print(violin)
  dev.off()
  plot1 <- FeatureScatter(sc_dataset, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(sc_dataset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  combplot <- CombinePlots(plots = list(plot1, plot2))
  png(paste0(results_save_dir, "/QC_scatter.png"))
  print(combplot)
  dev.off()
  sc_dataset <- subset(sc_dataset, subset = nFeature_RNA >= quantile(sc_dataset$nFeature_RNA, 0.05) &
                         percent.mt <= quantile(sc_dataset$percent.mt, 0.95))
  # Normalizing the data
  sc_dataset <- NormalizeData(sc_dataset, normalization.method = "LogNormalize", scale.factor = scale.factor)

  # Scaling the data and then run PCA
  sc_dataset <- FindVariableFeatures(sc_dataset, selection.method = "vst", nfeatures = nfeatures)
  sc_dataset <- ScaleData(sc_dataset, features = rownames(sc_dataset))
  sc_dataset <- RunPCA(sc_dataset, features = VariableFeatures(object = sc_dataset))

  # Cluster and Run Umap
  sc_dataset <- FindNeighbors(sc_dataset, dims = 1:10)
  sc_dataset <- FindClusters(sc_dataset, resolution = 0.5)
  sc_dataset <- RunUMAP(sc_dataset, dims = 1:10)
  umap <- DimPlot(sc_dataset, reduction = "umap")
  png(paste0(results_save_dir, "/umap.png"))
  print(umap)
  dev.off()
  # find marker genes for each cluster
  sc.markers <- FindAllMarkers(sc_dataset, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25) %>% group_by(cluster) %>% top_frac(select_frac, wt = avg_log2FC) %>% top_frac(-select_frac, wt = p_val)
  write.csv(sc.markers, paste0(results_save_dir,  "/sc_diff_expr.csv"))

  feature_map <- FeaturePlot(sc_dataset, features = sc.markers$gene[1:9])
  png(paste0(results_save_dir, "/feature_map.png"))
  print(feature_map)
  dev.off()

  saveRDS(sc_dataset, paste0(results_save_dir, "/", sc_name, ".rds"))

}

