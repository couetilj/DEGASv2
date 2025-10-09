# zscore normalization
normFunc <- function(x){return((x-mean(x, na.rm = T))/(sd(x, na.rm = T)+1e-3))}

# scaling from 0-1
scaleFunc <- function(x){return((x- min(x, na.rm = T)) /(max(x, na.rm = T)-min(x, na.rm = T)+1e-3))}

# Preprocess count data
normalizeScale <-function(X){
  return(t(apply(t(apply(as.matrix(t(X)),1,normFunc)),1,scaleFunc)))
}

preprocessCounts <- function(X){
  return(normalizeScale(1.5^log2(X+1)))
}

# Calculate UMAP coordinate
umap_coordinate <- function(count, metadata = NULL, min.cells = 3, min.features = 200, nfeatures = 2000, dims = 1:10, resolution = 0.5) {

  obj <- CreateSeuratObject(
    counts = as.matrix(count),
    meta.data = metadata,
    min.cells = min.cells,
    min.features = min.features
  )

  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, nfeatures = nfeatures)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj)
  obj <- FindNeighbors(obj, dims = dims)
  obj <- FindClusters(obj, resolution = resolution)
  obj <- RunUMAP(obj, dims = dims)

  umap_df <- Embeddings(obj, "umap") %>% as.data.frame()
  umap_df$cell <- rownames(umap_df)
  colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")

  return(umap_df)
}

# gene selection
select_genes <- function(scdata, sclab, patdata, phenotype, add_genes = NULL, bulk_hvg = TRUE, bulk_de = TRUE, sc_de = TRUE,
                         n_hvg = 250, n_bulk_de = 250, n_sc_de = 200, padj.thresh = 0.05, model_type = "category") {
  genes <- c()

  # Bulk HVG
  if (bulk_hvg) {
    common_genes <- rownames(patdata)
    gene_std_df <- data.frame(
      gene_names = common_genes,
      stdev = apply(patdata[common_genes, ], 1, sd, na.rm = TRUE)
    )
    gene_std_df <- gene_std_df[order(-gene_std_df$stdev), ]
    high_var_genes <- gene_std_df$gene_names[1:min(n_hvg, nrow(gene_std_df))]

    genes <- union(genes, as.character(high_var_genes))
  }

  # Bulk DE
  if (bulk_de && !is.null(phenotype)) {
    if (model_type == "survival") {
      if (!all(c("time","status") %in% colnames(phenotype))) {
        stop("For survival mode, phenotype must be a data.frame with columns 'time' and 'status'")
      }
    } else {
      phenotype <- as.vector(phenotype)
    }

    if (model_type == "survival") {
      surv_mid <- median(phenotype$time)
      patLab <- (phenotype$time < surv_mid) * phenotype$status + 1
    } else {
      patLab <- phenotype + 1
    }
    patLab <- as.factor(patLab)

    tryCatch({
      bulk_counts <- patdata
      bulk_counts[bulk_counts < 0] <- 0
      bulk_counts <- apply(bulk_counts, c(1, 2), as.integer)

      dds <- DESeqDataSetFromMatrix(
        countData = bulk_counts,
        colData   = data.frame(id = colnames(patdata), label = patLab),
        design    = ~label
      )
      dds <- DESeq(dds)
      res <- na.omit(results(dds))
      res <- res[order(res$padj), ]
      res <- res[res$padj < padj.thresh, ]

      high_diff_genes <- rownames(res)[1:min(n_bulk_de, nrow(res))]
      genes <- union(genes, as.character(high_diff_genes))
    }, error = function(e) {
      message("Bulk DE gene selection failed: ", e$message)
    })
  }

  # 3. SC/ST cluster DE
  if (sc_de && !is.null(scdata) && !is.null(sclab)) {
    obj <- CreateSeuratObject(counts = as.matrix(scdata), meta.data = sclab)
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj)
    obj <- ScaleData(obj)
    obj <- RunPCA(obj)
    obj <- FindNeighbors(obj, dims = 1:10)
    obj <- FindClusters(obj, resolution = 0.5)

    sc.markers <- FindAllMarkers(
      obj, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25
    )
    sc.markers <- na.omit(sc.markers)

    sc_genes <- sc.markers %>%
      arrange(p_val_adj, desc(avg_log2FC)) %>%
      head(n_sc_de) %>%
      pull(gene)

    genes <- union(genes, as.character(unique(sc_genes)))
  }

  # whether to add individual gene list
  if (!is.null(add_genes)) {
    add_genes <- as.character(add_genes)
    genes <- union(genes, intersect(add_genes, rownames(patdata)))
  }
  genes <- genes[!is.na(genes)]
  return(as.character(genes))
}

DEGAS_preprocessing <- function(
    scst_list, patdata, phenotype, sclab = NULL,
    bulk_hvg = TRUE, bulk_de = TRUE, sc_de = TRUE, add_genes = NULL,
    n_hvg = 250, n_bulk_de = 250, n_sc_de = 200,
    padj.thresh = 0.05, model_type = "category") {

  # common genes
  common_genes <- rownames(patdata)
  if (!is.list(scst_list) || inherits(scst_list, "data.frame")) {
    scst_list <- list(as.matrix(scst_list))
  }
  for (i in seq_along(scst_list)) {
    common_genes <- intersect(common_genes, rownames(scst_list[[i]]))
  }

  patdata <- patdata[common_genes, , drop = FALSE]

  scst_list <- lapply(scst_list, function(x) {
    x <- as.matrix(x)
    x[common_genes, , drop = FALSE]
  })

  # Gene selection

  gene_list <- select_genes(
    scdata      = scst_list[[1]],
    sclab       = sclab,
    patdata     = patdata,
    phenotype   = phenotype,
    add_genes   = add_genes,
    bulk_hvg    = bulk_hvg,
    bulk_de     = bulk_de,
    sc_de       = sc_de,
    n_hvg       = n_hvg,
    n_bulk_de   = n_bulk_de,
    n_sc_de     = n_sc_de,
    padj.thresh = padj.thresh,
    model_type = model_type
  )

  # Normalization
  norm_out <- normalize_counts_with_selected_genes(
    bulk_dataset = patdata,
    scst_list    = scst_list,
    gene_list    = gene_list
  )

  # clean phenotype
  if (model_type != "survival") {
    phenotype <- as.factor(phenotype)
    phenotype <- as.integer(phenotype) - 1
  }

  if (!is.null(sclab)) {
    sclab <- as.integer(as.factor(sclab)) - 1
    cat("sclab: ", unique(sclab), "\n")
  }

  return(list(
    patDat    = norm_out$patDat,
    phenotype = phenotype,
    scstDat   = norm_out$scstDat,
    scstName  = norm_out$scstName,
    sclab     = sclab
  ))
}



normalize_counts_with_selected_genes <- function(bulk_dataset, scst_list, gene_list) {
  # normalize bulk
  bulk_dataset <- bulk_dataset[gene_list, , drop = FALSE]
  patDat <- preprocessCounts(bulk_dataset)

  scst_expr_mat <- NULL
  scst_names <- c()

  # SC/ST datasets
  if (!is.list(scst_list)) {
    scst_list <- list(scst_list)
  }

  for (i in seq_along(scst_list)) {
    st_counts <- scst_list[[i]][gene_list, , drop = FALSE]
    st_counts <- preprocessCounts(as.matrix(st_counts))

    scst_expr_mat <- cbind(scst_expr_mat, st_counts)
    scst_names <- c(scst_names, rep(paste0("Dataset", i), ncol(st_counts)))
  }

  return(list(
    patDat   = patDat,
    scstDat  = scst_expr_mat,
    scstName = scst_names
  ))
}

