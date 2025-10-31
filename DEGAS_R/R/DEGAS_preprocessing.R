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
#' scRNA-seq UMAP calculation
#'
#' Using Seurat to calculate UMAP corrdinate for scRNA-seq data for visualization
#'
#' @param count A count matrix or data frame with scRNA-seq data (rows = genes, columns = cells).
#' @param metadata metadata for scRNA-seq, default = NULL
#' @param min.cells filtering cells, default = 3
#' @param min.features filtering genes, default = 200.
#' @param nfeatures The number of high variable genes, default = 2000
#' @param dims Dimensions of reduction to use as input
#' @param resolution Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.
#'
#' @return A UMAP coordinate data.frame
#' @export
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

#' Gene selection for DEGAS preprocessing
#'
#' Selects genes from bulk and single-cell data using high variability,
#' differential expression, or user-defined lists.
#'
#' @param scdata A scRNA-seq count matrix (genes × cells).
#' @param sclab Metadata for single-cell data.
#' @param patdata Bulk RNA-seq count matrix (genes × samples).
#' @param phenotype Phenotype vector or data.frame (depends on model_type).
#' @param add_genes Optional vector of user-specified genes.
#' @param bulk_hvg Logical. Use bulk HVG selection.
#' @param bulk_de Logical. Use bulk differential expression.
#' @param sc_de Logical. Use scRNA-seq differential expression.
#' @param n_hvg,n_bulk_de,n_sc_de Integers. Numbers of genes to select per method.
#' @param padj.thresh FDR threshold for DE.
#' @param model_type "category" or "survival".
#'
#' @return Character vector of selected gene names.
#' @export
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


#' DEGAS preprocessing
#'
#' Preprocess 4 necessary input data, scData, sclab, Patdata and phenotype with gene selection and normalization.
#'
#' @param scst_list Count matrix or data frame with scRNA-seq/spatial transcriptomics data (rows = genes, columns = cells).
#' @param sclab Metadata for scRNA-seq e.g. cell type, clusters (required), default = NULL
#' @param patdata Count Bulk RNA-seq data (genes × samples).
#' @param phenotype A phenotype vector or data frame describing sample labels as disease status.
#'   If \code{model_type = "category"}, this should be a categorical variable.
#'   If \code{model_type = "survival"}, it must include \code{time} and \code{status} columns.
#' @param bulk_hvg Whether user wants to use sample (bulk) level high variable genes as input in DEGAS. Default is \code{TRUE}.
#' @param bulk_de Whether user wants to use sample (bulk) differential expression genes for phenotype using DEseq2, commonly set Normal or Control as reference. Default is \code{TRUE}.
#' @param sc_de Whether user wants to use single cell differential expression genes for each cluster. Default is \code{TRUE}.
#' @param add_genes The user decides whether to add their own gene list by preference. If the user set bulk_hvg, bulk_de and sc_de as FALSE, this gene list is necessary. Default is \code{NULL}.
#' @param n_hvg Number of bulk high variable genes, default = 250.
#' @param n_bulk_de Number of bulk differential expression genes for phenotype, default = 250.
#' @param n_sc_de Number of genes for each single cell clusters, default = 200.
#' @param padj.thresh Threshold for differential expression analysis, default = 0.05
#' @param model_type Choose model type (category or survival), default = \code{category}.
#'
#' @return A list ready for DEGAS analysis, includes SC/ST RNA-seq data, SC metadata, bulk RNA-seq data, bulk sample phenotype, SC/ST datalist name if multiple datasets.
#' @export
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


#' Normalization
#'
#' This function selects genes from bulk and single-cell/spatial transcriptomics data
#' using high-variability, differential expression, or user-specified gene lists.
#'
#' @param scdata A count matrix or data frame containing single-cell RNA-seq
#'   expression data (rows = genes, columns = cells).
#' @param sclab A metadata data frame for the single-cell data (e.g., cluster or cell type labels).
#' @param patdata A count matrix or data frame containing bulk RNA-seq
#'   expression data (rows = genes, columns = samples).
#' @param phenotype A phenotype vector or data frame describing sample labels.
#'   If \code{model_type = "category"}, this should be a categorical variable.
#'   If \code{model_type = "survival"}, it must include \code{time} and \code{status} columns.
#' @param add_genes An optional character vector of user-specified genes to include.
#' @param bulk_hvg Logical. Whether to select bulk high-variable genes (HVG).
#'   Default is \code{TRUE}.
#' @param bulk_de Logical. Whether to include bulk differential expression (DE) genes.
#'   Default is \code{TRUE}.
#' @param sc_de Logical. Whether to include single-cell differential expression (DE) genes.
#'   Default is \code{TRUE}.
#' @param n_hvg Integer. Number of high-variable bulk genes to select. Default is \code{250}.
#' @param n_bulk_de Integer. Number of bulk differential expression genes to select. Default is \code{250}.
#' @param n_sc_de Integer. Number of single-cell differential expression genes to select. Default is \code{200}.
#' @param padj.thresh Numeric. Adjusted p-value threshold for differential expression filtering. Default is \code{0.05}.
#' @param model_type Character. Either \code{"category"} or \code{"survival"}. Determines phenotype handling.
#'
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

