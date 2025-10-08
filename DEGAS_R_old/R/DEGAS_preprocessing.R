# normalization functions from orig DEGAS
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

normalize_counts_with_selected_genes <- function(bulk_dataset, phenotype, st_list, gene_list, save_name, results_dir, st_lab_list = c()) {
  # create results save dir
  if (!file.exists(results_dir)) {
    dir.create(results_dir)
  }
  patDat <- preprocessCounts(bulk_dataset[gene_list,])
  st_expr_mat <- c()
  st_meta_mat <- c()
  st_index_mat <- c() # indicator of which st or sc the sample comes from
  st_names_list <- c() # list of st or sc names, correspond to st_index_mat, only useful when merging multiple st or sc datasets
  for (i in c(1:length(st_list))) {
    st_dataset <- st_list[[i]]
    st_counts <- GetAssayData(st_dataset, slot = "counts", assay = "RNA")
    st_counts <- preprocessCounts(as.matrix(st_counts[gene_list, ]))
    st_expr_mat <- rbind(st_expr_mat, st_counts)
    tryCatch(
      {
        st_coords <- GetTissueCoordinates(st_dataset)
        meta_data <- st_dataset@meta.data
        meta_data <- cbind(meta_data, st_coords)
        st_dataset@meta.data <- meta_data
      }, error = function(e) {
        cat("No Coordinate exist, save meta data")
      }
    )
    st_meta <- st_dataset@meta.data
    st_meta_mat <- rbind(st_meta_mat, st_meta)
    st_index_mat <- append(st_index_mat, rep(i - 1, dim(st_counts)[1]))
    st_name <- as.character(st_dataset@meta.data$orig.ident[1])
    st_names_list <- append(st_names_list, st_name)
  }
  # create sc or st label if exist
  if (length(st_lab_list) > 0) {
    if (is.numeric(st_lab_list) == FALSE) {
      strings_factor <- factor(st_lab_list)
      unique_strings <- levels(strings_factor)
      st_lab_map <- setNames(unique_strings, as.character((1:length(unique_strings)) - 1))
      st_lab_mat <- as.integer(strings_factor) - 1
    } else {
      if (min(st_lab_list) == 1) {
        st_lab_mat <- st_lab_list - 1
      } else {
        st_lab_mat <- st_lab_list
      }
      st_lab_map <- sort(unique(st_lab_mat))
    }
  } else {
    st_lab_mat <- st_index_mat # which slides the sample comes from
    st_lab_map <- st_names_list
  }

  save(patDat, phenotype, st_expr_mat, st_meta_mat, st_lab_mat, st_lab_map, st_index_mat, st_names_list, file = paste0(results_dir, "/", save_name, ".RData"))
}

DEGAS_preprocessing <- function(bulk_dataset, phenotype, st_list, results_dir, st_lab_list = c(), model_type = "survival", diff_expr_files = c(), n_high_var_genes = 250, n_from = 500, n_by = 50, n_until = 250, n_high_diff_genes = 250, padj.thresh = 0.05) {
  # step 0. create folder save dir
  if (!file.exists(results_dir)) {
    dir.create(results_dir)
  }
  if (! is.list(st_list)) {
    st_list <- list(st_list)
  }

  # Step 1. get common genes across patients and st datasets
  common_genes <- rownames(bulk_dataset)
  for (i in c(1:length(st_list))) {
    st_counts <- GetAssayData(st_list[[i]], slot = "counts", assay = "RNA")
    common_genes <- intersect(common_genes, rownames(st_counts))
  }
  # Step 2. get high var genes in patients
  gene_std_df <- data.frame(gene_names = common_genes, stdev = apply(bulk_dataset[common_genes, ], 1, sd))
  gene_std_df %<>% arrange(desc(stdev))
  high_var_gene <- gene_std_df$gene_names[1:n_high_var_genes]
  normalize_counts_with_selected_genes(bulk_dataset, phenotype, st_list, high_var_gene, "Pat_std", results_dir, st_lab_list)
  # Step 3. select high var sc or st genes
  if (length(diff_expr_files) > 0) {
    st_var_gene_list <- list()
    for (diff_expr_file in diff_expr_files) {
      high_var_st_gene <- read.csv(diff_expr_file)
      st_var_gene_list <- union(st_var_gene_list, intersect(common_genes, high_var_st_gene$gene))
    }
    st_var_gene_list <- unlist(st_var_gene_list, recursive = TRUE)
    while (length(st_var_gene_list) > n_until) {
      st_var_gene_list <- intersect(st_var_gene_list, gene_std_df$gene_names[1:n_from])
      n_from <- n_from - n_by
    }
    normalize_counts_with_selected_genes(bulk_dataset, phenotype, st_list, st_var_gene_list, "SDE", results_dir, st_lab_list)
  }
  # Step 4. select diff expr genes in patients
  if (model_type == "survival") {
    surv_mid <- median(phenotype$time)
    patLab <- (phenotype$time < surv_mid) * phenotype$status + 1
  } else {
    patLab <- phenotype + 1
  }
  patLab <- as.factor(patLab)
  # creat DEseq object (This step may takes longer...)
  tryCatch(
    {
      bulk_counts <- bulk_dataset[common_genes, ]
      bulk_counts[bulk_counts < 0] <- 0
      bulk_counts <- apply(bulk_counts, c(1, 2), as.integer)
      dds <- DESeqDataSetFromMatrix(countData = bulk_counts,
                                    colData = data.frame(id = colnames(bulk_dataset), label = patLab), design = ~label)
      dds <- DESeq(dds)
      results <- na.omit(results(dds))
      results <- results[order(results$padj), ]
      results <- results[results$padj < padj.thresh, ]
      high_diff_genes <- rownames(results)[1:min(n_high_diff_genes, length(results$padj))]
      normalize_counts_with_selected_genes(bulk_dataset, phenotype, st_list, high_diff_genes, "Pat_Diff", results_dir, st_lab_list)
    }, error = function(e) {
      cat("There are some error in finding DEG. We terminate this step.")
    }
  )
}
