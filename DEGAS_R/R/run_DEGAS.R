#' Disease associated or risk score prediction using DEGAS.
#'
#' Run DEGAS with different model type.
#'
#' @param data_list A list ready for DEGAS analysis, includes SC/ST RNA-seq data (rows = cells, columns = genes) after preprocessing, SC metadata, bulk RNA-seq data (rows = samples, columns = genes) after preprocessing, bulk sample phenotype, SC/ST datalist name if multiple datasets.
#' @param model_type String specifying which DEGAS Python model configuration to call.
#'   Supported options include `BlankClass`, `ClassClass`, `BlankCox`, `ClassCox`, `BlankBCE` and `ClassBCE`.
#'   Models with the `Blank*` prefix indicate that single-cell/spatial data have no labels and only the patient branch is trained.
#'   Models with the `Class*` prefix incorporate single-cell labels to supervise the high-resolution branch.
#'   The suffix `Class` corresponds to classification tasks, `Cox` to survival analysis (log-negative-likelihood or ranking loss), and `BCE` to binary cross-entropy output.
#'
#' @param data_name Dataset name. Used to generate the output directory and log file prefixes,
#' @param loss_type String specifying the loss function used for training the low-resolution (patient) branch.
#'   Typically, classification models use `cross_entropy`, while survival models can use `log_neg` or `rank_loss`.
#'   This value must match the expected loss type of the chosen `model_type`, otherwise the backend will raise an error.
#'
#' @param transfer_type String indicating the cross-domain alignment method.
#'   The backend currently supports `Wasserstein` (adversarial with gradient penalty) and `MMD` (maximum mean discrepancy).
#'   This parameter determines how the transfer alignment loss is computed.
#'
#' @param model_save_dir String specifying the output directory for model files.
#'   The function will automatically create the directory if it does not exist and will save configuration files, prediction results, and optionally embedding vectors for each submodel.
#'
#' @param lambda1 Numeric value, a weighting coefficient that controls the contribution
#'   of the high-resolution branch (single-cell/spatial supervision loss).
#'   Only applicable to models with the `Class` prefix. Default = 1.
#'
#' @param lambda2 Numeric value, a weighting coefficient that controls
#'   the loss strength for the low-resolution (patient) branch.
#'   This parameter applies to all model types. Default = 2.
#'
#' @param lambda3 Numeric value, a weighting coefficient that controls
#'   the contribution of the domain-alignment (transfer) loss.
#'   Works jointly with `transfer_type` to regulate the alignment between
#'   high- and low-resolution embedding distributions. Default = 3.
#'
#' @param tot_seeds Integer specifying the number of bagging iterations. Default = 10.
#'
#' @param tot_iters Integer specifying the total number of training iterations (epochs).
#'
#' @param extract_embs Logical value. If `TRUE`, saves the high- and low-resolution hidden embeddings
#'   as `.npy` files at each save interval for downstream visualization or analysis. Default = `FALSE`.
#'
#' @param random_feat Logical value. If `TRUE`, randomly selects a subset of gene features
#'   before each submodel training. Default = `FALSE`.
#'
#' @param random_perc Numeric value between 0 and 1 specifying the proportion of features to retain
#'   when `random_feat` is enabled default = 0.8. The same proportion is applied to both patient-level and single-cell/spatial expression matrices.
#'
#' @param early_stopping Logical value. A reserved flag for early stopping;
#'
#' @return A data frame where the column `hazard` stores the predicted score for each cell or spot. For classification models, this represents the
#'   predicted disease associated score, and for Cox models, it represents the estimated risk value. The output also includes the original `st_lab_mat`
#'   metadata to facilitate downstream analyses.
#'
#' @export
run_DEGAS_SCST <- function(data_list, model_type, data_name, loss_type, transfer_type, model_save_dir,
                           lambda1 = 1.0, lambda2 = 3.0, lambda3 = 3.0, tot_seeds = 10, tot_iters = 300, extract_embs = FALSE, random_feat = FALSE, random_perc = 0.8, early_stopping = FALSE) {
  # load required packages
  numpy <- import("numpy")
  DEGAS_python <- import("DEGAS_python")

  st_expr_mat <- data_list$scstDat
  patDat      <- data_list$patDat
  phenotype   <- data_list$phenotype
  st_lab_mat  <- data_list$sclab

  # transform data input numpy array format
  st_expr_mat <- numpy$array(r_to_py(st_expr_mat))
  st_lab_mat <- numpy$array(r_to_py(st_lab_mat))
  pat_expr_mat <- numpy$array(r_to_py(patDat))


  if (grepl("Cox", model_type) | grepl("BCE", model_type)) {
    phenotype <- as.matrix(phenotype)
    if (length(unique(phenotype[, 1])) == 2) {
      phenotype <- phenotype[, c(2, 1), drop = FALSE]
    }
  }


  pat_lab_mat <- numpy$array(r_to_py(phenotype))

  if (model_type == "BlankClass") {
    opt <- DEGAS_python$BlankClass_opt
  } else if (model_type == "ClassClass") {
    opt <- DEGAS_python$ClassClass_opt
  } else if (model_type == "BlankCox") {
    opt <- DEGAS_python$BlankCox_opt
  } else if (model_type == "ClassCox"){
    opt <- DEGAS_python$ClassCox_opt
  } else if (model_type == "BlankBCE") {
    opt <- DEGAS_python$BlankBCE_opt
  } else if (model_type == "ClassBCE") {
    opt <- DEGAS_python$ClassBCE_opt
  }
  opt$data_name <- data_name
  opt$loss_type <- loss_type
  opt$transfer_type <- transfer_type
  opt$lambda1 <- lambda1
  opt$lambda2 <- lambda2
  opt$lambda3 <- lambda3
  opt$tot_seeds <- as.integer(tot_seeds)
  opt$tot_iters <- as.integer(tot_iters)
  opt$save_dir <- model_save_dir
  opt$extract_embs <- extract_embs
  opt$random_feat <- random_feat
  opt$random_perc <- random_perc
  opt$early_stopping <- early_stopping
  opt$high_reso_output_shape <- n_st_classes

  if (grepl("Cox", model_type) | grepl("BCE", model_type)) {
    opt$low_reso_output_shape <- 1L
  } else {
    opt$low_reso_output_shape <- length(unique(as.vector(phenotype)))
  }

  if (!file.exists(model_save_dir)) {
    dir.create(model_save_dir)
  }

#  cat("unique patient labels:", sort(unique(as.vector(phenotype))), "\n")
#  cat("range patient labels:", range(as.vector(phenotype)), "\n")

#  cat("unique sc labels:", sort(unique(st_lab_mat)), "\n")
#  cat("n_st_classes:", n_st_classes, "\n")


  # Run the model
  degas_results <- DEGAS_python$bagging_all_results(opt, pat_expr_mat, pat_lab_mat, st_expr_mat, sc_lab_mat = st_lab_mat)

  # Join the results with meta information
  degas_results <- cbind(degas_results, st_lab_mat)

  return(degas_results)
}

plot_hidden_feat <- function(folder_path, phenotype, random_seed = 0, fold = -1, epoch_from = 50, epoch_to = 300, epoch_by = 50, dis_label = "AD", NC_label = "NC") {
  library("ggplot2")
  library("reticulate")
  numpy <- import("numpy")
  sklearn <- import("sklearn.decomposition")
  for (epoch in seq(epoch_from, epoch_to, epoch_by)) {
    sub_folder_path <- paste0(folder_path, "/fold_", fold, "_random_seed_", random_seed)
    sc_hidden <- numpy$load(paste0(sub_folder_path, "/high_reso_embs_epoch_", epoch, ".npy"))
    pat_hidden <- numpy$load(paste0(sub_folder_path, "/low_reso_embs_epoch_", epoch, ".npy"))
    # calculate PCA
    pca <- sklearn$PCA()
    sc_PCA <- pca$fit_transform(sc_hidden)
    pat_PCA <- pca$transform(pat_hidden)
    sc_df <- data.frame(PC1 = sc_PCA[, 1], PC2 = sc_PCA[, 2], group = "Single Cells", label = "Unknown")
    pat_df <- data.frame(PC1 = pat_PCA[, 1], PC2 = pat_PCA[, 2], group = "Patients", label = phenotype)
    cb_df <- rbind(sc_df, pat_df)
    p <- ggplot(cb_df, aes(x = PC1, y = PC2, color = group, shape = label, size = group)) +
      geom_point() +
      labs(x = "PC 1", y = "PC 2", title = paste0("PCA Plot of Hidden Features (Epoch = ", epoch, ")")) +
      scale_color_manual(values = c("Single Cells" = "gray", "Patients" = "darkred")) +
      scale_alpha_manual(values = c("Group 1" = 0.1, "Group 2" = 1)) +
      # scale_shape_manual(values = c(dis_label = 16, NC_label = 17, "Unknown" = 15)) +
      scale_size_manual(values = c("Single Cells" = 1, "Patients" = 2.5)) +
      theme_minimal()
    ggsave(file = paste0(folder_path, "/fold_", fold, "_random_seed_", random_seed, "_epoch_", epoch, ".pdf"), plot = p)
  }
}







