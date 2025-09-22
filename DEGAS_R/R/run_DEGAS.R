run_DEGAS_SCST1 <- function(data_list, model_type, data_name, loss_type, transfer_type, model_save_dir,
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
    if (length(unique(phenotype[,1])) == 2) {
      phenotype <- phenotype[, c(2, 1)]
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

  if (!file.exists(model_save_dir)) {
    dir.create(model_save_dir)
  }

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







