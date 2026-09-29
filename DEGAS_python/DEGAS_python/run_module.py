import os
import pandas as pd
import numpy as np
from glob import glob
from .models import *
from .datasets import load_datasets
import time
from tqdm import tqdm
from .oof import patient_group_folds, write_calibration_outputs

# random seed for reproducibility
import random
seed = 42
torch.manual_seed(seed)
torch.cuda.manual_seed(seed)
torch.cuda.manual_seed_all(seed)  # if you are using multi-GPU.
np.random.seed(seed)  # Numpy module.
random.seed(seed)  # Python random module.
torch.manual_seed(seed)
torch.backends.cudnn.benchmark = False
torch.backends.cudnn.deterministic = True

def run_model(opt, pat_expr_mat, pat_lab_mat, sc_expr_mat, sc_loc_mat = None, sc_lab_mat = None,
              pat_eval_expr_mat = None, pat_eval_lab_mat = None, pat_eval_ids = None,
              heldout_group = None, return_model = False):
    """Train one DEGAS model and write its usual outputs.

    With return_model=True also return the fitted model for direct prediction
    and SHAP: ``directory, model = run_model(...)``. The default still returns
    only the directory. Call model.set_evaluate_mode() before inference.
    """
    # define data loaders and model
    high_reso_loader, low_reso_loader = load_datasets("train", opt, pat_expr_mat, pat_lab_mat, sc_expr_mat, sc_loc_mat, sc_lab_mat)
    eval_expr = pat_expr_mat if pat_eval_expr_mat is None else pat_eval_expr_mat
    eval_label = pat_lab_mat if pat_eval_lab_mat is None else pat_eval_lab_mat
    high_reso_eval_loader, low_reso_eval_loader = load_datasets("eval", opt, eval_expr, eval_label, sc_expr_mat, sc_loc_mat, sc_lab_mat)

    # define the model
    first_item = next(iter(low_reso_eval_loader))
    opt["input_shape"] = first_item["data"].shape[1]
    # check optionals for the model
    if opt["high_reso_output_shape"] < 0: # not specified yet
        raise ValueError("Please specify the number of single cell or tissue types in options (high_reso_output_shape = your number).")
    model = load_models(opt)

    # train and evaluate the model
    epoch = 1
    for high_reso_data, low_reso_data in tqdm(zip(high_reso_loader, low_reso_loader)):
        # print("Training Epoch {}".format(epoch))

        model.set_input(high_reso_data, low_reso_data) # notice that when training, the data shape will be 1 X batch_size X feature size, we will squeeze it in our backend code
        model.optimize_parameters(epoch)
    
        if (epoch % opt["save_freq"] == 0):    
            if opt["is_save"]:
                print("Saving the model...")
                model.save_networks(epoch)
            
            # print("Evaluate the model...")
            model.set_evaluate_mode()
            if opt["graph_type"] is None:
                high_reso_results, high_reso_embs = model.linear_eval(high_reso_eval_loader, opt["extract_embs"])
                low_reso_results, low_reso_embs = model.linear_eval(low_reso_eval_loader, opt["extract_embs"])
            else:
                raise ValueError("Not Implemented yet")

            if opt["extract_embs"]:
                np.save(os.path.join(model.save_dir, "high_reso_embs_epoch_{}.npy".format(epoch)), high_reso_embs)
                np.save(os.path.join(model.save_dir, "low_reso_embs_epoch_{}.npy".format(epoch)), low_reso_embs)
            high_reso_results.to_csv(os.path.join(model.save_dir, "high_reso_results_epoch_{}.csv".format(epoch)))
            low_reso_results.to_csv(os.path.join(model.save_dir, "low_reso_results_epoch_{}.csv".format(epoch)))
            if pat_eval_ids is not None:
                local_pid = low_reso_results["pid"].astype(int).to_numpy()
                patient_oof = low_reso_results.copy()
                patient_oof["patient_id"] = np.asarray(pat_eval_ids).reshape(-1)[local_pid]
                labels = np.asarray(eval_label)
                if labels.ndim > 1 and labels.shape[1] > 1:
                    labels = np.argmax(labels, axis=1)
                patient_oof["label"] = labels.reshape(-1)[local_pid]
                patient_oof["heldout_group"] = str(heldout_group)
                patient_oof.rename(columns={"hazard": "raw_probability"}, inplace=True)
                patient_oof.to_csv(os.path.join(model.save_dir, "patient_oof_results_epoch_{}.csv".format(epoch)), index=False)
            # print("Back to Training phase")
            if opt["is_save"]:
                model.load_networks(epoch)
            model.set_train_mode()  
        epoch += 1   

    model.loss_rec.to_csv(os.path.join(model.save_dir, "losses.csv".format(epoch)))
    return (model.save_dir, model) if return_model else model.save_dir


def _aggregate_high_resolution(save_results_folder, opt):
    results_file_list = glob(os.path.join(save_results_folder, "*", "high_reso_results_epoch_{}.csv".format(opt["tot_iters"])))
    if not results_file_list:
        raise RuntimeError("No high-resolution DEGAS results were produced")
    results = [pd.read_csv(results_file, index_col = 0, header = 0) for results_file in results_file_list]
    for i, results_file in enumerate(results_file_list):
        meta_info = results_file.split("/")[-2].split("_")
        results[i]["fold"] = int(meta_info[1])
        results[i]["seed"] = int(meta_info[4])
    results = pd.concat(results, ignore_index = True)
    results.to_csv(os.path.join(save_results_folder, "summary.csv"))
    results_mean = results.groupby("index").mean().reset_index()
    results_mean.to_csv(os.path.join(save_results_folder, "summary_mean.csv"))
    return results_mean


def bagging_patient_oof(opt, pat_expr_mat, pat_lab_mat, sc_expr_mat,
                        pat_groups, pat_ids = None, sc_loc_mat = None, sc_lab_mat = None):
    labels = np.asarray(pat_lab_mat)
    flat_labels = np.argmax(labels, axis=1) if labels.ndim > 1 and labels.shape[1] > 1 else labels.reshape(-1)
    groups = np.asarray(pat_groups).reshape(-1)
    patient_ids = np.arange(len(flat_labels)) if pat_ids is None else np.asarray(pat_ids).reshape(-1)
    if len(patient_ids) != len(flat_labels):
        raise ValueError("pat_ids must align with patient labels")
    folds = patient_group_folds(flat_labels, groups)
    save_results_folder = None
    for fold_info in folds:
        opt["fold"] = int(fold_info["fold"])
        opt["patient_oof_group"] = str(fold_info["heldout_group"])
        train = fold_info["train_indices"]
        test = fold_info["test_indices"]
        for seed_value in range(opt["tot_seeds"]):
            opt["seed"] = seed_value
            print("Run patient OOF group {} fold {} submodel {}...".format(
                fold_info["heldout_group"], fold_info["fold"], seed_value))
            train_expr = pat_expr_mat[train]
            eval_expr = pat_expr_mat[test]
            high_expr = sc_expr_mat
            if "random_feat" in opt and opt["random_feat"] and "random_perc" in opt:
                np.random.seed(seed_value)
                num_select_feats = np.floor(pat_expr_mat.shape[1] * opt["random_perc"]).astype(int)
                select_feats = np.sort(np.random.choice(list(range(pat_expr_mat.shape[1])), num_select_feats, replace=False))
                train_expr = train_expr[:, select_feats]
                eval_expr = eval_expr[:, select_feats]
                high_expr = high_expr[:, select_feats]
            save_results_folder = run_model(
                opt,
                train_expr,
                labels[train],
                high_expr,
                sc_loc_mat,
                sc_lab_mat,
                pat_eval_expr_mat=eval_expr,
                pat_eval_lab_mat=labels[test],
                pat_eval_ids=patient_ids[test],
                heldout_group=fold_info["heldout_group"],
            )
    save_results_folder = os.path.dirname(save_results_folder)
    oof_files = glob(os.path.join(save_results_folder, "*", "patient_oof_results_epoch_{}.csv".format(opt["tot_iters"])))
    if len(oof_files) != len(folds) * opt["tot_seeds"]:
        raise RuntimeError("Patient OOF files are incomplete")
    submodels = pd.concat([pd.read_csv(path) for path in oof_files], ignore_index=True)
    submodels.to_csv(os.path.join(save_results_folder, "patient_oof_submodels.csv"), index=False)
    ensemble = submodels.groupby(
        ["patient_id", "label", "heldout_group"], as_index=False
    ).agg(
        raw_probability=("raw_probability", "mean"),
        raw_probability_sd=("raw_probability", "std"),
        n_submodels=("raw_probability", "size"),
    )
    if len(ensemble) != len(flat_labels) or not (ensemble.n_submodels == opt["tot_seeds"]).all():
        raise RuntimeError("Patient OOF ensemble coverage is incomplete")
    ensemble.to_csv(os.path.join(save_results_folder, "patient_oof_ensemble.csv"), index=False)
    write_calibration_outputs(
        save_results_folder,
        ensemble,
        penalty=float(opt.get("calibration_penalty", 1.0)),
        bins=int(opt.get("calibration_bins", 5)),
    )
    return _aggregate_high_resolution(save_results_folder, opt)


def bagging_all_results(opt, pat_expr_mat, pat_lab_mat, sc_expr_mat, sc_loc_mat = None, sc_lab_mat = None,
                        pat_groups = None, pat_ids = None):
    """
    sc_expr_mat: single cell or spatial transcriptomic data gene expression
    sc_loc_mat: spatial transcriptomic data (optional, for graph NN)
    sc_lab_mat: single cell labels (optional)
    pat_expr_mat: patient gene expression value
    pat_lab_mat: patient labels
    """
    if opt.get("patient_oof", False):
        if pat_groups is None:
            raise ValueError("patient_oof requires pat_groups")
        return bagging_patient_oof(
            opt, pat_expr_mat, pat_lab_mat, sc_expr_mat, pat_groups, pat_ids, sc_loc_mat, sc_lab_mat
        )
    if opt["tot_folds"] == 1:
        for seed in range(opt["tot_seeds"]):
            opt["seed"] = seed
            print("Run submodel {}...".format(seed))
            if "random_feat" in opt.keys() and opt["random_feat"] and "random_perc" in opt.keys():
                np.random.seed(opt["seed"])
                num_select_feats = np.floor(pat_expr_mat.shape[1] * opt["random_perc"]).astype(int)
                select_feats = np.sort(np.random.choice(list(range(pat_expr_mat.shape[1])), num_select_feats, replace = False))
                save_results_folder = run_model(opt, pat_expr_mat[:, select_feats], pat_lab_mat, sc_expr_mat[:, select_feats], sc_loc_mat, sc_lab_mat)
            else:
                save_results_folder = run_model(opt, pat_expr_mat, pat_lab_mat, sc_expr_mat, sc_loc_mat, sc_lab_mat)
    else:
        for fold in range(opt["tot_folds"]):
            opt["fold"] = fold
            for seed in range(opt["tot_seeds"]):
                opt["seed"] = seed
                print("Run fold {} submodel {}...".format(fold, seed))
                if "random_feat" in opt.keys() and opt["random_feat"] and "random_perc" in opt.keys():
                    np.random.seed(opt["seed"])
                    num_select_feats = np.floor(pat_expr_mat.shape[1] * opt["random_perc"]).astype(int)
                    select_feats = np.sort(np.random.choice(list(range(pat_expr_mat.shape[1])), num_select_feats, replace = False))
                    save_results_folder = run_model(opt, pat_expr_mat[:, select_feats], pat_lab_mat, sc_expr_mat[:, select_feats], sc_loc_mat, sc_lab_mat)
                else:
                    save_results_folder = run_model(opt, pat_expr_mat, pat_lab_mat, sc_expr_mat, sc_loc_mat, sc_lab_mat)
    print("Finish Run and Eval all models")
    print("Aggregate all results")
    # aggregate all results
    save_results_folder = os.path.dirname(save_results_folder) # get the parent folder which include all submodules
    return _aggregate_high_resolution(save_results_folder, opt)


    
    
