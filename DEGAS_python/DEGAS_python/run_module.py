import os
import pandas as pd
import numpy as np
from glob import glob
from .models import *
from .datasets import load_datasets
import time
from tqdm import tqdm

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

def run_model(opt, pat_expr_mat, pat_lab_mat, sc_expr_mat, sc_loc_mat = None, sc_lab_mat = None):
    # define data loaders and model
    high_reso_loader, low_reso_loader = load_datasets("train", opt, pat_expr_mat, pat_lab_mat, sc_expr_mat, sc_loc_mat, sc_lab_mat)
    high_reso_eval_loader, low_reso_eval_loader = load_datasets("eval", opt, pat_expr_mat, pat_lab_mat, sc_expr_mat, sc_loc_mat, sc_lab_mat)

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
            # print("Back to Training phase")
            if opt["is_save"]:
                model.load_networks(epoch)
            model.set_train_mode()  
        epoch += 1   

    model.loss_rec.to_csv(os.path.join(model.save_dir, "losses.csv".format(epoch)))
    return model.save_dir


def bagging_all_results(opt, pat_expr_mat, pat_lab_mat, sc_expr_mat, sc_loc_mat = None, sc_lab_mat = None):
    """
    sc_expr_mat: single cell or spatial transcriptomic data gene expression
    sc_loc_mat: spatial transcriptomic data (optional, for graph NN)
    sc_lab_mat: single cell labels (optional)
    pat_expr_mat: patient gene expression value
    pat_lab_mat: patient labels
    """
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
    results_file_list = glob(os.path.join(save_results_folder, "*", "high_reso_results_epoch_{}.csv".format(opt["tot_iters"])))
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


    
    
