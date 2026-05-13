from torch.utils.data import Dataset
import pandas as pd
import numpy as np
import os
from .tools import *

try:
    from sklearn.model_selection import GroupKFold
except ImportError:
    GroupKFold = None


class STSCDataset(Dataset):
    # load spatial transcriptomics datasets
    # the logic of this Dataset is different from the pytorch
    # in this STDataset, each batch samples is saved as a single sample in pytorch Dataset
    # which means, each sample is batch_size X feature_dim
    def __init__(self, st_expr_mat, sc_lab_mat = None,
                random_seed = 0, fold = 0, tot_folds = 1, tot_iters = 300, batch_size = 200,
                phase = "train", sample_method = "balance", groups = None):
        """
        random_seed (int): random seed for shuffling the dataset
        sample_method: which method to sample
        tot_iters: total number of training iteration
        batch_size (int): batch size
        phase (string): should be 'train' or 'eval'.
        groups (np.ndarray | None): optional per-sample group identifier (e.g. patient
            barcode). When provided and tot_folds > 1, fold splits use
            sklearn.GroupKFold so that no single group spans two folds — required to
            avoid intra-patient data leakage in scRNA-seq / spatial data.
        """
        super(STSCDataset, self).__init__()
        assert phase in ["train", "eval", "eval_all"]
        self.n_samples = st_expr_mat.shape[0]
        if sc_lab_mat is not None:
            if len(sc_lab_mat.shape) > 1:
                if sc_lab_mat.shape[1] > 1:
                    sc_lab_mat = np.argmax(sc_lab_mat, axis = 1)
                sc_lab_mat = sc_lab_mat.flatten()
            assert self.n_samples == sc_lab_mat.shape[0]
            self.scLab = sc_lab_mat
        else:
            self.scLab = np.zeros((self.n_samples)).astype(np.int8) # dummy value

        self.stDat = st_expr_mat
        self.st_idx = np.array(range(self.n_samples))

        self.phase = phase
        self.tot_iters = tot_iters
        self.random_seed = random_seed
        self.fold = fold
        self.tot_folds = tot_folds
        self.sample_method  = sample_method
        self.batch_size = batch_size
        self.groups = groups

        # ZL: add Nov, 21st, for multiple folds
        if self.fold >= 0 and self.phase == "train" and self.tot_folds > 1:
            if groups is not None:
                # Patient-grouped split — each group ends up entirely in one fold.
                # Use GroupKFold's test_idx as the fold's training subset, matching
                # the existing semantics ('each fold gets its own slice of the data,
                # disjoint from other folds').
                if GroupKFold is None:
                    raise ImportError("groups requires scikit-learn (GroupKFold)")
                groups_arr = np.asarray(groups)
                if len(groups_arr) != self.n_samples:
                    raise ValueError(
                        f"groups length {len(groups_arr)} != n_samples {self.n_samples}"
                    )
                n_unique = len(np.unique(groups_arr))
                if n_unique < self.tot_folds:
                    raise ValueError(
                        f"GroupKFold needs at least tot_folds={self.tot_folds} groups, "
                        f"got {n_unique}"
                    )
                gkf = GroupKFold(n_splits=self.tot_folds)
                splits = list(gkf.split(np.zeros((self.n_samples, 1)),
                                         y=self.scLab, groups=groups_arr))
                _, fold_indices = splits[fold]
                self.fold_indices = np.asarray(fold_indices)
            else:
                # Original random-shuffle behavior (preserved for backward compat).
                np.random.seed(0)
                indices = np.random.choice(range(self.n_samples), self.n_samples, replace = False)
                fold_size = self.n_samples // self.tot_folds
                self.fold_indices = indices[fold * fold_size : (fold + 1) * fold_size]
            self.scSubLab = self.scLab[self.fold_indices]


    def __len__(self):
        if self.phase == "train":
            return self.tot_iters
        else:
            return self.n_samples

    def get(self, index):
        assert self.phase == "train"
        if self.sample_method == "balance":
            if self.fold >= 0 and self.tot_folds > 1:
                sub_idx_select = balance_sampling(self.scSubLab, self.batch_size, self.random_seed + index)
                self.idx_select = self.fold_indices[sub_idx_select]
            else:
                self.idx_select = balance_sampling(self.scLab, self.batch_size, self.random_seed + index)
        else:
            raise ValueError("Other sampling method is not implemented yet!")

    def __getitem__(self, index):
        if self.phase == "train":
            self.get(index)
        else:
            self.idx_select = index

        return {"index": self.idx_select, "data": self.stDat[self.idx_select, :], "label": self.scLab[self.idx_select]}
#####################################################################################

class PatDataset(Dataset):
    def __init__(self, pat_expr_mat, pat_lab_mat, random_seed = 0, batch_size = 200, tot_iters = 300, phase = "train", model_type = "phenotype", sample_method = "balance"):
        super(PatDataset, self).__init__()
        assert phase in ["train", "eval", "eval_all"]

        self.patDat = pat_expr_mat
        self.n_samples = pat_expr_mat.shape[0]
        assert self.n_samples == pat_lab_mat.shape[0]
        if ("Cox" not in model_type):
            if len(pat_lab_mat.shape) > 1:
                if (pat_lab_mat.shape[1] > 1):
                    pat_lab_mat = np.argmax(np.array(pat_lab_mat), axis = 1)
                pat_lab_mat = pat_lab_mat.flatten()
            self.patLab = pat_lab_mat
        else:
            self.time = np.array(pat_lab_mat[:, 0]).astype(int).flatten()
            self.status = np.array(pat_lab_mat[:, 1]).flatten()
        self.pat_idx = np.array(range(self.n_samples))

        self.random_seed = random_seed
        self.batch_size = batch_size
        self.tot_iters = tot_iters
        self.phase = phase
        self.model_type = model_type
        self.sample_method = sample_method



    def __len__(self):
        if self.phase == "train":
            return self.tot_iters
        else:
            return self.n_samples

    def get(self, index):
        assert self.phase == "train"
        if ("Cox" in self.model_type):
            np.random.seed(self.random_seed + index)
            self.idx_select = np.random.choice(self.pat_idx, self.batch_size, replace = False)
        elif (self.sample_method == "balance") and ("Cox" not in self.model_type):
            self.idx_select = balance_sampling(self.patLab, self.batch_size, self.random_seed + index)
        else:
            raise ValueError("Other sampling method is not implemented yet!")


    def __getitem__(self, index):
        if self.phase == "train":
            self.get(index)
        else:
            self.idx_select = index
        if ("Cox" in self.model_type):
            return {"pid": self.pat_idx[self.idx_select], "data": self.patDat[self.idx_select, :], "time": self.time[self.idx_select], "status": self.status[self.idx_select]}
        else:
            return {"pid": self.pat_idx[self.idx_select], "data": self.patDat[self.idx_select, :], "label": self.patLab[self.idx_select]}
