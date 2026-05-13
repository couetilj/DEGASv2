from .dataset import *
# from .graph_dataset import *
# from torch_geometric.loader import DataLoader as GDataLoader
from torch.utils.data import DataLoader

def load_datasets(phase, opt, pat_expr_mat, pat_lab_mat, st_expr_mat, st_loc_mat = None, sc_lab_mat = None, pat_groups = None):
    """Build the high-resolution (SC/ST) and low-resolution (patient) DataLoaders.

    pat_groups : np.ndarray | None
        Per-sample patient identifier for the bulk training matrix. When provided
        and tot_folds > 1, the bulk fold split uses sklearn.GroupKFold so that
        matched T+N aliquots from the same patient stay in the same fold. None
        preserves the original behavior (every fold sees the full bulk set).
    """
    bs = 1 if phase == "train" else opt["batch_size"]
    # load single cell or st dataset
    if opt["graph_type"] is None:
        high_reso_dataset = DataLoader(STSCDataset(st_expr_mat, sc_lab_mat,
            random_seed = opt["seed"], fold = opt["fold"], tot_folds = opt["tot_folds"], tot_iters = opt["tot_iters"], batch_size = opt["batch_size"], phase = phase, sample_method = opt["sample_method"]),
            batch_size = bs, shuffle = False)
        low_reso_dataset = DataLoader(PatDataset(pat_expr_mat, pat_lab_mat, random_seed = opt["seed"], batch_size = opt["pat_batch_size"],
            tot_iters = opt["tot_iters"], phase = phase, model_type = opt["model_type"], sample_method = opt["sample_method"],
            fold = opt["fold"], tot_folds = opt["tot_folds"], groups = pat_groups),
            batch_size = bs, shuffle = False)
    else:
        print("Not implement yet")
        pass

    return high_reso_dataset, low_reso_dataset


