import numpy
import numpy as np
import scipy.stats
import scanpy as sc

from fisheye.utils import get_logger


log = get_logger(__name__)


def func_pcc(cell_phalanx_one, cell_phalanx_two):
    pcc_value = scipy.stats.pearsonr(cell_phalanx_one, cell_phalanx_two)
    pcc_score = abs(pcc_value[0])
    return pcc_score


def simi_pcc(adata, label_kcg, label_kc):
    sum_cell_score = 0
    for cell in range(adata.X.shape[0]):
        kc_phalanx = np.mean(adata.X[np.where(label_kc == label_kc[cell])[0], :], axis=0)
        kcg_phalanx = np.mean(adata.X[np.where(label_kcg == label_kcg[cell])[0], :], axis=0)
        sum_cell_score += func_pcc(kc_phalanx, kcg_phalanx)
    return sum_cell_score


def reading(path):
    log.info("Reading single cell expression matrix.")
    if path.endswith(".h5"):
        adata = sc.read_h5ad(path)
    else:
        adata = sc.read_text(path)
    adata.var_names_make_unique()
    return(adata)


def preprocessing(adata):
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    sc.pp.log1p(adata)
    return(adata)


def clustering(adata, n_pcs=19, alg="louvain"):
    if adata.var.shape[0] >= n_pcs:
        sc.tl.pca(adata, svd_solver='arpack')
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=19)
    else:
        sc.pp.neighbors(adata, n_neighbors=10)
    getattr(sc.tl, alg)(adata)
    clustered_cell = sc.get.obs_df(adata, keys=[alg]) 
    labels = np.array(clustered_cell[alg])
    return labels

