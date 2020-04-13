#!/E:/Python-3.7.0(pycharm)/python
# -*- coding:UTF-8 -*-
# author: wangcong
# time: 2020/4/10
import fire
import numpy as np
from anndata import AnnData
from sklearn.metrics.cluster import adjusted_rand_score
from pathos.pools import ProcessPool as Pool
from .utils import reading, preprocessing, clustering
from .utils import get_logger
from .compare import read_clustering_result
import scanpy as sc
import pandas as pd


log = get_logger(__name__)


def take_sub(adata, G):
    G_ = list(G)
    var = adata.var.loc[G_]
    indexes = [adata.var.index.get_loc(g) for g in G_]
    X = adata.X[:, indexes]
    sub_adata = AnnData(X, adata.obs, var)
    return sub_adata


def take_and_clustering(adata, G, c_alg):
    sub_adata = take_sub(adata, G)
    labels = clustering(sub_adata, alg=c_alg)
    return labels


def get_G(adata, G_path, clustering_alg, compare_method='t-test'):
    G_set = set()
    if G_path:
        with open(G_path) as f:
            for line in f:
                G_set.add(line.strip())
    else:
        sc.tl.rank_genes_groups(adata, clustering_alg, method=compare_method)
        m_markers = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(10)
        markers = pd.concat([m_markers.iloc[i] for i in range(m_markers.shape[0])]).reset_index(drop=True).unique()
        G_set = set(markers)
    return G_set


class FindMarkers():
    """
    Find marker gene set G with a greedy algorithm,
    which optimize the similarity (Adjusted Rand Index) between
    clustering result using all genes (labels_kc) and
    clustering result only using G (labels_kc(G)).

    This algorithm iteratively:
        1. For each gene g which not in G,
           calculate the incremental value(inc):
           inc = ARI(labels_kc(G + {g}), labels_kc)
                 - ARI(labels_kc(G), labels_kc)
        2. Add the gene g_max to G, which has largest inc.
        3. For each gene g which in G,
           calculate the incremental value(inc):
           inc = ARI(labels_kc(G), labels_kc)
                 - ARI(labels_kc(G - {g}), labels_kc)
        4. If there are g has negative inc value,
           remove the g with smallest inc from G.

    Reference:
    [1] Qian, X. et al. Probabilistic cell typing enables fine mapping
        of closely related cell types in situ. Nat Methods 17, 101–106 (2020).
    """
    def __init__(self, adata, G, cycles_times=100,
                 clustering_alg="louvain", genes=None,
                 n_workers=1):
        self.adata = adata
        self.cycles_times = cycles_times
        self.map_ = map if n_workers <= 1 else Pool(n_workers).map
        self.G = G
        self.clustering_alg = clustering_alg
        self.labels_kc = np.array(adata.obs[clustering_alg])
        self.genes = set(list(adata.var.index))

    def cal_score(self, adata, G, labels_kc):
        if not G:
            return 0
        else:
            labels_G = take_and_clustering(adata, G, self.clustering_alg)
            score_G = adjusted_rand_score(labels_G, labels_kc)
            return score_G

    def run(self):
        for idx in range(self.cycles_times):
            log.info(f"cycle {idx}:")
            self._break_cnt = 0
            self._run_reduce_loop()
            if self._break_cnt >= 1:
                log.info("Can't improve anymore.")
                break

    def _run_reduce_loop(self):
        """Find smallest negative inc gene, remove from G."""
        inc2name = {}  # score inc -> gene name
        genes, G, adata, labels_kc, map_, cal_score = (
            self.genes, self.G, self.adata, self.labels_kc,
            self.map_, self.cal_score)
        score_G = cal_score(adata, G, labels_kc)
        cal_score_ = lambda g: (g, cal_score(adata, G - {g}, labels_kc))
        for g, score_G_g in map_(cal_score_, G):
            inc = score_G - score_G_g
            inc2name[inc] = g
        min_inc = min(inc2name)
        if min_inc <= 0:
            min_g = inc2name[min_inc]
            G.remove(min_g)
            genes.add(min_g)
            log.info(f"Min dec: {min_inc}, remove gene: {min_g}")
        else:
            self._break_cnt += 1


def main(path, G_path=None, outpath=None, n_workers=1, clustering_alg="louvain"):
    """De nove gene panel design using greedy algorithm."""
    adata = reading(path)
    adata = preprocessing(adata)
    if clustering_alg not in adata.obs:
        clustering(adata, alg=clustering_alg)
        out_h5 = ".".join(path.split('.')[:-1]) + '.clust.h5'
        log.info(f"Write clustering result to anndata: {out_h5}")
        adata.write_h5ad(out_h5)
    G_set = get_G(adata, G_path, clustering_alg)
    find_marker = FindMarkers(adata, G_set, 100, clustering_alg=clustering_alg)
    find_marker.run()
    if outpath:
        log.info(f"Write result to path {outpath}")
        with open(outpath, "w") as f:
            for g in list(find_marker.G):
                f.write(g + "\n")


if __name__ == "__main__":
    fire.Fire(main)
