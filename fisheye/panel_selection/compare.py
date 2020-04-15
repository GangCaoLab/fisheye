import fire
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from .utils import reading, preprocessing, clustering
from .utils import get_logger


log = get_logger(__name__)


def read_clustering_result(path):
    res = []
    with open(path) as f:
        for line in f:
            res.append(line.strip())
    return res


def main(path, outpath=None, cluster_res=None,
         markers_per_group=30,
         clustering_alg="louvain", compare_method='t-test'):
    adata = reading(path)
    adata = preprocessing(adata)
    if cluster_res:
        clustered_cell = read_clustering_result(cluster_res)
        clustering_alg = "from_outside"
        adata.obs[clustering_alg] = clustered_cell
    else:
        clustered_cell = clustering(adata)
    log.debug(clustered_cell)
    sc.tl.rank_genes_groups(adata, clustering_alg, method=compare_method)
    sc.pl.rank_genes_groups(adata, n_genes=markers_per_group, sharey=False)
    outfig = f"marker_genes_{compare_method}.png"
    log.info(f"Write compare result figure to {outfig}")
    plt.savefig(outfig)
    # TODO:
    # 1. add score threshold option
    #    (1) output the number of markers per cluster
    # 2. Make output more verbose
    #    (1) sort genes by score
    #    (2) output a table, columns: [gene, mean_score, n_repeates]
    # 3. Calculate ARI of current markers
    m_markers = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(markers_per_group) 
    markers = pd.concat([m_markers.iloc[i] for i in range(m_markers.shape[0])]).reset_index(drop=True).unique() 
    log.info(f"marker genes: {markers}")
    if outpath:
        with open(outpath, "w") as f:
            for m in markers:
                f.write(m+"\n")


if __name__ == "__main__":
    fire.Fire(main)

