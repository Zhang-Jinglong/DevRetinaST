import os
from pathlib import Path

import cell2location
import numpy as np
import pandas as pd
import scanpy as sc
from cell2location.models import RegressionModel
from cell2location.utils.filtering import filter_genes

os.environ["CUDA_VISIBLE_DEVICES"] = "3"

proj_dir = "/data/home/zhangjl/Proj/Dev_Retina/Deconvolution/Cell2location/"

# Load data
st_count = sc.read_csv(Path(proj_dir + "retina_st_pub.csv")).T
sc_count = sc.read_csv(Path(proj_dir + "retina_sc.csv")).T
sc_count.obs_names = sc_count.obs_names.str.replace("\\.1", "-1")
sc_meta = pd.read_csv(proj_dir + "retina_sc_meta.csv", index_col=0)
st_meta = pd.read_csv(proj_dir + "retina_st_meta_pub.csv", index_col=0)
st_meta.index = "X" + st_meta.index

# Pre-processing
st_count.var["SYMBOL"] = st_count.var_names
st_count.obs["batch"] = st_meta.loc[st_count.obs_names, "batch"]

sc_count.var["SYMBOL"] = sc_count.var_names
sc_count.obs["cell_type"] = sc_meta.loc[sc_count.obs_names, "cell_type"]
sc_count.obs["batch"] = sc_meta.loc[sc_count.obs_names, "batch"]

selected = filter_genes(
    sc_count,
    cell_count_cutoff=15,
    cell_percentage_cutoff2=0.05,
    nonz_mean_cutoff=1.12
)
sc_count = sc_count[:, selected].copy()

# Estimation of reference cell type signatures
cell2location.models.RegressionModel.setup_anndata(
    adata=sc_count,
    batch_key="batch",
    labels_key="cell_type"
)
mod = RegressionModel(sc_count)
mod.train(max_epochs=250, use_gpu=True)

sc_count = mod.export_posterior(
    sc_count,
    sample_kwargs={
        "num_samples": 1000,
        "batch_size": 2500,
        "use_gpu": True
    }
)
mod.plot_QC()

if "means_per_cluster_mu_fg" in sc_count.varm.keys():
    inf_aver = sc_count.varm["means_per_cluster_mu_fg"][[
        f"means_per_cluster_mu_fg_{i}"
        for i in sc_count.uns["mod"]["factor_names"]
    ]].copy()
else:
    inf_aver = sc_count.var[[
        f"means_per_cluster_mu_fg_{i}"
        for i in sc_count.uns["mod"]["factor_names"]
    ]].copy()
inf_aver.columns = sc_count.uns["mod"]["factor_names"]

# Cell2location: spatial mapping
intersect = np.intersect1d(st_count.var_names, inf_aver.index)
st_count = st_count[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

cell2location.models.Cell2location.setup_anndata(
    adata=st_count,
    batch_key="batch"
)
mod = cell2location.models.Cell2location(
    st_count,
    cell_state_df=inf_aver,
    N_cells_per_location=20,
    detection_alpha=20
)
mod.train(
    max_epochs=30000,
    batch_size=None,
    train_size=1,
    use_gpu=True
)

st_count = mod.export_posterior(
    st_count,
    sample_kwargs={
        "num_samples": 1000,
        "batch_size": mod.adata.n_obs,
        "use_gpu": True
    }
)
results = st_count.obsm["q05_cell_abundance_w_sf"].T
results.index = results.index.str[23:]
results.to_csv(
    proj_dir + "cell2location_pub.csv"
)
