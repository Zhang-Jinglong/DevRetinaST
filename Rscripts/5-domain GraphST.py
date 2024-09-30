import os
import pathlib

import numpy as np
import pandas as pd
import scanpy as sc
import torch
from GraphST import GraphST
from GraphST.utils import clustering

os.environ["CUDA_VISIBLE_DEVICES"] = "4"
os.environ["R_HOME"] = "/ssdata/users/zhangjl/Apps/miniconda3/envs/graphst/lib/R"

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Load data
base_dir = pathlib.Path("/ssdata/users/zhangjl/Proj/Dev_Retina/outputs/Domain/GraphST")
data = sc.read_csv(base_dir / "data.csv").T
meta = pd.read_csv(base_dir / "meta.csv").set_index("spot")
coo = pd.read_csv(base_dir / "coo.csv", index_col=0)
data.obs["data"] = meta.loc[data.obs_names, "batch"]
data.obsm["spatial"] = coo.loc[data.obs_names, :].values

# Model
n_clusters = 9

model = GraphST.GraphST(data, dim_input=data.n_vars, device=device)
data = model.train()

clustering(data, n_clusters, radius=50, method="mclust", refinement=False)

emb = data.obsm["emb_pca"]
domain = np.array(data.obs["domain"])

# Save
np.savetxt(base_dir / "emb.csv", emb)
np.savetxt(base_dir / "domain.csv", domain)
