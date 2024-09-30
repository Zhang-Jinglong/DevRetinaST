import os
import pathlib

import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import torch

os.environ["CUDA_VISIBLE_DEVICES"] = "4"

scvi.settings.seed = 0
torch.set_float32_matmul_precision("high")

# Load data
base_dir = pathlib.Path("/ssdata/users/zhangjl/Proj/Dev_Retina/outputs/Domain/scVI")
data = sc.read_csv(base_dir / "data.csv").T
meta = pd.read_csv(base_dir / "meta.csv").set_index("spot")
data.obs["batch"] = meta.loc[data.obs_names, "batch"]

# Model
scvi.model.SCVI.setup_anndata(data, batch_key="batch")
model = scvi.model.SCVI(data, n_layers=2, n_latent=30, gene_likelihood="zinb")
model.train(max_epochs=500)

emb = model.get_latent_representation()

# Save
np.savetxt(base_dir / "emb.csv", emb)
