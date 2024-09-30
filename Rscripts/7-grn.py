import pathlib
import pickle
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import scripro
from tqdm import tqdm

warnings.filterwarnings("ignore")

base_dir = pathlib.Path("/ssdata/users/zhangjl/Proj/Dev_Retina/outputs/Development/SCRIPro")
if (base_dir / "scripro.pkl").exists():
    with open(base_dir / "scripro.pkl", "rb") as file:
        rna_seq_data = pickle.load(file)
else:
    # Load data
    data = sc.read_csv(base_dir / "data.csv").T
    meta = pd.read_csv(base_dir / "meta.csv").set_index("spot")
    data.obs["sample"] = meta.loc[data.obs_names, "sample"].values
    data.obs["domain"] = meta.loc[data.obs_names, "domain"].values.astype(str)
    data.obs["group"] = data.obs["sample"] + "_" + data.obs["domain"]
    count = data.obs["group"].value_counts()
    data = data[data.obs["group"].isin(count.index[count > 3]), :].copy()
    data.raw = data
    sc.pp.normalize_total(data, target_sum=1e4)
    sc.pp.log1p(data)

    # SCRIPro
    test_data = scripro.Ori_Data(
        adata=data, cluster_method="group",
        Cell_num=3, min_cell=2, cores=100
    )
    test_data.get_positive_marker_gene_parallel()

    rna_seq_data = scripro.SCRIPro_RNA(
        100, "hg38", test_data, assays=["Direct", "DNase", "H3K27ac"]
    )
    rna_seq_data.cal_ISD_cistrome()
    rna_seq_data.get_tf_score()

    with open(base_dir / "scripro.pkl", "wb") as file:
        pickle.dump(rna_seq_data, file)

# TF-target co-expression
if (base_dir / "scripro_tf_target.pkl").exists():
    with open(base_dir / "scripro_tf_target.pkl", "rb") as file:
        tf_target_dict = pickle.load(file)
else:
    reg_tf_list = pd.read_csv(base_dir / "reg_tf_list.csv")["x"]
    p_value = rna_seq_data.P_value_matrix
    p_value = p_value.loc[:, np.intersect1d(p_value.columns, reg_tf_list)]

    tf_ava_idx = (p_value < 0.05).sum(axis=0) > 0
    tf_ava = p_value.columns[tf_ava_idx].tolist()

    tf_target_dict = {}
    for tf_i in tqdm(tf_ava):
        score_i = rna_seq_data.get_tf_target(tf_i)
        tf_target_dict[tf_i] = score_i

    with open(base_dir / "scripro_tf_target.pkl", "wb") as file:
        pickle.dump(tf_target_dict, file)

# TF-target network
tf_target_strength_list = []
for tf_i in tf_target_dict:
    df_i = tf_target_dict[tf_i].copy()
    df_i_max = pd.DataFrame(df_i.max(axis=0)).reset_index()
    tf_target_strength_list.append(pd.DataFrame({
        "TF": tf_i, "Target": df_i_max.iloc[:, 0], "Strength": df_i_max.iloc[:, 1]
    }))

tf_target_strength = pd.concat(tf_target_strength_list)

tf_target_strength.to_csv(base_dir / "scripro_tf_target_strength.csv")
