import pathlib

import commot as ct
import numpy as np
import pandas as pd
import scanpy as sc

np.random.seed(3407)

# Load data
proj_dir = pathlib.Path("/Users/zhangjl/Documents/Proj/Dev_Retina/outputs/CCC/")

st_data = sc.read_csv(proj_dir / "COMMOT/retina_st.csv").T
st_meta = pd.read_csv(proj_dir / "COMMOT/retina_st_meta.csv", index_col=0)
st_data.obs = st_meta.loc[st_data.obs_names.to_numpy(), :]
st_data = st_data[st_data.obs["domain"] != "D9", :]

# COMMOT database
# ccc_db = ct.pp.ligand_receptor_database(
#     species="human", database="CellPhoneDB_v4.0", signaling_type=None
# )
# ccc_db.to_csv(proj_dir / "COMMOT/retina_db.csv")
ccc_db = pd.read_csv(
    proj_dir / "COMMOT/retina_db_filter.csv", index_col=0
)
ccc_db.columns = ["0", "1", "2", "3"]

# COMMOT
for sample_i in ["PCW9", "PCW12_early", "PCW12_late", "PCW14", "PCW16", "PCW17"]:
    st_data_sub = st_data[st_data.obs["sample"] == sample_i, :].copy()

    # Pre-processing
    sc.pp.normalize_total(st_data_sub, inplace=True)
    sc.pp.log1p(st_data_sub)

    st_data_sub.obsm["spatial"] = st_data_sub.obs.loc[:, ["coord.x", "coord.y"]].values

    # Spatial communication inference
    ccc_filtered = ct.pp.filter_lr_database(
        ccc_db, st_data_sub, min_cell_pct=0.05
    )

    ct.tl.spatial_communication(
        st_data_sub, database_name="CellPhoneDB_v4.0",
        df_ligrec=ccc_filtered, dis_thr=0.8, heteromeric=True, pathway_sum=True
    )

    st_data_sub.write_h5ad(
        proj_dir / ("COMMOT/commot_" + sample_i + ".h5ad")
    )

    # Overview of ligand-receptor pairs
    st_data_sub_lr = st_data_sub.obsm["commot-CellPhoneDB_v4.0-sum-sender"]
    st_data_sub_lr_sum = (
        st_data_sub_lr.sum(axis=0)
        .drop(["s-nan", "s-total-total"])
    )
    st_data_sub_lr_sum.index = st_data_sub_lr_sum.index.str[2:]

    st_data_sub_lr_sum.to_csv(
        proj_dir / ("COMMOT/commot_" + sample_i + "_lr.csv"),
        header=None
    )
