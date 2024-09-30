import pathlib

import commot as ct
import scanpy as sc

# Domain communication network (Fig. S15a)
proj_dir = pathlib.Path("/Users/zhangjl/Documents/Proj/Dev_Retina/outputs/")

sample_list = ["PCW9", "PCW12_early", "PCW12_late", "PCW14", "PCW16", "PCW17"]
for sample_i in sample_list:
    st_data_i = sc.read_h5ad(
        proj_dir / ("CCC/COMMOT/commot_" + sample_i + ".h5ad")
    )

    ct.tl.cluster_communication(
        st_data_i, database_name="CellPhoneDB_v4.0",
        pathway_name=None, lr_pair=None,
        clustering="domain", n_permutations=1000, random_seed=0
    )
    # noinspection PyTypeChecker
    ct.pl.plot_cluster_communication_network(
        st_data_i, uns_names=["commot_cluster-domain-CellPhoneDB_v4.0-total-total"],
        nx_node_pos=None, nx_bg_pos=False,
        quantile_cutoff=1.0, p_value_cutoff=0.05, nx_node_cmap="Light24",
        filename=str(proj_dir / ("Visualization/fig-S15a_" + sample_i + ".pdf"))
    )
