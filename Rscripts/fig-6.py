import pathlib

import commot as ct
import matplotlib.pyplot as plt
import scanpy as sc

# Spatial communication direction (Fig. 6a)
proj_dir = pathlib.Path("/Users/zhangjl/Documents/Proj/Dev_Retina/outputs/")

sample_list = ["PCW9", "PCW12_early", "PCW12_late", "PCW14", "PCW16", "PCW17"]
scale_list = [40, 40, 70, 50, 60, 60]
for sample_i, scale_i in zip(sample_list, scale_list):
    st_data_i = sc.read_h5ad(
        proj_dir / ("CCC/COMMOT/commot_" + sample_i + ".h5ad")
    )

    ct.tl.communication_direction(st_data_i, database_name="CellPhoneDB_v4.0", k=5)
    st_data_i.obsm["spatial"] = st_data_i.obs.loc[:, ["coord.x", "coord.y"]].values
    plt.figure(figsize=(3, 3))
    ct.pl.plot_cell_communication(
        st_data_i, database_name="CellPhoneDB_v4.0",
        plot_method="grid", grid_density=0.6, scale=scale_i, ndsize=30,
        background_legend=False, background="cluster", clustering="domain",
        cluster_cmap={
            "D1": "#A9170E", "D2": "#E25B5F", "D3": "#FEB7BB",
            "D4": "#DBCAE4", "D5": "#906C9C", "D6": "#74B5D6",
            "D7": "#718CB2", "D8": "#6862A7"
        }, arrow_color="#670009"
    )
    plt.savefig(
        proj_dir / ("Visualization/fig-6a_" + sample_i + ".pdf")
    )
