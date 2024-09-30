import json

import pandas as pd
from matplotlib.image import imread
from shapely.geometry import Point, Polygon

prefix_1 = ["GSM74414" + str(i) for i in range(42, 58)]
prefix_2 = ["A", "B", "C", "D"] * 4
prefix_3 = ["15817_8PCW"] * 4 + ["15870_10PCW"] * 4 + ["15685_11PCW"] * 4 + ["14749_13PCW"] * 4

meta_list = []
for i in range(16):
    prefix_cb = prefix_1[i] + "_" + prefix_3[i] + "_" + prefix_2[i]
    img = imread("../data/" + prefix_cb + "_tissue_lowres_image.png")
    with open("../data/" + prefix_cb + "_tissue_lowres_image.json", "r") as f:
        data = json.load(f)["shapes"]
    with open("../data/" + prefix_cb + "_scalefactors_json.json", "r") as f:
        scale = json.load(f)["tissue_lowres_scalef"]
    meta = pd.read_csv("../data/GSE234035_" + prefix_3[i] + "_md.csv", index_col=0)
    meta.index = prefix_3[i] + "_" + meta.index
    meta = meta.loc[meta["section"].str[0] == prefix_2[i], ["Image_X", "Image_Y"]] * scale
    meta["Image_Y"] = meta["Image_Y"]

    ref_1 = pd.DataFrame(data[0]["points"])
    ref_1[1] = ref_1[1]
    ref_2 = pd.DataFrame(data[1]["points"])
    ref_2[1] = ref_2[1]

    polygon = Polygon(ref_1).union(Polygon(ref_2))

    ct = []
    for p_i in range(meta.shape[0]):
        point = Point(meta.iloc[p_i, :].values)
        ct.append(polygon.contains(point))
    meta["retina"] = ct
    meta_list.append(meta)

    # plt.figure(figsize=(8, 8), dpi=300)
    # plt.imshow(img)
    # plt.scatter(x=meta["Image_X"], y=meta["Image_Y"], c=meta["retina"], s=1)
    # plt.title(prefix_3[i] + "_" + prefix_2[i])
    # plt.show()

meta_retina = pd.concat(meta_list)
meta_retina.to_csv("retina_area.csv")
