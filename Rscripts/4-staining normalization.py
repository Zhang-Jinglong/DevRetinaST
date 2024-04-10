import cv2
import matplotlib.pyplot as plt
import numpy as np
from pathml.preprocessing import StainNormalizationHE
from skimage import io

np.random.seed(3407)

for slice_i in ["PCW9_12", "PCW12", "PCW14", "PCW16", "PCW17"]:
    input_path = "/data/home/zhangjl/Proj/Dev_Retina/images/" + slice_i + ".tif"
    output_path = "/data/home/zhangjl/Proj/Dev_Retina/images/" + slice_i + "_norm.jpg"

    img = io.imread(input_path)

    # Image augmentation
    clahe = cv2.createCLAHE(
        clipLimit=2, tileGridSize=(8, 8)
    )
    img[:, :, 0] = clahe.apply(img[:, :, 0])
    img[:, :, 1] = clahe.apply(img[:, :, 1])
    img[:, :, 2] = clahe.apply(img[:, :, 2])

    # H&E normalization
    normalizer = StainNormalizationHE(
        target="hematoxylin", stain_estimation_method="macenko"
    )
    img_norm = normalizer.F(img)

    _, axes = plt.subplots(1, 2, figsize=(10, 5), dpi=300)
    axes[0].imshow(img)
    axes[1].imshow(img_norm)
    plt.show()

    # Adaptive threshold
    gray = cv2.cvtColor(img_norm, cv2.COLOR_RGB2GRAY)
    thresh = cv2.adaptiveThreshold(
        gray, 255,
        cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY,
        1025, 30
    )

    _, axes = plt.subplots(1, 2, figsize=(10, 5), dpi=300)
    axes[0].imshow(img)
    axes[1].imshow(thresh, "gray")
    plt.show()

    io.imsave(output_path, thresh, quality=7)
