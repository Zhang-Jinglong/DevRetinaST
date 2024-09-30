import pathlib

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import KFold, cross_val_score
from sklearn.svm import SVC

# Load data
base_dir = pathlib.Path("/ssdata/users/zhangjl/Proj/Dev_Retina/outputs/Domain/CrossVal")
x = pd.read_csv(base_dir / "expr.csv", index_col=0).T
y = pd.read_csv(base_dir / "domain.csv", index_col=0)

# Model
rf_clf = RandomForestClassifier(random_state=0)
lr_clf = LogisticRegression(random_state=0)
svm_clf = SVC(random_state=0)

# 10-fold
k_fold = KFold(n_splits=10, random_state=0, shuffle=True)

score_Harmony = pd.DataFrame({
    "Accuracy": np.concatenate([
        cross_val_score(rf_clf, x, y["domain.harmony"].values, cv=k_fold),
        cross_val_score(lr_clf, x, y["domain.harmony"].values, cv=k_fold),
        cross_val_score(svm_clf, x, y["domain.harmony"].values, cv=k_fold)
    ]),
    "Classifier": ["Random forest"] * 10 + ["Logistic regression"] * 10 + ["SVM"] * 10,
    "method": "Harmony"
})
print("Median ACC (Harmony): %.4f" % np.median(score_Harmony["Accuracy"]))

score_BASS = pd.DataFrame({
    "Accuracy": np.concatenate([
        cross_val_score(rf_clf, x, y["domain.bass"].values, cv=k_fold),
        cross_val_score(lr_clf, x, y["domain.bass"].values, cv=k_fold),
        cross_val_score(svm_clf, x, y["domain.bass"].values, cv=k_fold)
    ]),
    "Classifier": ["Random forest"] * 10 + ["Logistic regression"] * 10 + ["SVM"] * 10,
    "method": "BASS"
})
print("Median ACC (BASS): %.4f" % np.median(score_BASS["Accuracy"]))

score_GraphST = pd.DataFrame({
    "Accuracy": np.concatenate([
        cross_val_score(rf_clf, x, y["domain.graphst"].values, cv=k_fold),
        cross_val_score(lr_clf, x, y["domain.graphst"].values, cv=k_fold),
        cross_val_score(svm_clf, x, y["domain.graphst"].values, cv=k_fold)
    ]),
    "Classifier": ["Random forest"] * 10 + ["Logistic regression"] * 10 + ["SVM"] * 10,
    "method": "GraphST"
})
print("Median ACC (GraphST): %.4f" % np.median(score_GraphST["Accuracy"]))

score_CCA = pd.DataFrame({
    "Accuracy": np.concatenate([
        cross_val_score(rf_clf, x, y["domain.cca"].values, cv=k_fold),
        cross_val_score(lr_clf, x, y["domain.cca"].values, cv=k_fold),
        cross_val_score(svm_clf, x, y["domain.cca"].values, cv=k_fold)
    ]),
    "Classifier": ["Random forest"] * 10 + ["Logistic regression"] * 10 + ["SVM"] * 10,
    "method": "Seurat-CCA"
})
print("Median ACC (Seurat-CCA): %.4f" % np.median(score_CCA["Accuracy"]))

score_scVI = pd.DataFrame({
    "Accuracy": np.concatenate([
        cross_val_score(rf_clf, x, y["domain.scvi"].values, cv=k_fold),
        cross_val_score(lr_clf, x, y["domain.scvi"].values, cv=k_fold),
        cross_val_score(svm_clf, x, y["domain.scvi"].values, cv=k_fold)
    ]),
    "Classifier": ["Random forest"] * 10 + ["Logistic regression"] * 10 + ["SVM"] * 10,
    "method": "scVI"
})
print("Median ACC (scVI): %.4f" % np.median(score_scVI["Accuracy"]))

# Save
score_all = pd.concat([score_Harmony, score_BASS, score_GraphST, score_CCA, score_scVI])
score_all.to_csv(base_dir / "acc_score.csv")
