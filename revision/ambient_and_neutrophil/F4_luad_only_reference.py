"""F4 — Neutrophil label transfer using adenocarcinoma-only reference cells (label transfer restricted to adenocarcinoma reference cells).

Concern addressed: the Salcher reference atlas is from NSCLC (including squamous carcinoma); whether its neutrophil states transfer directly to pure LUAD.
Observed disease origin of 19,368 reference neutrophils: squamous 10,452 (54.0%), adenocarcinoma 7,937 (41.0%),
NSCLC NOS 808, normal 166 — squamous is indeed the majority; the concern is valid.

This script keeps only the 7,937 adenocarcinoma-derived reference cells and re-jointly trains with the 8,549 study neutrophils
scVI + scANVI for label transfer, then compares agreement with the original (full NSCLC reference) labels.
If the two label sets agree closely, neutrophil state definitions do not depend on histological subtype.

Reuse original pipeline parameters (step25c_scanvi_transfer.py):
  scVI  n_layers=2, n_latent=30, gene_likelihood="nb", max_epochs=200, early_stopping
  scANVI max_epochs=25, initialised from the trained scVI model
"""
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, os, warnings, time
warnings.filterwarnings("ignore")
import scvi, torch
from sklearn.metrics import cohen_kappa_score, adjusted_rand_score

OUT = "${PROJECT_ROOT}/results/sensitivity"
os.makedirs(OUT, exist_ok=True)
t0 = time.time()
def log(m): print(f"[{time.time()-t0:6.1f}s] {m}", flush=True)

log("Loading joint object ...")
j = ad.read_h5ad("${PROJECT_ROOT}/data/processed/luad_neutrophil_joint_scanvi.h5ad")
log(f"  {j.shape}  source: {j.obs['source'].value_counts().to_dict()}")

dis = pd.read_csv(f"{OUT}/salcher_neu_disease.csv").set_index("barcode")["disease"]
is_sal = j.obs["source"].astype(str) == "salcher"
d = pd.Series(index=j.obs_names, dtype=object)
d[is_sal.values] = dis.reindex(j.obs_names[is_sal.values]).values
keep = (~is_sal.values) | (d.values == "lung adenocarcinoma")
log(f"keep {keep.sum()} cells (this study {int((~is_sal.values).sum())} + adenocarcinoma reference"
    f"{int(((d.values=='lung adenocarcinoma')).sum())}）")

sub = j[keep].copy()
# Original labels (from full NSCLC reference transfer), for post-hoc comparison
orig = sub.obs["scanvi_predicted"].astype(str).values.copy() if "scanvi_predicted" in sub.obs else None
ann = ad.read_h5ad("${PROJECT_ROOT}/data/processed/luad_neutrophil_own_annotated.h5ad", backed="r")
own_lab = ann.obs["neu_subtype"].astype(str)

# scANVI requires true labels on reference cells and Unknown on query cells
ref_lab = sub.obs["transfer_label"].astype(str) if "transfer_label" in sub.obs else None
lab = np.array(["Unknown"] * sub.n_obs, dtype=object)
is_ref = (sub.obs["source"].astype(str) == "salcher").values
if ref_lab is not None:
    lab[is_ref] = ref_lab.values[is_ref]
else:
    raise SystemExit("Joint object lacks the reference label column; cannot retrain")
sub.obs["scanvi_label"] = pd.Categorical(lab)
sub.obs["batch_key"] = (sub.obs["source"].astype(str) + "-" + sub.obs["study"].astype(str)).values
log(f"reference label classes: {sorted(set(lab[is_ref]))}")

X = sub.layers["counts"] if "counts" in sub.layers else sub.X
sub.X = sp.csr_matrix(X).astype(np.float32)
scvi.settings.seed = 42
scvi.model.SCVI.setup_anndata(sub, batch_key="batch_key")
log("Training scVI ...")
m = scvi.model.SCVI(sub, n_layers=2, n_latent=30, gene_likelihood="nb")
m.train(max_epochs=200, early_stopping=True, early_stopping_patience=15,
        plan_kwargs={"lr": 1e-3}, check_val_every_n_epoch=1)
log("Training scANVI ...")
sm = scvi.model.SCANVI.from_scvi_model(m, adata=sub, labels_key="scanvi_label",
                                       unlabeled_category="Unknown")
sm.train(max_epochs=25, n_samples_per_label=100)
pred = sm.predict()
prob = sm.predict(soft=True)
sub.obs["new_label"] = pred
sub.obs["new_uncertainty"] = 1 - prob.max(axis=1).values
log("Training finished")

q = sub.obs.loc[~is_ref].copy()
q["own_prev"] = own_lab.reindex(q.index).values
res = q[["new_label", "own_prev", "new_uncertainty"]].dropna()
k = cohen_kappa_score(res.own_prev, res.new_label)
ari = adjusted_rand_score(res.own_prev, res.new_label)
agree = (res.own_prev == res.new_label).mean()
print(f"\n=== Agreement with original (full NSCLC reference) labels, n={len(res)} ===", flush=True)
print(f"exact match rate {agree*100:.1f}%   Cohen's kappa {k:.3f}   ARI {ari:.3f}", flush=True)
print(f"mean uncertainty of new labels {res.new_uncertainty.mean():.3f} (original pipeline 0.192)", flush=True)
ct = pd.crosstab(res.own_prev, res.new_label)
print("\nConfusion matrix (rows=original labels, columns=adenocarcinoma-reference new labels):", flush=True)
print(ct.to_string(), flush=True)
ct.to_csv(f"{OUT}/F4_label_confusion.csv")
res.to_csv(f"{OUT}/F4_labels.csv")
print(f"\nRetention rate by original label:", flush=True)
for s in ct.index:
    r = ct.loc[s]
    print(f"{s:<20} n={r.sum():>5}  retained {r.get(s,0)/r.sum()*100:5.1f}%"
          f"most often reassigned to {r.drop(s, errors='ignore').idxmax() if len(r.drop(s, errors='ignore'))>0 else '—'}", flush=True)
log(f"Results written to {OUT}/")
