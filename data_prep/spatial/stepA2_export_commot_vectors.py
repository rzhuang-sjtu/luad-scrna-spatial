"""Step A2: export COMMOT sender/receiver vector fields per section, append to CSVs."""
import os
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc

OUT = Path("${DATA_ROOT}/ST/results/r_data/per_section")

CONFIGS = [
    ("EMTAB13530", Path("${DATA_ROOT}/ST/results/step06_commot/section_h5ad")),
    ("Okamura",   Path("${DATA_ROOT}/ST/results/step09_okamura_validation/commot")),
]
PATHWAYS = ["OSM", "IL1"]

for cohort, sec_dir in CONFIGS:
    for h5 in sorted(sec_dir.glob("*.h5ad")):
        sample = h5.stem
        a = sc.read_h5ad(str(h5), backed="r")
        coords = a.obsm["spatial"]
        rows = pd.DataFrame({"section_barcode": a.obs_names})
        for p in PATHWAYS:
            for kind in ("sender", "receiver"):
                key = f"commot_{kind}_vf-cellchat-{p}"
                if key in a.obsm:
                    vf = np.asarray(a.obsm[key])
                    rows[f"vf_{kind[0]}_{p}_dx"] = vf[:, 0]
                    rows[f"vf_{kind[0]}_{p}_dy"] = vf[:, 1]
        a.file.close()

        csv_path = OUT / f"{cohort}__{sample}.csv"
        if not csv_path.exists():
            print(f"[skip] {csv_path.name} (CSV not found)")
            continue
        df = pd.read_csv(csv_path, index_col=0)
        # cohort obs_names look like "<barcode>-<sample>" → strip suffix
        suf = "-" + sample
        section_bc = [n[:-len(suf)] if n.endswith(suf) else n for n in df.index]
        rows = rows.set_index("section_barcode").reindex(section_bc)
        rows.index = df.index
        # remove existing vf_* columns from prior runs
        df = df.loc[:, ~df.columns.str.startswith("vf_")]
        df = df.join(rows)
        df.to_csv(csv_path)
        added = [c for c in rows.columns]
        print(f"  {cohort}/{sample}: added {len(added)} vf columns")
print("[done]")
