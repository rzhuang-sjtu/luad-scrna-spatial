"""Step C: copy each section's hires H&E image + scalefactors to r_data/he/<cohort>__<sample>/."""
import json, shutil, os
from pathlib import Path
import numpy as np, scanpy as sc

OUT = Path("${DATA_ROOT}/ST/results/r_data/he")
OUT.mkdir(parents=True, exist_ok=True)

# E-MTAB: spatial dirs already exist alongside h5 files
EMTAB_RAW = Path("${DATA_ROOT}/ST/E-MTAB-13530/E-MTAB-13530")
EMTAB_QC  = Path("${DATA_ROOT}/ST/results/step01_qc/section_h5ad")

# Okamura: already extracted to step09_okamura_validation/raw
OKA_RAW = Path("${DATA_ROOT}/ST/results/step09_okamura_validation/raw")
OKA_QC  = Path("${DATA_ROOT}/ST/results/step09_okamura_validation/section_h5ad")


def copy_one(cohort, sample, hires_src, sf_src):
    out_dir = OUT / f"{cohort}__{sample}"
    out_dir.mkdir(parents=True, exist_ok=True)
    if hires_src.exists():
        shutil.copy2(hires_src, out_dir / "tissue_hires_image.png")
    if sf_src.exists():
        shutil.copy2(sf_src, out_dir / "scalefactors_json.json")
    print(f"  {cohort}/{sample} -> {out_dir}")


# E-MTAB
for h5 in sorted(EMTAB_QC.glob("*.h5ad")):
    s = h5.stem
    sp_dir = EMTAB_RAW / f"{s}-spatial"
    copy_one("EMTAB13530", s, sp_dir / "tissue_hires_image.png", sp_dir / "scalefactors_json.json")

# Okamura
for h5 in sorted(OKA_QC.glob("*.h5ad")):
    s = h5.stem
    sp_dir = OKA_RAW / s / "spatial"
    copy_one("Okamura", s, sp_dir / "tissue_hires_image.png", sp_dir / "scalefactors_json.json")

print("[done]")
