"""Step 24c: GSE93157 NanoString 730-gene panel — skip MP scoring.

Rationale:
  - Panel is immune-focused (730 genes). MP top-50 overlap:
      MP1: 11/50, MP2: 5/50, MP3: 8/50, MP4: 9/50  → all < 15/50
  - Per spec, < 15/50 → skip ssGSEA (insufficient signature representation).
  - Still record metadata + overlap report for transparency.
"""
from __future__ import annotations
import os, sys, gzip
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
sys.path.insert(0, str(os.path.expanduser("~/luad/scripts")))
from pathlib import Path
import pandas as pd
from importlib import import_module
C = import_module("24_common")

GSE = "gse93157"
RAW = Path("${DATA_ROOT}/GSE93157/GSE93157_raw_data_values.txt.gz")


def main():
    C.log(f"=== {GSE.upper()} NanoString 730-gene (SKIP MP scoring) ===")
    # Read panel genes
    with gzip.open(RAW, "rt") as f:
        lines = [ln.rstrip("\n") for ln in f]
    hdr_idx = next(i for i, ln in enumerate(lines) if ln.startswith("ID_REF"))
    panel_genes = []
    for ln in lines[hdr_idx + 1:]:
        if ln.strip() == "":
            continue
        g = ln.split("\t")[0].strip()
        if g and not g.startswith("#"):
            panel_genes.append(g)
    panel = set(panel_genes)
    C.log(f"panel size: {len(panel)} genes")

    # MP overlap report
    mp_sigs = C.load_mp_signatures()
    rows = []
    for mp, gl in mp_sigs.items():
        ov = [g for g in gl if g in panel]
        rows.append({"MP": mp, "n_signature": len(gl),
                      "n_in_panel": len(ov),
                      "overlap_genes": ";".join(ov),
                      "overlap_pct": round(len(ov) / len(gl) * 100, 1)})
    rep = pd.DataFrame(rows)
    C.log("\nMP × NanoString panel overlap:")
    print(rep[["MP", "n_signature", "n_in_panel", "overlap_pct"]].to_string(index=False))

    # Sample-level metadata for completeness
    gse_obj = C.fetch_geo_soft("GSE93157")
    md = C.gsm_characteristics(gse_obj)
    md["cancer_source"] = [gse_obj.gsms[g].metadata.get("source_name_ch1", [""])[0]
                            for g in md["GSM"]]
    nsclc = md[md["cancer_source"].isin(
        ["LUNG NON-SQUAMOUS CANCER", "SQUAMOUS LUNG CANCER"])].copy()
    C.log(f"\nNSCLC subset: {len(nsclc)} (LUAD-like {sum(nsclc['cancer_source']=='LUNG NON-SQUAMOUS CANCER')}, "
          f"SQ {sum(nsclc['cancer_source']=='SQUAMOUS LUNG CANCER')})")
    C.log(f"best.resp distribution: {nsclc['best.resp'].value_counts().to_dict()}")

    # Define R/NR by RECIST: CR/PR = R; SD/PD = NR (alternate: include SD as NR strict)
    nsclc["response_group"] = nsclc["best.resp"].map(
        {"CR": "R", "PR": "R", "SD": "NR", "PD": "NR"}).fillna("NA")
    C.log(f"R/NR groups (NSCLC): {nsclc['response_group'].value_counts().to_dict()}")

    OUT = C.OUT; OUT.mkdir(parents=True, exist_ok=True)
    rep.to_csv(OUT / f"{GSE}_mp_panel_overlap.csv", index=False)
    nsclc[["GSM", "title", "cancer_source", "best.resp", "response_group",
           "pfs", "pfse", "drug", "age", "sex"]].to_csv(
        OUT / f"{GSE}_metadata_nsclc.csv", index=False)

    # Status note
    note = (f"# GSE93157 NanoString 730 — SKIPPED for MP scoring\n\n"
            f"Panel size: {len(panel)} genes; MP top-50 overlap < 15/50 for all 4 MPs.\n\n"
            + rep.to_markdown(index=False) + "\n\n"
            f"NSCLC samples available: {len(nsclc)}\n"
            f"  - LUNG NON-SQUAMOUS (≈LUAD): "
            f"{sum(nsclc['cancer_source']=='LUNG NON-SQUAMOUS CANCER')}\n"
            f"  - SQUAMOUS LUNG: "
            f"{sum(nsclc['cancer_source']=='SQUAMOUS LUNG CANCER')}\n\n"
            f"Best response distribution: "
            f"{nsclc['best.resp'].value_counts().to_dict()}\n")
    (OUT / f"{GSE}_skip_note.md").write_text(note, encoding="utf-8")
    C.log("written outputs")


if __name__ == "__main__":
    main()
