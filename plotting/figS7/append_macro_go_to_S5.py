"""Add the 7-macrophage GO BP enrichment data (Fig S6 c-i source) as a new
sheet inside Supplementary Table S5."""
from pathlib import Path
import pandas as pd

OUT = Path("${WORK_ROOT}/Supplementary_Tables/Supplementary_Table_S5_Myeloid_Subtype_DEGs.xlsx")
GO  = Path("${WORK_ROOT}/luad_figures/fig4/myeloid_go_enrichment.csv")

go = pd.read_csv(GO)
bp = go[go["gene_set"] == "GO_Biological_Process_2023"].copy()
macros = ["Macro_SPP1", "Macro_C1QC", "Macro_FCN1", "Macro_FOLR2",
          "Macro_MARCO", "Macro_general", "Macro_prolif"]
bp = bp[bp["subtype"].isin(macros)]
bp = bp[["subtype", "term", "overlap", "p_value", "adj_p_value",
        "odds_ratio", "combined_score", "genes"]]
bp["subtype"] = pd.Categorical(bp["subtype"], categories=macros, ordered=True)
bp = bp.sort_values(["subtype", "adj_p_value"]).reset_index(drop=True)

# Read existing S5, copy two existing sheets, then append new one
existing = pd.read_excel(OUT, sheet_name=None, header=None)
with pd.ExcelWriter(OUT, engine="openpyxl") as xw:
    for sheet, df in existing.items():
        df.to_excel(xw, sheet_name=sheet[:31], index=False, header=False)

    # New sheet: Macrophage GO BP enrichment
    title = ("Table S5c. GO Biological Process enrichment for the seven "
             "macrophage subtypes (Figure S6 c-i source).")
    subtitle = ("Source: Enrichr GO_Biological_Process_2023 on the per-subtype "
                "DE genes (logFC > 0.5 & adj p < 0.05). 10 top terms per "
                "subtype shown, ranked by combined_score.")
    rows = [[title] + [None]*7,
            [subtitle] + [None]*7,
            [None]*8,
            list(bp.columns)]
    rows.extend(bp.values.tolist())
    pd.DataFrame(rows).to_excel(xw, sheet_name="Macrophage GO BP",
                                 index=False, header=False)

print(f"appended Macrophage GO BP sheet -> {OUT}")
print(f"rows: {len(bp)} (7 subtypes x 10 BP terms)")
