"""Inspect MP1-MP5 signatures + Fig3 outputs to understand the semantics."""
import os, glob
import pandas as pd
sig = pd.read_csv("${PROJECT_ROOT}/results/step6_mp_signatures_top100.csv")
print("MPs in sig:", sorted(sig['MP'].unique()))
print("rank cols:", sig.columns.tolist())
for mp in sorted(sig['MP'].unique()):
    top = sig[(sig['MP']==mp) & (sig['rank']<=20)]['gene'].tolist()
    print(f"  {mp} top-20: {top}")

# Look for any semantic labels in fig3 dir
print("\n=== fig3 archive ===")
for f in sorted(glob.glob('${WORK_ROOT}/luad_figures/fig3/*')):
    print(f"  {os.path.basename(f)}  ({os.path.getsize(f)} bytes)")

# any md/summary
for f in glob.glob('${WORK_ROOT}/luad_figures/fig3/*.md'):
    print(f"\n--- {f} ---")
    print(open(f).read()[:2000])
