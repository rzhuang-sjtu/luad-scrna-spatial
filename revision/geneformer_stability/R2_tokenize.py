"""R2 — tokenize the rebuilt full-gene neutrophil inputs with the Geneformer V2 tokenizer.

Reads  results/fig8_geneformer/inputs_v2/{name}/data.h5ad
Writes results/fig8_geneformer/tokenized_v2/{name}/{name}.dataset

Prints the resulting token-length distribution so the fix can be verified against the
published run (mean 94.6 / median 89 / max 225 tokens on the HVG-truncated input).
"""
import sys
from pathlib import Path
import numpy as np
from geneformer import TranscriptomeTokenizer
from datasets import load_from_disk

INP_ROOT = Path("${PROJECT_ROOT}/results/fig8_geneformer/inputs_v2")
TOK_ROOT = Path("${PROJECT_ROOT}/results/fig8_geneformer/tokenized_v2")
TOK_ROOT.mkdir(parents=True, exist_ok=True)

NAMES = sys.argv[1:] or [
    "neu_osm_priming_to_low",
    "neu_osm_priming_to_low_half1",
    "neu_osm_priming_to_low_half2",
]

for n in NAMES:
    print(f"\n=== tokenize {n} ===")
    tk = TranscriptomeTokenizer(
        custom_attr_name_dict={"cell_state": "cell_state"},
        nproc=8,
        model_version="V2",
    )
    tk.tokenize_data(
        data_directory=str(INP_ROOT / n),
        output_directory=str(TOK_ROOT / n),
        output_prefix=n,
        file_format="h5ad",
    )
    ds = load_from_disk(str(TOK_ROOT / n / f"{n}.dataset"))
    L = np.array(ds["length"])
    st = np.array(ds["cell_state"])
    print(f"  cells={len(ds)}  sender={(st=='sender').sum()}  receiver={(st=='receiver').sum()}")
    print(f"  tokens/cell: mean={L.mean():.1f} median={int(np.median(L))} "
          f"min={int(L.min())} max={int(L.max())}")
    ls = L[st == "sender"]
    print(f"  sender tokens/cell: mean={ls.mean():.1f} median={int(np.median(ls))} "
          f"(published HVG-truncated run: mean 94.6 / median 89)")

print("\nDONE. Next: R3_perturb.py")
