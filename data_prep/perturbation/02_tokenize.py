"""Tokenize each transition input with Geneformer V2 TranscriptomeTokenizer."""
import sys
from pathlib import Path
from geneformer import TranscriptomeTokenizer

TRANSITIONS = ["macro_spp1_to_c1qc", "mal_mp3_to_mp1", "neu_osm_priming_to_low"]
INP_ROOT = Path("${PROJECT_ROOT}/results/fig8_geneformer/inputs")
TOK_ROOT = Path("${PROJECT_ROOT}/results/fig8_geneformer/tokenized")
TOK_ROOT.mkdir(parents=True, exist_ok=True)

for t in TRANSITIONS:
    print(f"\n=== tokenize {t} ===")
    tk = TranscriptomeTokenizer(
        custom_attr_name_dict={"cell_state": "cell_state"},
        nproc=8,
        model_version="V2",
    )
    tk.tokenize_data(
        data_directory=str(INP_ROOT / t),
        output_directory=str(TOK_ROOT / t),
        output_prefix=t,
        file_format="h5ad",
    )
print("\nDONE.")
