# Paths

Every script refers to data through four roots, written in the code as
placeholders.

| Placeholder | What it holds | Example |
|---|---|---|
| `${PROJECT_ROOT}` | this repository plus its `results/` tree | `/data/luad` |
| `${DATA_ROOT}` | raw downloads: GEO/ArrayExpress matrices, Visium sections | `/data/raw/LUAD` |
| `${WORK_ROOT}` | large intermediates and figure source tables | `/data/work` |
| `${HOME}` | user home; a few helper paths only | `/home/you` |

## Setting them

The placeholders are plain text inside ordinary string literals. Exporting
shell variables does **not** reach them — Python and R both read
`Path("${PROJECT_ROOT}/x")` as a directory named literally `${PROJECT_ROOT}`.
Substitute them once instead:

```bash
python config/set_paths.py \
    --project-root /data/luad \
    --data-root    /data/raw/LUAD \
    --work-root    /data/work \
    --apply
```

Without `--apply` the script only reports what it would change. It rewrites
files in place, so keep a pristine copy (or a git checkout) if you may want to
point the code at a different machine later.

The original scripts used absolute paths on the development machine. Those were
replaced by the placeholders above and the directory layout beneath each root is
unchanged, so substituting the four roots is sufficient.

## Two things the substitution does not cover

**`~/luad` and `Path.home() / "luad"`.** Some scripts reach the project tree
this way rather than through `${PROJECT_ROOT}` — 140 occurrences across 54 files.
They carry no user name, but they do assume the project sits at `~/luad`. Treat
`~/luad` as equivalent to `${PROJECT_ROOT}`: either place the tree there, or
substitute those as well:

```bash
grep -rl 'Path.home()\|~/luad' --include='*.py' --include='*.R' --include='*.sh' .
```

**`INFERCNV_ROOT`.** inferCNV was run on a rented compute node with its own
scratch directory. `revision/copy_number_and_cnv/G2_run_infercnv.R` and
`G3_launch_infercnv.sh` read `INFERCNV_ROOT`, falling back to
`${WORK_ROOT}/infercnv`. Set it if your scratch space lives elsewhere.

## Directories created under these roots

| Path | Written by | Read by |
|---|---|---|
| `${PROJECT_ROOT}/data/processed/*.h5ad` | `data_prep/atlas_build/` | most of `analysis/` |
| `${PROJECT_ROOT}/results/` | `analysis/`, `revision/` | `plotting/` |
| `${WORK_ROOT}/数据清洗/` | you, by copying the QC output here | `data_prep/atlas_build/01–03` |
| `${WORK_ROOT}/luad_figures/` | `analysis/`, `revision/` | `plotting/` |
| `${DATA_ROOT}/<accession>/` | you, when downloading | `data_prep/qc/` |

`${WORK_ROOT}/数据清洗` is a Chinese directory name meaning "data cleaning". It
is kept because the result tables reference it as a literal string; renaming it
would silently break filters that compare against the stored value.

## Where each input comes from

`config/data_provenance.md` lists every file the scripts read but do not
produce: which public source to download it from, and which of the derived
tables are distributed through the Zenodo deposit rather than regenerated
from this code.

## Data

No data are included in this repository. Every dataset used by the released
workflow is publicly accessible; accessions are listed in the README and cited
in the paper. Raw records for one spatial cohort (JGAS000613, JGAS000677) are
controlled-access and are neither required by the released workflow nor
redistributed here — see the README for what the workflow actually reads.
