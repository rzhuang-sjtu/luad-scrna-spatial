"""Check var format / gene-name harmonization between Salcher and own myeloid."""
import anndata as ad

s = ad.read_h5ad("${DATA_ROOT}/High-resolution/neutrophil_final.h5ad", backed="r")
m = ad.read_h5ad("${PROJECT_ROOT}/data/processed/luad_myeloid.h5ad", backed="r")

print("Salcher var.columns:", list(s.var.columns))
print("Salcher var.head():")
print(s.var.head().to_string())
print("\nSalcher var_names sample:", s.var_names[:10].tolist())

print("\n" + "=" * 60)
print("Myeloid var.columns:", list(m.var.columns))
print("Myeloid var.head():")
print(m.var.head().to_string())
print("\nMyeloid var_names sample:", m.var_names[:10].tolist())

# look for symbol column
for col in s.var.columns:
    if "symbol" in col.lower() or "name" in col.lower() or "feature" in col.lower():
        print(f"\nSalcher var['{col}'].head():")
        print(s.var[col].head(10).to_string())
        sv = set(s.var[col].astype(str))
        mv = set(m.var_names.astype(str))
        print(f"  intersection with own var_names: {len(sv & mv)}")

for col in m.var.columns:
    if "symbol" in col.lower() or "ensembl" in col.lower() or "ensg" in col.lower():
        print(f"\nMyeloid var['{col}'].head():")
        print(m.var[col].head(10).to_string())
