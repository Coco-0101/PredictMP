"""
Step 01: read the 49 GTEx v8 elastic-net model DBs, merge all SNP weights
into one wide table.
  -> weight.csv (gene, rsid, varID, ref_allele, eff_allele, <49 tissues>)
  -> Num_of_genes_tissue.csv
"""
import sqlite3
import pandas as pd
import re

from config import ELASTIC_NET_DIR, WEIGHT_FILE, GENE_NUM_FILE

# --- all 49 tissues ---
db_files = sorted(ELASTIC_NET_DIR.glob("en_*.db"))

connections = [sqlite3.connect(f) for f in db_files]

snp_weight_dfs = [pd.read_sql("SELECT * FROM weights", c) for c in connections]

# how many genes are available for each given tissue?
gene_nums = [
    pd.read_sql("SELECT count(*) FROM extra", c).iloc[0, 0]
    for c in connections
]

tissue_names = [re.sub(r"^en_(.*?)\.db$", r"\1", f.name) for f in db_files]
gene_num_df = pd.DataFrame({"tissue": tissue_names, "gene_number": gene_nums})
print(gene_num_df)

gene_num_df.to_csv(GENE_NUM_FILE, index=False)

print(f"Mean genes: {gene_num_df['gene_number'].mean():.0f}")  # ~5752
print(f"Range: {gene_num_df['gene_number'].min()} – {gene_num_df['gene_number'].max()}")

# check row counts per tissue
row_counts = [len(df) for df in snp_weight_dfs]
print(f"Max rows: {max(row_counts)}")  # 309155
print(f"Tissue with max rows: {tissue_names[row_counts.index(max(row_counts))]}")  # Nerve_Tibial

# rename the weight column to tissue name for each tissue dataframe
key_cols = ["gene", "rsid", "varID", "ref_allele", "eff_allele"] # key columns
for i, (df, tname) in enumerate(zip(snp_weight_dfs, tissue_names)):
    df.columns = key_cols + [tname]

# outer-join merge all tissues on key columns
merged = snp_weight_dfs[0]
for df in snp_weight_dfs[1:]:
    merged = merged.merge(df, on=key_cols, how="outer")

# 5318682 unique SNP-gene pairs across all tissues, with 49 tissue-specific weight columns + 5 key columns
print(merged.shape)  # (5318682, 54)

# output
merged.to_csv(WEIGHT_FILE, index=False)

for c in connections:
    c.close()
