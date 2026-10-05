"""
Step 04: convert gnomAD v4.1 per-chromosome rsid CSV to Parquet
(step 05 reads only the needed columns).
"""

import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import logging
from tqdm import tqdm

from config import GNOMAD_RSID_DIR as GNOMAD_DIR, GNOMAD_PARQUET_DIR as OUTPUT_DIR, CHROMS

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)


#  Per-chromosome conversion
def convert_chrom(chrom: str) -> str:
    csv_path     = GNOMAD_DIR  / f"gnomad_{chrom}_populations_rsid.csv"
    parquet_path = OUTPUT_DIR  / f"gnomad_{chrom}_populations_rsid.parquet"

    if parquet_path.exists():
        return f"{chrom}: already exists, skipped"

    df = pd.read_csv(csv_path, dtype={"CHROM": str}, low_memory=False)
    af_cols = [c for c in df.columns if c.startswith("AF_")]
    df[af_cols] = df[af_cols].apply(pd.to_numeric, errors="coerce")
    df.to_parquet(parquet_path, index=False)
    return f"{chrom}: done ({csv_path.stat().st_size / 1e9:.1f} GB → {parquet_path.stat().st_size / 1e9:.1f} GB)"



def run():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    logging.info(f"Converting {len(CHROMS)} chromosomes ({os.cpu_count()} workers)...")

    with ProcessPoolExecutor(max_workers=os.cpu_count()) as exe:
        futures = {exe.submit(convert_chrom, chrom): chrom for chrom in CHROMS}
        for fut in tqdm(as_completed(futures), total=len(CHROMS), desc="chromosomes"):
            logging.info(fut.result())

    logging.info("Done.")


# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    run()
