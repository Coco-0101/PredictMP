"""
Step 08: build PredictMP reference data (read by PredictMP_pkg/loader.py)

Step 1: Read SNP usage counts from SQLite databases
        -> Num_of_rs_used_x.csv, tissue_index.txt

Step 2: Compute per-tissue gene mean and SD for each population
        (NFE / EAS / AFR / AMI / AMR / ASJ / FIN / MID / SAS)
        -> result_table_{pop_lower}.csv  (one per population)

Step 3: Pack all population tissue data into parquet files (for Python package)
        -> PredictMP/data/tissue/{pop_lower}_x1.parquet ... {pop_lower}_x49.parquet
        -> PredictMP/data/result_table/result_table_{pop_lower}.parquet

Step 4: Pre-compute KS test (each population vs NFE) + NFE.{POP}.mean.diff
        -> PredictMP/data/ks/ks_{pop_lower}_nfe.parquet

Usage
-----
    python 08_build_package_data.py            # run all 4 steps
    python 08_build_package_data.py --step 1
    python 08_build_package_data.py --step 2
    python 08_build_package_data.py --step 3
    python 08_build_package_data.py --step 4
    python 08_build_package_data.py --step 1 2 3
"""

import argparse
import sqlite3
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats as sp_stats

from config import (
    DATA_DIR, ELASTIC_NET_DIR, PKG_TISSUE_DIR, PKG_RT_DIR, PKG_KS_DIR,
    POPULATIONS, REF_POP, pop_results_dir,
)

NON_REF_POPS = [p for p in POPULATIONS if p != REF_POP]
PRED_GLOB    = "*_predicted_expression.txt"
KS_THRESHOLD = 5e-8     # KS test pop vs NFE


# =============================================
# Step 1: SNP usage counts & tissue index
# =============================================

def step1():
    for _d in [PKG_TISSUE_DIR, PKG_RT_DIR, PKG_KS_DIR]:
        _d.mkdir(parents=True, exist_ok=True)

    db_files = sorted(ELASTIC_NET_DIR.glob("en_*.db"))
    n_tissues = len(db_files)
    print(f"Found {n_tissues} tissue databases")

    tis_name = [f.stem[3:] for f in db_files]

    tissue_table = pd.DataFrame({"index": range(1, n_tissues + 1), "tissue": tis_name})
    tissue_table.to_csv(DATA_DIR / "tissue_index.txt", index=False)
    print("Saved tissue_index.txt")

    gene_count_list = []
    for i, dbpath in enumerate(db_files):
        conn = sqlite3.connect(dbpath)
        weights = pd.read_sql("SELECT gene FROM weights", conn)
        conn.close()
        counts = (weights["gene"]
                  .value_counts()
                  .rename(tis_name[i])
                  .rename_axis("gene"))
        gene_count_list.append(counts)
        print(f"  [{i+1:02d}/{n_tissues}] {tis_name[i]}: {len(counts)} genes")

    snp_gene_used = pd.concat(gene_count_list, axis=1, join="outer").reset_index()
    snp_gene_used.to_csv(DATA_DIR / "Num_of_rs_used.csv", index=False)

    snp_gene_used_x = snp_gene_used.copy()
    snp_gene_used_x.columns = ["gene"] + [f"x{i+1}" for i in range(n_tissues)]
    snp_gene_used_x.to_csv(DATA_DIR / "Num_of_rs_used_x.csv", index=False)
    print("Step 1 done\n")

    return n_tissues, tis_name


# =============================================
# Step 2: Compute per-tissue mean and SD for each population
# =============================================

def step2():
    db_files  = sorted(ELASTIC_NET_DIR.glob("en_*.db"))
    n_tissues = len(db_files)

    def compute_result_table(pop: str) -> pd.DataFrame:
        pop_dir   = pop_results_dir(pop)
        pop_lower = pop.lower()
        pop_files = sorted(pop_dir.glob(PRED_GLOB))

        if len(pop_files) != n_tissues:
            raise ValueError(
                f"[{pop}] expected {n_tissues} .txt files, found {len(pop_files)}"
            )

        print(f"[{pop}] computing mean/sd ({len(pop_files)} tissues)...")

        ms_list = []
        for i, fpath in enumerate(pop_files, start=1):
            dt = pd.read_csv(fpath, sep="\t", low_memory=False)
            gene_cols = dt.columns[2:]

            ms = pd.DataFrame({
                "gene":        gene_cols,
                f"x{i}.mean": dt[gene_cols].mean().values,
                f"x{i}.sd":   dt[gene_cols].std(ddof=1).values,
            }).set_index("gene")
            ms_list.append(ms)

        rs_used = pd.read_csv(DATA_DIR / "Num_of_rs_used_x.csv", index_col="gene")
        for ms in ms_list:
            rs_used = rs_used.join(ms, how="outer")

        out_path = PKG_RT_DIR / f"result_table_{pop_lower}.csv"
        rs_used.reset_index().to_csv(out_path, index=False)
        print(f"  Saved {out_path.name}")
        return rs_used

    result_tables = {}
    for pop in POPULATIONS:
        result_tables[pop] = compute_result_table(pop)

    print("\nStep 2 done\n")
    return result_tables


# =============================================
# Step 3: Pack into parquet files
# =============================================

def step3(result_tables: dict | None = None):
    db_files  = sorted(ELASTIC_NET_DIR.glob("en_*.db"))
    n_tissues = len(db_files)

    if result_tables is None:
        print("[Step 3] result_tables not provided — loading from CSV files...")
        result_tables = {}
        for pop in POPULATIONS:
            pop_lower = pop.lower()
            csv_path  = PKG_RT_DIR / f"result_table_{pop_lower}.csv"
            if not csv_path.exists():
                raise FileNotFoundError(
                    f"Missing {csv_path}. Run Step 2 first."
                )
            result_tables[pop] = pd.read_csv(csv_path, index_col="gene")

    for _d in [PKG_TISSUE_DIR, PKG_RT_DIR, PKG_KS_DIR]:
        _d.mkdir(parents=True, exist_ok=True)

    print("[Step 3] Packing all populations to parquet...")

    for pop in POPULATIONS:
        pop_lower = pop.lower()
        pop_files = sorted(pop_results_dir(pop).glob(PRED_GLOB))

        if len(pop_files) != n_tissues:
            print(f"  WARNING: [{pop}] expected {n_tissues} files, found {len(pop_files)} — skipping")
            continue

        print(f"  [{pop}] packing {len(pop_files)} tissue files...")
        for i, fpath in enumerate(pop_files, start=1):
            df  = pd.read_csv(fpath, sep="\t", low_memory=False)
            out = PKG_TISSUE_DIR / f"{pop_lower}_x{i}.parquet"
            df.to_parquet(out, index=False, compression="snappy")

        result_table = result_tables[pop].reset_index()
        rt_out = PKG_RT_DIR / f"result_table_{pop_lower}.parquet"
        result_table.to_parquet(rt_out, index=False, compression="snappy")
        print(f"  [{pop}] result_table/{rt_out.name} done")

    total_mb = sum(p.stat().st_size for d in (PKG_TISSUE_DIR, PKG_RT_DIR, PKG_KS_DIR) for p in d.glob("*.parquet")) / 1e6
    print(f"\nTotal parquet size: {total_mb:.1f} MB")
    if total_mb > 60:
        print("WARNING: exceeds PyPI 60 MB recommended limit; consider splitting into a separate data package.")

    print("\nStep 3 done")
    print("Install: cd PredictMP_pkg && pip install -e .")


# =============================================
# Step 4: Pre-compute KS test (population vs NFE) + NFE.{POP}.mean.diff
# =============================================

def step4():
    db_files  = sorted(ELASTIC_NET_DIR.glob("en_*.db"))
    n_tissues = len(db_files)
    nfe_files = sorted(pop_results_dir(REF_POP).glob(PRED_GLOB))

    PKG_KS_DIR.mkdir(parents=True, exist_ok=True)

    print("[Step 4] Computing KS tests (each population vs NFE) ...")

    for pop in NON_REF_POPS:
        pop_lower = pop.lower()
        pop_files = sorted(pop_results_dir(pop).glob(PRED_GLOB))

        if len(pop_files) != n_tissues or len(nfe_files) != n_tissues:
            print(f"  WARNING: [{pop}] file count mismatch — skipping")
            continue

        print(f"  [{pop}] vs NFE ...")
        ks_frames = []

        for i, (pop_path, nfe_path) in enumerate(zip(pop_files, nfe_files), start=1):
            pop_dt = pd.read_csv(pop_path, sep="\t", low_memory=False)
            nfe_dt = pd.read_csv(nfe_path, sep="\t", low_memory=False)

            common_genes = sorted(set(pop_dt.columns[2:]) & set(nfe_dt.columns[2:]))
            if not common_genes:
                continue

            pop_mat = pop_dt[common_genes].to_numpy(dtype=float)
            nfe_mat = nfe_dt[common_genes].to_numpy(dtype=float)

            pop_means = pop_mat.mean(axis=0)
            nfe_means = nfe_mat.mean(axis=0)

            with np.errstate(divide="ignore", invalid="ignore"):
                mean_diffs = np.where(
                    np.abs(nfe_means) > 1e-15,
                    (pop_means - nfe_means) / nfe_means * 100,
                    np.nan,
                )

            ks_sigs = [
                "Sig." if sp_stats.ks_2samp(pop_mat[:, j], nfe_mat[:, j]).pvalue < KS_THRESHOLD
                else "Nonsig."
                for j in range(len(common_genes))
            ]

            tissue_ks = pd.DataFrame({
                "gene":              common_genes,
                f"x{i}.ks.sig":     ks_sigs,
                f"x{i}.mean.diff":  mean_diffs,
            })
            ks_frames.append(tissue_ks.set_index("gene"))

        if ks_frames:
            ks_df  = ks_frames[0].join(ks_frames[1:], how="outer")
            ks_out = PKG_KS_DIR / f"ks_{pop_lower}_nfe.parquet"
            ks_df.reset_index().to_parquet(ks_out, index=False, compression="snappy")
            print(f"    Saved ks/ks_{pop_lower}_nfe.parquet "
                  f"({ks_df.shape[0]} genes × {ks_df.shape[1]} cols)")

    print("\nStep 4 done")


# =============================================
# Entry point
# =============================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Step 08 — build PredictMP reference data (parquet files)"
    )
    parser.add_argument(
        "--step", nargs="+", type=int,
        choices=[1, 2, 3, 4],
        metavar="N",
        help="Step(s) to run (1 2 3 4). Omit to run all steps."
    )
    args = parser.parse_args()

    steps = set(args.step) if args.step else {1, 2, 3, 4}

    result_tables = None

    if 1 in steps:
        step1()

    if 2 in steps:
        result_tables = step2()

    if 3 in steps:
        step3(result_tables)

    if 4 in steps:
        step4()


# =============================================================================
# NOTES — C7 轉 Python 後需注意的事項
# =============================================================================
#
# 1. sd() / std() 的 ddof 差異
#    R sd() 預設 ddof=1（Bessel's correction）
#    numpy.std() 預設 ddof=0（population std）→ 會算出不同數值！
#    pandas DataFrame.std() 預設 ddof=1，與 R 一致。
#    本腳本已用 .std(ddof=1) 確保一致。
#
# 2. ecdf_fn 的精度
#    R 的 ecdf(x)(perc) = mean(x <= perc)
#    Python 的 np.searchsorted(sorted_ref, query, side="right") / n
#    兩者邏輯相同，結果一致（functions.py 已正確實作）。
#
# 3. Python package 多族群支援
#    rankcal() 現已支援 population 參數（nfe/eas/afr/ami/amr/asj/fin/mid/sas）。
#    parquet 檔命名規則：{pop_lower}_x{i}.parquet, result_table_{pop_lower}.parquet。
#
# 4. PyPI 套件大小限制
#    49 個 parquet 檔 + result_table，可能超過 PyPI 單一檔案 60MB 上限。
#    Step 3 執行後會自動印出總大小並警告。
#    若超過，考慮分拆成獨立 data package。
#
# 5. @lru_cache 記憶體
#    loader.py 用 @lru_cache 快取已載入的 tissue，
#    若一次跑所有 49 tissues 記憶體用量會累積，
#    長期服務部署時注意設置 maxsize 或定期 cache_clear()。
# =============================================================================
