"""
PredictAP — Percentile rank calculator for PrediXcan predictions
==================================================================

Compute ECDF-based percentile ranks of PrediXcan predicted gene expression
against a reference population (default: NFE).

Core functions:
  - ecdf_fn(reference, query)        ECDF percentile rank calculator
  - rankcal(data, tissue)            Main function; requires reference data loaded
  - build_reference_data(...)        One-time setup to prepare reference data
  - load_reference_data(...)         Load reference data from disk

Usage
-----
  python -c "
from PredictAP_code_python import rankcal, build_reference_data
import pandas as pd

# One-time setup:
build_reference_data(
    result_table_path='result_table.csv',
    tissue_results_dir='dosage/results/'
)

# Regular use:
data = pd.read_csv('your_predicted_expression.txt', sep='\t')
result = rankcal(data, tissue='x32')  # x32 = Lung
result.to_csv('rank_result.csv', index=False)
result.head(10)
  "
"""

import pandas as pd
import numpy as np
from pathlib import Path
import logging
import sys


# ─────────────────────────────────────────────
# Internal reference data loader
# ─────────────────────────────────────────────

_DATA_DIR = Path(__file__).resolve().parent / "data"

_tissue_cache: dict[str, pd.DataFrame] = {}
_result_table_cache: pd.DataFrame | None = None


def _load_result_table() -> pd.DataFrame:
    global _result_table_cache
    if _result_table_cache is None:
        path = _DATA_DIR / "result_table.csv"
        if not path.exists():
            raise FileNotFoundError(
                f"Reference result_table not found at {path}.\n"
                "Please run build_reference_data() first to generate it."
            )
        _result_table_cache = pd.read_csv(path)
    return _result_table_cache


def _load_tissue_ref(tissue: str) -> pd.DataFrame:
    """Lazy-load reference population data for one tissue (e.g. 'x32')."""
    if tissue not in _tissue_cache:
        path = _DATA_DIR / f"{tissue}.csv"
        if not path.exists():
            raise FileNotFoundError(
                f"Reference data for tissue '{tissue}' not found at {path}.\n"
                "Please run build_reference_data() first."
            )
        _tissue_cache[tissue] = pd.read_csv(path)
    return _tissue_cache[tissue]


# ─────────────────────────────────────────────
# Core functions
# ─────────────────────────────────────────────

def ecdf_fn(reference: pd.Series | np.ndarray,
            query: pd.Series | np.ndarray) -> np.ndarray:
    """
    ECDF-based percentile rank.

    For each value in `query`, compute the fraction of `reference` values
    that are less than or equal to it (i.e. empirical CDF evaluated at query).

    Parameters
    ----------
    reference : array-like  Reference population values (e.g. NFE expression)
    query     : array-like  Query subject values

    Returns
    -------
    np.ndarray of percentile ranks in [0, 1]
    """
    ref = np.sort(np.asarray(reference, dtype=float))
    q   = np.asarray(query, dtype=float)
    return np.searchsorted(ref, q, side="right") / len(ref)


def rankcal(data: pd.DataFrame, tissue: str) -> pd.DataFrame:
    """
    Compute ECDF percentile ranks for predicted gene expression.

    Compares each subject's PrediXcan-predicted expression against the
    European (NFE) reference population distribution for the specified tissue.

    Parameters
    ----------
    data   : pd.DataFrame
        Predicted expression table with columns:
        FID | IID | ENSG000... | ENSG000... | ...
        (output of PrediXcan --predict)
    tissue : str
        Tissue index, e.g. "x1" (Adipose_Subcutaneous) through "x49" (Whole_Blood).
        See tissue_index.csv for the full mapping.

    Returns
    -------
    pd.DataFrame with columns:
        FID, IID, gene, expression, percentile.rank,
        N.snps, mean, sd, ks.test.pval, European.Asian.mean.difference

    Examples
    --------
    >>> data = pd.read_csv("example_data.csv")
    >>> result = rankcal(data, tissue="x1")
    >>> result.head()
    """
    data         = data.copy()
    ref_data     = _load_tissue_ref(tissue)
    result_table = _load_result_table()

    # ── align to common genes ──────────────────────────────────────────────
    common_genes = data.columns[2:].intersection(ref_data.columns[2:])
    data_sub = data[["FID", "IID"] + list(common_genes)].reset_index(drop=True)
    ref_sub  = ref_data[["FID", "IID"] + list(common_genes)].reset_index(drop=True)

    # ── compute percentile rank for every gene ─────────────────────────────
    rank_values = {
        gene: ecdf_fn(ref_sub[gene].values, data_sub[gene].values)
        for gene in common_genes
    }
    rank_df = pd.DataFrame(rank_values)
    rank_df.insert(0, "IID", data_sub["IID"])
    rank_df.insert(0, "FID", data_sub["FID"])

    # ── long format: expression ────────────────────────────────────────────
    long_expr = data_sub.melt(
        id_vars=["FID", "IID"],
        var_name="gene",
        value_name="expression"
    )

    # ── long format: percentile rank ───────────────────────────────────────
    long_rank = rank_df.melt(
        id_vars=["FID", "IID"],
        var_name="gene",
        value_name="percentile.rank"
    )

    # ── merge expression + rank ────────────────────────────────────────────
    long_merged = long_expr.merge(long_rank, on=["FID", "IID", "gene"])

    # ── pull metadata from result_table for this tissue ────────────────────
    meta_cols = (
        ["gene", tissue]
        + [c for c in result_table.columns if c.startswith(f"{tissue}.")]
    )
    meta = (
        result_table[meta_cols]
        .dropna(subset=[tissue])
    )

    # ── final merge & column rename ────────────────────────────────────────
    result = meta.merge(long_merged, on="gene")

    col_order = (
        ["FID", "IID", "gene", "expression", "percentile.rank"]
        + [c for c in meta_cols if c != "gene"]
    )
    result = result[col_order]

    result = result.rename(columns={
        tissue:                    "N.snps",
        f"{tissue}.mean":          "mean",
        f"{tissue}.sd":            "sd",
        f"{tissue}.ks.test":       "ks.test.pval",
        f"{tissue}.nfe.eas.diff":  "European.Asian.mean.difference",
    })

    result = result.sort_values(["FID", "IID"]).reset_index(drop=True)
    return result


# ─────────────────────────────────────────────
# Reference data builder (run once)
# ─────────────────────────────────────────────

def build_reference_data(
    result_table_path: str | Path,
    tissue_results_dir: str | Path,
    output_dir: str | Path | None = None,
) -> None:
    """
    Build and save reference data files required by rankcal().

    This is the Python equivalent of the use_data() call in the R package —
    it packages the NFE reference population predictions and result_table
    into the data/ directory so rankcal() can load them at runtime.

    Parameters
    ----------
    result_table_path  : path to result_table.csv (SNP counts + mean/sd per tissue)
    tissue_results_dir : directory containing x1.txt ... x49.txt
                         (NFE PrediXcan prediction outputs)
    output_dir         : where to save output CSVs (default: package data/ dir)

    Run once after generating reference data:
    >>> build_reference_data(
    ...     result_table_path="result_table_v2.csv",
    ...     tissue_results_dir="dosage/results/",
    ... )
    """
    out_dir = Path(output_dir) if output_dir else _DATA_DIR
    out_dir.mkdir(parents=True, exist_ok=True)

    # save result_table
    rt = pd.read_csv(result_table_path)
    rt.to_csv(out_dir / "result_table.csv", index=False)
    logging.info(f"✓ Saved result_table → {out_dir / 'result_table.csv'}")

    # save x1 … x49
    tissue_dir = Path(tissue_results_dir)
    txt_files  = sorted(tissue_dir.glob("*.txt"))
    if not txt_files:
        raise FileNotFoundError(f"No .txt files found in {tissue_dir}")

    for idx, path in enumerate(txt_files, start=1):
        df = pd.read_csv(path, sep="\t")
        out_path = out_dir / f"x{idx}.csv"
        df.to_csv(out_path, index=False)
        if idx % 10 == 1 or idx == len(txt_files):
            logging.info(f"  Saved x{idx} / {len(txt_files)} ({path.name})")

    # save tissue index
    tissue_names = [p.stem for p in txt_files]
    tissue_index = pd.DataFrame({
        "index":  range(1, len(txt_files) + 1),
        "tissue": tissue_names,
    })
    tissue_index.to_csv(out_dir / "tissue_index.csv", index=False)
    logging.info(f"✓ Saved tissue_index → {out_dir / 'tissue_index.csv'}")
    logging.info("✓ build_reference_data complete.")


def load_reference_data(tissue_results_dir: str | Path, result_table_path: str | Path) -> None:
    """
    Convenience function to build reference data (shorter alias for build_reference_data).
    """
    build_reference_data(result_table_path, tissue_results_dir)


# ─────────────────────────────────────────────
# Command-line interface
# ─────────────────────────────────────────────

def main():
    """Command-line entry point."""
    import argparse

    parser = argparse.ArgumentParser(
        description="PredictAP — compute ECDF percentile ranks for PrediXcan predictions",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # One-time: build reference data
  python PredictAP_code_python.py --build-ref \\
    --result-table result_table.csv \\
    --tissue-dir dosage/results/

  # Regular use: compute percentile ranks
  python PredictAP_code_python.py \\
    --input predicted_expression.txt \\
    --tissue x32 \\
    --output rank_result.csv
        """
    )

    parser.add_argument("--build-ref", action="store_true",
                        help="Build reference data (one-time setup)")
    parser.add_argument("--result-table", type=str,
                        help="Path to result_table.csv (for --build-ref)")
    parser.add_argument("--tissue-dir", type=str,
                        help="Directory containing x1.txt...x49.txt (for --build-ref)")
    parser.add_argument("--input", type=str,
                        help="Input PrediXcan predicted expression file")
    parser.add_argument("--tissue", type=str, default="x1",
                        help="Tissue index (x1...x49, default: x1)")
    parser.add_argument("--output", type=str, default=None,
                        help="Output CSV file (default: print to stdout)")
    parser.add_argument("--quiet", action="store_true",
                        help="Suppress log messages")

    args = parser.parse_args()

    # Configure logging
    if not args.quiet:
        logging.basicConfig(
            level=logging.INFO,
            format="%(levelname)s: %(message)s"
        )

    # Build reference data mode
    if args.build_ref:
        if not args.result_table or not args.tissue_dir:
            parser.error("--build-ref requires --result-table and --tissue-dir")
        build_reference_data(args.result_table, args.tissue_dir)
        return

    # Compute percentile ranks mode
    if not args.input:
        parser.error("--input is required (or use --build-ref for setup)")

    logging.info(f"Loading input: {args.input}")
    data = pd.read_csv(args.input, sep="\t")
    logging.info(f"  {data.shape[0]} subjects × {data.shape[1]} columns")

    logging.info(f"Computing percentile ranks (tissue={args.tissue})...")
    result = rankcal(data, tissue=args.tissue)
    logging.info(f"  Result: {result.shape[0]} rows × {result.shape[1]} columns")

    if args.output:
        result.to_csv(args.output, index=False)
        logging.info(f"✓ Saved → {args.output}")
    else:
        print(result.to_csv(index=False))


if __name__ == "__main__":
    main()