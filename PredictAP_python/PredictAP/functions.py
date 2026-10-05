"""
functions.py
------------
Core functions for PredictAP.

Output columns of rankcal() match the paper (Chan et al. 2024, Briefings in
Bioinformatics 25(6) bbae549):

    FID, IID, gene, expression, percentile.rank,
    N.snps, mean, sd, ks.test, NFE.EAS.mean.diff
"""

import numpy as np
import pandas as pd
from .loader import load_tissue, load_result_table, load_ks_table


def ecdf_fn(reference: np.ndarray | pd.Series,
            query:     np.ndarray | pd.Series) -> np.ndarray:
    """
    Compute the empirical CDF percentile rank.

    Equivalent to R's ecdf(reference)(query): for each value in *query*,
    returns the fraction of *reference* values <= it.

    Parameters
    ----------
    reference : array-like
        Reference distribution for one gene.
    query : array-like
        Values to rank against the reference.

    Returns
    -------
    np.ndarray
        Percentile ranks in [0, 1].  Multiply by 100 for percentage.

    Examples
    --------
    >>> ecdf_fn([1, 2, 3, 4, 5], [3])
    array([0.6])
    """
    reference  = np.asarray(reference, dtype=float)
    query      = np.asarray(query,     dtype=float)
    ref_sorted = np.sort(reference)
    return np.searchsorted(ref_sorted, query, side="right") / len(ref_sorted)


def rankcal(data: pd.DataFrame,
            tissue:     str = "x1",
            population: str = "nfe") -> pd.DataFrame:
    """
    Compute gene expression percentile ranks — full PredictAP output.

    Implements the algorithm described in Chan et al. (2024):
      1. Align user's PrediXcan output with the chosen reference population.
      2. Compute ECDF-based percentile rank per subject × gene.
      3. Attach reference metadata: N.snps, mean, sd.
      4. Attach pre-computed {population} vs NFE comparison: ks.test, NFE.{POP}.mean.diff
         (absent when population="nfe").

    Parameters
    ----------
    data : pd.DataFrame
        PrediXcan predicted expression for query subjects.
        - Column 1 : FID  (family ID)
        - Column 2 : IID  (individual ID)
        - Columns 3+: ENSG gene IDs
    tissue : str, optional
        Tissue index, e.g. ``"x1"`` (Adipose_Subcutaneous) ... ``"x49"``
        (Whole_Blood).  Default ``"x1"``.
    population : str, optional
        Reference population (lowercase).  One of:
        ``"nfe"``, ``"eas"``, ``"afr"``, ``"ami"``, ``"amr"``,
        ``"asj"``, ``"fin"``, ``"mid"``, ``"sas"``.
        Default ``"nfe"``.

    Returns
    -------
    pd.DataFrame
        Long-format table — one row per subject × gene — with columns:

        ==================  =============================================
        FID                 Family ID from the user's data
        IID                 Individual ID from the user's data
        gene                Ensembl gene ID
        expression          Predicted gene expression (PrediXcan output)
        percentile.rank     Position within the reference distribution
                            (0–100; ~50 = typical; near 0/100 = outlier)
        N.snps              Variants used in predicting this gene
        mean                Reference population mean for this gene
        sd                  Reference population SD for this gene
        ks.test             {population} vs NFE KS test: ``"Sig."`` / ``"Nonsig."``
                            (threshold 5×10⁻⁸; pre-computed at build time)
        NFE.{POP}.mean.diff (POP_mean − NFE_mean) / NFE_mean × 100 (%)
                            Negative → POP expression lower than NFE
        ==================  =============================================
    """
    data = pd.DataFrame(data).copy()
    id_vars = [data.columns[0], data.columns[1]]          # FID, IID

    # ── 1. Load reference ────────────────────────────────────────────────────
    ref_data     = load_tissue(tissue, population)
    result_table = load_result_table(population)

    query_genes  = data.columns[2:].tolist()
    ref_genes    = set(ref_data.columns[2:])
    common_genes = [g for g in query_genes if g in ref_genes]

    if not common_genes:
        raise ValueError(
            f"No gene columns overlap between data and reference "
            f"(tissue={tissue!r}, population={population!r}). "
            "Check that your data uses Ensembl gene IDs (ENSG...)."
        )

    data_sub = data[id_vars + common_genes].reset_index(drop=True)
    ref_sub  = ref_data[["FID", "IID"] + common_genes].reset_index(drop=True)

    # ── 2. Percentile ranks (ECDF) ───────────────────────────────────────────
    rank_matrix = pd.DataFrame(index=data_sub.index,
                               columns=common_genes, dtype=float)
    for gene in common_genes:
        ranks = ecdf_fn(ref_sub[gene].values, data_sub[gene].values)
        rank_matrix[gene] = np.round(ranks * 100, decimals=1)

    long_expr = data_sub.melt(id_vars=id_vars,
                              var_name="gene", value_name="expression")
    long_rank = pd.concat([data_sub[id_vars], rank_matrix], axis=1).melt(
        id_vars=id_vars, var_name="gene", value_name="percentile.rank"
    )
    long_merged = long_expr.merge(long_rank, on=id_vars + ["gene"])

    # ── 3. N.snps / mean / sd from result_table ──────────────────────────────
    meta_cols = (["gene", tissue] +
                 [c for c in result_table.columns if c.startswith(f"{tissue}.")])
    meta   = result_table[meta_cols].dropna(subset=[tissue])
    result = meta.merge(long_merged, on="gene")

    col_order = id_vars + ["gene", "expression", "percentile.rank"] + \
                [c for c in meta_cols if c != "gene"]
    result = result[[c for c in col_order if c in result.columns]]
    result = result.rename(columns={
        tissue:           "N.snps",
        f"{tissue}.mean": "mean",
        f"{tissue}.sd":   "sd",
    })

    # ── 4. KS test (population vs NFE) + NFE.{POP}.mean.diff ────────────────
    # NFE is the baseline; comparing NFE to itself is skipped.
    mean_diff_col = f"NFE.{population.upper()}.mean.diff"
    try:
        ks_table = load_ks_table(population)   # None when population=="nfe"
        if ks_table is not None:
            ks_cols = [c for c in ks_table.columns if c.startswith(f"{tissue}.")]
            if ks_cols:
                ks_sub = ks_table[["gene"] + ks_cols].copy()
                result = result.merge(ks_sub, on="gene", how="left")
                result = result.rename(columns={
                    f"{tissue}.ks.sig":   "ks.test",
                    f"{tissue}.mean.diff": mean_diff_col,
                })
                if mean_diff_col in result.columns:
                    result[mean_diff_col] = (
                        result[mean_diff_col].round(1).astype(str) + "%"
                    )
    except FileNotFoundError:
        import warnings
        warnings.warn(
            f"KS table not found for population={population!r}. "
            "Columns 'ks.test' and 'NFE.{POP}.mean.diff' will be absent. "
            "Run code/08_build_package_data.py --step 4 to generate the KS parquet files.",
            UserWarning, stacklevel=2
        )

    return result.sort_values(id_vars).reset_index(drop=True)
