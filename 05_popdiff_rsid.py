"""
Step 05: aligns PrediXcan SNP alleles against gnomAD v4.1 multi-population frequency
data (9 populations: AFR, AMI, AMR, ASJ, EAS, FIN, MID, NFE, SAS), then
identifies population-differentiated loci.

Allele alignment (4 passes per chromosome, run in parallel):
  I.   Direct match        (rsid + ref + eff identical)
  II.  Flip ref/eff        → AF = 1 - AF
  III. Complement strand   (A↔T, C↔G), AF unchanged
  IV.  Complement + flip   → AF = 1 - AF

  gnomAD v4 provides AF directly; no AC computation needed.
  AF is flipped only when eff_allele ends up in the REF slot (Parts II, IV).

Population differentiation:
  Unpooled two-proportion Z-test, each population vs NFE (reference).
  diff_index = 1 if ANY comparison has p < PVALUE_THRESH AND |Δfreq| > FREQ_DIFF_THRESH.

Output:
  gnomad_rsid_info_with_popdiff.csv  — AN + freq + p_value_*_vs_nfe + diff_index
  gnomad_rsid_withoutinfo.csv        — SNPs not found in gnomAD v4
"""

import pandas as pd
import numpy as np
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import logging
from scipy import stats
from tqdm import tqdm


# ── Config ────────────────────────────────────────────────────────────────────

from config import (
    POPULATIONS, REF_POP, CHROMS,
    WEIGHT_FILE, GNOMAD_PARQUET_DIR as GNOMAD_DIR,
    POPDIFF_FILE, UNMATCHED_FILE,
)

PVALUE_THRESH    = 1e-8     # pop-diff Z-test
FREQ_DIFF_THRESH = 0.05     # |Δfreq| vs NFE

# Complement translation table for alleles
_COMP = str.maketrans("ACGT", "TGCA")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)


# ── Data loading ──────────────────────────────────────────────────────────────

def load_weights() -> pd.DataFrame:
    """Return unique (chrom, rsid, ref_allele, eff_allele) from weight.csv."""
    logging.info("Loading PrediXcan weights...")
    comb = pd.read_csv(WEIGHT_FILE, usecols=["rsid", "varID", "ref_allele", "eff_allele"])
    uni = (
        comb.drop_duplicates(subset=["rsid", "ref_allele", "eff_allele"])
            .copy()
            .reset_index(drop=True)
    )
    # Extract chromosome from varID (format: chr20_51637621_T_A_b38)
    uni["chrom"] = uni["varID"].str.split("_").str[0]
    uni = uni.drop(columns=["varID"])
    logging.info(f"Unique SNPs: {len(uni):,}")
    return uni


# ── Per-chromosome allele alignment ──────────────────────────────────────────

def _flip_af(df: pd.DataFrame, af_cols: list) -> pd.DataFrame:
    """Flip allele frequencies: AF(eff_allele) = 1 - AF(ref_allele)."""
    df = df.copy()
    for af in af_cols:
        df[af] = 1.0 - df[af]
    return df

# rem(remaining): intermediate DataFrame of unmatched SNPs at each pass, used to track progress and avoid redundant merges
def align_alleles_chrom(args: tuple) -> pd.DataFrame:
    """
    For one chromosome, align PrediXcan alleles against gnomAD using 4 passes.
    Returns matched rows with original (uni) allele labels preserved.
    """
    chrom, uni_df, gnomad_dir = args

    # gnomad_file = Path(gnomad_dir) / f"gnomad_{chrom}_populations_rsid.csv"
    gnomad_file = Path(gnomad_dir) / f"gnomad_{chrom}_populations_rsid.parquet"
    if not gnomad_file.exists():
        logging.warning(f"Missing: {gnomad_file.name}")
        return pd.DataFrame()

    # Filter uni to this chromosome only
    uni = uni_df[uni_df["chrom"] == chrom][["rsid", "ref_allele", "eff_allele"]].copy()
    if uni.empty:
        return pd.DataFrame()

    an_cols = [f"AN_{p}" for p in POPULATIONS]
    af_cols = [f"AF_{p}" for p in POPULATIONS]
    meta    = ["rsid", "ref_allele", "eff_allele", "CHROM", "POS"]

    # gnomad = pd.read_csv(gnomad_file, dtype={"CHROM": str})
    # gnomad = gnomad.rename(columns={"ID": "rsid", "REF": "ref_allele", "ALT": "eff_allele"})
    # gnomad = gnomad[meta + an_cols + af_cols]
    gnomad = pd.read_parquet(gnomad_file, columns=["ID", "REF", "ALT", "CHROM", "POS"] + an_cols + af_cols)
    gnomad = gnomad.rename(columns={"ID": "rsid", "REF": "ref_allele", "ALT": "eff_allele"})

    results  = []
    matched  = set()

    # ── Part I: direct match ─────────────────────────────────────────────────
    m1 = uni.merge(gnomad, on=["rsid", "ref_allele", "eff_allele"], how="inner")
    results.append(m1)
    matched |= set(m1["rsid"])

    # ── Part II: flip ref/eff ────────────────────────────────────────────────
    rem = uni[~uni["rsid"].isin(matched)]
    gnomad_ii = gnomad.rename(columns={"ref_allele": "eff_allele", "eff_allele": "ref_allele"})
    m2 = rem.merge(gnomad_ii, on=["rsid", "ref_allele", "eff_allele"], how="inner")
    if not m2.empty:
        m2 = _flip_af(m2, af_cols)
    results.append(m2)
    matched |= set(m2["rsid"])

    # ── Part III: complement strand ──────────────────────────────────────────
    rem = uni[~uni["rsid"].isin(matched)].copy()
    rem_comp = rem.assign(
        ref_allele=rem["ref_allele"].str.translate(_COMP),
        eff_allele=rem["eff_allele"].str.translate(_COMP),
    )
    m3_meta = rem_comp.merge(gnomad, on=["rsid", "ref_allele", "eff_allele"], how="inner")
    if not m3_meta.empty:
        # Restore original allele labels
        m3 = (
            rem[rem["rsid"].isin(m3_meta["rsid"])]
            .merge(m3_meta.drop(columns=["ref_allele", "eff_allele"]), on="rsid")
        )
    else:
        m3 = pd.DataFrame()
    results.append(m3)
    matched |= set(m3_meta["rsid"]) if not m3_meta.empty else set()

    # ── Part IV: complement + flip ───────────────────────────────────────────
    rem = uni[~uni["rsid"].isin(matched)].copy()
    rem_comp_flip = rem.assign(
        ref_allele=rem["eff_allele"].str.translate(_COMP),  # comp(eff) → ref slot
        eff_allele=rem["ref_allele"].str.translate(_COMP),  # comp(ref) → eff slot
    )
    m4_meta = rem_comp_flip.merge(gnomad, on=["rsid", "ref_allele", "eff_allele"], how="inner")
    if not m4_meta.empty:
        m4_meta = _flip_af(m4_meta, af_cols)
        m4 = (
            rem[rem["rsid"].isin(m4_meta["rsid"])]
            .merge(m4_meta.drop(columns=["ref_allele", "eff_allele"]), on="rsid")
        )
    else:
        m4 = pd.DataFrame()
    results.append(m4)

    combined = pd.concat([r for r in results if not r.empty], ignore_index=True)
    combined = combined.drop_duplicates(subset=["rsid", "ref_allele", "eff_allele"])
    return combined


# ── Population differentiation ────────────────────────────────────────────────

def compute_pop_diff(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute allele frequencies for all populations and run unpooled Z-tests
    (each population vs NFE).  diff_index = 1 if ANY comparison meets both
    thresholds (p < PVALUE_THRESH and |Δfreq| > FREQ_DIFF_THRESH).
    """
    logging.info("Computing allele frequencies...")
    df = df.copy()

    # Use gnomAD AF directly (already flipped for Part II/IV during alignment)
    for pop in POPULATIONS:
        df[f"freq_{pop}"] = df[f"AF_{pop}"]

    # Drop SNPs where any population has no samples (AN == 0)
    an_cols   = [f"AN_{p}" for p in POPULATIONS]
    freq_cols = [f"freq_{p}" for p in POPULATIONS]
    before = len(df)
    df = df[df[an_cols].gt(0).all(axis=1)].copy()
    removed = before - len(df)
    if removed:
        logging.info(f"Removed {removed:,} SNPs with AN=0 in at least one population")

    logging.info("Running pairwise Z-tests vs NFE...")
    other_pops = [p for p in POPULATIONS if p != REF_POP]

    ref_freq = df[f"freq_{REF_POP}"].values
    ref_an   = df[f"AN_{REF_POP}"].values

    for pop in other_pops:
        f2   = df[f"freq_{pop}"].values
        n2   = df[f"AN_{pop}"].values
        se   = np.sqrt(
            ref_freq * (1 - ref_freq) / ref_an
            + f2 * (1 - f2) / n2
        )
        with np.errstate(divide="ignore", invalid="ignore"):
            Z = np.abs((ref_freq - f2) / se)
        p_value = 2 * stats.norm.sf(Z)
        p_value = np.where(np.isfinite(p_value), p_value, np.nan)

        df[f"p_value_{pop}_vs_{REF_POP}"] = p_value

        sig = (p_value < PVALUE_THRESH) & (np.abs(ref_freq - f2) > FREQ_DIFF_THRESH)
        df[f"diff_index_{pop}"] = sig.astype(int).astype(str)

    return df


# ── Orchestration ─────────────────────────────────────────────────────────────

def run():
    uni = load_weights()
    args = [(chrom, uni, str(GNOMAD_DIR)) for chrom in CHROMS]

    logging.info(f"Aligning alleles across {len(CHROMS)} chromosomes "
                 f"({os.cpu_count()} workers)...")
    chunk_results = []
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as exe:
        for result in tqdm(
            exe.map(align_alleles_chrom, args),
            total=len(CHROMS),
            desc="chromosomes",
        ):
            chunk_results.append(result)

    all_snps = pd.concat(
        [r for r in chunk_results if not r.empty], ignore_index=True
    )
    all_snps = all_snps.drop_duplicates(subset=["rsid", "ref_allele", "eff_allele"])

    n_uni     = len(uni)
    n_matched = len(all_snps)
    logging.info(
        f"Matched: {n_matched:,} / {n_uni:,} ({n_matched / n_uni * 100:.2f}%)"
    )

    # Save unmatched SNPs
    unmatched = uni[~uni["rsid"].isin(all_snps["rsid"])].drop(columns=["chrom"])
    unmatched.to_csv(UNMATCHED_FILE, index=False)
    logging.info(f"Unmatched: {len(unmatched):,} → {UNMATCHED_FILE.name}")

    # Population differentiation
    final = compute_pop_diff(all_snps)

    other_pops    = [p for p in POPULATIONS if p != REF_POP]
    an_cols       = [f"AN_{p}" for p in POPULATIONS]
    freq_cols     = [f"freq_{p}" for p in POPULATIONS]
    p_value_cols  = [f"p_value_{p}_vs_{REF_POP}" for p in other_pops]
    diff_idx_cols = [f"diff_index_{p}" for p in other_pops]
    col_order = (
        ["rsid", "ref_allele", "eff_allele", "CHROM", "POS"]
        + an_cols
        + freq_cols
        + p_value_cols
        + diff_idx_cols
    )
    final = final[col_order]

    final.to_csv(POPDIFF_FILE, index=False)
    logging.info(f"Saved: {POPDIFF_FILE.name}")

    for pop in other_pops:
        n_diff = (final[f"diff_index_{pop}"] == "1").sum()
        logging.info(
            f"diff_index_{pop}=1: {n_diff:,} / {len(final):,} "
            f"({n_diff / len(final) * 100:.2f}%)"
        )


# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    run()
