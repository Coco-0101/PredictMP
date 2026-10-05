"""
Shared paths and populations for the PrediXcan / PredictAP pipeline.
Every numbered script does `from config import ...`.

Root directory defaults to two levels above this file
(/mnt/data_40T/betty0101/Project); override with env var PREDIXCAN_ROOT.
Step-specific parameters (thresholds, N, seed) live at the top of each script.
"""

import os
from pathlib import Path

# ── Paths ────────────────────────────────────────────────────────────────────
CODE_DIR      = Path(__file__).resolve().parent
PROJECT_ROOT  = Path(os.environ.get("PREDIXCAN_ROOT", CODE_DIR.parent.parent)).resolve()
PREDIXCAN_DIR = PROJECT_ROOT / "PredictMP"
DATA_DIR      = PREDIXCAN_DIR / "data"
RESULTS_DIR   = PREDIXCAN_DIR / "results"

ELASTIC_NET_DIR = PROJECT_ROOT / "elastic_net_models"              # en_{tissue}.db

# Raw gnomAD v4.1 joint VCFs (external, read-only)
GNOMAD_RAW_DIR = Path("/mnt/data_35T/kinomoto/gnomad_joint_freq")

# Step 02 (gnomAD processing) outputs
GNOMAD_INFO_DIR      = DATA_DIR / "gnomad_joint_v4.1_info"                  # all variants
GNOMAD_RSID_DIR      = DATA_DIR / "gnomad_joint_v4.1_info_with_rsid"        # ID != "."
GNOMAD_PARQUET_DIR   = DATA_DIR / "gnomad_joint_v4.1_info_with_rsid_parquet"
GNOMAD_V2_FILE       = DATA_DIR / "gnomad_v2_allchr.csv"                    # QC comparison only

# Step 01
WEIGHT_FILE     = DATA_DIR / "weight.csv"
GENE_NUM_FILE   = DATA_DIR / "Num_of_genes_tissue.csv"

# Step 05
POPDIFF_FILE    = DATA_DIR / "gnomad_rsid_info_with_popdiff.csv"
UNMATCHED_FILE  = DATA_DIR / "gnomad_rsid_withoutinfo.csv"

# Step 06 / 07
PHENO_FILE      = DATA_DIR / "phenotype_file500.txt"
DOSAGE_DIR      = DATA_DIR / "dosage"                                       # dosage_{POP}/{tissue}/

# Step 08 (PredictAP package data; loader.py reads the same dirs)
PKG_TISSUE_DIR  = DATA_DIR / "tissue"
PKG_RT_DIR      = DATA_DIR / "result_table"
PKG_KS_DIR      = DATA_DIR / "ks"

# ── Populations ──────────────────────────────────────────────────────────────
POPULATIONS = ["afr", "ami", "amr", "asj", "eas", "fin", "mid", "nfe", "sas"]
REF_POP     = "nfe"

CHROMS = [f"chr{i}" for i in range(1, 23)]


def sample_file(pop: str, n: int) -> Path:
    """PrediXcan --samples file for one population, e.g. sample_file500_eas.txt"""
    return DATA_DIR / f"sample_file{n}_{pop.lower()}.txt"


def pop_results_dir(pop: str) -> Path:
    """PrediXcan predicted expression output, e.g. results/EAS/"""
    return RESULTS_DIR / pop.upper()


if __name__ == "__main__":
    for k, v in dict(globals()).items():
        if k.isupper():
            print(f"{k:18s} {v}")
