"""
Step 09: run PredictAP on a PrediXcan output file.

Workflow (Figure 1 & Figure 3 in paper):
  Step 1  Obtain predicted gene expression from PrediXcan (GTEx v8 elastic-net)
          for subjects of a target ancestry (e.g. EAS).
  Step 2  Feed the PrediXcan output into PredictAP (rankcal).
          Each gene's predicted expression is ranked against 500 reference
          values derived from the target ancestry's allele frequencies (gnomAD).
  Step 3  Evaluate output:
          - percentile.rank : position within the reference distribution
                              (~50 = typical; near 0/100 = outlier)
          - N.snps          : variants used to predict this gene
          - mean / sd       : reference population statistics
          - ks.test         : EAS vs NFE KS test significance (Sig./Nonsig.)
                              threshold = 5×10⁻⁸ (pre-computed at build time)
          - NFE.EAS.mean.diff: (EAS_mean − NFE_mean) / NFE_mean × 100 (%)
                              negative → EAS expression lower than NFE

Usage
-----
    python 09_run_predictap.py [--tissue x32] [--population eas]
                          [--input path/to/predicted_expression.txt]
                          [--output path/to/result.csv]
"""

import argparse
import sys
import pandas as pd
from pathlib import Path

from config import POPULATIONS, pop_results_dir


# ── tissue index  ───────────────────────────────────────────
TISSUE_INDEX = {
    "x1": "Adipose_Subcutaneous",     "x2": "Adipose_Visceral_Omentum",
    "x3": "Adrenal_Gland",            "x4": "Artery_Aorta",
    "x5": "Artery_Coronary",          "x6": "Artery_Tibial",
    "x7": "Brain_Amygdala",           "x8": "Brain_Anterior_cingulate_cortex_BA24",
    "x9": "Brain_Caudate_basal_ganglia","x10": "Brain_Cerebellar_Hemisphere",
    "x11": "Brain_Cerebellum",         "x12": "Brain_Cortex",
    "x13": "Brain_Frontal_Cortex_BA9", "x14": "Brain_Hippocampus",
    "x15": "Brain_Hypothalamus",       "x16": "Brain_Nucleus_accumbens_basal_ganglia",
    "x17": "Brain_Putamen_basal_ganglia","x18": "Brain_Spinal_cord_cervical_c-1",
    "x19": "Brain_Substantia_nigra",   "x20": "Breast_Mammary_Tissue",
    "x21": "Cells_Cultured_fibroblasts","x22": "Cells_EBV-transformed_lymphocytes",
    "x23": "Colon_Sigmoid",            "x24": "Colon_Transverse",
    "x25": "Esophagus_Gastroesophageal_Junction",
    "x26": "Esophagus_Mucosa",         "x27": "Esophagus_Muscularis",
    "x28": "Heart_Atrial_Appendage",   "x29": "Heart_Left_Ventricle",
    "x30": "Kidney_Cortex",            "x31": "Liver",
    "x32": "Lung",                     "x33": "Minor_Salivary_Gland",
    "x34": "Muscle_Skeletal",          "x35": "Nerve_Tibial",
    "x36": "Ovary",                    "x37": "Pancreas",
    "x38": "Pituitary",                "x39": "Prostate",
    "x40": "Skin_Not_Sun_Exposed_Suprapubic",
    "x41": "Skin_Sun_Exposed_Lower_leg",
    "x42": "Small_Intestine_Terminal_Ileum",
    "x43": "Spleen",                   "x44": "Stomach",
    "x45": "Testis",                   "x46": "Thyroid",
    "x47": "Uterus",                   "x48": "Vagina",
    "x49": "Whole_Blood",
}


def run(input_path: Path,
        tissue:     str,
        population: str,
        output_path: Path | None) -> pd.DataFrame:
    """Execute the full PredictAP pipeline and return the result table."""

    # ── Step 1: load PrediXcan output ─────────────────────────────────────────
    print(f"[Step 1] Reading PrediXcan output: {input_path.name}")
    data = pd.read_csv(input_path, sep="\t", low_memory=False)
    n_subjects = data.shape[0]
    n_genes    = data.shape[1] - 2          # exclude FID, IID
    print(f"         {n_subjects} subjects  ×  {n_genes} genes")

    # ── Step 2: run PredictAP ─────────────────────────────────────────────────
    tissue_name = TISSUE_INDEX.get(tissue, tissue)
    print(f"[Step 2] rankcal(tissue={tissue!r} [{tissue_name}], "
          f"population={population!r})")

    try:
        from PredictAP import rankcal
    except ImportError:
        sys.exit(
            "ERROR: PredictAP package not found.\n"
            "       Install with:  cd PredictAP_python && pip install -e ."
        )

    result = rankcal(data, tissue=tissue, population=population)

    # ── Step 3: display & save ────────────────────────────────────────────────
    print(f"\n[Step 3] Result shape: {result.shape[0]} rows × {result.shape[1]} cols")
    print(f"         Columns: {list(result.columns)}\n")
    print(result.head(10).to_string(index=False))

    # Summary statistics
    print(f"\n── Percentile rank summary (population={population!r}) ──")
    print(result["percentile.rank"].describe().round(2))

    if "ks.test" in result.columns:
        sig_genes   = result[result["ks.test"] == "Sig."]["gene"].nunique()
        total_genes = result["gene"].nunique()
        print(f"\n── KS test ({population.upper()} vs NFE, threshold=5e-8) ──")
        print(f"   Significant genes: {sig_genes} / {total_genes} "
              f"({sig_genes / total_genes * 100:.1f}%)")

    mean_diff_col = f"NFE.{population.upper()}.mean.diff"
    if mean_diff_col in result.columns:
        print(f"\n── {mean_diff_col} (sample) ──")
        print(result[["gene", mean_diff_col]]
              .drop_duplicates("gene")
              .head(5)
              .to_string(index=False))

    if output_path is not None:
        result.to_csv(output_path, index=False)
        print(f"\nSaved: {output_path}")

    return result


def main():
    parser = argparse.ArgumentParser(
        description="PredictAP — evaluate PrediXcan predictions against a "
                    "population reference (Chan et al. 2024)"
    )
    parser.add_argument(
        "--input", type=Path,
        default=pop_results_dir("eas") / "Adipose_Subcutaneous_predicted_expression.txt",
        help="PrediXcan predicted expression file (.txt, tab-separated)"
    )
    parser.add_argument(
        "--tissue", type=str, default="x1",
        choices=list(TISSUE_INDEX.keys()),
        help="Tissue index (default: x1 = Adipose_Subcutaneous)"
    )
    parser.add_argument(
        "--population", type=str, default="eas",
        choices=POPULATIONS,
        help="Reference population (default: eas)"
    )
    parser.add_argument(
        "--output", type=Path, default=None,
        help="Output CSV path (optional; prints to stdout if omitted)"
    )
    args = parser.parse_args()

    run(
        input_path=args.input,
        tissue=args.tissue,
        population=args.population,
        output_path=args.output,
    )


if __name__ == "__main__":
    main()
