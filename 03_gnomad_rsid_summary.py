"""
Step 03 (optional QC): gnomAD v4.1 rsid coverage per chromosome, and
comparison with gnomAD v2.
"""

import subprocess

from config import GNOMAD_INFO_DIR, GNOMAD_RSID_DIR, GNOMAD_V2_FILE


def line_count(path):
    """Return number of data rows (excluding header) in a CSV file."""
    result = subprocess.run(["wc", "-l", path], capture_output=True, text=True)
    total_lines = int(result.stdout.strip().split()[0])
    return total_lines - 1  # subtract header


def main():
    rows = []
    total_orig = 0
    total_rsid = 0

    for chrom in range(1, 23):
        orig_file = GNOMAD_INFO_DIR / f"gnomad_chr{chrom}_populations.csv"
        rsid_file = GNOMAD_RSID_DIR / f"gnomad_chr{chrom}_populations_rsid.csv"

        n_orig = line_count(orig_file)
        n_rsid = line_count(rsid_file)
        ratio  = n_rsid / n_orig if n_orig > 0 else 0.0

        rows.append((chrom, n_orig, n_rsid, ratio))
        total_orig += n_orig
        total_rsid += n_rsid
        print(f"chr{chrom}: total={n_orig:>12,}  with_rsid={n_rsid:>12,}  ratio={ratio:.4f}")

    total_ratio = total_rsid / total_orig if total_orig > 0 else 0.0
    print("-" * 70)
    print(f"ALL:   total={total_orig:>12,}  with_rsid={total_rsid:>12,}  ratio={total_ratio:.4f}")

    # --- comparison with gnomAD v2 ---
    n_v2 = line_count(GNOMAD_V2_FILE)
    print()
    print("=== Comparison with gnomAD v2 (gnomad_v2_allchr.csv) ===")
    print(f"gnomAD v4.1 with rsid (chr1-22): {total_rsid:>12,}")
    print(f"gnomAD v2 all chr:               {n_v2:>12,}")
    print(f"v4.1_rsid / v2:                  {total_rsid / n_v2:.4f}")
    print(f"v2 / v4.1_rsid:                  {n_v2 / total_rsid:.4f}")
    diff = total_rsid - n_v2
    print(f"Difference (v4.1_rsid - v2):     {diff:>12,}")


if __name__ == "__main__":
    main()
