"""
Step 06: generate sample_file for PrediXcan (NO HEADER)
  -> sample_file500_{pop}.txt   (one per population)

Parameters
----------
output_path : str or Path
    Output file path
populations : list
    e.g. ["EAS", "AFR"]
n_per_pop : int or list
    number of samples per population
cohort_name : str
    e.g. TWB
start_index : int
    starting index for sample naming
"""

from pathlib import Path

def generate_sample_file(
    output_path,
    populations,
    n_per_pop,
    cohort_name="SIM",
    start_index=1
):
    """
    Generate sample_file for PrediXcan (NO HEADER)
    """

    # 確保資料夾存在
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # 如果是單一數字 → 每族群同樣數量
    if isinstance(n_per_pop, int):
        n_per_pop = [n_per_pop] * len(populations)

    assert len(populations) == len(n_per_pop), "populations and n_per_pop length mismatch"

    current_id = start_index

    with open(output_path, "w") as f:
        for pop, n in zip(populations, n_per_pop):
            for _ in range(n):
                sample_id = f"{pop}{current_id}"

                # FID IID cohort ancestry sex
                line = f"{sample_id}\t{sample_id}\t{cohort_name}\t{pop}\t0\n"
                f.write(line)

                current_id += 1

    print(f"Sample file generated → {output_path}")


N_SAMPLES = 500     # simulated subjects per population; must match 07

if __name__ == "__main__":
    import argparse
    from config import POPULATIONS, sample_file

    ap = argparse.ArgumentParser(description="Step 06: write PrediXcan sample files")
    ap.add_argument("--pop", nargs="+", default=POPULATIONS, choices=POPULATIONS,
                    help="populations (default: all)")
    args = ap.parse_args()

    for pop in args.pop:
        generate_sample_file(
            output_path=sample_file(pop, N_SAMPLES),
            populations=[pop.upper()],
            n_per_pop=N_SAMPLES,
        )
