"""
Generate sample_file for PrediXcan (NO HEADER)

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
                line = f"{sample_id}\t{sample_id}\tSIM\t{pop}\t0\n"
                f.write(line)

                current_id += 1

    print(f"Sample file generated → {output_path}")


if __name__ == "__main__":
    BASE_DIR = Path(__file__).resolve().parent
    DATA_DIR = BASE_DIR.parent / "data"

    OUTPUT_DIR = DATA_DIR 
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # ===== 參數設定 =====
    # populations = ["EAS", "AFR"]
    # n_per_pop = 250   # 250 + 250 = 500
    populations = ["AFR"]
    n_per_pop = 500   # 250 + 250 = 500

    output_file = OUTPUT_DIR / "sample_file500_AFR.txt"

    # ===== 執行 =====
    generate_sample_file(
        output_path=output_file,
        populations=populations,
        n_per_pop=n_per_pop,
        cohort_name="TWB"
    )