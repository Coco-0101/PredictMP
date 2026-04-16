import pandas as pd
import numpy as np
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import shutil
import logging
from tqdm import tqdm
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR.parent / "data"
OUTPUT_DIR = BASE_DIR.parent / "data" /"dosage_all_tissues_AFR"
RESULTS_DIR = BASE_DIR.parent / "results"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

N_SAMPLES = 500
np.random.seed(123)

# ===============================
# Logging
# ===============================
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)

# ===============================
# Step 1: Load data
# ===============================
def load_data():
    logging.info("Loading data...")
    comb = pd.read_csv(DATA_DIR / "weight.csv")
    gnomad = pd.read_csv(DATA_DIR / "predixcan_gnomad_rsid_freq.csv")
    return comb, gnomad


# ===============================
# Step 2: 建 tissue models
# ===============================
def build_tissue_models(comb):
    logging.info("Building tissue models...")
    tissue_models = {}

    for tis in comb.columns.tolist()[5:]:
        df = comb.dropna(subset=[tis])[
            ['gene','rsid','varID','ref_allele','eff_allele', tis]
        ]
        tissue_models[tis] = df

    return tissue_models


# ===============================
# Step 3: merge gnomad
# ===============================
def merge_with_gnomad(args):
    tis, model, gnomad = args

    merged = pd.merge(
        model,
        gnomad,
        on=['rsid','ref_allele','eff_allele'],
        how='inner'
    )

    # merged = merged[['chr','rsid','POS','ref_allele','eff_allele','freq_eas']] # for EAS frequency
    merged = merged[['chr','rsid','POS','ref_allele','eff_allele','freq_afr']] # for AFR frequency


    return tis, merged


# ===============================
# Step 4: 模擬 dosage（超快）
# ===============================
def simulate_dosage(args):
    tis, df = args

    n_snps = len(df)
    # p = df['freq_eas'].values.reshape(-1,1) # for EAS frequency
    p = df['freq_afr'].values.reshape(-1,1) # for AFR frequency


    prob0 = (1-p)**2
    prob1 = 2*p*(1-p)

    rand = np.random.rand(n_snps, N_SAMPLES)

    dosage = np.zeros((n_snps, N_SAMPLES), dtype=np.int8)
    dosage[rand > prob0] = 1
    dosage[rand > (prob0 + prob1)] = 2

    dosage_df = pd.concat(
        [df.reset_index(drop=True),
         pd.DataFrame(dosage)],
        axis=1
    )

    dosage_df['chr'] = dosage_df['chr'].apply(
        lambda x: f"chr{x}" if not str(x).startswith("chr") else x
    )

    dosage_df = dosage_df.sort_values('chr')

    dosage_df = dosage_df.drop_duplicates(
        subset=['rsid','ref_allele','eff_allele']
    )

    return tis, dosage_df


# ===============================
# Step 5: 寫檔
# ===============================
def write_tissue(args):
    tis, df = args

    TISSUE_DIR = OUTPUT_DIR / tis
    TISSUE_DIR.mkdir(parents=True, exist_ok=True)

    outputs = []

    for chrom, group in df.groupby('chr'):
        output_file = TISSUE_DIR / f"{tis}_{chrom}.dosage.txt.gz"

        group.to_csv(
            output_file,
            sep='\t',
            index=False,
            header=False,
            compression={'method':'gzip','compresslevel':1}
        )

        outputs.append(output_file.name)

    sample_file = DATA_DIR / "sample_file500_AFR.txt"

    if sample_file.exists():
        dest = TISSUE_DIR / "sample_file500_AFR.txt"
        shutil.copy(sample_file, dest)

    return tis, outputs



# ===============================
# 主：產生 dosage（全部 tissues）
# ===============================
def generate_dosages_all_tissues():

    comb, gnomad = load_data()
    tissue_models = build_tissue_models(comb)

    logging.info("Merging with gnomad...")
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as exe:
        merged_results = list(tqdm(
            exe.map(
                merge_with_gnomad,
                [(k,v,gnomad) for k,v in tissue_models.items()]
            ),
            total=len(tissue_models)
        ))

    merged_dict = dict(merged_results)

    logging.info("Simulating dosage...")
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as exe:
        dosage_results = list(tqdm(
            exe.map(simulate_dosage, merged_dict.items()),
            total=len(merged_dict)
        ))

    dosage_dict = dict(dosage_results)

    logging.info("Writing files...")
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as exe:
        write_results = list(tqdm(
            exe.map(write_tissue, dosage_dict.items()),
            total=len(dosage_dict)
        ))

    for tis, files in write_results:
        logging.info(f"{tis}: {len(files)} files written")


# ===============================
# PrediXcan 單 tissue
# ===============================
def run_predixcan_single(tissue):

    cmd = [
        "python",
        "./PrediXcan.py",
        "--predict",
        "--weights", f"../../elastic_net_models/en_{tissue}.db",
        "--dosages", f"../data/dosage_all_tissues_AFR/{tissue}",
        "--dosages_prefix", f"{tissue}_chr",
        "--samples", "./sample_file500_AFR.txt",
        "--pheno", "../data/phenotype_file500.txt",
        "--output_prefix", f"../results/{tissue}_AFR"
    ]

    subprocess.run(cmd, check=True, cwd=BASE_DIR)

    return tissue



# PrediXcan
def run_predixcan_all_tissues():

    comb = pd.read_csv(DATA_DIR / "weight.csv")
    tissues = comb.columns.tolist()[5:]

    logging.info("Running PrediXcan...")

    with ProcessPoolExecutor(max_workers=os.cpu_count()) as exe:
        list(tqdm(
            exe.map(run_predixcan_single, tissues),
            total=len(tissues)
        ))

    logging.info("PrediXcan finished")


# plot (buggy, 待修正)
def plot_density(file_path, output_path):
    pred = pd.read_csv(file_path, sep='\t')

    genes = pred.columns[1:7]

    plt.figure(figsize=(10,8))

    for i, gene in enumerate(genes,1):
        plt.subplot(3,2,i)
        data = pred[gene].dropna()

        density = gaussian_kde(data)
        x = np.linspace(data.min(), data.max(), 1000)

        plt.plot(x, density(x))
        plt.title(gene)

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


# ===============================
def analyze_results_all():

    logging.info("Analyzing results...")

    for file in RESULTS_DIR.glob("*_predicted_expression.txt"):
        output = RESULTS_DIR / f"{file.stem}.png"
        plot_density(file, output)


# ===============================
if __name__ == "__main__":
    # generate_dosages_all_tissues()
    # run_predixcan_all_tissues()
    analyze_results_all()