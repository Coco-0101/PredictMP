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

    # subprocess.run(cmd, check=True, cwd=BASE_DIR)
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=BASE_DIR
    )

    print(f"\n===== {tissue} =====")
    print(result.stdout)

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

        data = pd.to_numeric(pred[gene], errors='coerce').dropna()

        if len(data) < 2:
            continue

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
    generate_dosages_all_tissues()
    run_predixcan_all_tissues()
    analyze_results_all()


# Note
"""
- found SNPs
1. Adipose_Subcutaneous: 214410/214645(99.89%)
2. Adipose_Visceral_Omentum: 178519/178760(99.87%)
3. Adrenal_Gland: 132813/132957(99.89%)
4. Artery_Aorta: 196533/196728(99.90%)
5. Artery_Coronary: 112951/113130(99.84%)
6. Artery_Tibial: 221667/221893(99.90%)
7. Brain_Amygdala: 86694/86805(99.87%)
8. Brain_Anterior_cingulate_cortex_BA24: 106173/106328(99.85%)
9. Brain_Caudate_basal_ganglia:144574/144730(99.89%)
10. Brain_Cerebellar_Hemisphere: 173469/173681(99.88%)
11. Brain_Cerebellum: 196506/196723(99.89%)
12. Brain_Cortex: 157413/157608(99.88%)
13. Brain_Frontal_Cortex_BA9: 131926/132061(99.90%)
14. Brain_Hippocampus: 107355/107480(99.88%)
15. Brain_Hypothalamus: 105704/105843(99.87%)
16. Brain_Nucleus_accumbens_basal_ganglia: 136818/137024(99.85%)
17. Brain_Putamen_basal_ganglia: 129590/129768(99.86%)
18. Brain_Spinal_cord_cervical_c1: 104023/104122(99.90%)
19. Brain_Substantia_nigra: 83639/83791(99.82%)
20. Breast_Mammary_Tissue: 157680/157911(99.85%)
21. Cells_Cultured_fibroblasts: 236066/236279(99.91%)
22. Cells_EBV_transformed_lymphocytes: 93125/93235(99.88%)
23. Colon_Sigmoid: 159493/159704(99.87%)
24. Colon_Transverse: 159328/159506(99.89%)
25. Esophagus_Gastroesophageal_Junction: 164262/164531(99.84%)
26. Esophagus_Mucosa: 213769/214022(99.88%)
27. Esophagus_Muscularis: 207254/207468(99.90%)
28. Heart_Atrial_Appendage: 170548/170801(99.85%)
29. Heart_Left_Ventricle: 150650/150861(99.86%)
30. Kidney_Cortex: 60627/60744(99.81%)
31. Liver: 105248/105390(99.87%)
32. Lung: 193929/194164(99.88%)
33. Minor_Salivary_Gland: 91606/91751(99.84%)
34. Muscle_Skeletal: 186348/186517(99.91%)
35. Nerve_Tibial: 259192/259465(99.89%)
36. Ovary: 106982/107125(99.87%)
37. Pancreas: 158620/158802(99.89%)
38. Pituitary: 151811/152008(99.87%)
39. Prostate: 118909/119076(99.86%)
40. Skin_Not_Sun_Exposed_Suprapubic: 216036/216250(99.90%)
41. Skin_Sun_Exposed_Lower_leg: 229886/230117(99.90%)
42. Small_Intestine_Terminal_Ileum: 107157/107287(99.88%)
43. Spleen: 161155/161387(99.86%)
44. Stomach: 133176/133323(99.89%)
45. Testis: 264459/264789(99.88%)
46. Thyroid: 243654/243941(99.88%)
47. Uterus: 80099/80231(99.84%)
48. Vagina: 77732/77872(99.82%)
49. Whole_Blood: 176832/177016(99.90%)
"""