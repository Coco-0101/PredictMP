import pandas as pd
import numpy as np
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import os
import shutil

# ✅ 平行寫檔 function（移到外面）
def write_group(args):
    chrom, group, output_base = args
    output_path = output_base / f"{chrom}_500.dosage.txt.gz"
    
    group.to_csv(
        output_path, 
        sep='\t', 
        index=False, 
        header=False, 
        compression={'method': 'gzip', 'compresslevel': 1}  # ✅ 加速
    )
    return output_path.name


def generate_dosages():

    BASE_DIR = Path(__file__).resolve().parent
    DATA_DIR = BASE_DIR.parent / "data"
    OUTPUT_BASE = (BASE_DIR.parent / "weight" / "package_setting" / "package_test" / "phenotype_file" / "N500")
    OUTPUT_BASE.mkdir(parents=True, exist_ok=True)

    RESULTS_DIR = BASE_DIR.parent / "results"
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    snp_weight = pd.read_csv(DATA_DIR / "weight.csv")
    gnomad_snp_freq = pd.read_csv(DATA_DIR / "predixcan_gnomad_rsid_freq.csv")

    # Filter and Merge (Whole Blood model example)
    # Delete rows with NA in Whole_Blood and select relevant columns
    blood_model = snp_weight.dropna(subset=['Whole_Blood'])[
    ['gene','rsid','varID','ref_allele','eff_allele','Whole_Blood']]

    # Merge on variant identifiers
    blood_data = pd.merge(
        blood_model, 
        gnomad_snp_freq, 
        on=['rsid', 'ref_allele', 'eff_allele'], # insure matching on these key columns
        how='inner' # only keep SNPs present in both datasets
    ).sort_values(by=['rsid', 'ref_allele', 'eff_allele'])

    # print(blood_data.columns)
    """
    Index(
    ['gene', 'rsid', 'varID', 'ref_allele', 'eff_allele', 'Whole_Blood',
     'chr', 'POS', 'AC_nfe', 'AN_nfe', 'AC_afr', 'AN_afr', 'AC_eas',
     'AN_eas', 'freq_nfe', 'freq_afr', 'freq_eas', 'pval_ne', 'diff_index'],
      dtype='str')
    """

    # Select necessary columns for PrediXcan dosage format
    cols_to_keep = ['chr', 'rsid', 'POS', 'ref_allele', 'eff_allele', 'freq_eas']
    blood_data = blood_data[cols_to_keep].copy()
    print(blood_data.columns) # Check column names and order

    # Hardy–Weinberg Equilibrium (HWE) Dosage Generation
    # High-performance vectorized generation for 500 samples
    print("Simulating dosage matrix via HWE...")
    np.random.seed(123)
    n_snps = len(blood_data)
    n_samples = 500
    p = blood_data['freq_eas'].values.reshape(-1, 1) # for EAS allele frequency

    # Calculate genotype probabilities
    prob_0 = (1 - p)**2 # P(AA) = (1-p)^2
    prob_1 = 2 * p * (1 - p) # P(Aa) = 2p(1-p)

    # Generate random matrix and assign genotypes (0, 1, 2)
    rand_mat = np.random.rand(n_snps, n_samples)
    dosage = np.zeros((n_snps, n_samples), dtype=int)
    dosage[rand_mat > prob_0] = 1 # Assign heterozygous (1) where random number exceeds P(AA)
    dosage[rand_mat > (prob_0 + prob_1)] = 2 # P(aa) where random number exceeds P(AA) + P(Aa)

    # Combine dosage with variant info
    dosage_df = pd.DataFrame(dosage, index=blood_data.index)
    dosage_df = pd.concat([blood_data, dosage_df], axis=1) 
    # Ensure chromosome format (e.g., 1 -> chr1)
    dosage_df['chr'] = dosage_df['chr'].apply(lambda x: f"chr{x}" if not str(x).startswith('chr') else x)
    dosage_df = dosage_df.sort_values('chr')
    dosage_df = dosage_df.drop_duplicates(subset=['rsid', 'ref_allele', 'eff_allele']) # Drop duplicates to prevent PrediXcan errors
    # print(dosage_df.head())
    
    # Export by Chromosome with Gzip compression
    print("Writing GZipped dosage files...")
    # turn in to list of tuples for multiprocessing
    groups = [(chrom, group, OUTPUT_BASE) 
              for chrom, group in dosage_df.groupby('chr')]
    # multiprocessing
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
        results = list(executor.map(write_group, groups))
    for name in results:
        print(f"Created: {name}")

    # Copy sample file
    sample_file = BASE_DIR.parent / "data" / "sample_file500.txt"
    print(f"Checking for sample file at: {sample_file}")
    if not sample_file.exists():
        print("Sample file NOT found, skipping copy.")
    else:
        try:
            dest = OUTPUT_BASE / "sample_file500.txt"
            shutil.copy(sample_file, dest)
            print(f"Sample file copied successfully → {dest}")
        except Exception as e:
            print(f"Error copying sample file: {e}")

def run_predixcan():
    """
    Instructions for running PrediXcan.py
    """
    # Define paths for the command
    BASE_DIR = Path(__file__).resolve().parent
    
    # cmd = [
    #     "./PrediXcan.py",
    #     "--predict",
    #     "--weights", "./data/dosage/weight_db/en_Whole_Blood.db",
    #     "--dosages", "../weight/package_setting/package_test/phenotype_file/N500",
    #     "--samples", "./sample_file500.txt",
    #     "--pheno", "./data/dosage/Whole_Blood/phenotype_file500.txt",
    #     "--output_prefix", "results/N500"
    # ]

    cmd = [
    "python",
    "./PrediXcan.py",
    "--predict",
    "--weights", "../data/dosage/weight_db/en_Whole_Blood.db",
    "--dosages", "../weight/package_setting/package_test/phenotype_file/N500",
    "--samples", "./sample_file500.txt",
    "--pheno", "../data/dosage/Whole_Blood/phenotype_file500.txt",
    "--output_prefix", "../results/N500"
    ]
    
    subprocess.run(cmd, check=True, cwd=BASE_DIR)
    # print("\n[INFO] Run the following command in your terminal:")
    # print(" ".join(cmd))

if __name__ == "__main__":
    generate_dosages()
    run_predixcan()