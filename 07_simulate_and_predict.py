"""
Step 07: for one population, simulate N_SAMPLES genotypes from gnomAD allele
frequencies (HWE) for every tissue model, then run PrediXcan --predict.
  -> data/dosage/dosage_{POP}/{tissue}/{tissue}_chr*.dosage.txt.gz
  -> results/{POP}/{tissue}_predicted_expression.txt

    python 07_simulate_and_predict.py --pop eas
    python 07_simulate_and_predict.py --pop all
"""
import pandas as pd
import numpy as np
import subprocess
from concurrent.futures import ProcessPoolExecutor
import os
import shutil
import logging
import sys
from tqdm import tqdm
from functools import partial

from config import (
    CODE_DIR, POPULATIONS, ELASTIC_NET_DIR,
    WEIGHT_FILE, POPDIFF_FILE, PHENO_FILE, DOSAGE_DIR,
    sample_file, pop_results_dir,
)

N_SAMPLES   = 500    # simulated subjects per population; must match 06
RANDOM_SEED = 123    # each (population, tissue) gets its own stream, see simulate_dosage

# Worker processes read the population from the environment; main() sets it.
POPULATION = os.environ.get("PREDIXCAN_POP", "").lower()
FREQ_COL    = f"freq_{POPULATION}"
SAMPLE_FILE = sample_file(POPULATION, N_SAMPLES) if POPULATION else None
OUTPUT_DIR  = DOSAGE_DIR / f"dosage_{POPULATION.upper()}"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)


# Load data
def load_data():
    logging.info("Loading data...")
    comb   = pd.read_csv(WEIGHT_FILE)
    gnomad = pd.read_csv(POPDIFF_FILE)
    gnomad = gnomad.rename(columns={"CHROM": "chr"})
    return comb, gnomad


# Build tissue models
def build_tissue_models(comb):
    logging.info("Building tissue models...")
    tissue_models = {}
    for tis in comb.columns.tolist()[5:]:
        df = comb.dropna(subset=[tis])[
            ['gene', 'rsid', 'varID', 'ref_allele', 'eff_allele', tis]
        ]
        tissue_models[tis] = df
    return tissue_models


# Merge gnomad
def merge_with_gnomad(args):
    tis, model, gnomad = args
    merged = pd.merge(
        model, gnomad,
        on=['rsid', 'ref_allele', 'eff_allele'],
        how='inner'
    )
    merged = merged[['chr', 'rsid', 'POS', 'ref_allele', 'eff_allele', FREQ_COL]]
    return tis, merged


# Simulate dosage
def simulate_dosage(args):
    tis, df, tis_idx = args
    # Seed per (population, tissue): independent across tissues/populations and
    # reproducible regardless of which worker process runs the job.
    rng = np.random.default_rng([RANDOM_SEED, POPULATIONS.index(POPULATION), tis_idx])
    n_snps = len(df)
    p = df[FREQ_COL].values.reshape(-1, 1)

    prob0 = (1 - p) ** 2
    prob1 = 2 * p * (1 - p)

    rand   = rng.random((n_snps, N_SAMPLES))
    dosage = np.zeros((n_snps, N_SAMPLES), dtype=np.int8)
    dosage[rand > prob0]             = 1
    dosage[rand > (prob0 + prob1)]   = 2

    dosage_df = pd.concat(
        [df.reset_index(drop=True), pd.DataFrame(dosage)],
        axis=1
    )
    dosage_df['chr'] = dosage_df['chr'].apply(
        lambda x: f"chr{x}" if not str(x).startswith("chr") else x
    )
    dosage_df = dosage_df.sort_values('chr')
    dosage_df = dosage_df.drop_duplicates(subset=['rsid', 'ref_allele', 'eff_allele'])
    return tis, dosage_df



# Write files
def write_tissue(args):
    tis, df = args
    TISSUE_DIR = OUTPUT_DIR / tis
    TISSUE_DIR.mkdir(parents=True, exist_ok=True)

    outputs = []
    for chrom, group in df.groupby('chr'):
        output_file = TISSUE_DIR / f"{tis}_{chrom}.dosage.txt.gz"
        if output_file.exists():
            outputs.append(output_file.name)
            continue
        group.to_csv(
            output_file, sep='\t', index=False, header=False,
            compression={'method': 'gzip', 'compresslevel': 1}
        )
        outputs.append(output_file.name)

    # PrediXcan.py looks for --samples inside the dosage directory
    shutil.copy(SAMPLE_FILE, TISSUE_DIR / SAMPLE_FILE.name)

    return tis, outputs



# Generate dosage
def _dosage_exists(tis):
    tis_dir = OUTPUT_DIR / tis
    return tis_dir.is_dir() and any(tis_dir.glob("*.dosage.txt.gz"))


def generate_dosages_all_tissues():
    comb, gnomad = load_data()
    tissue_models = build_tissue_models(comb)

    pending = {k: v for k, v in tissue_models.items() if not _dosage_exists(k)}
    n_skip = len(tissue_models) - len(pending)
    if n_skip:
        logging.info(f"Skipping {n_skip} tissues with existing dosage files")
    if not pending:
        logging.info("All dosage files already exist, nothing to do")
        return

    logging.info("Merging with gnomad...")
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as exe:
        merged_results = list(tqdm(
            exe.map(merge_with_gnomad, [(k, v, gnomad) for k, v in pending.items()]),
            total=len(pending)
        ))
    merged_dict = dict(merged_results)
    tissue_index = {tis: i for i, tis in enumerate(tissue_models, start=1)}   # weight.csv column order

    logging.info("Simulating dosage...")
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as exe:
        dosage_results = list(tqdm(
            exe.map(simulate_dosage, [(tis, df, tissue_index[tis]) for tis, df in merged_dict.items()]),
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



# PrediXcan single tissue
def run_predixcan_single(tissue, population=POPULATION):
    pop_dir = pop_results_dir(population)
    pop_dir.mkdir(parents=True, exist_ok=True)

    predict_file = pop_dir / f"{tissue}_predicted_expression.txt"
    if predict_file.exists():
        logging.info(f"Skipping {tissue} — result already exists")
        return tissue

    cmd = [
        sys.executable, str(CODE_DIR / "PrediXcan.py"),
        "--predict",
        "--weights",        str(ELASTIC_NET_DIR / f"en_{tissue}.db"),
        "--dosages",        str(OUTPUT_DIR / tissue),
        "--dosages_prefix", f"{tissue}_chr",
        "--samples",        SAMPLE_FILE.name,
        "--pheno",          str(PHENO_FILE),
        "--output_prefix",  str(pop_dir / f"{tissue}")
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    print(f"\n===== {population} | {tissue} =====")
    print(result.stdout)
    return tissue



# PrediXcan all tissues
def run_predixcan_all_tissues(population=POPULATION):
    tissues = pd.read_csv(WEIGHT_FILE, nrows=0).columns.tolist()[5:]

    logging.info(f"Running PrediXcan for {population}...")
    func = partial(run_predixcan_single, population=population)

    with ProcessPoolExecutor(max_workers=os.cpu_count()) as exe:
        list(tqdm(exe.map(func, tissues), total=len(tissues)))

    logging.info("PrediXcan finished")


def main():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pop", nargs="+", required=True, choices=POPULATIONS + ["all"])
    pops = ap.parse_args().pop
    pops = POPULATIONS if "all" in pops else pops

    for pop in pops:
        sf = sample_file(pop, N_SAMPLES)
        if not sf.exists():
            sys.exit(f"Missing {sf} — run 06_generate_sample_file.py --pop {pop}")
        n_rows = sum(1 for line in open(sf) if line.strip())
        if n_rows != N_SAMPLES:
            sys.exit(f"{sf} has {n_rows} subjects but N_SAMPLES={N_SAMPLES} — rerun 06 with the same N")
        print(f"\n{'='*50}\nRunning population: {pop.upper()}\n{'='*50}")
        env = os.environ | {"PREDIXCAN_POP": pop}
        subprocess.run([sys.executable, __file__, "--worker"], env=env, check=True)


if __name__ == "__main__":
    if sys.argv[1:] == ["--worker"]:
        logging.info(f"Pipeline:  POPULATION={POPULATION}  |  FREQ_COL={FREQ_COL}  |  N_SAMPLES={N_SAMPLES}")
        generate_dosages_all_tissues()
        run_predixcan_all_tissues()
    else:
        main()


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
