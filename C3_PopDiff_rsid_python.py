import dask.dataframe as dd
import numpy as np
from scipy.stats import norm
import os
from dask.distributed import Client
import logging
import glob
import subprocess

"""
nporc  # cpu cores
ls -lh PATH/TO/FILE  # file size

partition = cpu cores * (2 ~ 4)
blocksize = file size ÷ partition

e.g. file size: 2.4GB, cpu cores: 8
2.4GB ÷ 24 ≈ 100MB
--> blocksize ≈ 100MB ~ 128MB
"""

# Merge partitioned CSV files into a single CSV with header
def merge_dask_csv(output_dir, output_file):
    files = sorted(glob.glob(os.path.join(output_dir, "*")))

    if not files:
        raise ValueError(f"No files found in {output_dir}")

    print(f"Found {len(files)} partition files")

    first_file = files[0]

    # write header
    subprocess.run(
        f"head -n 1 {first_file} > {output_file}",
        shell=True,
        check=True
    )

    # merge all files (skip header) using tail and cat
    files_str = " ".join(files)
    subprocess.run(
        f"tail -n +2 -q {files_str} >> {output_file}",
        shell=True,
        check=True
    )

    print(f"Merged file saved to: {output_file}")



# Linux sort command to sort CSV by rsid (first column)
def sort_csv(input_file, output_file):
    print(f"Sorting file: {input_file}")

    subprocess.run(
        f"(head -n 1 {input_file} && tail -n +2 {input_file} | sort -t, -k1,1) > {output_file}",
        shell=True,
        check=True
    )

    print(f"Sorted file saved to: {output_file}")

def run_analysis():
    logging.getLogger("distributed").setLevel(logging.ERROR) 
    client = Client(n_workers=55, threads_per_worker=1, memory_limit='10GB')
    print(f"Dask Dashboard: {client.dashboard_link}")

    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    DATA_DIR = os.path.abspath(os.path.join(BASE_DIR, "../data"))

    # snp_weight: PrediXcan weight file
    snp_weight = dd.read_csv(
        os.path.join(DATA_DIR, "weight.csv"),
        blocksize="512MB",
        dtype={'rsid': 'string', 'ref_allele': 'string', 'eff_allele': 'string'}
    )
    # gnomad_data: gnomAD dataset
    gnomad_data = dd.read_csv(
        os.path.join(DATA_DIR, "allchr_gnomad_info.csv"),
        blocksize="512MB",
        dtype={'rsid': 'string', 'A1': 'string', 'A2': 'string'}
    )

    # Remove duplicate SNPs by varID and keep rsid and allele definitions
    uni_snp_weight = snp_weight.drop_duplicates(subset=['varID'])[
        ['rsid', 'ref_allele', 'eff_allele']
    ].persist()

    gnomad_data = gnomad_data.rename(columns={
        "A1": "ref_allele",
        "A2": "eff_allele"
    }).persist()

    # DNA complement mapping (strand flip)
    comp_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    def flip_series(s):
        return s.replace(comp_map)

    # swap: ref <-> effect allele
    gnomad_swap = gnomad_data.rename(
        columns={'ref_allele': 'eff_allele', 'eff_allele': 'ref_allele'}
    )
    # Adjust allele counts after swapping
    for p in ['nfe', 'afr', 'eas']:
        gnomad_swap[f'AC_{p}'] = gnomad_swap[f'AN_{p}'] - gnomad_swap[f'AC_{p}']
    gnomad_swap = gnomad_swap.persist()

    # flip: strand complement
    uni_flip = uni_snp_weight.copy()
    uni_flip['ref_allele'] = uni_flip['ref_allele'].map_partitions(flip_series)
    uni_flip['eff_allele'] = uni_flip['eff_allele'].map_partitions(flip_series)
    uni_flip = uni_flip.persist() 

    # Matching (4 cases)
    # Case 1: direct match
    m1 = uni_snp_weight.merge(gnomad_data, on=['rsid','ref_allele','eff_allele'])
    m1['match_type'] = 1

    # Case 2: ref/eff swapped
    m2 = uni_snp_weight.merge(gnomad_swap, on=['rsid','ref_allele','eff_allele'])
    m2['match_type'] = 2

    # Case 3: strand flipped
    # flip the alleles back to original after merging to maintain consistency in downstream analysis
    m3 = uni_flip.merge(gnomad_data, on=['rsid','ref_allele','eff_allele'])
    m3['ref_allele'] = m3['ref_allele'].map_partitions(flip_series)
    m3['eff_allele'] = m3['eff_allele'].map_partitions(flip_series)
    m3['match_type'] = 3

    # Case 4: flipped + swapped
    m4 = uni_flip.merge(gnomad_swap, on=['rsid','ref_allele','eff_allele'])
    m4['ref_allele'] = m4['ref_allele'].map_partitions(flip_series)
    m4['eff_allele'] = m4['eff_allele'].map_partitions(flip_series)
    m4['match_type'] = 4



    # merge all matches and keep best match (priority: 1 > 2 > 3 > 4)
    merged_all = dd.concat([m1, m2, m3, m4]).persist()

    # Keep best match (priority: 1 > 2 > 3 > 4)
    merged = merged_all.drop_duplicates(subset=['rsid'], keep='first').persist()

    # Population diff
    def analyze_population_diff(df):
        df = df[df['AN_eas'] > 0].copy()
        df['freq_nfe'] = df['AC_nfe'] / df['AN_nfe']
        df['freq_afr'] = df['AC_afr'] / df['AN_afr']
        df['freq_eas'] = df['AC_eas'] / df['AN_eas']

        var_nfe = (df['freq_nfe'] * (1 - df['freq_nfe'])) / df['AN_nfe']
        var_eas = (df['freq_eas'] * (1 - df['freq_eas'])) / df['AN_eas']
        z = (df['freq_nfe'] - df['freq_eas']).abs() / np.sqrt(var_nfe + var_eas + 1e-20)

        df['pval_ne'] = norm.sf(z) * 2
        mask = (df['pval_ne'] < 1e-8) & ((df['freq_nfe'] - df['freq_eas']).abs() > 0.05)
        df['diff_index'] = mask.astype(int)
        return df

    popdiff_df = merged.map_partitions(analyze_population_diff) 

    # --------------------------------------------------------------
    print("Writing output...")

    cols = [
        "rsid", "ref_allele", "eff_allele",
        "chr", "POS",
        "AC_nfe", "AN_nfe",
        "AC_afr", "AN_afr",
        "AC_eas", "AN_eas",
        "freq_nfe", "freq_afr", "freq_eas",
        "pval_ne", "diff_index"
    ]
    popdiff_df = popdiff_df[cols]
    popdiff_df = popdiff_df.repartition(partition_size="200MB")
    popdiff_dir = os.path.join(DATA_DIR, "popdiff_output")
    popdiff_df.to_csv(popdiff_dir, index=False)
    popdiff_file = os.path.join(DATA_DIR, "predixcan_gnomad_rsid_freq.csv")
    merge_dask_csv(popdiff_dir, popdiff_file)
    # sorted_file = os.path.join(DATA_DIR, "predixcan_gnomad_rsid_freq.csv")
    # sort_csv(output_file, sorted_file)

    # Missing SNPs
    print("Finding missing SNPs...")
    missing = uni_snp_weight.merge(
        merged[['rsid']].drop_duplicates(),
        on='rsid',
        how='left',
        indicator=True
    )
    missing = missing[missing['_merge'] == 'left_only'].drop(columns=['_merge'])
    missing = missing.repartition(partition_size="200MB")
    missing_dir = os.path.join(DATA_DIR, "missing_output")
    missing.to_csv(missing_dir, index=False)
    missing_file = os.path.join(DATA_DIR, "predixcan_gnomad_rsid_withoutinfo.csv")
    merge_dask_csv(missing_dir, missing_file)
    # sorted_file = os.path.join(DATA_DIR, "predixcan_gnomad_rsid_withoutinfo.csv")
    # sort_csv(missing_file, sorted_file)

    print("Done!")
    client.close()

if __name__ == '__main__':
    run_analysis()