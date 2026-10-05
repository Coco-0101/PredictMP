#!/usr/bin/env bash
# Step 02: dbSNP (GRCh38) -> rsid annotation of gnomAD v4.1 joint -> per-chr population CSV
#
#   conda activate bcf_env        # bcftools / bgzip / tabix
#   bash 02_gnomad_process.sh     # every stage is skipped if its output exists
#
# Stages a-d come from the 2026-05 notes / VCF header records.
# Stage e (bcftools query) was not recorded; it is reconstructed from the CSV column order.
set -euo pipefail

CODE_DIR=$(cd "$(dirname "$0")" && pwd)
ROOT=${PREDIXCAN_ROOT:-$(dirname "$(dirname "$CODE_DIR")")}
DATA_DIR=$ROOT/PredictMP/data
GNOMAD_RAW_DIR=/mnt/data_35T/kinomoto/gnomad_joint_freq
DBSNP_RAW=$HOME/GCF_000001405.40                     # dbSNP b156 VCF (RefSeq contig names)
THREADS=${THREADS:-8}
POPS="afr ami amr asj eas fin mid nfe sas"

DBSNP_RENAMED=$DATA_DIR/GCF_000001405.40_chr1_22_autosome.vcf
DBSNP_FIXED=$DATA_DIR/GCF_000001405.40_chr1_22_fixed.vcf.gz
DBSNP_SPLIT=$DATA_DIR/GCF_000001405.40_chr1_22_autosome_split.vcf.gz
HEADER=$DATA_DIR/all_chr_header.txt
ANNOT_DIR=$DATA_DIR/gnomad_joint_v4.1_annotated
INFO_DIR=$DATA_DIR/gnomad_joint_v4.1_info
RSID_DIR=$DATA_DIR/gnomad_joint_v4.1_info_with_rsid
mkdir -p "$ANNOT_DIR" "$INFO_DIR" "$RSID_DIR"

log() { echo "[$(date '+%F %T')] $*"; }

# a. RefSeq contig (NC_0000xx) -> chr1..chr22, keep autosomes only
if [[ ! -s $DBSNP_RENAMED ]]; then
  log "a. rename dbSNP contigs"
  awk 'BEGIN{OFS="\t"
    split("NC_000001.11 NC_000002.12 NC_000003.12 NC_000004.12 NC_000005.10 NC_000006.12 NC_000007.14 NC_000008.11 NC_000009.12 NC_000010.11 NC_000011.10 NC_000012.12 NC_000013.11 NC_000014.9 NC_000015.10 NC_000016.10 NC_000017.11 NC_000018.10 NC_000019.10 NC_000020.11 NC_000021.9 NC_000022.11", nc, " ")
    for (i = 1; i <= 22; i++) map[nc[i]] = "chr" i }
    /^#/ {print; next}
    ($1 in map) {$1 = map[$1]; print}' "$DBSNP_RAW" > "$DBSNP_RENAMED"
fi

# b. contig lines were not renamed in the header -> add ##contig=<ID=chrN>
if [[ ! -s $DBSNP_FIXED ]]; then
  log "b. add contig header"
  for i in $(seq 1 22); do echo "##contig=<ID=chr$i>"; done > "$HEADER"
  bcftools annotate -h "$HEADER" "$DBSNP_RENAMED" -Oz -o "$DBSNP_FIXED" --threads "$THREADS"
fi

# c. dbSNP stores multi-allelic ALT on one line, gnomAD is biallelic -> split
if [[ ! -s $DBSNP_SPLIT.tbi ]]; then
  log "c. split multi-allelic dbSNP records"
  bcftools norm -m -any "$DBSNP_FIXED" -Oz -o "$DBSNP_SPLIT" --threads "$THREADS"
  tabix -p vcf "$DBSNP_SPLIT"
fi

per_chrom() {
  local chr=$1
  local raw=$GNOMAD_RAW_DIR/gnomad.joint.v4.1.sites.$chr.vcf.bgz
  local annot=$ANNOT_DIR/gnomad.joint.v4.1.sites.${chr}_rsid_annotated.vcf.gz
  local info=$INFO_DIR/gnomad_${chr}_populations.csv
  local rsid=$RSID_DIR/gnomad_${chr}_populations_rsid.csv

  # d. copy dbSNP ID into gnomAD (matched on CHROM, POS, REF, ALT)
  [[ -s $annot ]] || bcftools annotate -a "$DBSNP_SPLIT" -c ID "$raw" -Oz -o "$annot"

  # e. (reconstructed) extract AC/AN/AF per population
  if [[ ! -s $info ]]; then
    local fmt="%CHROM,%POS,%ID,%REF,%ALT" hdr="CHROM,POS,ID,REF,ALT"
    for p in $POPS; do
      fmt+=",%INFO/AC_joint_$p,%INFO/AN_joint_$p,%INFO/AF_joint_$p"
      hdr+=",AC_$p,AN_$p,AF_$p"
    done
    { echo "$hdr"; bcftools query -f "$fmt\n" "$annot"; } > "$info"
  fi

  # f. keep variants with an rsid
  [[ -s $rsid ]] || awk -F, 'NR == 1 || $3 != "."' "$info" > "$rsid"
  echo "$chr done"
}
export -f per_chrom
export GNOMAD_RAW_DIR ANNOT_DIR INFO_DIR RSID_DIR DBSNP_SPLIT POPS

log "d-f. per-chromosome annotate / query / filter ($THREADS parallel)"
seq 1 22 | sed 's/^/chr/' | xargs -P "$THREADS" -I{} bash -c 'per_chrom {}'
log "done"
