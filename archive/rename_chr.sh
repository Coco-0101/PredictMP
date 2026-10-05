pv ~/GCF_000001405.40 | awk '
BEGIN{
    OFS="\t"

    map["NC_000001.11"]="chr1"
    map["NC_000002.12"]="chr2"
    map["NC_000003.12"]="chr3"
    map["NC_000004.12"]="chr4"
    map["NC_000005.10"]="chr5"
    map["NC_000006.12"]="chr6"
    map["NC_000007.14"]="chr7"
    map["NC_000008.11"]="chr8"
    map["NC_000009.12"]="chr9"
    map["NC_000010.11"]="chr10"
    map["NC_000011.10"]="chr11"
    map["NC_000012.12"]="chr12"
    map["NC_000013.11"]="chr13"
    map["NC_000014.9"]="chr14"
    map["NC_000015.10"]="chr15"
    map["NC_000016.10"]="chr16"
    map["NC_000017.11"]="chr17"
    map["NC_000018.10"]="chr18"
    map["NC_000019.10"]="chr19"
    map["NC_000020.11"]="chr20"
    map["NC_000021.9"]="chr21"
    map["NC_000022.11"]="chr22"

}

/^#/ {
    print
    next
}

{
    if ($1 in map) {
        $1 = map[$1]
    }

    print
}
' > ~/Predixcan_materials/Predixcan/data/GCF_000001405.40_chr1_22_rename.vcf