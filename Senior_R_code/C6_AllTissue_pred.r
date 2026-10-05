library(data.table)
setwd("~/Predixcan_materials/Predixcan/code")
# comb <- fread("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/weight.csv")
# pgno_rsinfo <- fread("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/predixcan_gnomad_rsid_freq.csv")
comb <- fread("../data/weight.csv")
pgno_rsinfo <- fread("../data/predixcan_gnomad_rsid_freq.csv")

#把name有-的去掉
# names(comb)[c(24, 28)] <- c("Brain_Spinal_cord_cervical_c1", "Cells_EBV_transformed_lymphocytes")
names(comb)[names(comb) == "Brain_Spinal_cord_cervical_c-1"] <- "Brain_Spinal_cord_cervical_c1"
names(comb)[names(comb) == "Cells_EBV-transformed_lymphocytes"] <- "Cells_EBV_transformed_lymphocytes"
tis_name <- names(comb)
for (i in 6:length(tis_name)) {
  a <- tis_name[i]
  assign(paste0(tis_name[i],"_model"), comb[!is.na(comb[,eval(parse(text=tis_name[i]))]),.(gene, rsid, varID, ref_allele, eff_allele, eval(parse(text = tis_name[i])))])
}
#checkk <- comb[!is.na(comb$Whole_Blood),c(1:5,54)]
setkeyv(pgno_rsinfo, c("rsid", "ref_allele", "eff_allele"))

##將每個tissue與database資料合併
temp_name1 <- paste0(tis_name[-c(1:5)],"_model")
for (j in 1:length(temp_name1)) {
  setkeyv(eval(parse(text = temp_name1[j])), c("rsid", "ref_allele", "eff_allele"))
  dt <- merge.data.table(eval(parse(text = temp_name1[j])), pgno_rsinfo, sort = T)
  assign(paste0(tis_name[j+5], "_data"), dt[,c(7,1,8,2,3,17)])
}
sum(is.na(Whole_Blood_data$freq_eas)) #check

##模擬dosage
library(parallel)
library(data.table)

temp_name2 <- paste0(tis_name[-c(1:5)], "_data")

cores <- detectCores() - 2   # 使用CPU核心

results <- mclapply(seq_along(temp_name2), function(i){

  dat <- get(temp_name2[i])   # 取得資料
  nnrow <- nrow(dat)

  dosage <- matrix(NA, nnrow, 500)

  for (j in seq_len(nnrow)) {

    p <- dat$freq_eas[j]

    dosage[j, ] <- sample(
      c(0,1,2),
      500,
      prob = c((1-p)^2, 2*p*(1-p), p^2),
      replace = TRUE
    )

  }

  dosage_data <- data.table(dat, dosage)

  dosage_data$chr <- paste0("chr", dosage_data$chr)

  setorder(dosage_data, chr)

  dosage_data <- dosage_data[
    !duplicated(dosage_data[,c("rsid","ref_allele","eff_allele")])
  ]

  list(
    name = paste0(tis_name[i+5], "_dosage"),
    data = dosage_data
  )

}, mc.cores = cores)

# 把結果 assign 回 workspace
for(res in results){
  assign(res$name, res$data)
}
"""
temp_name2 <- paste0(tis_name[-c(1:5)], "_data")
for (i in 1:length(temp_name2)) {
  nnrow <- nrow(eval(parse(text = temp_name2[i])))
  dosage <- matrix(NA,nrow(eval(parse(text = temp_name2[i]))),500, byrow = T)
  for (j in 1:nnrow) {
    dosage[j,] <- sample(c(0,1,2), 500, prob = c((1-eval(parse(text = temp_name2[i]))$freq_eas[j])^2,
     2*eval(parse(text = temp_name2[i]))$freq_eas[j]*(1-eval(parse(text = temp_name2[i]))$freq_eas[j]),
      eval(parse(text = temp_name2[i]))$freq_eas[j]^2), replace = T)
  }
  dosage_data <- data.table(eval(parse(text = temp_name2[i])), dosage)
  dosage_data$chr <- paste0("chr",dosage_data$chr)
  dosage_data <- dosage_data[order(chr)]
  assign(paste0(tis_name[i+5], "_dosage"), dosage_data[!duplicated(dosage_data[,c("rsid", "ref_allele", "eff_allele")]),])
}
"""


##依chr輸出
library(data.table)

temp_name3 <- paste0(tis_name[-c(1:5)], "_dosage")
output_folder <- "../data/dosage_betty/"

for (j in seq_along(temp_name3)) {

  dat <- get(temp_name3[j])

  tissue_folder <- paste0(output_folder, tis_name[j+5])

  # 如果資料夾不存在就建立
  if (!dir.exists(tissue_folder)) {
    dir.create(tissue_folder, recursive = TRUE)
  }

  chr_group <- split(dat, dat$chr)

  for (chr_name in names(chr_group)) {

    output_file <- paste0(
      tissue_folder, "/",
      tis_name[j+5], "_", chr_name, ".dosage.txt"
    )

    fwrite(
      chr_group[[chr_name]],
      output_file,
      col.names = FALSE,
      row.names = FALSE,
      quote = FALSE,
      sep = "\t"
    )
  }
}
"""
temp_name3 <- paste0(tis_name[-c(1:5)], "_dosage")
output_folder <- "../data/dosage_betty/"
for (j in 1:49) {
  chr_group <- split(eval(parse(text = temp_name3[j])), eval(parse(text = temp_name3[j]))$chr)
  for (i in names(chr_group)) {
    output_file <- paste0(output_folder, tis_name[j+5], "/", tis_name[j+5], "_", i, ".dosage.txt")
    fwrite(chr_group[[i]], output_file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
  }
}
"""

##輸出完整檔以備份
output_folder2 <- "../data/dosage_betty/"
for (j in 1:49) {
  output_file <- paste0(output_folder2, tis_name[j+5], ".dosage.txt")
  fwrite(eval(parse(text = temp_name3[j])), output_file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
}

##linux壓縮成gz檔
# cd /Users/chgsh14414/Desktop/Mac/Predixcan/統整/weight/package_setting/package_test/
# find /Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/dosage -type f -name "*chr*" -exec gzip {} \;
cd ../weight/package_setting/package_test/
cd ~/Predixcan_materials/Predixcan/code
# find ../data/dosage_betty -type f -name "*chr*" -exec gzip {} \;
find ../data/dosage_betty -type f -name "*chr*" | xargs pigz


##predict
#conda create --name pp python=2.7
source activate myenv
#conda install numpy
conda deactivate


##PrediXcan for all 49 tissues
# cd /Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/dosage
cd  ~/Predixcan_materials/Predixcan/data/dosage_betty

# Build tissue list
tissues="Adipose_Subcutaneous
Adipose_Visceral_Omentum
Adrenal_Gland
Artery_Aorta
Artery_Coronary
Artery_Tibial
Brain_Amygdala
Brain_Anterior_cingulate_cortex_BA24
Brain_Caudate_basal_ganglia
Brain_Cerebellar_Hemisphere
Brain_Cerebellum
Brain_Cortex
Brain_Frontal_Cortex_BA9
Brain_Hippocampus
Brain_Hypothalamus
Brain_Nucleus_accumbens_basal_ganglia
Brain_Putamen_basal_ganglia
Brain_Spinal_cord_cervical_c1
Brain_Substantia_nigra
Breast_Mammary_Tissue
Cells_Cultured_fibroblasts
Cells_EBV_transformed_lymphocytes
Colon_Sigmoid
Colon_Transverse
Esophagus_Gastroesophageal_Junction
Esophagus_Mucosa
Esophagus_Muscularis
Heart_Atrial_Appendage
Heart_Left_Ventricle
Kidney_Cortex
Liver
Lung
Minor_Salivary_Gland
Muscle_Skeletal
Nerve_Tibial
Ovary
Pancreas
Pituitary
Prostate
Skin_Not_Sun_Exposed_Suprapubic
Skin_Sun_Exposed_Lower_leg
Small_Intestine_Terminal_Ileum
Spleen
Stomach
Testis
Thyroid
Uterus
Vagina
Whole_Blood"


echo "$tissues" | parallel -j 16 '
./PrediXcan.py \
--predict \
--weights ~/Predixcan_materials/elastic_net_models/en_{}.db \
--dosages {} \
--dosages_prefix {}_chr \
--samples ~/Predixcan_materials/Predixcan/sample_file500.txt \
--pheno ~/Predixcan_materials/Predixcan/phenotype_file500.txt \
--output_prefix results/{}
'

"""
for tissue in Adipose_Subcutaneous Adipose_Visceral_Omentum Adrenal_Gland Artery_Aorta Artery_Coronary Artery_Tibial Brain_Amygdala Brain_Anterior_cingulate_cortex_BA24 Brain_Caudate_basal_ganglia Brain_Cerebellar_Hemisphere Brain_Cerebellum Brain_Cortex Brain_Frontal_Cortex_BA9 Brain_Hippocampus Brain_Hypothalamus Brain_Nucleus_accumbens_basal_ganglia Brain_Putamen_basal_ganglia Brain_Spinal_cord_cervical_c1 Brain_Substantia_nigra Breast_Mammary_Tissue Cells_Cultured_fibroblasts Cells_EBV_transformed_lymphocytes Colon_Sigmoid Colon_Transverse Esophagus_Gastroesophageal_Junction Esophagus_Mucosa Esophagus_Muscularis Heart_Atrial_Appendage Heart_Left_Ventricle Kidney_Cortex Liver Lung Minor_Salivary_Gland Muscle_Skeletal Nerve_Tibial Ovary Pancreas Pituitary Prostate Skin_Not_Sun_Exposed_Suprapubic Skin_Sun_Exposed_Lower_leg Small_Intestine_Terminal_Ileum Spleen Stomach Testis Thyroid Uterus Vagina Whole_Blood
do 
# ./PrediXcan.py --predict --weights weight_db/en_"$tissue".db --dosages "$tissue" --dosages_prefix "$tissue"_chr --samples sample_file500.txt --pheno phenotype_file500.txt --output_prefix results/"$tissue"
./PrediXcan.py \
--predict \
--weights ~/Predixcan_materials/elastic_net_models/en_"$tissue".db \
--dosages "$tissue" \
--dosages_prefix "$tissue"_chr \
--samples ~/Predixcan_materials/Predixcan/sample_file500.txt \
--pheno ~/Predixcan_materials/Predixcan/phenotype_file500.txt \
--output_prefix results/"$tissue"
done
"""



##rename prediction file
for f in *.txt; do mv -- "$f" "nfe_$f"; done


################ *↓Note↓* ################
 -found SNPs
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