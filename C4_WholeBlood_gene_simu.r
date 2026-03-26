################ *↓模擬1000個Whole blood基因表現值↓* ################
setwd("C:\\Predixcan_materials\\Predixcan\\code")
## read data
library(data.table)
comb <- fread("../data/weight.csv")
pgno_rsinfo <- fread("../data/predixcan_gnomad_rsid_freq.csv")


##以Whole blood model舉例
blood_model <- comb[!is.na(comb$Whole_Blood),c(1:5,54)]
setkeyv(blood_model, c("rsid", "ref_allele", "eff_allele"))
setkeyv(pgno_rsinfo, c("rsid", "ref_allele", "eff_allele"))
blood_data <- merge.data.table(blood_model, pgno_rsinfo, sort = T)
blood_data <- blood_data[,c(7,1,8,2,3,17)]
sum(is.na(blood_data$freq_eas)) #check

# dosage <- sample(c(0,1,2), 10, prob = c((1-blood_data$freq_eas)^2, 2*blood_data$freq_eas*(1-blood_data$freq_eas), blood_data$freq_eas^2), replace = T)
#### simulate genotype dosage for N=1000 individuals
N <- 1000
dosage <- matrix(NA, nrow(blood_data), N)
for (i in seq_len(nrow(blood_data))) {
  p <- blood_data$freq_eas[i]
  dosage[i, ] <- sample(
    c(0, 1, 2),
    N,
    replace = TRUE,
    prob = c((1 - p)^2, 2 * p * (1 - p), p^2)
  )
}
####

dosage <- matrix(NA,nrow(blood_data),1000, byrow = T)
for (i in 1:nrow(blood_data)) {
  dosage[i,] <- sample(c(0,1,2), 1000, prob = c((1-blood_data$freq_eas[i])^2, 2*blood_data$freq_eas[i]*(1-blood_data$freq_eas[i]), blood_data$freq_eas[i]^2), replace = T)
}

dosage_data <- data.table(blood_data, dosage)
dosage_data$chr <- paste0("chr",dosage_data$chr)
dosage_data <- dosage_data[order(chr)]
dosage_data <- dosage_data[!duplicated(dosage_data[,c("rsid", "ref_allele", "eff_allele")]),]

chr_group <- split(dosage_data, dosage_data$chr)

output_folder <- "../weight/package_setting/package_test/"
for (i in names(chr_group)) {
  output_file <- paste0(output_folder, i, ".dosage.txt")
  fwrite(chr_group[[i]], output_file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
}

######################################
##linux壓縮成gz檔
cd ../weight/package_setting/package_test/
#linux
for i in {1..22}
do 
gzip chr${i}.dosage.txt
done

mkdir phenotype_file//N1000
mv chr*.dosage.txt.gz phenotype_file

##predict
conda create --name pp python=2.7
source activate myenv
conda install numpy
conda deactivate

mkdir results
# ./PrediXcan.py --predict --weights en_Whole_Blood.db --dosages phenotype_file/N1000 --samples sample_file1000.txt --pheno phenotype_file1000.txt --output_prefix results/
./PrediXcan.py --predict --weights ./data/dosage/weight_db/en_Whole_Blood.db --dosages ./weight/package_setting/package_test/phenotype_file --samples ./weight/package_setting/package_test/phenotype_file/sample_file1000.txt --pheno ./weight/package_setting/package_test/phenotype_file/phenotype_file1000.txt --output_prefix results/

./PrediXcan.py \
  --predict \
  --weights ./data/dosage/weight_db/en_Whole_Blood.db \
  --dosages ./weight/package_setting/package_test/phenotype_file \
  --samples sample_file1000.txt \
  --pheno phenotype_file1000.txt \
  --output_prefix results/N1000


#######################################

##查看預測是否有常態趨勢
# pred <- fread("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/weight/package_setting/package_test/results/_predicted_expression.txt")
pred <- fread("../weight/package_setting/package_test/results/_predicted_expression.txt")
par(mfrow = c(3, 2))
plot(density(pred$ENSG00000145244.11))
plot(density(pred$ENSG00000000457.13))
plot(density(pred$ENSG00000001561.6))
plot(density(pred$ENSG00000145020.15))
plot(density(pred$ENSG00000144635.8))
plot(density(pred$ENSG00000143740.14))

View(blood_model[gene=="ENSG00000144635.8",])


################ *↓探討模擬N=500 or 2000是否有差異↓* ################
##500
setwd("~/Predixcan_materials/Predixcan/code")
library(data.table)
comb <- fread("../data/weight.csv")
pgno_rsinfo <- fread("../data/predixcan_gnomad_rsid_freq.csv")


##以Whole blood model舉例
blood_model <- comb[!is.na(comb$Whole_Blood),c(1:5,54)]
setkeyv(blood_model, c("rsid", "ref_allele", "eff_allele"))
setkeyv(pgno_rsinfo, c("rsid", "ref_allele", "eff_allele"))
blood_data <- merge.data.table(blood_model, pgno_rsinfo, sort = T)


names(blood_data) # list
blood_data <- blood_data[,c(7,1,8,2,3,17)]
sum(is.na(blood_data$freq_eas)) #check
#--------------------------------fix below--------------------------------#
#Hardy–Weinberg equilibrium
# dosage <- sample(c(0,1,2), 10, prob = c((1-blood_data$freq_eas)^2, 2*blood_data$freq_eas*(1-blood_data$freq_eas), blood_data$freq_eas^2), replace = T)
# dosage <- matrix(NA,nrow(blood_data),500, byrow = T)
# for (i in 1:nrow(blood_data)) {
#   dosage[i,] <- sample(c(0,1,2), 500, prob = c((1-blood_data$freq_eas[i])^2, 2*blood_data$freq_eas[i]*(1-blood_data$freq_eas[i]), blood_data$freq_eas[i]^2), replace = T)
# }
set.seed(123)

## 建立 dosage matrix
dosage <- matrix(NA, nrow(blood_data), 500)

for (i in 1:nrow(blood_data)) {

  p <- blood_data$freq_eas[i]

  dosage[i,] <- sample(
    c(0,1,2),
    size = 500,
    replace = TRUE,
    prob = c((1-p)^2, 2*p*(1-p), p^2)
  )
}

#--------------------------------fix above--------------------------------#

## merge data
dosage_data <- data.table(blood_data, dosage)
dosage_data$chr <- paste0("chr",dosage_data$chr)
dosage_data <- dosage_data[order(chr)]
dosage_data <- dosage_data[!duplicated(dosage_data[,c("rsid", "ref_allele", "eff_allele")]),]
# head(dosage_data)
# dim(dosage_data)
chr_group <- split(dosage_data, dosage_data$chr)

# output_folder <- "/Users/chgsh14414/Desktop/Mac/Predixcan/統整/weight/package_setting/package_test/"

path <- "../weight/package_setting/package_test"
if (!dir.exists(path)) {
  dir.create(path, recursive = TRUE)
}
output_folder <- "../weight/package_setting/package_test/"  
for (i in names(chr_group)) {
  output_file <- paste0(output_folder, i, "_500.dosage.txt")
  fwrite(chr_group[[i]], output_file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
}

########################
##linux壓縮成gz檔
# cd /Users/chgsh14414/Desktop/Mac/Predixcan/統整/weight/package_setting/package_test/
cd ../weight/package_setting/package_test/
for i in {1..22}
do 
gzip chr${i}_500.dosage.txt
done

mkdir phenotype_file/N500
mv chr*.dosage.txt.gz phenotype_file/N500

##predict
conda create --name pp python=2.7
source activate myenv
conda install numpy
conda deactivate

mkdir results
#./PrediXcan.py --predict --weights en_Whole_Blood.db --dosages phenotype_file --samples sample_file500.txt --pheno phenotype_file500.txt --output_prefix results/
### P.S --samples's path: dosages_path + samples_path
cp ./sample_file500.txt ./weight/package_setting/package_test/phenotype_file/N500/

./PrediXcan.py \
  --predict \
  --weights ./data/dosage/weight_db/en_Whole_Blood.db \
  --dosages ./weight/package_setting/package_test/phenotype_file/N500 \
  --samples ./sample_file500.txt \
  --pheno ./data/dosage/Whole_Blood/phenotype_file500.txt \
  --output_prefix results/N500

###########################

##查看預測是否有常態趨勢
# pred <- fread("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/weight/package_setting/package_test/results/_predicted_expression.txt")
pred <- fread("../results/N500_predicted_expression.txt")
pdf("../results/density_plot.pdf", width=10, height=8)

par(mfrow=c(3,2))

plot(density(pred$ENSG00000145244.11))
plot(density(pred$ENSG00000000457.13))
plot(density(pred$ENSG00000001561.6))
plot(density(pred$ENSG00000145020.15))
plot(density(pred$ENSG00000144635.8))
plot(density(pred$ENSG00000143740.14))

dev.off()
#linux command to convert pdf to png
pdftoppm -png density_plot.pdf density_plot

##########################################################
## 2000
library(data.table)
# comb <- fread("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/weight.csv")
# pgno_rsinfo <- fread("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/predixcan_gnomad_rsid_info.csv")
comb <- fread("../data/weight.csv")
pgno_rsinfo <- fread("../data/predixcan_gnomad_rsid_info.csv")

##以Whole blood model舉例
blood_model <- comb[!is.na(comb$Whole_Blood),c(1:5,54)]
setkeyv(blood_model, c("rsid", "ref_allele", "eff_allele"))
setkeyv(pgno_rsinfo, c("rsid", "ref_allele", "eff_allele"))
blood_data <- merge.data.table(blood_model, pgno_rsinfo, sort = T)

blood_data <- blood_data[,c(7,1,8,2,3,17)]
sum(is.na(blood_data$freq_eas)) #check
dosage <- sample(c(0,1,2), 10, prob = c((1-blood_data$freq_eas)^2, 2*blood_data$freq_eas*(1-blood_data$freq_eas), blood_data$freq_eas^2), replace = T)

dosage <- matrix(NA,nrow(blood_data),2000, byrow = T)
for (i in 1:nrow(blood_data)) {
  dosage[i,] <- sample(c(0,1,2), 2000, prob = c((1-blood_data$freq_eas[i])^2, 2*blood_data$freq_eas[i]*(1-blood_data$freq_eas[i]), blood_data$freq_eas[i]^2), replace = T)
}

dosage_data <- data.table(blood_data, dosage)
dosage_data$chr <- paste0("chr",dosage_data$chr)
dosage_data <- dosage_data[order(chr)]
dosage_data <- dosage_data[!duplicated(dosage_data[,c("rsid", "ref_allele", "eff_allele")]),]

chr_group <- split(dosage_data, dosage_data$chr)
# output_folder <- "/Users/chgsh14414/Desktop/Mac/Predixcan/統整/weight/package_setting/package_test/"
output_folder <- "../weight/package_setting/package_test/"
for (i in names(chr_group)) {
  output_file <- paste0(output_folder, i, "_2000.dosage.txt")
  fwrite(chr_group[[i]], output_file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
}

#############################
##linux壓縮成gz檔
# cd /Users/chgsh14414/Desktop/Mac/Predixcan/統整/weight/package_setting/package_test/
cd ../weight/package_setting/package_test/
#linux
for i in {1..22}
do 
gzip chr${i}_2000.dosage.txt
done

mv chr*.dosage.txt.gz phenotype_file

##pedict
conda create --name pp python=2.7
source activate myenv
conda install numpy
conda deactivate

./PrediXcan.py --predict --weights en_Whole_Blood.db --dosages phenotype_file/N2000 --samples sample_file2000.txt --pheno phenotype_file2000.txt --output_prefix results/
#############################

##查看預測是否有常態趨勢
# pred <- fread("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/weight/package_setting/package_test/results/_predicted_expression.txt")
pred <- fread("../weight/package_setting/package_test/results/_predicted_expression.txt")
par(mfrow = c(3, 2))
plot(density(pred$ENSG00000145244.11))
plot(density(pred$ENSG00000000457.13))
plot(density(pred$ENSG00000001561.6))
plot(density(pred$ENSG00000145020.15))
plot(density(pred$ENSG00000144635.8))
plot(density(pred$ENSG00000143740.14))

###read N=500, 1000, 2000
library(data.table)
# pred500 <- fread("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/weight/package_setting/package_test/results/N500_predicted_expression.txt")
# pred1000 <- fread("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/weight/package_setting/package_test/results/N1000_predicted_expression.txt")
# pred2000 <- fread("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/weight/package_setting/package_test/results/N2000_predicted_expression.txt")
pred500 <- fread("../results/N500_predicted_expression.txt")
pred1000 <- fread("../results/N1000_predicted_expression.txt")
pred2000 <- fread("../results/N2000_predicted_expression.txt")

gene_list <- colnames(pred500)[-c(1,2)]
p_500_1000 <- c()
for (i in 1:7252) {
  p_500_1000[i] <- ks.test(pred500[,eval(parse(text = gene_list[i]))], pred1000[,eval(parse(text = gene_list[i]))])$p.value
}
summary(p_500_1000)
sum(p_500_1000 < 0.05/7252) #2

p_1000_2000 <- c()
for (i in 1:7252) {
  p_1000_2000[i] <- ks.test(pred1000[,eval(parse(text = gene_list[i]))], pred2000[,eval(parse(text = gene_list[i]))])$p.value
}
summary(p_1000_2000)
sum(p_1000_2000 < 0.05/7252) #0

################ *↑探討模擬N=1000 or 2000是否有差異↑* ################
