###############################
## PrediXcan needed files    ##
###############################

################ *↓資料前處理↓* ################
###Status:已將資料進行imputation，並取出Info score>0.75的位點
###資料包含：(1)GSE33355 lung normal tissue (N=61) ; (2)GSE26853 gastric cancer blood tissue (N=95); 合計：156 samples
##資料位置：/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/imputed/INFO0.75

#Plink
#dosage 先行準備檔案：traw & frq file
lung_data_path=/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/imputed/INFO0.75
dosage_data_path=/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/dosage
cd /Users/chgsh14414/Downloads/plink_mac_20220402
#traw
for i in {1..22}
do 
./plink --bfile $lung_data_path/all_0.75.dose -chr $i --recode A-transpose --out $dosage_data_path/lung_chr$i.dosage
done
#freq
for i in {1..22}
do 
./plink --bfile $lung_data_path/all_0.75.dose --chr ${i} --freq --out $dosage_data_path/lung_chr$i.frq --allow-extra-chr
done

################ *↑資料前處理↑* ################

################ *↓整理PrediXcan所需要的資料↓* ################
#R
library(data.table)
library(parallel)
library("readr")

### (1) Sample file
data_fam = fread("/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/imputed/INFO0.75/all_0.75.dose.fam", header = FALSE)
sample_file = data.table("FID"=data_fam[[1]],
                         "IID"=data_fam[[1]],
                         "POP"="TWB",
                         "Super_POP"="EAS",
                         "SEX"=data_fam[[5]])
write.table(sample_file, file = "/Users/chgsh14414/Desktop/Mac/Predixcan/lung/sample_file.txt", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

### (2) Phenotype file
phenotype_file = data.table("FID"=data_fam[[1]],
                            "IID"=data_fam[[1]],
                            "POP"="TWB",
                            "SEX"=data_fam[[5]],
                            "phenotype"=data_fam[[6]])

write.table(phenotype_file, file = "/Users/chgsh14414/Desktop/Mac/Predixcan/lung/phenotype_file.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)

################ *↓所需資料 - dosage分為使用lung & whole blood tissue預測的↓* ################
##由於資料variant為GRCh37版本，需轉換為PrediXcan所使用的rsid，基於資料量大，先將該tissue會使用到的rsid取出再進行轉換

### (3-1) Dosage file - for Lung tissue
##這裡透過先前已有gnomAD的rsid表，其版本也為GRCh37做為對照表
rs_info <- fread("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/predixcan_gnomad_rsid_info.csv")
comb <- fread("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/weight.csv")
which(colnames(comb)=="Lung")
comb <- comb[,c(1:5,38)]
comb <- comb[complete.cases(comb),]
uni_comb <- comb[!duplicated(comb$varID),c(2,4,5)]
lung_snp <- uni_comb$rsid
lung_rs_info <- subset(rs_info, rsid %in% lung_snp)
colnames(lung_rs_info)[4] <- "CHR"
#fwrite(lung_rs_info, file = paste0("/Users/chgsh14414/Desktop/Mac/Predixcan/lung/lung_rs_info.txt"))
lung_rs_info <- fread("/Users/chgsh14414/Desktop/Mac/Predixcan/lung/lung_rs_info.txt")
setkeyv(lung_rs_info, c("CHR", "POS"))

##製作並輸出dosage file
for (i in seq(22)) {
  print(i)
  tped = fread(paste0("/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/dosage/lung_chr", i,".dosage.traw"), header = TRUE)
  freq = fread(paste0("/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/dosage/lung_chr", i,".frq.frq"), header = TRUE)
  freq[,"NCHROBS":=NULL]
  freq[,"CHR":=paste0("chr", CHR)]
  setkeyv(tped, c("CHR", "POS"))
  tped2 = merge.data.table(lung_rs_info, tped,sort = T)
  tped2 = tped2[,-c(1,4:11)]
  setkeyv(freq, c("SNP"))
  setkeyv(tped2, c("SNP"))
  freq2 = merge.data.table(freq, tped2,sort = T)
  if(identical(tped2[["SNP"]], freq2[["SNP"]])){
    output = data.table(freq2[,c("CHR","rsid","POS","A1", "A2", "MAF")], freq2[,-seq(1:10), with = FALSE])
    rm(tped, freq, tped2,freq2)
    fwrite(output, file = paste0("/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/dosage/final/chr", i, ".dosage.txt"),
           col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
  }else{
    stop("Wrong SNP order!")
  }
} 
i=1

### (3-2) Dosage file - for whole blood tissue
##這裡透過先前已有gnomAD的rsid表，其版本也為GRCh37做為對照表
rs_info <- fread("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/predixcan_gnomad_rsid_info.csv")
comb <- fread("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/weight.csv")
which(colnames(comb)=="Whole_Blood")
comb <- comb[,c(1:5,54)]
comb <- comb[complete.cases(comb),]
uni_comb <- comb[!duplicated(comb$varID),c(2,4,5)]
blood_snp <- uni_comb$rsid
blood_rs_info <- subset(rs_info, rsid %in% blood_snp)
colnames(blood_rs_info)[4] <- "CHR"
#fwrite(blood_rs_info, file = paste0("/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/blood_tissue/blood_rs_info.txt"))
blood_rs_info <- fread("/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/blood_tissue/blood_rs_info.txt")
setkeyv(blood_rs_info, c("CHR", "POS"))

##製作並輸出dosage file
for (i in seq(22)) {
  print(i)
  tped = fread(paste0("/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/dosage/lung_chr", i,".dosage.traw"), header = TRUE)
  freq = fread(paste0("/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/dosage/lung_chr", i,".frq.frq"), header = TRUE)
  freq[,"NCHROBS":=NULL]
  freq[,"CHR":=paste0("chr", CHR)]
  setkeyv(tped, c("CHR", "POS"))
  tped2 = merge.data.table(blood_rs_info, tped,sort = T)
  tped2 = tped2[,-c(1,4:11)]
  setkeyv(freq, c("SNP"))
  setkeyv(tped2, c("SNP"))
  freq2 = merge.data.table(freq, tped2,sort = T)
  if(identical(tped2[["SNP"]], freq2[["SNP"]])){
    output = data.table(freq2[,c("CHR","rsid","POS","A1", "A2", "MAF")], freq2[,-seq(1:10), with = FALSE])
    rm(tped, freq, tped2,freq2)
    fwrite(output, file = paste0("/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/blood_tissue/phenotype_file/chr", i, ".dosage.txt"),
           col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
  }else{
    stop("Wrong SNP order!")
  }
} 
#i=1


### (4) Compress file 
##for lung
dosage_final_path=/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/dosage/final
#linux
for i in {1..22}
do 
gzip $dosage_final_path/chr${i}.dosage.txt
done

mv $dosage_final_path/chr*.dosage.txt.gz /Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/phenotype_file

##for whole blood
dosage_final_path=/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/blood_tissue/phenotype_file
#linux
for i in {1..22}
do 
gzip $dosage_final_path/chr${i}.dosage.txt
done


############################
## PrediXcan prediction   ##
############################
###這邊的預測用到的樣本是156 samples一起預測，後續再分開兩個GSE datasets
################ *↓Lung tissue↓* ################
cd /Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關
### (5-1) Construct environment
#conda create --name mypy27 python=2.7
source activate py27
#conda install numpy
source deactivate myenv

### (5-2) Prediction
./PrediXcan.py --predict \
--weights en_Lung.db \
--dosages phenotype_file \--samples sample_file.txt \--pheno phenotype_file.txt \
--output_prefix results/

##結果：134227/194164(69.13%) snps found
################ *↑Lung tissue↑* ################

################ *↓Whole Blood tissue↓* ################
cd /Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/blood_tissue
### (6-1) Construct environment
source activate py27
source deactivate myenv

### (6-2) Prediction
./PrediXcan.py --predict \
--weights en_Whole_Blood.db \
--dosages phenotype_file \--samples sample_file.txt \--pheno phenotype_file.txt \
--output_prefix results/

##結果：120755/177016(68.22%) snps found
################ *↑Whole Blood tissue↑* ################


########################################
## Result v.s. Reference comparsion   ##
########################################

################ *↓Lung tissue↓* ################
#R
library(data.table)
library(dplyr)
library(broom)
### (7-1) 首先計算data實際在每個gene用到多少個SNP預測
#取出Lung tissue
comb <- fread("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/weight.csv")
lung_rs_info <- fread("/Users/chgsh14414/Desktop/Mac/Predixcan/lung/lung_rs_info.txt")
setkeyv(lung_rs_info, c("CHR", "POS"))
which(colnames(comb)=="Lung")
comb <- comb[,c(1:5,38)]
comb <- comb[complete.cases(comb),]

##計算每個用到的SNP個數
n_used <- comb %>% group_by(gene) %>% do(tidy(nrow(.)))
lung_match_snp <- data.frame()
for (i in seq(22)) {
  print(i)
  tped = fread(paste0("/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/dosage/lung_chr", i,".dosage.traw"), header = TRUE)
  #freq = fread(paste0("/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/dosage/lung_chr", i,".frq.frq"), header = TRUE)
  #freq[,"NCHROBS":=NULL]
  #freq[,"CHR":=paste0("chr", CHR)]
  setkeyv(tped, c("CHR", "POS"))
  tped2 = merge.data.table(lung_rs_info, tped,sort = T)
  tped2 = tped2[,-c(1,4:11)]
  lung_match_snp <- rbind(lung_match_snp, tped2[,2])
} 
lung_match_snp <- lung_match_snp[!duplicated(lung_match_snp$rsid),]
fwrite(lung_match_snp, "/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/variant個數計算/lung_match_snp.txt")
data_list <- lung_match_snp$rsid
data_num_snp <- comb %>% group_by(gene) %>% do(tidy(sum(ifelse(.$rsid %in% data_list,1,0))))

##合併資料跟ref的
names(n_used)[2] <-"ref"
names(data_num_snp)[2] <-"data"
setkeyv(n_used, "gene")
setkeyv(data_num_snp, "gene")
lung_snp_merge <- merge(data_num_snp, n_used, by = "gene")
lung_snp_merge$prop <- lung_snp_merge$data/lung_snp_merge$ref
summary(lung_snp_merge$prop)
fwrite(lung_snp_merge, "/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/variant個數計算/lung_snp_merge.txt")

### (7-2) 整理計算
lung <- fread("/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/results/_predicted_expression.txt")
rank <- fread("/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/results/rank_result.csv")
lung_snp_merge <- fread("/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/variant個數計算/lung_snp_merge.txt")

##排除GSE26853 gastric cancer的人
rank <- subset(rank, substr(rank$IID,1,5)=="GSM82")
lung <- subset(lung, substr(lung$IID,1,5)=="GSM82")

ref_mean <- rank[,c(3,7,8)]
ref_mean <- ref_mean[!duplicated(ref_mean$gene),]
colnames(ref_mean)[2] <- "ref"
#ref_mean$lower_bound <- ref_mean$ref-2*ref_mean$sd
#ref_mean$upper_bound <- ref_mean$ref+2*ref_mean$sd
#ref_mean <- data.table(ref_mean[,c(1,2,4,5)])

##計算data mean &sd
ref_mean$data_mean <- apply(lung[,-c(1:2)], 2, mean)
ref_mean$data_sd <- apply(lung[,-c(1:2)], 2, sd)
ref_mean$abs_diff <- abs(ref_mean$ref-ref_mean$data_mean)
ref_mean$prop_diff <- (ref_mean$data_mean-ref_mean$ref)/ref_mean$ref
summary(ref_mean$abs_diff)
summary(ref_mean$prop_diff)
summary(ref_mean$data_sd)
summary(ref_mean$pop_diff[-which(ref_mean$pop_diff==Inf)])

##
names(ref_mean)[2] <- "ref.mean"
#lung_snp_merge <- data.table(lung_snp_merge)
setkeyv(ref_mean2, "gene")
setkeyv(lung_snp_merge, "gene")
lung_eval = merge.data.table(ref_mean, lung_snp_merge,sort = T)

##計算每個gene有多少人符合reference的範圍
rank_eval= rank %>% group_by(gene) %>% do(tidy(sum(ifelse(.$percentile.rank %in% c(0,100),0,1))))
rank_eval$index <- rank_eval$x/61
rank_eval <- data.table(rank_eval)
setkeyv(rank_eval, c("gene"))
setkeyv(lung_eval, c("gene"))
lung_rank_eval = merge.data.table(lung_eval, rank_eval,sort = T)
colnames(lung_rank_eval) <- c("Gene","Ref.Mean","Ref.Sd","Data.Mean","Data.Sd","Abs.Difference","Ratio.Differece","Data.SNP.Num","PrediXcan.SNP.Num","SNP.Ratio","GE.Match.Ref.Num","GE.Match.Ref.Ratio")
fwrite(lung_rank_eval, "/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/results/lung_rank_eval.csv")
lung_rank_eval <- fread("/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/results/lung_rank_eval.csv")

################ *↓Whole Blood tissue↓* ################
#R
library(data.table)
library(dplyr)
library(broom)
### (7-1) 首先計算data實際在每個gene用到多少個SNP預測
#取出Lung tissue
comb <- fread("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/weight.csv")
blood_rs_info <- fread("/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/blood_tissue/blood_rs_info.txt")
setkeyv(blood_rs_info, c("CHR", "POS"))
which(colnames(comb)=="Whole_Blood")
comb <- comb[,c(1:5,54)]
comb <- comb[complete.cases(comb),]

##計算每個用到的SNP個數
n_used <- comb %>% group_by(gene) %>% do(tidy(nrow(.)))
blood_match_snp <- data.frame()
for (i in seq(22)) {
  print(i)
  tped = fread(paste0("/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/dosage/lung_chr", i,".dosage.traw"), header = TRUE)
  #freq = fread(paste0("/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/dosage/lung_chr", i,".frq.frq"), header = TRUE)
  #freq[,"NCHROBS":=NULL]
  #freq[,"CHR":=paste0("chr", CHR)]
  setkeyv(tped, c("CHR", "POS"))
  tped2 = merge.data.table(blood_rs_info, tped,sort = T)
  tped2 = tped2[,-c(1,4:11)]
  blood_match_snp <- rbind(blood_match_snp, tped2[,2])
} 
blood_match_snp <- blood_match_snp[!duplicated(blood_match_snp$rsid),]
fwrite(blood_match_snp, "/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/variant個數計算/blood_match_snp.txt")
data_list <- blood_match_snp$rsid
data_num_snp <- comb %>% group_by(gene) %>% do(tidy(sum(ifelse(.$rsid %in% data_list,1,0))))

##合併資料跟ref的
data_num_snp <- data.table(data_num_snp)
n_used <- data.table(n_used)
names(n_used)[2] <-"ref"
names(data_num_snp)[2] <-"data"
setkeyv(n_used, "gene")
setkeyv(data_num_snp, "gene")
blood_snp_merge <- merge(data_num_snp, n_used, by = "gene")
blood_snp_merge$prop <- blood_snp_merge$data/blood_snp_merge$ref
summary(blood_snp_merge$prop)
fwrite(blood_snp_merge, "/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/variant個數計算/blood_snp_merge.txt")

### (7-2) 整理計算
library(PredictAP)
blood <- fread("/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/blood_tissue/results/_predicted_expression.txt")
rank <- rankcal(blood, "x49")
#fwrite(rank, "/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/blood_tissue/results/rank_blood_result.csv")
rank <- fread("/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/blood_tissue/results/rank_blood_result.csv")
blood_snp_merge <- fread("/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/variant個數計算/blood_snp_merge.txt")

##排除GSE33356 Lung cancer的人
rank <- subset(rank, substr(rank$IID,1,5)=="GSM66")
blood <- subset(blood, substr(blood$IID,1,5)=="GSM66")

ref_mean <- rank[,c(3,7,8)]
ref_mean <- ref_mean[!duplicated(ref_mean$gene),]
colnames(ref_mean)[2] <- "ref"
#ref_mean$lower_bound <- ref_mean$ref-2*ref_mean$sd
#ref_mean$upper_bound <- ref_mean$ref+2*ref_mean$sd
#ref_mean <- data.table(ref_mean[,c(1,2,4,5)])

##計算data mean &sd
ref_mean$data_mean <- apply(blood[,-c(1:2)], 2, mean)
ref_mean$data_sd <- apply(blood[,-c(1:2)], 2, sd)
ref_mean$abs_diff <- abs(ref_mean$ref-ref_mean$data_mean)
ref_mean$prop_diff <- (ref_mean$data_mean-ref_mean$ref)/ref_mean$ref
summary(ref_mean$abs_diff)
summary(ref_mean$prop_diff)
summary(ref_mean$data_sd)
summary(ref_mean$pop_diff[-which(ref_mean$pop_diff==Inf)])

##
names(ref_mean)[2] <- "ref.mean"
#blood_snp_merge <- data.table(blood_snp_merge)
setkeyv(ref_mean, "gene")
setkeyv(blood_snp_merge, "gene")
blood_eval = merge.data.table(ref_mean, blood_snp_merge,sort = T)

##計算每個gene有多少人符合reference的範圍
rank_eval= rank %>% group_by(gene) %>% do(tidy(sum(ifelse(.$percentile.rank %in% c(0,100),0,1))))
rank_eval$index <- rank_eval$x/95
rank_eval <- data.table(rank_eval)
setkeyv(rank_eval, c("gene"))
setkeyv(blood_eval, c("gene"))
blood_rank_eval = merge.data.table(blood_eval, rank_eval,sort = T)
colnames(blood_rank_eval) <- c("Gene","Ref.Mean","Ref.Sd","Data.Mean","Data.Sd","Abs.Difference","Ratio.Differece","Data.SNP.Num","PrediXcan.SNP.Num","SNP.Ratio","GE.Match.Ref.Num","GE.Match.Ref.Ratio")
fwrite(blood_rank_eval, "/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/blood_tissue/results/blood_rank_eval.csv")
blood_rank_eval <- fread("/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/blood_tissue/results/blood_rank_eval.csv")

################ *↑blood tissue↑* ################