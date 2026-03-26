#############################################
## Extract the unique rsid from weight.csv ##
#############################################
##R
library(data.table)
comb <- fread("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/weight.csv")
comb <- comb[,c(1:5)]
uni_comb <- comb[!duplicated(comb$varID),2] #1,657,598
fwrite(uni_comb, "/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/gnomAD_use_rslist.txt",
       col.names = F, row.names = F, quote = F)


###################################################
## Variants data downloaded from gnomAD database ##
###################################################
##server
gunzip -c gnomad.genomes.r2.1.1.sites.vcf.bgz > gnomad.genomes.r2.1.1.sites.vcf
sed -i '1,597d' gnomad.genomes.r2.1.1.sites.vcf
sed -i '1,2d' gnomad.genomes.r2.1.1.sites.vcf

##整理vcf檔，將所需的rsid撈出來
awk -v FS=" " 'FNR==NR {a[$1]; next}; $3 in a' gnomAD_use_rslist.txt gnomad.genomes.r2.1.1.sites.vcf > extr_gnomad.genomes.r2.1.1.sites.vcf
#extr_gnomad.genomes.r2.1.1.sites.vcf: 1,720,048 (比rs_list多有可能是因爲同樣rsid eff_allele不一樣)
sed 's/;/\n/g' extr_gnomad.genomes.r2.1.1.sites.vcf > extr2_gnomad.genomes.r2.1.1.sites.vcf
sed 's/\t/\n/g' extr2_gnomad.genomes.r2.1.1.sites.vcf > extr3_gnomad.genomes.r2.1.1.sites.vcf
#檢查rsid數量: 1,720,048 SNPs
awk -v FS=" " '{print $1,$2,$3,$4,$5}' extr_gnomad.genomes.r2.1.1.sites.vcf > rs_id.txt

##額外計算原始database variants個數
awk -v FS=" " '{print $3}' gnomad.genomes.r2.1.1.sites.vcf > all_rs_id.txt


##R
library(data.table)
myCon = file(description = "/home/chgsh14414/workspace/gnomad_data/extr3_gnomad.genomes.r2.1.1.sites.vcf", open="r", blocking = TRUE)

#############
close(myCon)
rm(myCon)
#############

##找尋特定rsid所在的行數
grep -n "rs1077232" extr3_gnomad.genomes.r2.1.1.sites.vcf | cut -d: -f1

count <- 1720048 ##total snp number

##############不一定要用到############
count <- 1720048 - 200000  ##殘餘
count <- 400000 ##一次負荷量
nskip <- 226068988 ##第二次起始位置
, skip = nskip
##############不一定要用到############

rs_id <- c()
AC_afr <- c()
AN_afr <- c()
AC_eas <- c()
AN_eas <- c()
AC_nfe <- c()
AN_nfe <- c()
t1 <- proc.time()
for (i in 1:count) {
  AC_nef_temp <- 0
  AN_nef_temp <- 0
  repeat{
    pl = readLines(myCon, n = 1) # Read one line from the connection.
    if(identical(substr(pl,1,2), "rs")){
      #print(pl)
      rs_id[i] <- pl
    }
    if(identical(substr(pl,1,11), "AC_nfe_seu=")){
      #print(pl)
      AC_nef_temp <- AC_nef_temp + as.numeric(substr(pl,12,nchar(pl)))
    }
    if(identical(substr(pl,1,11), "AN_nfe_seu=")){
      #print(pl)
      AN_nef_temp <- AN_nef_temp + as.numeric(substr(pl,12,nchar(pl)))
    }
    if(identical(substr(pl,1,7), "AC_afr=")){
      #print(pl)
      AC_afr[i] <- as.numeric(substr(pl,8,nchar(pl)))
    }
    if(identical(substr(pl,1,7), "AN_afr=")){
      #print(pl)
      AN_afr[i] <- as.numeric(substr(pl,8,nchar(pl)))
    }
    if(identical(substr(pl,1,11), "AC_nfe_onf=")){
      #print(pl)
      AC_nef_temp <- AC_nef_temp + as.numeric(substr(pl,12,nchar(pl)))
    }
    if(identical(substr(pl,1,11), "AN_nfe_onf=")){
      #print(pl)
      AN_nef_temp <- AN_nef_temp + as.numeric(substr(pl,12,nchar(pl)))
    }
    if(identical(substr(pl,1,7), "AC_eas=")){
      #print(pl)
      AC_eas[i] <- as.numeric(substr(pl,8,nchar(pl)))
    }
    if(identical(substr(pl,1,7), "AN_eas=")){
      #print(pl)
      AN_eas[i] <- as.numeric(substr(pl,8,nchar(pl)))
    }
    if(identical(substr(pl,1,11), "AC_nfe_nwe=")){
      #print(pl)
      AC_nef_temp <- AC_nef_temp + as.numeric(substr(pl,12,nchar(pl)))
    }
    if(identical(substr(pl,1,11), "AN_nfe_nwe=")){
      #print(pl)
      AN_nef_temp <- AN_nef_temp + as.numeric(substr(pl,12,nchar(pl)))
    }
    if(identical(substr(pl,1,11), "AC_nfe_est=")){
      #print(pl)
      AC_nef_temp <- AC_nef_temp + as.numeric(substr(pl,12,nchar(pl)))
    }
    if(identical(substr(pl,1,11), "AN_nfe_est=")){
      #print(pl)
      AN_nef_temp <- AN_nef_temp + as.numeric(substr(pl,12,nchar(pl)))
    }
    if(identical(substr(pl,1,3), "vep")){break} # Otherwise, print and repeat next iteration.
    AC_nfe[i] <- AC_nef_temp
    AN_nfe[i] <- AN_nef_temp
  }
}
tt <- data.frame(rs_id, AC_nfe, AN_nfe, AC_afr, AN_afr, AC_eas, AN_eas)
t2 <- proc.time()
t <- t2 - t1
print(paste0("執行時間", t[3][[1]]/60, "分")) #全跑大約花17.5個小時

fwrite(tt, "/home/chgsh14414/workspace/gnomad_data/result/allchr_result.csv", row.names = F)

##合併加上位點位置與ref/eff alleles資訊

rs_info <- fread("/home/chgsh14414/workspace/gnomad_data/rs_id.txt", header = F)
all_chr_info <- cbind(rs_info, tt)
#check rsid對應無誤
sum(all_chr_info$V3!=all_chr_info$rs_id)
#整理colname format
all_chr_info <- all_chr_info[,-6]
colnames(all_chr_info)[1:5] <- c("chr","POS","rsid","A1","A2") #Note:此處POS為b37版本
fwrite(all_chr_info, "/home/chgsh14414/workspace/gnomad_data/result/allchr_gnomad_info.csv")

##下載至本機
scp chgsh14414@140.112.117.141://home/chgsh14414/workspace/gnomad_data/result/allchr_gnomad_info.csv /Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/
