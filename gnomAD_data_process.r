awk -v FS=" " '{print $3}' /Users/chgsh14414/Downloads/test.txt > /Users/chgsh14414/Downloads/rs_test.txt
awk -v FS=";" '{print $79}' /Users/chgsh14414/Downloads/test_genomes_chr21.txt
awk 'BEGIN {FS=";"}''{/AC_afr=/; print $FNR}' /Users/chgsh14414/Downloads/test_genomes_chr21.txt
awk -v FS=";" '{/AC_afr=/; print $NF}' /Users/chgsh14414/Downloads/test_genomes_chr21.txt
echo 'AC_afr=0' | awk '{print substr($0,7,4)}'

##snp
#check
awk -v FS=" " '{print $3}' /Users/chgsh14414/Downloads/test_genomes_chr21.txt
awk -v FS=" " 'FNR==NR {a[$1]; next}; $3 in a' /Users/chgsh14414/rs_list_0221.txt /Users/chgsh14414/Downloads/test_genomes_chr21.txt > /Users/chgsh14414/Downloads/extract_test.txt
#run
sed -i '1,598d' /home/chgsh14414/workspace/gnomad_data/gnomad.genomes.r2.1.1.sites.21.vcf
awk -v FS=" " 'FNR==NR {a[$1]; next}; $3 in a' rs_list_0221.txt gnomad.genomes.r2.1.1.sites.21.vcf > extr_gnomad.genomes.r2.1.1.sites.21.vcf

##換行調整
sed 's/;/\n/g' extr_gnomad.genomes.r2.1.1.sites.21.vcf > extrv2_gnomad.genomes.r2.1.1.sites.21.vcf





gunzip -c /Users/chgsh14414/Downloads/gnomad.genomes.r2.1.1.sites.21.vcf.bgz > /Users/chgsh14414/Downloads/gnomad.genomes.r2.1.1.sites.21.vcf

awk 'NR>=1 && NR <=500' /Users/chgsh14414/Downloads/gnomad.genomes.r2.1.1.sites.21.vcf > /Users/chgsh14414/Downloads/name_genomes_chr21.txt
awk 'NR>=500 && NR <=1500' /Users/chgsh14414/Downloads/gnomad.genomes.r2.1.1.sites.21.vcf > /Users/chgsh14414/Downloads/test_genomes_chr21.txt

#################################
##先處理vcf把要的variant撈出來
./plink --vcf /Users/chgsh14414/Downloads/gnomad.genomes.r2.1.1.sites.21.vcf.bgz --allow-no-samples --keep /Users/chgsh14414/rs_list_0221.txt --make-bed --out /Users/chgsh14414/try.vcf.gz
./plink --vcf /Users/chgsh14414/Downloads/gnomad.genomes.r2.1.1.sites.21.vcf.gz --allow-no-samples --make-bed --out /Users/chgsh14414/try00
./bcftools sort /Users/chgsh14414/Downloads/gnomad.genomes.r2.1.1.sites.21.vcf -Oz -o /Users/chgsh14414/Downloads/gnomad.genomes.r2.1.1.sites.21.vcf.gz

./plink --vcf /Users/chgsh14414/Downloads/gnomad.exomes.r2.1.1.sites.21.vcf.bgz --allow-no-samples --extract /Users/chgsh14414/rs_list_0221.txt --make-bed --out /Users/chgsh14414/try.vcf.gz
./plink --vcf /Users/chgsh14414/Downloads/gnomad.exomes.r2.1.1.sites.21.vcf.bgz --extract /Users/chgsh14414/rs_list_0221.txt --make-bed --out /Users/chgsh14414/try.vcf.gz
##################################

##########server
gunzip -c gnomad.genomes.r2.1.1.sites.21.vcf.bgz > gnomad.genomes.r2.1.1.sites.21.vcf
gunzip -c gnomad.genomes.r2.1.1.sites.vcf.bgz > gnomad.genomes.r2.1.1.sites.vcf

awk -v FS=";" '{NR<50; print $79}' gnomad.genomes.r2.1.1.sites.21.vcf

sed 's/;/\n/g' /Users/chgsh14414/Downloads/test_genomes_chr21.txt > /Users/chgsh14414/Downloads/test_genomes_chr21_v2.txt

sed -i '1,598d' /home/chgsh14414/workspace/gnomad_data/gnomad.genomes.r2.1.1.sites.21.vcf
sed -i '1,2d' /home/chgsh14414/workspace/gnomad_data/gnomad.genomes.r2.1.1.sites.21.vcf

##total chrosome
sed -i '1,597d' gnomad.genomes.r2.1.1.sites.vcf
sed -i '1,2d' gnomad.genomes.r2.1.1.sites.vcf
awk -v FS=" " 'FNR==NR {a[$1]; next}; $3 in a' rs_list_0221.txt gnomad.genomes.r2.1.1.sites.vcf > extr_gnomad.genomes.r2.1.1.sites.vcf
sed 's/;/\n/g' extr_gnomad.genomes.r2.1.1.sites.vcf > extr2_gnomad.genomes.r2.1.1.sites.vcf
sed 's/\t/\n/g' extr2_gnomad.genomes.r2.1.1.sites.vcf > extr3_gnomad.genomes.r2.1.1.sites.vcf


##extracr rsid
awk -v FS=" " '{print $3}' extr_gnomad.genomes.r2.1.1.sites.vcf > rs_id.txt
awk -v FS=" " '{print $1,$2,$3,$4,$5}' extr_gnomad.genomes.r2.1.1.sites.vcf > rs_id.txt



##R part
library(data.table)
close(myCon)
rm(myCon)
myCon = file(description = "/home/chgsh14414/workspace/gnomad_data/extr3_gnomad.genomes.r2.1.1.sites.vcf", open="r", blocking = TRUE)
myCon = file(description = "/home/chgsh14414/workspace/gnomad_data/check_file3.txt", open="r", blocking = TRUE)

#########################不太需要用到#########################
count <- 0
repeat{
  pl = readLines(myCon, n = 1) # Read one line from the connection.
  if(identical(pl, character(0))){break} # If the line is empty, exit.
  if(identical(substr(pl,1,2), "21")){
    count <- count + 1
  } # Otherwise, print and repeat next iteration.
}
#########################不太需要用到#########################

##找尋特定rsid所在的行數
grep -n "rs3748594" extr2_gnomad.genomes.r2.1.1.sites.vcf | cut -d: -f1

count <- 269228 ##total snp number
count <- 269228 - 240000  ##殘餘
count <- 80000 ##一次負荷量
nskip <- 44655449 ##第二次起始位置
, skip = nskip

#close(myCon)
#rm(myCon)
#myCon = file(description = "/Users/chgsh14414/Downloads/test_genomes_chr21_v2.txt", open="r", blocking = TRUE)

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
print(paste0("執行時間", t[3][[1]]/60, "分"))
write(tt, "/home/chgsh14414/workspace/gnomad_data/allchr_result_sub1.csv", row.names = F)
fwrite(tt, "/home/chgsh14414/workspace/gnomad_data/allchr_result_sub1.csv", row.names = F)
fwrite(tt, "/home/chgsh14414/workspace/gnomad_data/allchr_result_sub2.csv", row.names = F)
fwrite(tt, "/home/chgsh14414/workspace/gnomad_data/allchr_result_sub3.csv", row.names = F)
fwrite(tt, "/home/chgsh14414/workspace/gnomad_data/allchr_result_sub4.csv", row.names = F)



write.table(tt, "/home/chgsh14414/workspace/gnomad_data/allchr_result.txt")


##merge snp id & frequency
library(data.table)
##frequency
sub1 <- fread("/home/chgsh14414/workspace/gnomad_data/allchr_result_sub1.csv", sep = ",")
sub2 <- fread("/home/chgsh14414/workspace/gnomad_data/allchr_result_sub2.csv", sep = ",")
sub3 <- fread("/home/chgsh14414/workspace/gnomad_data/allchr_result_sub3.csv", sep = ",")
sub4 <- fread("/home/chgsh14414/workspace/gnomad_data/allchr_result_sub4.csv", sep = ",")
mm <- rbind(sub1, sub2, sub3, sub4)
fwrite(mm, "/home/chgsh14414/workspace/gnomad_data/allchr_result.csv", row.names = F)

##rsid
rsid <- fread("/home/chgsh14414/workspace/gnomad_data/rs_id.txt", header = F)
colnames(rsid) <- c("chr","POS","rsid","A1","A2")
f_merge <- cbind(rsid, mm)
fwrite(f_merge, "/home/chgsh14414/workspace/gnomad_data/allchr_result_rsid.csv", row.names = F)


#########################不太需要用到 確認資料#########################
#linux
head -n 1 extr_gnomad.genomes.r2.1.1.sites.vcf > check_file.txt
sed 's/;/\n/g' check_file.txt > check_file2.txt
sed 's/\t/\n/g' check_file2.txt > check_file3.txt
#########################不太需要用到 確認資料#########################


