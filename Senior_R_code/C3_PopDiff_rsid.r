################ 先再次整合PrediXcan有用到ref/eff alleles位點 ################

setwd("~/Predixcan_materials/Predixcan/code")

library(data.table)
comb <- fread("../data/weight.csv")
data <- fread("../data/allchr_gnomad_info.csv")

uni_comb <- comb[!duplicated(comb$varID),c(2,4,5)] # 1,657,598

# *** part I ***
setkeyv(uni_comb, c("rsid", "ref_allele", "eff_allele"))
colnames(data)[4:5] <- c("ref_allele", "eff_allele")
setkeyv(data, c("rsid", "ref_allele", "eff_allele"))
data2 <- merge.data.table(uni_comb, data,sort = T)
nrow(uni_comb) - nrow(data2) # 7355

# *** part II ***
##翻轉eff/ref
colnames(data)[4:5] <- c("eff_allele", "ref_allele")
data3 <- merge.data.table(uni_comb, data,sort = T)
data3$AC_nfe <- data3$AN_nfe - data3$AC_nfe
data3$AC_afr <- data3$AN_afr - data3$AC_afr
data3$AC_eas <- data3$AN_eas - data3$AC_eas
#colnames(data3)[4:5] <- c("ref_allele", "eff_allele")
nrow(uni_comb) - nrow(data2) - nrow(data3) # 2011

# *** part III ***
##翻轉ATCG
colnames(data)[4:5] <- c("ref_allele", "eff_allele")
trun_data <- merge.data.table(uni_comb, data,sort = T, all.x = T)
trun_data <- trun_data[is.na(trun_data$AC_nfe),]
sum(trun_data$rsid %in% data3$rsid) # 5344
trun_data <- trun_data[!(trun_data$rsid %in% data3$rsid),1:3]
trun_data$eff_allele <- ifelse(trun_data$eff_allele == "A", "T",
                               ifelse(trun_data$eff_allele == "C", "G",
                               ifelse(trun_data$eff_allele == "T", "A", "C")))

trun_data$ref_allele <- ifelse(trun_data$ref_allele == "A", "T",
                               ifelse(trun_data$ref_allele == "C", "G",
                                      ifelse(trun_data$ref_allele == "T", "A", "C")))
setkeyv(trun_data, c("rsid", "ref_allele", "eff_allele"))
data4 <- merge.data.table(data, trun_data,sort = T)
#最後合併成結果時先轉回來
data4$eff_allele <- ifelse(data4$eff_allele == "A", "T",
                               ifelse(data4$eff_allele == "C", "G",
                                      ifelse(data4$eff_allele == "T", "A", "C")))

data4$ref_allele <- ifelse(data4$ref_allele == "A", "T",
                               ifelse(data4$ref_allele == "C", "G",
                                      ifelse(data4$ref_allele == "T", "A", "C")))
nrow(uni_comb) - nrow(data2) - nrow(data3) - nrow(data4) # 1101

# *** part IV ***
##翻轉ATCG & 翻轉eff/ref
colnames(trun_data)[2:3] <- c("eff_allele", "ref_allele")
data5 <- merge.data.table(trun_data, data, sort = T)
data5$AC_nfe <- data5$AN_nfe - data5$AC_nfe
data5$AC_afr <- data5$AN_afr - data5$AC_afr
data5$AC_eas <- data5$AN_eas - data5$AC_eas
data5$eff_allele <- ifelse(data5$eff_allele == "A", "T",
                           ifelse(data5$eff_allele == "C", "G",
                                  ifelse(data5$eff_allele == "T", "A", "C")))

data5$ref_allele <- ifelse(data5$ref_allele == "A", "T",
                           ifelse(data5$ref_allele == "C", "G",
                                  ifelse(data5$ref_allele == "T", "A", "C")))
colnames(data5)[2:3] <- c("ref_allele", "eff_allele")
nrow(uni_comb) - nrow(data2) - nrow(data3) - nrow(data4) - nrow(data5) # 946

##check
View(uni_comb[rsid == "rs10",])
View(uni_comb[rsid == "rs1000319",])
View(uni_comb[rsid == "rs1003007",])
View(uni_comb[rsid == "rs1008660",])
### View in console
# cat("rs10\n")
# print(uni_comb[rsid == "rs10",])
# cat("rs1000319\n")
# print(uni_comb[rsid == "rs1000319",])
# cat("rs1003007\n")
# print(uni_comb[rsid == "rs1003007",])
# cat("rs1008660\n")
# print(uni_comb[rsid == "rs1008660",])


# *** Final part ***
##merge all data
predixcan_rsid_info <- rbind(data2, data3, data4, data5)
nrow(uni_comb) - nrow(predixcan_rsid_info) ##剩餘946
fwrite(predixcan_rsid_info, "../data/predixcan_gnomad_rsid_info.csv")
##gnomAD 沒有的rsid list
aa <- uni_comb[!(rsid %in% predixcan_rsid_info$rsid),]
fwrite(aa, "../data/predixcan_gnomad_rsid_withoutinfo.csv")

################ *↓計算有族群差異顯著位點↓* ################
library(data.table)
pgno_rsinfo <- fread("../data/predixcan_gnomad_rsid_info.csv")
pgno_rsinfo$freq_nfe <- pgno_rsinfo$AC_nfe/pgno_rsinfo$AN_nfe
pgno_rsinfo$freq_afr <- pgno_rsinfo$AC_afr/pgno_rsinfo$AN_afr
pgno_rsinfo$freq_eas <- pgno_rsinfo$AC_eas/pgno_rsinfo$AN_eas

##刪掉分母為0的位點(代表沒有樣本) n=48
sum(is.na(pgno_rsinfo$freq_eas)) ##48
pgno_rsinfo <- pgno_rsinfo[!is.na(pgno_rsinfo$freq_eas),]

##nfe v.s. eas
p <- (pgno_rsinfo$AC_nfe+pgno_rsinfo$AC_eas)/(pgno_rsinfo$AN_nfe+pgno_rsinfo$AN_eas)
# unpooled Z-test / pooled proportion Z-test: choose one of way
## unpooled Z-test
Z_ne <- abs((pgno_rsinfo$freq_nfe - pgno_rsinfo$freq_eas)/sqrt((pgno_rsinfo$freq_nfe*(1-pgno_rsinfo$freq_nfe)/pgno_rsinfo$AN_nfe)+(pgno_rsinfo$freq_eas*(1-pgno_rsinfo$freq_eas)/pgno_rsinfo$AN_eas)))
## pooled proportion Z-test
# Z_ne <- abs((pgno_rsinfo$freq_nfe - pgno_rsinfo$freq_eas)/sqrt(p*(1-p)*(1/pgno_rsinfo$AN_nfe+1/pgno_rsinfo$AN_eas)))
pgno_rsinfo$pval_ne <- pnorm(Z_ne, lower.tail = F)*2
sum(is.na(pgno_rsinfo$pval_ne))
View(pgno_rsinfo[is.na(pgno_rsinfo$pval_ne),])

### View in console
# cat("Entries with NA p-values:\n")
# print(pgno_rsinfo[is.na(pgno_rsinfo$pval_ne),])

##再多設差距要0.05以上
pgno_rsinfo$diff_index <- ifelse(pgno_rsinfo$pval_ne < 10^-8 & abs(pgno_rsinfo$freq_nfe - pgno_rsinfo$freq_eas) > 0.05, "1", "0")
pgno_rsinfo$diff_index[which(is.na(pgno_rsinfo$pval_ne))] <- "0"
table(pgno_rsinfo$diff_index) #0:498442 ; 1:1158162
prop.table(table(pgno_rsinfo$diff_index))

################ *↑計算有族群差異顯著位點↑* ################

# fix
setcolorder(
  pgno_rsinfo,
  c(
    "rsid", "ref_allele", "eff_allele",
    "chr", "POS",
    "AC_nfe", "AN_nfe",
    "AC_afr", "AN_afr",
    "AC_eas", "AN_eas",
    "freq_nfe", "freq_afr", "freq_eas",
    "pval_ne", "diff_index"
  )
)

fwrite(
  pgno_rsinfo,
  "../data/predixcan_gnomad_rsid_info_with_popdiff.csv"
)
