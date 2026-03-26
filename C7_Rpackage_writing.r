##R package連結設定
git config user.name "chgsh14414"
git config --global user.email "chgsh14414@gmail.com"


#############################################
## Step1: 整理Package所需要的檔案           ##
#############################################
##先整理出各個tissue使用到rsid列表
#R
#test code
setwd("~/Predixcan_materials/Predixcan/code")
library("RSQLite")
sqlite <- dbDriver("SQLite")
# dbname <- "/Users/chgsh14414/Desktop/Mac/Predixcan/elastic_net_models/en_Whole_Blood.db"
dbname <- "../../elastic_net_models/en_Whole_Blood.db"

db = dbConnect(sqlite,dbname)
## list tables
dbListTables(db)
dbListFields(db, "weights")
query <- function(...) dbGetQuery(db, ...)
query('select * from weights limit 50')
query('select * from MAF limit 50')
weight_database <- query('select * from weights')
weight_database <- data.table(weight_database)
query('select count(*) from extra') #7572 genes

##for all 49 tissues
# setwd("D:\\Python_Program\\PredictAP_web\\Predixcan_materials\\elastic_net_models")
setwd("~/Predixcan_materials/Predixcan/elastic_net_models")

ped_temp <- list.files(pattern=".db")
# ped_temp2 <- paste0("D:\\Python_Program\\PredictAP_web\\Predixcan_materials\\elastic_net_models\\", ped_temp)
ped_temp2 <- paste0("~/Predixcan_materials/Predixcan/elastic_net_models", ped_temp)

for (i in 1:49) {
  assign(paste0("db",i) , dbConnect(sqlite,ped_temp2[i]))
}

pedfiles <- paste0("db", seq(1:49))
for (i in 1:length(pedfiles)) {
  query <- function(...) dbGetQuery(eval(parse(text = pedfiles[i])), ...)
  assign(paste0("snp_sub",i) , data.table(query('select * from weights')))
}

db_name <- paste0("snp_sub", seq(1:49))
##how many snps were used for each gene?
library(tidyverse)
tis_name <- c()
for (k in 1:length(pedfiles)) {
  tis_name[k] <- substr(ped_temp[k], regexpr("_", ped_temp[k])+1,regexpr(".db", ped_temp[k])-1)
}

g1 <- rename(data.table(table(snp_sub1$gene)), Adipose_Subcutaneous=N)
g2 <- rename(data.table(table(snp_sub2$gene)), Adipose_Visceral_Omentum=N)
setkeyv(g1, "V1")
setkeyv(g2, "V1")
snp_gene_used <- merge.data.table(g1, g2, all = T)
for (i in 3:49) {
  g <- data.table(table(eval(parse(text = db_name[i]))$gene))
  print(nrow(g)) ##check
  colnames(g)[2] <- tis_name[i]
  setkeyv(g, "V1")
  snp_gene_used <- merge.data.table(snp_gene_used, g, all = T)
}
snp_gene_used <- rename(snp_gene_used, gene=V1)
# fwrite(snp_gene_used, "/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/Num_of_rs_used.csv")
fwrite(snp_gene_used, "../data/Num_of_rs_used.csv")
colnames(snp_gene_used)[-1] <- paste0("x",seq(1:49))
# fwrite(snp_gene_used, "/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/Num_of_rs_used_x.csv")
fwrite(snp_gene_used, "../data/Num_of_rs_used_x.csv")
##製作tissue索引
tissue_table <- data.table(index=seq(1:49), tissue=tis_name)
# fwrite(tissue_table, "/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/tissue_index.txt")
fwrite(tissue_table, "../data/tissue_index.txt")



##計算每個tissue表的基因的mean及sd
library(data.table)
# setwd("D:\\Python_Program\\PredictAP_web\\Predixcan_materials\\統整\\PrediXcan統整完檔案\\dosage\\nfe\\results")
setwd("D:\\Python_Program\\PredictAP_web\\Predixcan_materials\\統整\\PrediXcan統整完檔案\\dosage\\nfe\\results")

ped_temp <- list.files(pattern=".txt")
# ped_temp2 <- paste0("D:\\Python_Program\\PredictAP_web\\Predixcan_materials\\統整\\PrediXcan統整完檔案\\dosage\\nfe\\results\\", ped_temp)
ped_temp2 <- paste0("D:\\Python_Program\\PredictAP_web\\Predixcan_materials\\統整\\PrediXcan統整完檔案\\dosage\\nfe\\results\\", ped_temp)

for (i in 1:49) {
  assign(substr(ped_temp[i], 1,regexpr("_pred", ped_temp[i])-1), fread(ped_temp2[i]))
}
tis_name <- substr(ped_temp, 1,regexpr("_pred", ped_temp)-1)
for (i in 1:49) {
  aa <- matrix(nrow = ncol(eval(parse(text = tis_name[i])))-2, ncol = 3, byrow = F)
  aa[,1] <- colnames(eval(parse(text = tis_name[i])))[-c(1:2)]
  aa[,2] <- as.numeric(apply(eval(parse(text = tis_name[i]))[,-c(1:2)], 2, mean))
  aa[,3] <- as.numeric(apply(eval(parse(text = tis_name[i]))[,-c(1:2)], 2, sd))
  colnames(aa) <- c("gene", paste0("x",i,".mean"), paste0("x",i,".sd"))
  assign(paste0("x",i,".ms"), data.table(aa))
}
# rs_used <- fread("D:\\Python_Program\\PredictAP_web\\Predixcan_materials\\統整\\PrediXcan統整完檔案\\Num_of_rs_used_x.csv")
rs_used <- fread("D:\\Python_Program\\PredictAP_web\\Predixcan_materials\\統整\\PrediXcan統整完檔案\\Num_of_rs_used_x.csv")
setkeyv(rs_used, "gene")
temp_nam <- paste0("x",1:49,".ms")
for (j in 1:49) {
  dt <- setkeyv(eval(parse(text = temp_nam[j])), c("gene"))
  rs_used <- merge.data.table(rs_used, dt, sort = T, all = T)
}

# fwrite(rs_used, "D:\\Python_Program\\PredictAP_web\\Predixcan_materials\\統整\\PrediXcan統整完檔案\\dosage\\nfe\\result_table.csv")
fwrite(rs_used, "D:\\Python_Program\\PredictAP_web\\Predixcan_materials\\統整\\PrediXcan統整完檔案\\dosage\\nfe\\result_table.csv")


##
rankcal <- function(data, tissue, nsubject){
  rankm <- matrix(NA, nsubject, ncol(tissue), byrow = F)
  for (i in 3:ncol(tissue)) {
    rankm[,i] <- ecdf_fun(tissue[[i]], data[[i]])
  }
  rankm[,1:2] <- data[,c(FID,IID)]
  rankm <- data.table(rankm)
  colnames(rankm) <- colnames(data)
  #return(rankm)
  l_data <- pivot_longer(data, cols = starts_with("ENSG"), names_to = "gene", values_to = "expression")
  l_rank <- pivot_longer(rankm, cols = starts_with("ENSG"), names_to = "gene", values_to = "rank")
  l_mer <- merge(l_data, l_rank, by = c("FID", "IID","gene"))
  res_por <- result_table %>% select(gene, starts_with("x1")) %>% merge(tt, by = "gene")
  result <- res_por %>% select(FID, IID, gene, x1, expression, rank, x1.mean, x1.sd)%>%
    rename(N.snps=x1, mean=x1.mean, sd=x1.sd) %>% arrange(FID, IID)
  return(result)
}

##前置git設定
cd /Users/chgsh14414/Desktop/Mac/Predixcan/predictEAS
git config user.name "Han-Ching Chan"
git config user.email "chgsh14414@gmail.com"

git config --global user.name "Han-Ching Chan"
git config --global user.email "chgsh14414@gmail.com"
git config --global --list

git remote -v
git branch -vv
git remote show origin
git config pull.rebase false


hint:   git config pull.rebase false  # merge
hint:   git config pull.rebase true   # rebase
hint:   git config pull.ff only       # fast-forward only
hint: You can replace "git config" with "git config --global" to set a default
hint: preference for all repositories. You can also pass --rebase, --no-rebase,
hint: or --ff-only on the command line to override the configured default per
hint: invocation.
fatal: Need to specify how to reconcile divergent branches.

##R writing
library(devtools)
library(tidyverse)
use_git()

use_r("ecdf_fn")
ecdf_fn <- function(x,perc) ecdf(x)(perc)
use_r("rankcal")
rankcal <- function(data, tissue, nsubject, rank.tis){
  rankm <- matrix(NA, nsubject, ncol(tissue), byrow = F)
  for (i in 3:ncol(tissue)) {
    rankm[,i] <- ecdf_fn(tissue[[i]], data[[i]])
  }
  rankm[,1:2] <- data[,c(FID,IID)]
  rankm <- data.table(rankm)
  colnames(rankm) <- colnames(data)
  #return(rankm)
  l_data <- pivot_longer(data, cols = starts_with("ENSG"), names_to = "gene", values_to = "expression")
  l_rank <- pivot_longer(rankm, cols = starts_with("ENSG"), names_to = "gene", values_to = "rank")
  l_mer <- merge(l_data, l_rank, by = c("FID", "IID","gene"))
  res_por <- result_table %>% select(gene, {{rank.tis}}, starts_with(paste0(rank.tis,"."))) %>% 
    filter(!is.na(eval(parse(text = rank.tis)))) %>% merge(l_mer, by = "gene")
  result <- res_por %>% select(FID, IID, gene, expression, rank, starts_with(rank.tis))%>% arrange(FID, IID)
  colnames(result)[c(4,7,8)] <- c("N.snps", "mean", "sd")
  return(result)
}
 
check <- aa(ss, tissue = x2, nsubject = 10, rank.tis = "x2")
# ss <- fread("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/example_data2.csv")
ss <- fread("../data/example_data2.csv")


library(data.table)
# result_table <- fread("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/result_table.csv")
# setwd("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/dosage/results")
result_table <- fread("../data/result_table.csv")
setwd("../data/dosage/results")
ped_temp <- list.files(pattern=".txt")
# ped_temp2 <- paste0("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/dosage/results/", ped_temp)
ped_temp2 <- paste0("../data/dosage/results/", ped_temp)
for (i in 1:49) {
  assign(paste0("x",i), fread(ped_temp2[i]))
}
use_data(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,result_table,
         x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,
         x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,
         x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,
         x41,x42,x43,x44,x45,x46,x47,x48,x49,  compress = T, overwrite = T)
