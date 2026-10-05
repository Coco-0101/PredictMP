#############################################
## Step1: 整理Package所需要的檔案           ##
#############################################
##先整理出各個tissue使用到rsid列表
#R
#test code
setwd("~/Predixcan_materials/Predixcan/code")
library("RSQLite")
library(data.table)
sqlite <- dbDriver("SQLite")
dbname <- "../../elastic_net_models/en_Whole_Blood.db"

db = dbConnect(sqlite,dbname)
## list tables
dbListTables(db)
dbListFields(db, "weights")
query <- function(...) dbGetQuery(db, ...)
query('select * from weights limit 50')
# query('select * from MAF limit 50') # no MAF table
weight_database <- query('select * from weights')
weight_database <- data.table(weight_database)
query('select count(*) from extra') #7572 genes

##for all 49 tissues
# setwd("D:\\Python_Program\\PredictAP_web\\Predixcan_materials\\elastic_net_models")
setwd("~/Predixcan_materials/elastic_net_models")

ped_temp <- list.files(pattern=".db")
# ped_temp2 <- paste0("D:\\Python_Program\\PredictAP_web\\Predixcan_materials\\elastic_net_models\\", ped_temp)
ped_temp2 <- paste0("~/Predixcan_materials/elastic_net_models/", ped_temp)

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
setwd("~/Predixcan_materials/Predixcan/code")
fwrite(snp_gene_used, "../data/Num_of_rs_used_betty.csv")
colnames(snp_gene_used)[-1] <- paste0("x",seq(1:49))
# fwrite(snp_gene_used, "/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/Num_of_rs_used_x.csv")
fwrite(snp_gene_used, "../data/Num_of_rs_used_x_betty.csv")
##製作tissue索引
tissue_table <- data.table(index=seq(1:49), tissue=tis_name)
# fwrite(tissue_table, "/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/tissue_index.txt")
fwrite(tissue_table, "../data/tissue_index_betty.txt")



##計算每個tissue表的基因的mean及sd
library(data.table)
# setwd("D:\\Python_Program\\PredictAP_web\\Predixcan_materials\\統整\\PrediXcan統整完檔案\\dosage\\nfe\\results")
setwd("~/Predixcan_materials/Predixcan/results/NFE")

ped_temp <- list.files(pattern=".txt")
# ped_temp2 <- paste0("D:\\Python_Program\\PredictAP_web\\Predixcan_materials\\統整\\PrediXcan統整完檔案\\dosage\\nfe\\results\\", ped_temp)
ped_temp2 <- paste0("~/Predixcan_materials/Predixcan/results/NFE/", ped_temp)

for (i in 1:49) {
  assign(substr(ped_temp[i], 1,regexpr("_pred", ped_temp[i])-1), fread(ped_temp2[i]))
}
# tis_name <- substr(ped_temp, 1,regexpr("_pred", ped_temp)-1)



# 1. 乾淨地提取 tissue 名稱（直接移除字尾）
tis_name <- gsub("_predicted_expression\\.txt$", "", ped_temp)

tis_name # 確認結果

# 2. 用 list 讀取所有資料，並直接命名
tissue_list <- list()
for (i in 1:49) {
  tissue_list[[tis_name[i]]] <- fread(ped_temp2[i])
}


# 3. 計算每個 tissue 表的基因 mean 及 sd
for (i in 1:49) {
  # 直接從 list 裡面取出對應的 data.table
  current_dt <- tissue_list[[tis_name[i]]]
  
  # 建立矩陣
  aa <- matrix(nrow = ncol(current_dt) - 2, ncol = 3, byrow = FALSE)
  
  # 填入基因名稱、Mean、SD
  aa[, 1] <- colnames(current_dt)[-c(1:2)]
  aa[, 2] <- as.numeric(apply(current_dt[, -c(1:2), with = FALSE], 2, mean))
  aa[, 3] <- as.numeric(apply(current_dt[, -c(1:2), with = FALSE], 2, sd))
  
  colnames(aa) <- c("gene", paste0("x", i, ".mean"), paste0("x", i, ".sd"))
  
  # 將結果存成 x1.ms, x2.ms ... (保留你原本後續程式需要的變數命名)
  assign(paste0("x", i, ".ms"), data.table(aa))
}

# for (i in 1:49) {
#   aa <- matrix(nrow = ncol(eval(parse(text = tis_name[i])))-2, ncol = 3, byrow = F)
#   aa[,1] <- colnames(eval(parse(text = tis_name[i])))[-c(1:2)]
#   aa[,2] <- as.numeric(apply(eval(parse(text = tis_name[i]))[,-c(1:2)], 2, mean))
#   aa[,3] <- as.numeric(apply(eval(parse(text = tis_name[i]))[,-c(1:2)], 2, sd))
#   colnames(aa) <- c("gene", paste0("x",i,".mean"), paste0("x",i,".sd"))
#   assign(paste0("x",i,".ms"), data.table(aa))
# }

# rs_used <- fread("D:\\Python_Program\\PredictAP_web\\Predixcan_materials\\統整\\PrediXcan統整完檔案\\Num_of_rs_used_x.csv")
setwd("~/Predixcan_materials/Predixcan/code")
rs_used <- fread("../data/Num_of_rs_used_x_betty.csv")
setkeyv(rs_used, "gene")
temp_nam <- paste0("x",1:49,".ms")
for (j in 1:49) {
  dt <- setkeyv(eval(parse(text = temp_nam[j])), c("gene"))
  rs_used <- merge.data.table(rs_used, dt, sort = T, all = T)
}

# fwrite(rs_used, "D:\\Python_Program\\PredictAP_web\\Predixcan_materials\\統整\\PrediXcan統整完檔案\\dosage\\nfe\\result_table.csv")
fwrite(rs_used, "../data/result_table_betty.csv")


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


##R writing
library(devtools)
library(tidyverse)
use_git()
use_r("ecdf_fn") # 建立一個新的 R 函式檔案，名稱為 ecdf_fn.R
ecdf_fn <- function(x,perc) ecdf(x)(perc) # 在 ecdf_fn.R 中撰寫 ecdf_fn 函式
use_r("rankcal")
rankcal <- function(data, population = "afr", tissue = "x1") {
  data <- as.data.table(data)

  tissue_varname <- paste0(population, "_", tissue)
  ref_data <- eval(parse(text = tissue_varname))

  result_table_name <- paste0(population, "_result_table")
  result_table <- eval(parse(text = result_table_name))

  # 只保留有共同欄位的資料
  common_cols <- intersect(colnames(data)[-c(1, 2)], colnames(ref_data)[-c(1, 2)])
  final_cols <- c("FID", "IID", common_cols)
  data <- data[, ..final_cols]
  ref_data <- ref_data[, ..final_cols]

  # 開始計算百分位
  nsubject <- nrow(data)
  nc <- ncol(ref_data)
  rankm <- matrix(NA, nsubject, nc, byrow = FALSE)
  for (i in 3:nc) {
    rankm[, i] <- ecdf_fn(ref_data[[i]], data[[i]])
  }
  rankm[, 1] <- data[[1]]
  rankm[, 2] <- data[[2]]
  rankm <- data.table(rankm)
  colnames(rankm) <- colnames(data)

  l_data <- pivot_longer(data, cols = tidyselect::starts_with("ENSG"),
                         names_to = "gene", values_to = "expression")
  l_rank <- pivot_longer(rankm, cols = tidyselect::starts_with("ENSG"),
                         names_to = "gene", values_to = "percentile.rank")
  l_mer <- merge(l_data, l_rank, by = c("FID", "IID", "gene"))

  # 取得對應的 result_table 中資料
  res_por <- result_table %>%
    select(gene, all_of(tissue), tidyselect::starts_with(paste0(tissue, "."))) %>%
    filter(!is.na(eval(parse(text = tissue)))) %>%
    merge(l_mer, by = "gene")

  result <- res_por %>%
    select(FID, IID, gene, expression, percentile.rank, tidyselect::starts_with(tissue)) %>%
    arrange(FID, IID)

  colnames(result)[c(6, 7, 8)] <- c("N.snps", "mean", "sd")

  return(result)
}
 
ss <- fread("../data/example_data2.csv")


library(data.table)
setwd("~/Predixcan_materials/Predixcan/code")

## 讀入各族群資料，計算 result_table (mean, sd)
for (pop in c("EAS", "AFR")) {
  pop_lower <- tolower(pop)
  pop_files <- sort(list.files(paste0("../results/", pop, "/"), pattern="\\.txt$", full.names=TRUE))

  # 讀入資料，命名為 eas_x1~eas_x49 / afr_x1~afr_x49
  for (i in 1:49) {
    assign(paste0(pop_lower, "_x", i), fread(pop_files[i]))
  }

  # 計算每個 tissue 的 mean 及 sd
  for (i in 1:49) {
    current_dt <- eval(parse(text = paste0(pop_lower, "_x", i)))
    aa <- matrix(nrow = ncol(current_dt) - 2, ncol = 3, byrow = FALSE)
    aa[, 1] <- colnames(current_dt)[-c(1:2)]
    aa[, 2] <- as.numeric(apply(current_dt[, -c(1:2), with = FALSE], 2, mean))
    aa[, 3] <- as.numeric(apply(current_dt[, -c(1:2), with = FALSE], 2, sd))
    colnames(aa) <- c("gene", paste0("x", i, ".mean"), paste0("x", i, ".sd"))
    assign(paste0("x", i, ".ms_", pop_lower), data.table(aa))
  }

  # Merge with N.snps
  pop_rs <- fread("../data/Num_of_rs_used_x_betty.csv")
  setkeyv(pop_rs, "gene")
  for (j in 1:49) {
    dt <- setkeyv(eval(parse(text = paste0("x", j, ".ms_", pop_lower))), "gene")
    pop_rs <- merge.data.table(pop_rs, dt, sort = TRUE, all = TRUE)
  }

  assign(paste0(pop_lower, "_result_table"), pop_rs)
  fwrite(pop_rs, paste0("../data/result_table_", pop_lower, "_betty.csv"))
}


usethis::create_package("/mnt/data_40T/betty0101/Predixcan_materials/PredictAP.test")


use_data(eas_x1,  eas_x2,  eas_x3,  eas_x4,  eas_x5,  eas_x6,  eas_x7,  eas_x8,  eas_x9,  eas_x10,
         eas_x11, eas_x12, eas_x13, eas_x14, eas_x15, eas_x16, eas_x17, eas_x18, eas_x19, eas_x20,
         eas_x21, eas_x22, eas_x23, eas_x24, eas_x25, eas_x26, eas_x27, eas_x28, eas_x29, eas_x30,
         eas_x31, eas_x32, eas_x33, eas_x34, eas_x35, eas_x36, eas_x37, eas_x38, eas_x39, eas_x40,
         eas_x41, eas_x42, eas_x43, eas_x44, eas_x45, eas_x46, eas_x47, eas_x48, eas_x49,
         eas_result_table,
         afr_x1,  afr_x2,  afr_x3,  afr_x4,  afr_x5,  afr_x6,  afr_x7,  afr_x8,  afr_x9,  afr_x10,
         afr_x11, afr_x12, afr_x13, afr_x14, afr_x15, afr_x16, afr_x17, afr_x18, afr_x19, afr_x20,
         afr_x21, afr_x22, afr_x23, afr_x24, afr_x25, afr_x26, afr_x27, afr_x28, afr_x29, afr_x30,
         afr_x31, afr_x32, afr_x33, afr_x34, afr_x35, afr_x36, afr_x37, afr_x38, afr_x39, afr_x40,
         afr_x41, afr_x42, afr_x43, afr_x44, afr_x45, afr_x46, afr_x47, afr_x48, afr_x49,
         afr_result_table,
         compress = TRUE, overwrite = TRUE)


# test code
devtools::load_all("/mnt/data_40T/betty0101/Predixcan_materials/PredictAP.test")
rankcal(ss, population = "eas", tissue = "x1")
