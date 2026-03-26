library(devtools)
use_git()
use_package("dplyr")
use_package("tidyr")
use_package("tidyverse", type = "depends")
use_package("data.table")
use_package("tidyselect")
use_package("magrittr")
use_r("ecdf_fn")
load_all()
check()
use_gpl_license()
##open ecdf_fn.R and insert roxygen skeleton. (?help)
document()
?ecdf_fn()

use_testthat()
use_test("ecdf_fn")
test()

##put data in package
library(data.table)
result_table <- fread("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/result_table.csv")
setwd("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/dosage/results")
ped_temp <- list.files(pattern=".txt")
ped_temp2 <- paste0("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/dosage/results/", ped_temp)
for (i in 1:49) {
  assign(paste0("x",i), fread(ped_temp2[i]))
}
setwd("/Users/chgsh14414/Desktop/Mac/Predixcan/predixcanEAS")
use_data(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,result_table,
         x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,
         x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,
         x31,x32,x33,x34,x35,x36,x37,x38,x39,x40,
         x41,x42,x43,x44,x45,x46,x47,x48,x49, internal = T, compress = "xz", overwrite = T)
load_all()
dim(x1)

##document data
use_r("data")


##main function
use_r("rankcal")
load_all()
##check if work
x1_ex <- fread("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/example_data.csv")
x2_ex <- fread("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/example_data2.csv")
rank_result <- rankcal(x2_ex, x2, 10, "x2")
rank_result <- rankcal(x1_ex, tissue = x1, nsubject = 10, rank.tis = "x1")

##add text in help
document()
?ecdf_fn()
?rankcal()

##create the data example in ?rankcal
x1_example <- fread("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/example_data.csv")
use_data(x1_example)
check()
build()

##solve lazydata load



#gitcreds::gitcreds_set()
#use_github(protocol = "https")
#use_github()

##2023//8/8 modify function
test2 <- function(data, tissue){
  ref_data <- eval(parse(text = tissue))
  nsubject <- nrow(data)
  rank.tis <- tissue
  nc <- ncol(ref_data)
  rankm <- matrix(NA, nsubject, nc, byrow = F)
  for (i in 3:nc) {
    rankm[,i] <- ecdf_fn(ref_data[[i]], data[[i]])
  }
  rankm[,1:2] <- data[,c(FID,IID)]
  rankm <- data.table(rankm)
  colnames(rankm) <- colnames(data)
  #return(rankm)
  l_data <- pivot_longer(data, cols = tidyselect::starts_with("ENSG"), names_to = "gene", values_to = "expression")
  l_rank <- pivot_longer(rankm, cols = tidyselect::starts_with("ENSG"), names_to = "gene", values_to = "rank")
  l_mer <- merge(l_data, l_rank, by = c("FID", "IID","gene"))
  res_por <- result_table %>% select(gene, all_of(tissue), tidyselect::starts_with(paste0(rank.tis,"."))) %>%
    filter(!is.na(eval(parse(text = tissue)))) %>% merge(l_mer, by = "gene")
  result <- res_por %>% select(FID, IID, gene, expression, rank, tidyselect::starts_with(rank.tis))%>% arrange(FID, IID)
  colnames(result)[c(4,7,8)] <- c("N.snps", "mean", "sd")
  return(result)
}
rank_result2 <- test2(x1_ex, tissue = "x1")

##retest
document()
load_all()
rank_result <- rankcal(x2_ex,"x2")
ex1 <- read.csv("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/example_data.csv")
ex1 <- as.data.table(ex1)
ee <- rankcal(ex1, "x1")
substitute((parse(text = tissue)))

##test code
x1_ex <- fread("example_data.csv")
x2_ex <- fread("example_data2.csv")
rank_result <- rankcal(x1_ex,"x1")
rank_result2 <- rankcal(x2_ex,"x2")
x1_ex <- as.data.table(x1_ex)

##2023/8/10 modify function
test2 <- function(data, tissue){
  data <- as.data.table(data)
  ref_data <- eval(parse(text = tissue))
  nsubject <- nrow(data)
  rank.tis <- tissue
  nc <- ncol(ref_data)
  rankm <- matrix(NA, nsubject, nc, byrow = F)
  for (i in 3:nc) {
    rankm[,i] <- ecdf_fn(ref_data[[i]], data[[i]])
  }
  rankm[,1:2] <- data[,c(FID,IID)]
  rankm <- data.table(rankm)
  colnames(rankm) <- colnames(data)
  #return(rankm)
  l_data <- pivot_longer(data, cols = tidyselect::starts_with("ENSG"), names_to = "gene", values_to = "expression")
  l_rank <- pivot_longer(rankm, cols = tidyselect::starts_with("ENSG"), names_to = "gene", values_to = "rank")
  l_mer <- merge(l_data, l_rank, by = c("FID", "IID","gene"))
  res_por <- result_table %>% select(gene, all_of(tissue), tidyselect::starts_with(paste0(rank.tis,"."))) %>%
    filter(!is.na(eval(parse(text = tissue)))) %>% merge(l_mer, by = "gene")
  result <- res_por %>% select(FID, IID, gene, expression, rank, tidyselect::starts_with(rank.tis))%>% arrange(FID, IID)
  colnames(result)[c(6,7,8)] <- c("N.snps", "mean", "sd")
  return(result)
}
ex2 <- x1_ex[,1:500]
r_ex2 <- test2(ex2, "x1")
rank_result <- rankcal(x1_ex,"x1")
##修改package version number
use_version()

##修改rancal函數example data
gene.name <- colnames(x1_ex)
use_data(gene.name)

#data(gene.name)
document()
load_all()
check()
build()
?rankcal()
rank_result <- rankcal(x1_ex, "x1")

use_r("data")

##將rank改成percentile
use_version()
load_all()
p_fn <- function(x,perc){round(ecdf(x)(perc)*100,1)}
x <- c(1:22)
x <- c(1,3,5,7,9,12)
y <- 10
p_fn(x,y)
round(p_fn(x,y),1)
p_fn(c(1:500), c(500,5,22))
fwrite(rank_result, "/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/rank_result.csv")
fwrite(rank_result2, "/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/rank_result2.csv")

