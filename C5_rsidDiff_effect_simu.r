################ *↓將全部的weight整理為一條向量↓* ################
# setwd("~/Predixcan_materials/Predixcan/code")
library("RSQLite")
library(data.table)
sqlite <- dbDriver("SQLite")
setwd("~/Predixcan_materials/elastic_net_models")
ped_temp <- list.files(pattern=".db")
# ped_temp2 <- paste0("/Users/chgsh14414/Desktop/Mac/Predixcan/elastic_net_models/", ped_temp)
ped_temp2 <- paste0("~/Predixcan_materials/elastic_net_models/", ped_temp)

for (i in 1:49) {
  assign(paste0("db",i) , dbConnect(sqlite,ped_temp2[i]))
}

db_name <- paste0("snp_sub", seq(1:49))
t_weigt <- c()
pedfiles <- paste0("db", seq(1:49))
for (i in 1:length(pedfiles)) {
  query <- function(...) dbGetQuery(eval(parse(text = pedfiles[i])), ...)
  assign(paste0("snp_sub",i) , data.table(query('select * from weights')))
  ww <- eval(parse(text = db_name[i]))$weight
  t_weigt <- c(t_weigt, ww)
}
t_weight <- data.table(abs(t_weigt))
colnames(t_weight) <- "value"
# fwrite(t_weight, "/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/t_weight.csv")
setwd("~/Predixcan_materials/Predixcan/code")
fwrite(t_weight, "../data/t_weight.csv")

##確認權重數量沒錯 
pp <- paste0("snp_sub",seq(1:49))
b <- 0
for (i in 1:49) {
  assign(paste0("db",i) , dbConnect(sqlite,ped_temp2[i]))
  a <- nrow(eval(parse(text = pp[i])))
  b= b+a
}
b
dbDisconnect(con)
##得全部權重數量：8558894
################ *↑將全部的weight整理為一條向量↑* ################

################ *↓Simulation↓* ################
library(data.table)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
# t_weight <- fread("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/t_weight.csv" )
t_weight <- fread("../data/t_weight.csv" )

quant <- quantile(t_weight$value, seq(0,1,0.2))
low.weight <- t_weight[abs(value) < quant[2]]
med.weight <- t_weight[abs(value) > quant[3] & abs(value) < quant[4]]
high.weight <- t_weight[abs(value) > quant[5]]

##抓取每個tissue預測一個基因最多會用到幾個snp，以方便決定模擬使用的snp數量
# comb <- fread("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/weight.csv")
comb <- fread("../data/weight.csv")
colnames(comb)[24] <- "Brain_Spinal_cord_cervical_c1"
colnames(comb)[28] <- "Cells_EBV_transformed_lymphocytes"

tis_name <- names(comb)[-c(1:5)]
minall <- c()
maxall <- c()
meanall <- c()
for (i in 1:49) {
  aa <- comb %>% select(gene, starts_with(paste0(tis_name[i]))) %>% filter(!is.na(eval(parse(text = tis_name[i]))))
  rscal <- table(aa[,1])
  minall[i] <- min(rscal)
  maxall[i] <- max(rscal)
  meanall[i] <- mean(rscal)
}
aa <- comb %>% select(gene, starts_with(paste0(tis_name[15]))) %>% filter(!is.na(eval(parse(text = tis_name[15]))))

maxall
meanall
minall
max(maxall)
summary(maxall)
#最終決定以1, 30, 200進行模擬
#最大是430! (Brain_Hippocampus的ENSG00000238083.7基因)

##模擬開始
library(data.table)
library(tidyr)
library(ggplot2)
library(ggpubr)
# t_weight <- fread("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/t_weight.csv" )
t_weight <- fread("../data/t_weight.csv" )

quant <- quantile(t_weight$value, seq(0,1,0.2))
low.weight <- t_weight[abs(value) < quant[2]]
med.weight <- t_weight[abs(value) > quant[3] & abs(value) < quant[4]]
high.weight <- t_weight[abs(value) > quant[5]]
#high.weight[1:139240] <- rnorm(139240, 10, 2)
summary(high.weight)

################ *↓設定有差異snp的weight權重高中低的情況↓* ################
##function
## diff_p:用來預測的SNP有多少比例是有族群差異的 (0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40)
## snp_n: 使用多少個SNP進行預測 (1 30 200)
## maf_p: Major allele frequency數值設定 (0.55 0.65 0.75 0.85 0.95)
## loop_n:模擬次數 (500)

#1 - high 有族群差異的SNP其權重很高
simul.fn.high <- function(diff_p, snp_n, maf_p, loop_n){
  pred.value <- c()
  for (i in 1:loop_n) {
    dosage <- sample(c(0,1,2), snp_n, prob = c(maf_p^2, 2*maf_p*(1-maf_p), (1-maf_p)^2), replace = T)
    effect.weight <- sample(high.weight$value, ceiling(snp_n*diff_p), replace = F)
    sample.weight <- sample(low.weight$value, floor(snp_n*(1-diff_p)), replace = F)
    merge.weight <- c(effect.weight, sample.weight)
    pred <- sum(dosage*merge.weight)
    pred.value[i] <- pred
  }
  #print(pred.value)
  return(pred.value)
}

#2 - low 有族群差異的SNP其權重很低
simul.fn.low <- function(diff_p, snp_n, maf_p, loop_n){
  pred.value <- c()
  for (i in 1:loop_n) {
    dosage <- sample(c(0,1,2), snp_n, prob = c(maf_p^2, 2*maf_p*(1-maf_p), (1-maf_p)^2), replace = T)
    sample.weight <- sample(low.weight$value, snp_n, replace = F)
    merge.weight <- sample.weight
    pred <- sum(dosage*merge.weight)
    pred.value[i] <- pred
  }
  #print(pred.value)
  return(pred.value)
}

#3 - medium 有族群差異的SNP其權重介於中間
simul.fn.med <- function(diff_p, snp_n, maf_p, loop_n){
  pred.value <- c()
  for (i in 1:loop_n) {
    dosage <- sample(c(0,1,2), snp_n, prob = c(maf_p^2, 2*maf_p*(1-maf_p), (1-maf_p)^2), replace = T)
    effect.weight <- sample(med.weight$value, ceiling(snp_n*diff_p), replace = F)
    sample.weight <- sample(low.weight$value, floor(snp_n*(1-diff_p)), replace = F)
    merge.weight <- c(effect.weight, sample.weight)
    pred <- sum(dosage*merge.weight)
    pred.value[i] <- pred
  }
  #print(pred.value)
  return(pred.value)
}
################ *↑設定有差異snp的weight權重高中低的情況↑* ################

################ *↓一次跑全部set的function↓* ################

mf_fn_mean <- function(freq, simtime){
  set.seed(230)
  ##generate plot
  ##1 - high, 1 snp
  diff.prop.set <- c(seq(0, 0.4, 0.05))
  value <- c()
  for (j in 1:length(diff.prop.set)) {
    print(paste("running", j))
    expron.value <- simul.fn.high(eval(parse(text = diff.prop.set[j])), 1, freq, simtime)
    value <- c(value, expron.value)
  }
  
  diffprop <- paste0(rep(seq(0, 0.4, 0.05), each = simtime))
  num <- rep(seq(1,simtime), 9)
  pp_pred <- data.frame(num, diffprop, value)
  pp_pred_short <- pp_pred %>% pivot_wider(names_from = diffprop, values_from = value)
  colnames(pp_pred_short)[-1] <- paste0("prop_", rep(seq(0, 0.4, 0.05)))
  pp_sh_no0 <- pp_pred_short[pp_pred_short$prop_0!=0,] ##計算關係將0排除
  abs_change1 <- apply(pp_sh_no0[,-c(1:2)], 2,
                       function(x) mean(abs(x - pp_sh_no0$prop_0)))
  mean_change1 <- apply(pp_sh_no0[,-c(1:2)], 2,
                        function(x) mean(abs(x - pp_sh_no0$prop_0)/pp_sh_no0$prop_0))
  
  plot1 <- ggplot(pp_pred, aes(x =value, group = diffprop, color = diffprop))+
    stat_ecdf(linewidth=1) +
    theme_bw(base_size = 28) +
    theme(legend.position ="top") +
    xlab("Gene expression") +
    ylab("High \n ECDF") +
    theme(legend.title=element_blank(), axis.text.x = element_text(angle=45, hjust=1)) +
    scale_x_continuous(n.breaks = 6)
  
  ##2 - high, 30 snp
  value <- c()
  for (j in 1:length(diff.prop.set)) {
    expron.value <- simul.fn.high(eval(parse(text = diff.prop.set[j])), 30, freq, simtime)
    value <- c(value, expron.value)
  }
  
  diffprop <- paste0(rep(seq(0, 0.4, 0.05), each = simtime))
  num <- rep(seq(1,simtime), 9)
  pp_pred <- data.frame(num, diffprop, value)
  pp_pred_short <- pp_pred %>% pivot_wider(names_from = diffprop, values_from = value)
  colnames(pp_pred_short)[-1] <- paste0("prop_", rep(seq(0, 0.4, 0.05)))
  pp_sh_no0 <- pp_pred_short[pp_pred_short$prop_0!=0,]
  abs_change2 <- apply(pp_sh_no0[,-c(1:2)], 2,
                       function(x) mean(abs(x - pp_sh_no0$prop_0)))
  mean_change2 <- apply(pp_sh_no0[,-c(1:2)], 2,
                        function(x) mean(abs(x - pp_sh_no0$prop_0)/pp_sh_no0$prop_0))
  
  
  plot2 <- ggplot(pp_pred, aes(x =value, group = diffprop, color = diffprop))+
    stat_ecdf(linewidth=1) +
    theme_bw(base_size = 28) +
    theme(legend.position ="top") +
    xlab("Gene expression") +
    ylab(" ") +
    theme(legend.title=element_blank(), axis.text.x = element_text(angle=45, hjust=1)) +
    scale_x_continuous(n.breaks = 6)
  
  ##3 - high, 200 snp
  value <- c()
  for (j in 1:length(diff.prop.set)) {
    expron.value <- simul.fn.high(eval(parse(text = diff.prop.set[j])), 200, freq, simtime)
    value <- c(value, expron.value)
  }
  
  diffprop <- paste0(rep(seq(0, 0.4, 0.05), each = simtime))
  num <- rep(seq(1,simtime), 9)
  pp_pred <- data.frame(num, diffprop, value)
  pp_pred_short <- pp_pred %>% pivot_wider(names_from = diffprop, values_from = value)
  colnames(pp_pred_short)[-1] <- paste0("prop_", rep(seq(0, 0.4, 0.05)))
  pp_sh_no0 <- pp_pred_short[pp_pred_short$prop_0!=0,]
  abs_change3 <- apply(pp_sh_no0[,-c(1:2)], 2,
                       function(x) mean(abs(x - pp_sh_no0$prop_0)))
  mean_change3 <- apply(pp_sh_no0[,-c(1:2)], 2,
                        function(x) mean(abs(x - pp_sh_no0$prop_0)/pp_sh_no0$prop_0))
  mean_change3
  
  plot3 <- ggplot(pp_pred, aes(x =value, group = diffprop, color = diffprop))+
    stat_ecdf(linewidth=1) +
    theme_bw(base_size = 28) +
    theme(legend.position ="top") +
    xlab("Gene expression") +
    ylab(" ") +
    theme(legend.title=element_blank(), axis.text.x = element_text(angle=45, hjust=1)) +
    scale_x_continuous(n.breaks = 6)
  
  ##4 - medium, 1 snp
  value <- c()
  for (j in 1:length(diff.prop.set)) {
    expron.value <- simul.fn.med(eval(parse(text = diff.prop.set[j])), 1, freq, simtime)
    value <- c(value, expron.value)
  }
  
  diffprop <- paste0(rep(seq(0, 0.4, 0.05), each = simtime))
  num <- rep(seq(1,simtime), 9)
  pp_pred <- data.frame(num, diffprop, value)
  pp_pred_short <- pp_pred %>% pivot_wider(names_from = diffprop, values_from = value)
  colnames(pp_pred_short)[-1] <- paste0("prop_", rep(seq(0, 0.4, 0.05)))
  pp_sh_no0 <- pp_pred_short[pp_pred_short$prop_0!=0,]
  abs_change4 <- apply(pp_sh_no0[,-c(1:2)], 2,
                       function(x) mean(abs(x - pp_sh_no0$prop_0)))
  mean_change4 <- apply(pp_sh_no0[,-c(1:2)], 2,
                        function(x) mean(abs(x - pp_sh_no0$prop_0)/pp_sh_no0$prop_0))
  mean_change4
  
  plot4 <- ggplot(pp_pred, aes(x =value, group = diffprop, color = diffprop))+
    stat_ecdf(linewidth=1) +
    theme_bw(base_size = 28) +
    theme(legend.position ="top") +
    xlab("Gene expression") +
    ylab("Medium \nECDF") +
    theme(legend.title=element_blank(), axis.text.x = element_text(angle=45, hjust=1)) +
    scale_x_continuous(n.breaks = 6, limits = c(0, max(pp_pred$value)))
  
  ##5 - medium, 30 snp
  value <- c()
  for (j in 1:length(diff.prop.set)) {
    expron.value <- simul.fn.med(eval(parse(text = diff.prop.set[j])), 30, freq, simtime)
    value <- c(value, expron.value)
  }
  
  diffprop <- paste0(rep(seq(0, 0.4, 0.05), each = simtime))
  num <- rep(seq(1,simtime), 9)
  pp_pred <- data.frame(num, diffprop, value)
  pp_pred_short <- pp_pred %>% pivot_wider(names_from = diffprop, values_from = value)
  colnames(pp_pred_short)[-1] <- paste0("prop_", rep(seq(0, 0.4, 0.05)))
  pp_sh_no0 <- pp_pred_short[pp_pred_short$prop_0!=0,]
  abs_change5 <- apply(pp_sh_no0[,-c(1:2)], 2,
                       function(x) mean(abs(x - pp_sh_no0$prop_0)))
  mean_change5 <- apply(pp_sh_no0[,-c(1:2)], 2,
                        function(x) mean(abs(x - pp_sh_no0$prop_0)/pp_sh_no0$prop_0))
  mean_change5
  plot5 <- ggplot(pp_pred, aes(x =value, group = diffprop, color = diffprop))+
    stat_ecdf(linewidth=1) +
    theme_bw(base_size = 28) +
    theme(legend.position ="top") +
    xlab("Gene expression") +
    ylab(" ") +
    theme(legend.title=element_blank(), axis.text.x = element_text(angle=45, hjust=1)) +
    scale_x_continuous(n.breaks = 6)
  
  ##6 - medium, 200 snp
  value <- c()
  for (j in 1:length(diff.prop.set)) {
    expron.value <- simul.fn.med(eval(parse(text = diff.prop.set[j])), 200, freq, simtime)
    value <- c(value, expron.value)
  }
  
  diffprop <- paste0(rep(seq(0, 0.4, 0.05), each = simtime))
  num <- rep(seq(1,simtime), 9)
  pp_pred <- data.frame(num, diffprop, value)
  pp_pred_short <- pp_pred %>% pivot_wider(names_from = diffprop, values_from = value)
  colnames(pp_pred_short)[-1] <- paste0("prop_", rep(seq(0, 0.4, 0.05)))
  pp_sh_no0 <- pp_pred_short[pp_pred_short$prop_0!=0,]
  abs_change6 <- apply(pp_sh_no0[,-c(1:2)], 2,
                       function(x) mean(abs(x - pp_sh_no0$prop_0)))
  mean_change6 <- apply(pp_sh_no0[,-c(1:2)], 2,
                        function(x) mean(abs(x - pp_sh_no0$prop_0)/pp_sh_no0$prop_0))
  mean_change6
  plot6 <- ggplot(pp_pred, aes(x =value, group = diffprop, color = diffprop))+
    stat_ecdf(linewidth=1) +
    theme_bw(base_size = 28) +
    theme(legend.position ="top") +
    xlab("Gene expression") +
    ylab(" ") +
    theme(legend.title=element_blank(), axis.text.x = element_text(angle=45, hjust=1)) +
    scale_x_continuous(n.breaks = 6)
  
  ##7 - low, 1 snp
  value <- c()
  for (j in 1:length(diff.prop.set)) {
    expron.value <- simul.fn.low(eval(parse(text = diff.prop.set[j])), 1, freq, simtime)
    value <- c(value, expron.value)
  }
  
  diffprop <- paste0(rep(seq(0, 0.4, 0.05), each = simtime))
  num <- rep(seq(1,simtime), 9)
  pp_pred <- data.frame(num, diffprop, value)
  pp_pred_short <- pp_pred %>% pivot_wider(names_from = diffprop, values_from = value)
  colnames(pp_pred_short)[-1] <- paste0("prop_", rep(seq(0, 0.4, 0.05)))
  pp_sh_no0 <- pp_pred_short[pp_pred_short$prop_0!=0,]
  abs_change7 <- apply(pp_sh_no0[,-c(1:2)], 2,
                       function(x) mean(abs(x - pp_sh_no0$prop_0)))
  mean_change7 <- apply(pp_sh_no0[,-c(1:2)], 2,
                        function(x) mean(abs(x - pp_sh_no0$prop_0)/pp_sh_no0$prop_0))
  mean_change7
  plot7 <- ggplot(pp_pred, aes(x =value, group = diffprop, color = diffprop))+
    stat_ecdf(linewidth=1) +
    theme_bw(base_size = 28) +
    theme(legend.position ="top") +
    xlab("1 variant \n Gene expression") +
    ylab("Low \n ECDF") +
    theme(legend.title=element_blank(), axis.text.x = element_text(angle=45, hjust=1)) +
    scale_x_continuous(n.breaks = 6)
  
  ##8 - low, 30 snp
  value <- c()
  for (j in 1:length(diff.prop.set)) {
    expron.value <- simul.fn.low(eval(parse(text = diff.prop.set[j])), 30, freq, simtime)
    value <- c(value, expron.value)
  }
  
  diffprop <- paste0(rep(seq(0, 0.4, 0.05), each = simtime))
  num <- rep(seq(1,simtime), 9)
  pp_pred <- data.frame(num, diffprop, value)
  pp_pred_short <- pp_pred %>% pivot_wider(names_from = diffprop, values_from = value)
  colnames(pp_pred_short)[-1] <- paste0("prop_", rep(seq(0, 0.4, 0.05)))
  pp_sh_no0 <- pp_pred_short[pp_pred_short$prop_0!=0,]
  abs_change8 <- apply(pp_sh_no0[,-c(1:2)], 2,
                       function(x) mean(abs(x - pp_sh_no0$prop_0)))
  mean_change8 <- apply(pp_sh_no0[,-c(1:2)], 2,
                        function(x) mean(abs(x - pp_sh_no0$prop_0)/pp_sh_no0$prop_0))
  mean_change8
  
  plot8 <- ggplot(pp_pred, aes(x =value, group = diffprop, color = diffprop))+
    stat_ecdf(linewidth=1) +
    theme_bw(base_size = 28) +
    theme(legend.position ="top") +
    xlab("30 variants \n Gene expression") +
    ylab(" ") +
    theme(legend.title=element_blank(), axis.text.x = element_text(angle=45, hjust=1)) +
    scale_x_continuous(n.breaks = 6)
  
  ##9 - low, 200 snp
  value <- c()
  for (j in 1:length(diff.prop.set)) {
    expron.value <- simul.fn.low(eval(parse(text = diff.prop.set[j])), 200, freq, simtime)
    value <- c(value, expron.value)
  }
  
  diffprop <- paste0(rep(seq(0, 0.4, 0.05), each = simtime))
  num <- rep(seq(1,simtime), 9)
  pp_pred <- data.frame(num, diffprop, value)
  pp_pred_short <- pp_pred %>% pivot_wider(names_from = diffprop, values_from = value)
  colnames(pp_pred_short)[-1] <- paste0("prop_", rep(seq(0, 0.4, 0.05)))
  pp_sh_no0 <- pp_pred_short[pp_pred_short$prop_0!=0,]
  abs_change9 <- apply(pp_sh_no0[,-c(1:2)], 2,
                       function(x) mean(abs(x - pp_sh_no0$prop_0)))
  mean_change9 <- apply(pp_sh_no0[,-c(1:2)], 2,
                        function(x) mean(abs(x - pp_sh_no0$prop_0)/pp_sh_no0$prop_0))
  
  change_table <- data.frame("high_1snp_abs"=abs_change1, "high_1snp_prop"=mean_change1,
                             "high_30snp_abs"=abs_change2, "high_30snp_prop"=mean_change2,
                             "high_200snp_abs"=abs_change3, "high_200snp_prop"=mean_change3,
                             "med_1snp_abs"=abs_change4, "med_1snp_prop"=mean_change4,
                             "med_30snp_abs"=abs_change5, "med_30snp_prop"=mean_change5,
                             "med_200snp_abs"=abs_change6, "med_200snp_prop"=mean_change6,
                             "low_1snp_abs"=abs_change7,"low_1snp_prop"=mean_change7,
                             "low_30snp_abs"=abs_change8, "low_30snp_prop"=mean_change8,
                             "low_200snp_abs"=abs_change9, "low_200snp_prop"=mean_change9)
  #return(change_table)
  plot9 <- ggplot(pp_pred, aes(x =value, group = diffprop, color = diffprop))+
    stat_ecdf(linewidth=1) +
    theme_bw(base_size = 28) +
    xlab("200 variants \n Gene expression") +
    ylab(" ") + 
    theme(legend.position ="top", legend.title=element_blank(), axis.text.x = element_text(angle=45, hjust=1)) +
    scale_x_continuous(n.breaks = 6)
  filename <- paste0("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/simu/mtcars_all_",freq,".png")
  all_plot <- ggarrange(plot1+rremove("legend"), plot2+rremove("legend"), 
                        plot3+rremove("legend"), plot4+rremove("legend"),
                        plot5+rremove("legend"), plot6+rremove("legend"),
                        plot7+rremove("legend"), plot8+rremove("legend"),
                        plot9+rremove("legend"), label.x = "try", common.legend = T,
                        nrow = 3, ncol = 3)
  #ggsave(filename, dpi = 300, width = 15, height = 10)
  return(list("value" = change_table, "plot" = all_plot))
}



##run all set
maf_freq_set <- seq(0.55, 0.95, 0.1)
# for (ff in 1:length(maf_freq_set)) {
#   result <- mf_fn_mean(maf_freq_set[ff],500)
#   result.mean <- result[[1]]
#   result.plot <- result[[2]]
#   #filename <- paste0("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/weight/simu_result/result_", maf_freq_set[ff], ".csv")
#   # filename2 <- paste0("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/simu/result_", maf_freq_set[ff], "_mean_v2.csv")
#   filename2 <- paste0("../data/simu_betty/result_", maf_freq_set[ff], "_mean_v2.csv")

#   write.csv(result.mean, filename2, row.names = T)
#   # plot.filename <- paste0("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/simu/thu_mtcars_all_",maf_freq_set[ff],"_v2.png")
#   plot.filename <- paste0("../data/simu_betty/thu_mtcars_all_",maf_freq_set[ff],"_v2.png")

#   ggsave(plot.filename, result.plot, dpi = 300, width = 15, height = 10)
# }

# Parallel processing
library(parallel)
cores <- detectCores() - 2
results <- mclapply(maf_freq_set, function(freq){

  result <- mf_fn_mean(freq, 500)

  result.mean <- result[[1]]
  result.plot <- result[[2]]

  filename2 <- paste0("../data/simu_betty/result_", freq, "_mean_v2.csv")
  write.csv(result.mean, filename2, row.names = TRUE)

  plot.filename <- paste0("../data/simu_betty/thu_mtcars_all_", freq, "_v2.png")
  ggsave(plot.filename, result.plot, dpi = 300, width = 15, height = 10)

  return(NULL)

}, mc.cores = cores)


################ *↓增加跑430最大值的結果↓* ################
################ *↓設定有差異snp的weight權重高中低的情況↓* ################
##function
## diff_p:用來預測的SNP有多少比例是有族群差異的 (0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40)
## snp_n: 使用多少個SNP進行預測 430
## maf_p: Major allele frequency數值設定 (0.55 0.65 0.75 0.85 0.95)
## loop_n:模擬次數 (500)

#1 - high 有族群差異的SNP其權重很高
simul.fn.high <- function(diff_p, snp_n, maf_p, loop_n){
  pred.value <- c()
  for (i in 1:loop_n) {
    dosage <- sample(c(0,1,2), snp_n, prob = c(maf_p^2, 2*maf_p*(1-maf_p), (1-maf_p)^2), replace = T)
    effect.weight <- sample(high.weight$value, ceiling(snp_n*diff_p), replace = F)
    sample.weight <- sample(low.weight$value, floor(snp_n*(1-diff_p)), replace = F)
    merge.weight <- c(effect.weight, sample.weight)
    pred <- sum(dosage*merge.weight)
    pred.value[i] <- pred
  }
  #print(pred.value)
  return(pred.value)
}

#2 - low 有族群差異的SNP其權重很低
simul.fn.low <- function(diff_p, snp_n, maf_p, loop_n){
  pred.value <- c()
  for (i in 1:loop_n) {
    dosage <- sample(c(0,1,2), snp_n, prob = c(maf_p^2, 2*maf_p*(1-maf_p), (1-maf_p)^2), replace = T)
    sample.weight <- sample(low.weight$value, snp_n, replace = F)
    merge.weight <- sample.weight
    pred <- sum(dosage*merge.weight)
    pred.value[i] <- pred
  }
  #print(pred.value)
  return(pred.value)
}

#3 - medium 有族群差異的SNP其權重介於中間
simul.fn.med <- function(diff_p, snp_n, maf_p, loop_n){
  pred.value <- c()
  for (i in 1:loop_n) {
    dosage <- sample(c(0,1,2), snp_n, prob = c(maf_p^2, 2*maf_p*(1-maf_p), (1-maf_p)^2), replace = T)
    effect.weight <- sample(med.weight$value, ceiling(snp_n*diff_p), replace = F)
    sample.weight <- sample(low.weight$value, floor(snp_n*(1-diff_p)), replace = F)
    merge.weight <- c(effect.weight, sample.weight)
    pred <- sum(dosage*merge.weight)
    pred.value[i] <- pred
  }
  #print(pred.value)
  return(pred.value)
}
################ *↓一次跑全部set的function↓* ################
################ *↓一次跑全部set的function↓* ################
set.seed(230)
run_sim <- function(maf, type){

  diff.prop.set <- seq(0,0.4,0.05)

  value_list <- lapply(diff.prop.set, function(p){

    if(type=="high"){
      simul.fn.high(p,430,maf,500)
    }else if(type=="med"){
      simul.fn.med(p,430,maf,500)
    }else{
      simul.fn.low(p,430,maf,500)
    }

  })

  value <- unlist(value_list)

  diffprop <- rep(diff.prop.set, each=500)
  num <- rep(1:500,9)

  pp_pred <- data.frame(num,diffprop,value)

  pp_pred_short <- tidyr::pivot_wider(pp_pred,
                     names_from=diffprop,
                     values_from=value)

  colnames(pp_pred_short)[-1] <- paste0("prop_",diff.prop.set)

  pp_sh_no0 <- pp_pred_short[pp_pred_short$prop_0!=0,]

  abs_change <- apply(pp_sh_no0[,-c(1:2)],2,
                      function(x) mean(abs(x-pp_sh_no0$prop_0)))

  mean_change <- apply(pp_sh_no0[,-c(1:2)],2,
                       function(x) mean(abs(x-pp_sh_no0$prop_0)/pp_sh_no0$prop_0))

  plot <- ggplot(pp_pred,aes(x=value,group=diffprop,color=diffprop))+
    stat_ecdf(linewidth=1)+
    theme_bw(base_size=18)

  list(abs=abs_change,mean=mean_change,plot=plot)
}

maf_set <- c(0.95,0.85,0.75,0.65,0.55)

task_grid <- expand.grid(
  maf = maf_set,
  type = c("high","med","low")
)

results <- mclapply(1:nrow(task_grid), function(i){

  maf  <- task_grid$maf[i]
  type <- task_grid$type[i]

  run_sim(maf,type)

}, mc.cores = cores)

"""
##generate plot
##1 - MAF=0.05,high, 430 snp
diff.prop.set <- c(seq(0, 0.4, 0.05))
value <- c()
for (j in 1:length(diff.prop.set)) {
  expron.value <- simul.fn.high(eval(parse(text = diff.prop.set[j])), 430, 0.95, 500)
  value <- c(value, expron.value)
}

diffprop <- paste0(rep(seq(0, 0.4, 0.05), each = 500))
num <- rep(seq(1,500), 9)
pp_pred <- data.frame(num, diffprop, value)
pp_pred_short <- pp_pred %>% pivot_wider(names_from = diffprop, values_from = value)
colnames(pp_pred_short)[-1] <- paste0("prop_", rep(seq(0, 0.4, 0.05)))
pp_sh_no0 <- pp_pred_short[pp_pred_short$prop_0!=0,] ##計算關係將0排除
abs_change1 <- apply(pp_sh_no0[,-c(1:2)], 2,
                     function(x) mean(abs(x - pp_sh_no0$prop_0)))
mean_change1 <- apply(pp_sh_no0[,-c(1:2)], 2,
                      function(x) mean(abs(x - pp_sh_no0$prop_0)/pp_sh_no0$prop_0))

plot1 <- ggplot(pp_pred, aes(x =value, group = diffprop, color = diffprop))+
  stat_ecdf(linewidth=1) +
  theme_bw(base_size = 18) +
  theme(legend.position ="top") +
  xlab("Gene expression") +
  ylab("High \n ECDF") +
  theme(legend.title=element_blank(), axis.text.x = element_text(angle=45, hjust=1)) +
  scale_x_continuous(n.breaks = 6)

##2 - MAF=0.15,high, 430 snp
diff.prop.set <- c(seq(0, 0.4, 0.05))
value <- c()
for (j in 1:length(diff.prop.set)) {
  expron.value <- simul.fn.high(eval(parse(text = diff.prop.set[j])), 430, 0.85, 500)
  value <- c(value, expron.value)
}

diffprop <- paste0(rep(seq(0, 0.4, 0.05), each = 500))
num <- rep(seq(1,500), 9)
pp_pred <- data.frame(num, diffprop, value)
pp_pred_short <- pp_pred %>% pivot_wider(names_from = diffprop, values_from = value)
colnames(pp_pred_short)[-1] <- paste0("prop_", rep(seq(0, 0.4, 0.05)))
pp_sh_no0 <- pp_pred_short[pp_pred_short$prop_0!=0,] ##計算關係將0排除
abs_change2 <- apply(pp_sh_no0[,-c(1:2)], 2,
                     function(x) mean(abs(x - pp_sh_no0$prop_0)))
mean_change2 <- apply(pp_sh_no0[,-c(1:2)], 2,
                      function(x) mean(abs(x - pp_sh_no0$prop_0)/pp_sh_no0$prop_0))

plot2 <- ggplot(pp_pred, aes(x =value, group = diffprop, color = diffprop))+
  stat_ecdf(linewidth=1) +
  theme_bw(base_size = 18) +
  theme(legend.position ="top") +
  xlab("Gene expression") +
  ylab(" ") + 
  theme(legend.title=element_blank(), axis.text.x = element_text(angle=45, hjust=1)) +
  scale_x_continuous(n.breaks = 6)

##3 - MAF=0.25,high, 430 snp
diff.prop.set <- c(seq(0, 0.4, 0.05))
value <- c()
for (j in 1:length(diff.prop.set)) {
  expron.value <- simul.fn.high(eval(parse(text = diff.prop.set[j])), 430, 0.75, 500)
  value <- c(value, expron.value)
}

diffprop <- paste0(rep(seq(0, 0.4, 0.05), each = 500))
num <- rep(seq(1,500), 9)
pp_pred <- data.frame(num, diffprop, value)
pp_pred_short <- pp_pred %>% pivot_wider(names_from = diffprop, values_from = value)
colnames(pp_pred_short)[-1] <- paste0("prop_", rep(seq(0, 0.4, 0.05)))
pp_sh_no0 <- pp_pred_short[pp_pred_short$prop_0!=0,] ##計算關係將0排除
abs_change3 <- apply(pp_sh_no0[,-c(1:2)], 2,
                     function(x) mean(abs(x - pp_sh_no0$prop_0)))
mean_change3 <- apply(pp_sh_no0[,-c(1:2)], 2,
                      function(x) mean(abs(x - pp_sh_no0$prop_0)/pp_sh_no0$prop_0))

plot3 <- ggplot(pp_pred, aes(x =value, group = diffprop, color = diffprop))+
  stat_ecdf(linewidth=1) +
  theme_bw(base_size = 18) +
  theme(legend.position ="top") +
  xlab("Gene expression") +
  ylab(" ") + 
  theme(legend.title=element_blank(), axis.text.x = element_text(angle=45, hjust=1)) +
  scale_x_continuous(n.breaks = 6)

##4 - MAF=0.35,high, 430 snp
diff.prop.set <- c(seq(0, 0.4, 0.05))
value <- c()
for (j in 1:length(diff.prop.set)) {
  expron.value <- simul.fn.high(eval(parse(text = diff.prop.set[j])), 430, 0.65, 500)
  value <- c(value, expron.value)
}

diffprop <- paste0(rep(seq(0, 0.4, 0.05), each = 500))
num <- rep(seq(1,500), 9)
pp_pred <- data.frame(num, diffprop, value)
pp_pred_short <- pp_pred %>% pivot_wider(names_from = diffprop, values_from = value)
colnames(pp_pred_short)[-1] <- paste0("prop_", rep(seq(0, 0.4, 0.05)))
pp_sh_no0 <- pp_pred_short[pp_pred_short$prop_0!=0,] ##計算關係將0排除
abs_change4 <- apply(pp_sh_no0[,-c(1:2)], 2,
                     function(x) mean(abs(x - pp_sh_no0$prop_0)))
mean_change4 <- apply(pp_sh_no0[,-c(1:2)], 2,
                      function(x) mean(abs(x - pp_sh_no0$prop_0)/pp_sh_no0$prop_0))

plot4 <- ggplot(pp_pred, aes(x =value, group = diffprop, color = diffprop))+
  stat_ecdf(linewidth=1) +
  theme_bw(base_size = 18) +
  theme(legend.position ="top") +
  xlab("Gene expression") +
  ylab(" ") + 
  theme(legend.title=element_blank(), axis.text.x = element_text(angle=45, hjust=1)) +
  scale_x_continuous(n.breaks = 6)

##5 - MAF=0.45,high, 430 snp
diff.prop.set <- c(seq(0, 0.4, 0.05))
value <- c()
for (j in 1:length(diff.prop.set)) {
  expron.value <- simul.fn.high(eval(parse(text = diff.prop.set[j])), 430, 0.55, 500)
  value <- c(value, expron.value)
}

diffprop <- paste0(rep(seq(0, 0.4, 0.05), each = 500))
num <- rep(seq(1,500), 9)
pp_pred <- data.frame(num, diffprop, value)
pp_pred_short <- pp_pred %>% pivot_wider(names_from = diffprop, values_from = value)
colnames(pp_pred_short)[-1] <- paste0("prop_", rep(seq(0, 0.4, 0.05)))
pp_sh_no0 <- pp_pred_short[pp_pred_short$prop_0!=0,] ##計算關係將0排除
abs_change5 <- apply(pp_sh_no0[,-c(1:2)], 2,
                     function(x) mean(abs(x - pp_sh_no0$prop_0)))
mean_change5 <- apply(pp_sh_no0[,-c(1:2)], 2,
                      function(x) mean(abs(x - pp_sh_no0$prop_0)/pp_sh_no0$prop_0))

plot5 <- ggplot(pp_pred, aes(x =value, group = diffprop, color = diffprop))+
  stat_ecdf(linewidth=1) +
  theme_bw(base_size = 18) +
  theme(legend.position ="top") +
  xlab("Gene expression") +
  ylab(" ") + 
  theme(legend.title=element_blank(), axis.text.x = element_text(angle=45, hjust=1)) +
  scale_x_continuous(n.breaks = 6)

##6 - MAF=0.05,med, 430 snp
diff.prop.set <- c(seq(0, 0.4, 0.05))
value <- c()
for (j in 1:length(diff.prop.set)) {
  expron.value <- simul.fn.med(eval(parse(text = diff.prop.set[j])), 430, 0.95, 500)
  value <- c(value, expron.value)
}

diffprop <- paste0(rep(seq(0, 0.4, 0.05), each = 500))
num <- rep(seq(1,500), 9)
pp_pred <- data.frame(num, diffprop, value)
pp_pred_short <- pp_pred %>% pivot_wider(names_from = diffprop, values_from = value)
colnames(pp_pred_short)[-1] <- paste0("prop_", rep(seq(0, 0.4, 0.05)))
pp_sh_no0 <- pp_pred_short[pp_pred_short$prop_0!=0,] ##計算關係將0排除
abs_change6 <- apply(pp_sh_no0[,-c(1:2)], 2,
                     function(x) mean(abs(x - pp_sh_no0$prop_0)))
mean_change6 <- apply(pp_sh_no0[,-c(1:2)], 2,
                      function(x) mean(abs(x - pp_sh_no0$prop_0)/pp_sh_no0$prop_0))

plot6 <- ggplot(pp_pred, aes(x =value, group = diffprop, color = diffprop))+
  stat_ecdf(linewidth=1) +
  theme_bw(base_size = 18) +
  theme(legend.position ="top") +
  xlab("Gene expression") +
  ylab("Medium \n ECDF") +
  theme(legend.title=element_blank(), axis.text.x = element_text(angle=45, hjust=1)) +
  scale_x_continuous(n.breaks = 6)

##7 - MAF=0.15,med, 430 snp
diff.prop.set <- c(seq(0, 0.4, 0.05))
value <- c()
for (j in 1:length(diff.prop.set)) {
  expron.value <- simul.fn.med(eval(parse(text = diff.prop.set[j])), 430, 0.85, 500)
  value <- c(value, expron.value)
}

diffprop <- paste0(rep(seq(0, 0.4, 0.05), each = 500))
num <- rep(seq(1,500), 9)
pp_pred <- data.frame(num, diffprop, value)
pp_pred_short <- pp_pred %>% pivot_wider(names_from = diffprop, values_from = value)
colnames(pp_pred_short)[-1] <- paste0("prop_", rep(seq(0, 0.4, 0.05)))
pp_sh_no0 <- pp_pred_short[pp_pred_short$prop_0!=0,] ##計算關係將0排除
abs_change7 <- apply(pp_sh_no0[,-c(1:2)], 2,
                     function(x) mean(abs(x - pp_sh_no0$prop_0)))
mean_change7 <- apply(pp_sh_no0[,-c(1:2)], 2,
                      function(x) mean(abs(x - pp_sh_no0$prop_0)/pp_sh_no0$prop_0))

plot7 <- ggplot(pp_pred, aes(x =value, group = diffprop, color = diffprop))+
  stat_ecdf(linewidth=1) +
  theme_bw(base_size = 18) +
  theme(legend.position ="top") +
  xlab("Gene expression") +
  ylab(" ") + 
  theme(legend.title=element_blank(), axis.text.x = element_text(angle=45, hjust=1)) +
  scale_x_continuous(n.breaks = 6)

##8 - MAF=0.25,med, 430 snp
diff.prop.set <- c(seq(0, 0.4, 0.05))
value <- c()
for (j in 1:length(diff.prop.set)) {
  expron.value <- simul.fn.med(eval(parse(text = diff.prop.set[j])), 430, 0.75, 500)
  value <- c(value, expron.value)
}

diffprop <- paste0(rep(seq(0, 0.4, 0.05), each = 500))
num <- rep(seq(1,500), 9)
pp_pred <- data.frame(num, diffprop, value)
pp_pred_short <- pp_pred %>% pivot_wider(names_from = diffprop, values_from = value)
colnames(pp_pred_short)[-1] <- paste0("prop_", rep(seq(0, 0.4, 0.05)))
pp_sh_no0 <- pp_pred_short[pp_pred_short$prop_0!=0,] ##計算關係將0排除
abs_change8 <- apply(pp_sh_no0[,-c(1:2)], 2,
                     function(x) mean(abs(x - pp_sh_no0$prop_0)))
mean_change8 <- apply(pp_sh_no0[,-c(1:2)], 2,
                      function(x) mean(abs(x - pp_sh_no0$prop_0)/pp_sh_no0$prop_0))

plot8 <- ggplot(pp_pred, aes(x =value, group = diffprop, color = diffprop))+
  stat_ecdf(linewidth=1) +
  theme_bw(base_size = 18) +
  theme(legend.position ="top") +
  xlab("Gene expression") +
  ylab(" ") + 
  theme(legend.title=element_blank(), axis.text.x = element_text(angle=45, hjust=1)) +
  scale_x_continuous(n.breaks = 6)

##9 - MAF=0.35,med, 430 snp
diff.prop.set <- c(seq(0, 0.4, 0.05))
value <- c()
for (j in 1:length(diff.prop.set)) {
  expron.value <- simul.fn.med(eval(parse(text = diff.prop.set[j])), 430, 0.65, 500)
  value <- c(value, expron.value)
}

diffprop <- paste0(rep(seq(0, 0.4, 0.05), each = 500))
num <- rep(seq(1,500), 9)
pp_pred <- data.frame(num, diffprop, value)
pp_pred_short <- pp_pred %>% pivot_wider(names_from = diffprop, values_from = value)
colnames(pp_pred_short)[-1] <- paste0("prop_", rep(seq(0, 0.4, 0.05)))
pp_sh_no0 <- pp_pred_short[pp_pred_short$prop_0!=0,] ##計算關係將0排除
abs_change9 <- apply(pp_sh_no0[,-c(1:2)], 2,
                     function(x) mean(abs(x - pp_sh_no0$prop_0)))
mean_change9 <- apply(pp_sh_no0[,-c(1:2)], 2,
                      function(x) mean(abs(x - pp_sh_no0$prop_0)/pp_sh_no0$prop_0))

plot9 <- ggplot(pp_pred, aes(x =value, group = diffprop, color = diffprop))+
  stat_ecdf(linewidth=1) +
  theme_bw(base_size = 18) +
  theme(legend.position ="top") +
  xlab("Gene expression") +
  ylab(" ") + 
  theme(legend.title=element_blank(), axis.text.x = element_text(angle=45, hjust=1)) +
  scale_x_continuous(n.breaks = 6)

##10 - MAF=0.45,med, 430 snp
diff.prop.set <- c(seq(0, 0.4, 0.05))
value <- c()
for (j in 1:length(diff.prop.set)) {
  expron.value <- simul.fn.med(eval(parse(text = diff.prop.set[j])), 430, 0.55, 500)
  value <- c(value, expron.value)
}

diffprop <- paste0(rep(seq(0, 0.4, 0.05), each = 500))
num <- rep(seq(1,500), 9)
pp_pred <- data.frame(num, diffprop, value)
pp_pred_short <- pp_pred %>% pivot_wider(names_from = diffprop, values_from = value)
colnames(pp_pred_short)[-1] <- paste0("prop_", rep(seq(0, 0.4, 0.05)))
pp_sh_no0 <- pp_pred_short[pp_pred_short$prop_0!=0,] ##計算關係將0排除
abs_change10 <- apply(pp_sh_no0[,-c(1:2)], 2,
                      function(x) mean(abs(x - pp_sh_no0$prop_0)))
mean_change10 <- apply(pp_sh_no0[,-c(1:2)], 2,
                       function(x) mean(abs(x - pp_sh_no0$prop_0)/pp_sh_no0$prop_0))

plot10 <- ggplot(pp_pred, aes(x =value, group = diffprop, color = diffprop))+
  stat_ecdf(linewidth=1) +
  theme_bw(base_size = 18) +
  theme(legend.position ="top") +
  xlab("Gene expression") +
  ylab(" ") + 
  theme(legend.title=element_blank(), axis.text.x = element_text(angle=45, hjust=1)) +
  scale_x_continuous(n.breaks = 6)

##11 - MAF=0.05,low, 430 snp
diff.prop.set <- c(seq(0, 0.4, 0.05))
value <- c()
for (j in 1:length(diff.prop.set)) {
  expron.value <- simul.fn.low(eval(parse(text = diff.prop.set[j])), 430, 0.95, 500)
  value <- c(value, expron.value)
}

diffprop <- paste0(rep(seq(0, 0.4, 0.05), each = 500))
num <- rep(seq(1,500), 9)
pp_pred <- data.frame(num, diffprop, value)
pp_pred_short <- pp_pred %>% pivot_wider(names_from = diffprop, values_from = value)
colnames(pp_pred_short)[-1] <- paste0("prop_", rep(seq(0, 0.4, 0.05)))
pp_sh_no0 <- pp_pred_short[pp_pred_short$prop_0!=0,] ##計算關係將0排除
abs_change11 <- apply(pp_sh_no0[,-c(1:2)], 2,
                      function(x) mean(abs(x - pp_sh_no0$prop_0)))
mean_change11 <- apply(pp_sh_no0[,-c(1:2)], 2,
                       function(x) mean(abs(x - pp_sh_no0$prop_0)/pp_sh_no0$prop_0))

plot11 <- ggplot(pp_pred, aes(x =value, group = diffprop, color = diffprop))+
  stat_ecdf(linewidth=1) +
  theme_bw(base_size = 18) +
  theme(legend.position ="top") +
  xlab("Gene expression \n 0.05") +
  ylab("Low \n ECDF") +
  theme(legend.title=element_blank(), axis.text.x = element_text(angle=45, hjust=1)) +
  scale_x_continuous(n.breaks = 6)

##12 - MAF=0.15,low, 430 snp
diff.prop.set <- c(seq(0, 0.4, 0.05))
value <- c()
for (j in 1:length(diff.prop.set)) {
  expron.value <- simul.fn.low(eval(parse(text = diff.prop.set[j])), 430, 0.85, 500)
  value <- c(value, expron.value)
}

diffprop <- paste0(rep(seq(0, 0.4, 0.05), each = 500))
num <- rep(seq(1,500), 9)
pp_pred <- data.frame(num, diffprop, value)
pp_pred_short <- pp_pred %>% pivot_wider(names_from = diffprop, values_from = value)
colnames(pp_pred_short)[-1] <- paste0("prop_", rep(seq(0, 0.4, 0.05)))
pp_sh_no0 <- pp_pred_short[pp_pred_short$prop_0!=0,] ##計算關係將0排除
abs_change12 <- apply(pp_sh_no0[,-c(1:2)], 2,
                      function(x) mean(abs(x - pp_sh_no0$prop_0)))
mean_change12 <- apply(pp_sh_no0[,-c(1:2)], 2,
                       function(x) mean(abs(x - pp_sh_no0$prop_0)/pp_sh_no0$prop_0))

plot12 <- ggplot(pp_pred, aes(x =value, group = diffprop, color = diffprop))+
  stat_ecdf(linewidth=1) +
  theme_bw(base_size = 18) +
  theme(legend.position ="top") +
  xlab("Gene expression \n 0.15") +
  ylab(" ") + 
  theme(legend.title=element_blank(), axis.text.x = element_text(angle=45, hjust=1)) +
  scale_x_continuous(n.breaks = 6)

##13 - MAF=0.25,low, 430 snp
diff.prop.set <- c(seq(0, 0.4, 0.05))
value <- c()
for (j in 1:length(diff.prop.set)) {
  expron.value <- simul.fn.low(eval(parse(text = diff.prop.set[j])), 430, 0.75, 500)
  value <- c(value, expron.value)
}

diffprop <- paste0(rep(seq(0, 0.4, 0.05), each = 500))
num <- rep(seq(1,500), 9)
pp_pred <- data.frame(num, diffprop, value)
pp_pred_short <- pp_pred %>% pivot_wider(names_from = diffprop, values_from = value)
colnames(pp_pred_short)[-1] <- paste0("prop_", rep(seq(0, 0.4, 0.05)))
pp_sh_no0 <- pp_pred_short[pp_pred_short$prop_0!=0,] ##計算關係將0排除
abs_change13 <- apply(pp_sh_no0[,-c(1:2)], 2,
                      function(x) mean(abs(x - pp_sh_no0$prop_0)))
mean_change13 <- apply(pp_sh_no0[,-c(1:2)], 2,
                       function(x) mean(abs(x - pp_sh_no0$prop_0)/pp_sh_no0$prop_0))

plot13 <- ggplot(pp_pred, aes(x =value, group = diffprop, color = diffprop))+
  stat_ecdf(linewidth=1) +
  theme_bw(base_size = 18) +
  theme(legend.position ="top") +
  xlab("Gene expression \n 0.25") +
  ylab(" ") + 
  theme(legend.title=element_blank(), axis.text.x = element_text(angle=45, hjust=1)) +
  scale_x_continuous(n.breaks = 6)

##14 - MAF=0.35,low, 430 snp
diff.prop.set <- c(seq(0, 0.4, 0.05))
value <- c()
for (j in 1:length(diff.prop.set)) {
  expron.value <- simul.fn.low(eval(parse(text = diff.prop.set[j])), 430, 0.65, 500)
  value <- c(value, expron.value)
}

diffprop <- paste0(rep(seq(0, 0.4, 0.05), each = 500))
num <- rep(seq(1,500), 9)
pp_pred <- data.frame(num, diffprop, value)
pp_pred_short <- pp_pred %>% pivot_wider(names_from = diffprop, values_from = value)
colnames(pp_pred_short)[-1] <- paste0("prop_", rep(seq(0, 0.4, 0.05)))
pp_sh_no0 <- pp_pred_short[pp_pred_short$prop_0!=0,] ##計算關係將0排除
abs_change14 <- apply(pp_sh_no0[,-c(1:2)], 2,
                      function(x) mean(abs(x - pp_sh_no0$prop_0)))
mean_change14 <- apply(pp_sh_no0[,-c(1:2)], 2,
                       function(x) mean(abs(x - pp_sh_no0$prop_0)/pp_sh_no0$prop_0))

plot14 <- ggplot(pp_pred, aes(x =value, group = diffprop, color = diffprop))+
  stat_ecdf(linewidth=1) +
  theme_bw(base_size = 18) +
  theme(legend.position ="top") +
  xlab("Gene expression \n 0.35") +
  ylab(" ") + 
  theme(legend.title=element_blank(), axis.text.x = element_text(angle=45, hjust=1)) +
  scale_x_continuous(n.breaks = 6)

##15 - MAF=0.45,low, 430 snp
diff.prop.set <- c(seq(0, 0.4, 0.05))
value <- c()
for (j in 1:length(diff.prop.set)) {
  expron.value <- simul.fn.low(eval(parse(text = diff.prop.set[j])), 430, 0.55, 500)
  value <- c(value, expron.value)
}

diffprop <- paste0(rep(seq(0, 0.4, 0.05), each = 500))
num <- rep(seq(1,500), 9)
pp_pred <- data.frame(num, diffprop, value)
pp_pred_short <- pp_pred %>% pivot_wider(names_from = diffprop, values_from = value)
colnames(pp_pred_short)[-1] <- paste0("prop_", rep(seq(0, 0.4, 0.05)))
pp_sh_no0 <- pp_pred_short[pp_pred_short$prop_0!=0,] ##計算關係將0排除
abs_change15 <- apply(pp_sh_no0[,-c(1:2)], 2,
                      function(x) mean(abs(x - pp_sh_no0$prop_0)))
mean_change15 <- apply(pp_sh_no0[,-c(1:2)], 2,
                       function(x) mean(abs(x - pp_sh_no0$prop_0)/pp_sh_no0$prop_0))

plot15 <- ggplot(pp_pred, aes(x =value, group = diffprop, color = diffprop))+
  stat_ecdf(linewidth=1) +
  theme_bw(base_size = 18) +
  theme(legend.position ="top") +
  xlab("Gene expression  \n 0.45") +
  ylab(" ") + 
  theme(legend.title=element_blank(), axis.text.x = element_text(angle=45, hjust=1)) +
  scale_x_continuous(n.breaks = 6)

"""

#change_table <- data.frame("high_430snp_abs"=abs_change1, "high_430snp_prop"=mean_change1,
#"med_430snp_abs"=abs_change4, "med_430snp_prop"=mean_change4,
#"low_430snp_abs"=abs_change7,"low_430snp_prop"=mean_change7)
#filename <- paste0("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/simu/mtcars_all_",freq,".png")

"""
all_plot <- ggarrange(plot1+rremove("legend"), plot2+rremove("legend"), plot3+rremove("legend"), plot4+rremove("legend"), plot5+rremove("legend"),
                      plot6+rremove("legend"), plot7+rremove("legend"), plot8+rremove("legend"), plot9+rremove("legend"), plot10+rremove("legend"),
                      plot11+rremove("legend"), plot12+rremove("legend"), plot13+rremove("legend"), plot14+rremove("legend"), plot15+rremove("legend"),
                      label.x = "try", common.legend = T,
                      nrow = 3, ncol = 5)
# plot.filename <- paste0("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/simu/thu_mtcars_max430.png")
plot.filename <- paste0("../data/simu_betty/thu_mtcars_max430.png")

ggsave(plot.filename, all_plot, dpi = 300, width = 15, height = 10)
"""
plots <- lapply(results, function(x) x$plot)
plots <- lapply(plots, function(p) p + rremove("legend"))

all_plot <- ggarrange(
  plotlist = plots,
  common.legend = TRUE,
  nrow = 3,
  ncol = 5
)

plot.filename <- "../data/simu_betty/thu_mtcars_max430.png"

ggsave(
  plot.filename,
  all_plot,
  dpi = 300,
  width = 15,
  height = 10
)