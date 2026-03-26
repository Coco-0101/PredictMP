#??local?覡?w???ɮסA?]?ɮ׸??j?ݵ??ݤ@?q?ɶ?
library(PredictAP)
library(data.table)
x1_ex <- fread("/Users/chgsh14414/Desktop/Mac/Predixcan/統整/PrediXcan統整完檔案/example_data.csv")
x2_ex <- fread("example_data2.csv")
rank_result <- rankcal(x1_ex,"x1")
rank_result2 <- rankcal(x2_ex,"x2")

?rankcal

#????package
remove.packages("PredictAP")

lung <- fread("/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/results/_predicted_expression.txt")
rank_result <- rankcal(lung,"x32")
?rankcal
fwrite(rank_result, "/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/results/rank_result.csv")
rank <- fread("/Volumes/TOSHIBA_EXT/研究類備存文件/Lung_data相關/results/rank_result.csv")

