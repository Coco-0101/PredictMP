## C1_PrediXcan_code
### 統整PrediXcan所有權重model的資訊
- 權重模型版本：elastic_net_models
- All mapped to the GRCh38
- 49 tissues
- range of gene prediction: (1,642, 10,012)
- range of # of used SNPs: (63,654, 309,155)
- total # of used SNPs after merging all tissues: 5,318,682 SNPs (contain duplicated SNPs due to diff gene prediction)

## C2_gnomAD_code
### 整理從gnomad網站下載之ALL CHR檔案，從中擷取出PrediXcan所使用到的1,657,598(已去重複)SNPs
- 版本：gnomad.genomes.r.2.1.1.site.vcf.bgz
- 原始共有261,942,336 個SNPs資訊
- All mapped to the GRCh37/hg19 reference sequence. (PrediXcan variant位置資訊是b38版本，撈取時使用rsid，須額外注意可能會有重複，差別在ref/eff allele)
- 初步撈完為1,720,048 SNPs: allchr_result.csv
- 合併POS,ref/eff alleles資訊:allchr_gnomad_info.csv
- 確認哪些variant在EAS population沒有資訊：
- 計算有顯著族群差異位點

## C3_PopDiff_rsid_code: 合併整理gnomAD & PrediXcan rsid資訊，計算allele freq. 及 population test
- All mapped to the GRCh37/hg19 reference sequence.
- 欲得到PrediXcan共1,657,598位點資訊，最終得到1,656,652個(缺946個)
- 計算allele frequency部分，排除分母為0之48個位點(代表該位點沒收集到樣本)
- 最終有freq. info的有1,656,604個，其中有57個位點計算出來的值為0，直接列為不顯著
- 以p-value < 10*e-8為臨界值，#0:478545(28.9%) ; 1:1178059(71.1%)
- 以p-value < 10*e-8為臨界值 & 兩差值 > 0.05，#0:497862(30%) ; 1:1158742(70%)

## C4_WholeBlood_gene_simu
### 以gnomAD database的EAS population來模擬預測基因表現
- 分別模擬N=500, 1000, 2000的gene prediction是否有差異(使用Kolmogorov-Smirnov test)
- 多重檢定以Bonferroni做校正，500 v.s. 1000有0個有顯著差異，1000 v.s. 2000則是0個
- 結論：用500個進行模擬就夠惹

## C5_rsidDiff_effect_simu
### 以gnomAD database的EAS population來模擬預測全部49個tissues的基因表現量
- 合併全部的權重：8558894個權重

## C6_AllTissue_pred
### 以gnomAD database的EAS population來模擬預測全部49個tissues的基因表現量
- 製作packge會用的x1-49大表
- 模擬500個樣本(結論來自C4)

## C7_Rpackage_writing
### 準備所需檔案&R code撰寫
- 整理預測基因所需SNP個數, mean及sd, across all 49 tissues
            
## C8_extra_nfe_pred
### 以gnomAD database的NFE population來模擬預測全部49個tissues的基因表現量
- 與EAS相同，模擬500個樣本
- 檔案位置：Predixcan/統整/PrediXcan統整完檔案/dosage/nfe/results           
          
## C9_simu_revision
### 模擬修正

## C10_revision_lung_PrediXcan
### 應reviewer要求加上real data應用情形