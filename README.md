# PrediXcan 參考族群 → PredictMP

用 gnomAD v4.1 各族群的 allele frequency 模擬 500 人的基因型，丟進 PrediXcan（GTEx v8 elastic-net，49 個組織）得到各族群的「參考預測表現量分佈」。
再打包成 `PredictMP` 套件：使用者的 PrediXcan 結果可以對參考分佈算 percentile rank，並附上 KS test（vs NFE）結果。

PredictMP 是 PredictAP（Chan et al. 2024, Brief Bioinform 25(6) bbae549，只有 EAS 參考族群）的多族群版本。

## 目錄

```
Project/
├── elastic_net_models/          en_{tissue}.db（49 個組織的權重模型）
└── PredictMP/
    ├── data/                    所有中間檔與套件資料
    ├── results/{POP}/           PrediXcan 預測結果（步驟 07 產生）
    └── code/
        ├── config.py            ← 所有路徑與族群都在這裡（參數在各程式開頭）
        ├── 01 … 09              主流程（照編號跑）
        ├── PrediXcan.py         PrediXcan 原始工具（07 會呼叫）
        ├── PredictMP_pkg/    PredictMP Python 套件
        ├── analysis/            側分析（不在主流程上）
        ├── archive/             舊版程式，保留備查
        └── Senior_R_code/       學長姐的 R 版本
```

## 設定：`config.py`

`config.py` 只放**路徑和族群**，所有程式都 `from config import ...`，**不要在個別程式裡寫死路徑**。

- 根目錄預設是 `code/` 往上兩層（`/mnt/data_40T/betty0101/Project`），可用環境變數 `PREDIXCAN_ROOT` 覆蓋。
- 族群：`POPULATIONS = ["afr", "ami", "amr", "asj", "eas", "fin", "mid", "nfe", "sas"]`，參考族群 `REF_POP = "nfe"`。
- `python config.py` 會印出目前所有設定。

分析參數寫在各程式開頭，直接改：

| 參數 | 程式 | 預設 |
|---|---|---|
| `PVALUE_THRESH`、`FREQ_DIFF_THRESH` | 05 | `1e-8`、`0.05` |
| `N_SAMPLES` | 06 **和** 07（兩邊要一樣，07 開跑前會檢查 sample file 人數） | `500` |
| `RANDOM_SEED` | 07 | `123` |
| `KS_THRESHOLD` | 08 | `5e-8` |

`02_gnomad_process.sh` 是 shell，不能 import config，所以同樣的路徑寫在檔案開頭的變數裡；改路徑時兩邊都要改。

## 環境

```bash
conda activate p_env                      # Python 3.14: pandas, numpy, scipy, pyarrow, tqdm, matplotlib, seaborn
pip install -e PredictMP_pkg           # 搬過目錄後要重裝一次
conda activate bcf_env                    # 只有步驟 02 需要（bcftools / bgzip / tabix）
```

## 主流程

| 步驟 | 指令 | 輸入 → 輸出（在 `data/`，除非另外註明） |
|---|---|---|
| 01 | `python 01_extract_weight.py` | `elastic_net_models/en_*.db` → `weight.csv`、`Num_of_genes_tissue.csv` |
| 02 | `bash 02_gnomad_process.sh` | dbSNP `~/GCF_000001405.40` + gnomAD v4.1 joint VCF → `gnomad_joint_v4.1_info/`、`gnomad_joint_v4.1_info_with_rsid/` |
| 03 | `python 03_gnomad_rsid_summary.py` | （選用 QC）各 chr 有 rsid 的比例、與 gnomAD v2 比較 |
| 04 | `python 04_convert_csv_to_parquet.py` | `…_with_rsid/*.csv` → `…_with_rsid_parquet/*.parquet` |
| 05 | `python 05_popdiff_rsid.py` | `weight.csv` + parquet → `gnomad_rsid_info_with_popdiff.csv`、`gnomad_rsid_withoutinfo.csv` |
| 06 | `python 06_generate_sample_file.py [--pop eas …]` | → `sample_file{N}_{pop}.txt`（預設 9 個族群全做，N 預設 500） |
| 07 | `python 07_simulate_and_predict.py --pop eas` 或 `--pop all` | → `data/dosage/dosage_{POP}/{tissue}/`、`results/{POP}/{tissue}_predicted_expression.txt` |
| 08 | `python 08_build_package_data.py [--step 1 2 3 4]` | `results/*/` → `tissue/{pop}_x{i}.parquet`、`result_table/result_table_{pop}.parquet`、`ks/ks_{pop}_nfe.parquet` |
| 09 | `python 09_run_predictmp.py --input … --tissue x32 --population eas --output …` | 使用者的 PrediXcan 結果 → percentile rank 表 |

每一步都可以單獨重跑；02、04、07 已存在的輸出會自動跳過。

### 02 gnomAD 處理細節

| 子步驟 | 做什麼 | 為什麼 |
|---|---|---|
| a | dbSNP contig `NC_0000xx` → `chr1..22`，只留常染色體 | gnomAD 用 `chrN` 命名 |
| b | `bcftools annotate -h` 補 `##contig` header | 改名後 header 裡的 contig 沒跟著改 |
| c | `bcftools norm -m -any` + tabix | dbSNP 一行多個 ALT，gnomAD 是一行一個 |
| d | `bcftools annotate -a dbSNP -c ID`（每條 chr） | 用 CHROM/POS/REF/ALT 把 rsid 填進 gnomAD |
| e | `bcftools query` 取 `AC/AN/AF_joint_{pop}` | 產生 `gnomad_chrN_populations.csv` |
| f | 留 `ID != "."` | 產生 `gnomad_chrN_populations_rsid.csv` |

a–d 的指令出自 2026-05 的筆記和 VCF header 裡的 `##bcftools_annotateCommand`。**e 當時沒有留下紀錄**，是依照現有 CSV 的欄位順序重建的；如果要重跑，請先拿一條 chr 的結果和現有檔案比對行數。

### 05 族群差異

PrediXcan 的 SNP 和 gnomAD 對 allele 時分四輪：直接對上 → ref/eff 對調（AF = 1−AF）→ 互補股 → 互補股加對調。
接著對每個族群都和 NFE 做 two-proportion Z-test。`diff_index_{pop} = 1` 的條件是 p < 1e-8 **且** |Δfreq| > 0.05。

> 門檻說明：Chan et al. 2024 論文寫的是 5×10⁻⁸，但學姊的原始程式（`Senior_R_code/C3_PopDiff_rsid.r` 第 123 行）實際用 `10^-8`，論文的「約 70% SNP 有族群差異」也是用 10⁻⁸ 算的。本流程沿用 10⁻⁸。`diff_index` 只用於統計，不影響 07 的模擬和網站的參考資料。

### 08 KS test

每個族群對 NFE、每個組織每個基因做兩樣本 KS test（scipy `ks_2samp`），p < 5×10⁻⁸ 為 `Sig.`，與論文及學姊套件（`result_table_v2.csv`）一致。
兩族群模擬值完全相同時 p = 1，判為 `Nonsig.`（學姊的 R 程式把 p 為 NA 的情況算成 `Sig.`，論文沒有說明；新資料只有 2 個組織×基因組合屬於這種情況）。
論文 Table S3 另外用「p < 1 ÷ 該組織基因數」做了補充分析（學姊的 `kstest_sign_num.csv`），沒有放進套件，本流程也沒有實作。

### 07 模擬

在 HWE 假設下，以 gnomAD 的 `freq_{pop}` 抽 0/1/2 dosage，每個族群 500 人，然後用 `PrediXcan.py --predict` 跑完 49 個組織。
`--pop all` 會依序跑 9 個族群，每個族群各開一個子程序（族群透過環境變數 `PREDIXCAN_POP` 傳給 worker）。

亂數：每個（族群, 組織）用自己的 seed `[RANDOM_SEED, 族群編號, 組織編號]`，所以組織之間、族群之間的模擬互相獨立，而且不管哪個 worker 跑都能重現。（學長姐的 R 版沒有設 seed；2026-10 之前的 Python 版每個 worker 都從 `seed(123)` 開始，所有組織和族群共用同一串亂數，舊結果放在 `PredictMP/old_seed123/`。）

## PredictMP 套件

```python
import pandas as pd
from PredictMP import rankcal

data = pd.read_csv("Lung_predicted_expression.txt", sep="\t")   # FID, IID, ENSG...
res  = rankcal(data, tissue="x32", population="eas")            # x32 = Lung
```

輸出欄位：`FID, IID, gene, expression, percentile.rank (0–100), N.snps, mean, sd`。
`population` 不是 `nfe` 時，會再多兩欄：`ks.test`（Sig./Nonsig.，門檻 5e-8）和 `NFE.{POP}.mean.diff`（%）。

- 參考資料從 `PredictMP/data/{tissue,result_table,ks}/` 讀取，可用 `PREDICTMP_DATA_DIR` 改位置。
- 組織代碼是 `x1`–`x49`，依 GTEx 組織名稱字母排序（`x1`=Adipose_Subcutaneous，`x32`=Lung，`x49`=Whole_Blood），完整對照見 `09_run_predictmp.py` 的 `TISSUE_INDEX`。
- `PredictMP_pkg/PredictMP/data/` 裡有 1.3 GB 舊格式的 NFE parquet，loader 已經不讀它，打包時也已排除。確認不需要後可以刪掉。

### 在其他主機使用

套件只有程式碼，參考資料（約 12 GB）要另外下載，再用 `PREDICTMP_DATA_DIR` 指過去。

```bash
# 1. 從 GitHub 安裝（private repo 用 ssh 網址：git+ssh://git@github.com/Coco-0101/PredictMP.git#subdirectory=PredictMP_pkg）
pip install "git+https://github.com/Coco-0101/PredictMP.git#subdirectory=PredictMP_pkg"
#    固定版本：在 .git 後面加 @<tag 或 commit>，例如 PredictMP.git@v1.1.0#subdirectory=PredictMP_pkg

# 2. 下載參考資料（rclone 設定見 PredictMP_GeneExprPredict-Web 的 README）
rclone copy gdrive:PredictMP_GeneExprPredict-Web/reference_data ~/predictmp_reference -P

# 3. 告訴套件資料在哪（寫進 ~/.bashrc 就不用每次設定）
export PREDICTMP_DATA_DIR=~/predictmp_reference

python -c "from PredictMP import rankcal; print('ok')"
```

更新：`pip install --force-reinstall --no-deps "git+https://..."`。

## 側分析（`analysis/`）

| 程式 | 內容 |
|---|---|
| `sample_size_check.py` | Whole_Blood 模擬 N=500/1000/2000 的分佈是否不同（結論：500 就夠）。讀的是舊的 gnomAD v2 檔 `predixcan_gnomad_rsid_freq.csv`，族群由檔案開頭的 `POP` 指定；模擬的 dosage 放在 `data/sample_size_check/N500/` |
| `weight_effect_simulation.py` | 依 weight 大小、MAF 分組，看族群差異 SNP 對預測值的影響 |

## 已知問題

- `data/sample_file500.txt`（NTUH/TWB/EAS 標籤）是學長姐的舊模板，舊版除了 AFR 以外的族群都用它；現在每個族群用自己的 `sample_file500_{pop}.txt`（ID 如 `NFE1`，cohort 標 `SIM`）。
- `archive/C6_AllTissue_pred_python.py` 是寫死 AFR 的舊版 07，只留著備查。
