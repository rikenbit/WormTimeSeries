source("src/functions_WTS4_Noisy_Label_TF.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
# path NeuronActivity Data
args_ReadData <- args[1]
# sample number サンプル番号の指定
args_sample <- args[2]
# マージデータの細胞名取得
# Cluster
args_merged_cls <- args[3]
# 試行回数
args_test_number <- args[4]
# save ALL_Ann_df
args_output <- args[5]

# #### test args####
# # path NeuronActivity Data
# args_ReadData <- c("data/normalize_1/ReadData_2.RData")
# # sample number サンプル番号の指定
# args_sample <- c("2")
# # マージデータの細胞名取得
# # Cluster
# args_merged_cls <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/Merged_cls/k_Number_9.RData")
# # 試行回数
# args_test_number <- c("1")
# # save ALL_Ann_df
# args_output <- c("data/n1_noise_1/ReadData_2.RData")

#### load ReadData####
load(args_ReadData)
# 元データがディレクトリごとではなく，ファイル名で各サンプルがわかれている対応
eval(parse(text=paste0("ReadData <- ReadData_",args_sample)))

#### 数字列のデータフレーム####
grep("^[0-9]", colnames(ReadData), invert=FALSE) %>%
  ReadData[,.] -> not_Ann_df

#### 細胞名列のデータフレーム####
grep("^[0-9]", colnames(ReadData), invert=TRUE) %>% 
  ReadData[,.] -> Ann_df

#### 和集合の細胞名を取得####
load(args_merged_cls)
as.data.frame(merged_cls) %>% 
  rownames_to_column("cell_type") %>% 
  .$cell_type -> all_celltype

#### 使われていない細胞名を取得####
setdiff(all_celltype, 
        colnames(Ann_df)
) -> not_Used_celltype

#### ランダムに割り当てる値を抽出####
# 乱数固定
set.seed(args_test_number)
# 重複を許して、数字列から値（6000行）を抽出
sample(colnames(not_Ann_df), 
       length(not_Used_celltype),
       replace = TRUE
) -> Ramdom_celltype

#### insert random celltype####
Labeled_Ann_df <- not_Ann_df[, Ramdom_celltype]
colnames(Labeled_Ann_df) <- not_Used_celltype

#### bind data.frame####
cbind.data.frame(Ann_df, Labeled_Ann_df) -> ALL_Ann_df

#### save####
ReadData <- ALL_Ann_df
save(ReadData, file=args_output)

#### sample code####
# 赤の細胞達（n=20とする）
A <- matrix(runif(20*6000), nrow=20, ncol=6000)

# 緑（数字名の細胞達、n=30とする）にあてがうためのデータをランダムに生成
random.sampling <- function(A, n){
  # 各行(細胞)x各列のランダムサンプリング
  out <- unlist(
    # 各行(各細胞)だけ、割り当てる細胞数分実行
    lapply(seq(n), 
           function(x){
             # MARGIN = 2,各列（各TimeFrame）に対して関数適用。
             apply(A, 2, 
                   # replace=FALSEなので重複なしでランダムサンプリング。
                   # 緑（数字名の細胞達）の「ある1TimeFrame」の「ある細胞」のデータを取り出す
                   function(xx){sample(xx, 1)}
             )
           }
    )
  )
  # 数値ベクトルの形になっているので、行列の形式に当てはめる
  dim(out) <- c(n, length(out)/n)
  out
}

B <- random.sampling(A, n=30)