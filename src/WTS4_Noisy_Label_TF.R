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

#### 時系列random.sampling####
not_Ann_df %>%
  as.matrix() %>%
  t() -> A
# 乱数固定
set.seed(args_test_number)
B <- random.sampling(A, n=length(not_Used_celltype))
# random.sampling関数が行列の想定が逆なのでt()で転値
B %>%
  t() %>%
  as.data.frame() -> Labeled_TF
colnames(Labeled_TF) <- colnames(Labeled_Ann_df)

#### bind data.frame####
# cbind.data.frame(Ann_df, Labeled_Ann_df) -> ALL_Ann_df
cbind.data.frame(Ann_df, Labeled_TF) -> ALL_Ann_df

#### save####
ReadData <- ALL_Ann_df
save(ReadData, file=args_output)