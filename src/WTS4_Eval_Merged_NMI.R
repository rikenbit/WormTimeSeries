source("src/functions_WTS4_Eval_Merged_NMI.R")

#### args setting####
#### test args####

#### samplecode aricode####
BiocManager::install("aricode",update = FALSE)
library(aricode)
# input dataframe
# 5column 150row
data(iris)
# input クラスタリング結果
# cl <- cutree(hclust(dist(iris[,-5])), 4)
iris[,-5] %>% 
  dist() %>%  # データ行列to距離行列
  hclust() %>%  #クラスタリング tree作成
  cutree(., 4) -> cl #クラスタ数を設定し、各行にクラスタリングナンバーを割り当て
# 150細胞のクラスタリング結果

# input
# cl: 150細胞のクラスタリング結果
# cl: WTS4_Eval_behavior.Rのcluster
# classes: 各細胞のラベル <- 行動ラベル
# classes: WTS4_Eval_behavior.Rのclasses
NMI(cl, iris$Species)