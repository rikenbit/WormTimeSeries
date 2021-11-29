#library
##################################################
library(tidyverse)
# ARI
library(mclust)
##################################################

####################################################################
############# 外的評価（クラスラベルとどれだけ一致しているか） #############
####################################################################
# purity
ClusterPurity <- function(clusters, classes) {
  sum(apply(table(classes, clusters), 2, max)) / length(clusters)
}
# Fmeasure
Fmeasure <- function(cluster, label){
  # クロス集計
  ctbl <- table(cluster, label)
  # 相当たりRecall
  R <- ctbl / colSums(ctbl)
  # 相当たりPrecision
  P <- ctbl / rowSums(ctbl)
  # 相当たりF-measure
  F <- 2 * R * P / (R + P)
  # 0補正
  F[is.nan(F)] <- 0
  # 重み
  w <- apply(ctbl, 2, sum) / sum(ctbl)
  # 全体のF-measure
  sum(w * apply(F, 1, max))
}