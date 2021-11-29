#library
##################################################
library(tidyverse)
# ARI
library(mclust)
##################################################

####################################################################
######### 内的評価（クラスタリング自体がどの程度うまくいっているか、 #########
############### クラスタ内の距離、クラスタ間の距離を利用）#################
####################################################################
# library(clValid) # connectivityの計算
# options(rgl.useNULL=TRUE) #https://scrapbox.io/Open-BioInfo-yamaken/TSclust（R）のインストールエラー
# library(clusterSim) # Pseudo-Fの計算

# # Pseudo-F値（大きいほどよい）
# .PseudoF <- function(data, cluster){
#     index.G1(data, cluster)
# }
# 
# # 結合度（kNNの失敗に対するペナルティの総和、小さいほどよい）
# .Connectivity <- function(data, cluster){
#     # 距離行列
#     Dist <- dist(data, method="euclidean")
#     # 結合度
#     connectivity(Dist, cluster)
# }

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