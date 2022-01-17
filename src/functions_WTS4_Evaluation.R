#library
##################################################
library(tidyverse)
library(clValid) # connectivityの計算
options(rgl.useNULL=TRUE) #https://scrapbox.io/Open-BioInfo-yamaken/TSclust（R）のインストールエラー
library(clusterSim) # Pseudo-Fの計算
library(FastKNN) # kNNの計算
##################################################

####################################################################
######### 内的評価（クラスタリング自体がどの程度うまくいっているか、 #########
############### クラスタ内の距離、クラスタ間の距離を利用）#################
####################################################################

# Pseudo-F値（大きいほどよい）
.PseudoF <- function(data, cluster){
    index.G1(data, cluster)
}

# 結合度（kNNの失敗に対するペナルティの総和、小さいほどよい）
.Connectivity <- function(data, cluster){
    # 距離行列
    Dist <- dist(data, method="euclidean")
    # 結合度
    connectivity(Dist, cluster)
}

# kNN（0〜1の値、大きいほどクラスタリングがうまくいっている）
# https://datachemeng.com/post-4657/
# .kNN <- function(data, cluster, k=1){
.kNN <- function(data, cluster, k){
  stopifnot(k <= nrow(data)-1)
  # kNNによる隣接行列
  dist_mat <- as.matrix(dist(data, method = "euclidean",
                             upper = TRUE, diag=TRUE))
  nrst <- lapply(1:nrow(dist_mat), function(i)
    k.nearest.neighbors(i, dist_mat, k = k))
  A <- matrix(0, nrow=nrow(data), ncol=nrow(data))
  for(i in seq_len(nrow(data))){
    A[i, nrst[[i]]] <- 1
  }
  # Cluster Label → Indicator Matrix
  H <- matrix(0, nrow=nrow(data), ncol=length(unique(cluster)))
  for(i in seq_len(nrow(data))){
    H[i, cluster[i]] <- 1
  }
  sum(H %*% t(H) * A) / (nrow(data) * k)
}