#library
##################################################
library(tidyverse)
library(clValid) # connectivityの計算
options(rgl.useNULL=TRUE) #https://scrapbox.io/Open-BioInfo-yamaken/TSclust（R）のインストールエラー
library(clusterSim) # Pseudo-Fの計算
library(FastKNN) # kNNの計算
library(openxlsx) # stimAfterのexcel読み込み
##################################################

####################################################################
######### 内的評価（クラスタリング自体がどの程度うまくいっているか、 #########
############### クラスタ内の距離、クラスタ間の距離を利用）#################
####################################################################

# Pseudo-F値（大きいほどよい）
# https://www.rdocumentation.org/packages/clusterSim/versions/0.49-2/topics/index.G1
.PseudoF <- function(data, cluster){
    index.G1(data, cluster)
}
.PseudoF_list = function(x) {
    data <- sample_data_list[[x]]
    cluster <- df_cls_list[[x]]$Clusters
    return_object <- .PseudoF(data, cluster)
    return(return_object)
}

# 結合度（kNNの失敗に対するペナルティの総和、小さいほどよい）
# https://rdrr.io/cran/clValid/man/connectivity.html
.Connectivity <- function(data, cluster){
    # 距離行列
    Dist <- dist(data, method="euclidean")
    # 結合度
    connectivity(Dist, cluster)
}
.Connectivity_list = function(x) {
    data <- sample_data_list[[x]]
    cluster <- df_cls_list[[x]]$Clusters
    return_object <- .Connectivity(data, cluster)
    return(return_object)
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
  #### fix クラスタ値の欠損####
  # 全てのサンプルを使うわけでないので、数字の細胞の削除後に使われるクラスタ数が減るケースがある
  # unique(cluster)だと、使わないクラスターがある場合は処理がうまくいかない
  # H <- matrix(0, nrow=nrow(data), ncol=length(unique(cluster)))
  H <- matrix(0, nrow=nrow(data), ncol=length(1:max(cluster)))
  ############################
  for(i in seq_len(nrow(data))){
    H[i, cluster[i]] <- 1
  }
  sum(H %*% t(H) * A) / (nrow(data) * k)
}

.kNN_list = function(x) {
    data <- sample_data_list[[x]]
    cluster <- df_cls_list[[x]]$Clusters
    k <- k
    return_object <- .kNN(data, cluster, k)
    return(return_object)
}

######### 時系列方向のtrim#########
.ReadData_stimAfter = function(x,y) {
    #### load stim timing extra####
    stimtimng_sheet <- read.xlsx(y,
                                 sheet = "Sheet1",
                                 rowNames = FALSE,
                                 colNames =TRUE)
    stimtimng_sheet %>% 
        dplyr::select(sample_number = 1, 
                      stim_first = 7,
                      ) %>% 
            filter(sample_number == args_sample) %>% 
                .$stim_first %>% 
                    trunc() -> stimtimng #切り捨て
                    # ceiling() -> stimtimng #切り上げ
    return_object <- x[stimtimng:nrow(x),]
    return(return_object)
}

.ReadData_all = function(x) {
    return_object <- x
    return(return_object)
}