#library
##################################################
library(tidyverse)
library(clValid) # connectivityの計算
options(rgl.useNULL=TRUE) #https://scrapbox.io/Open-BioInfo-yamaken/TSclust（R）のインストールエラー
library(clusterSim) # Pseudo-Fの計算
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