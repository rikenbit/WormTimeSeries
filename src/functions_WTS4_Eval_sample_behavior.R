#library
##################################################
library(tidyverse)
library(mclust) # ARI
library(openxlsx) # read.xlsx
library(aricode) # NMI
##################################################

####################################################################
############# 外的評価（クラスラベルとどれだけ一致しているか） #############
####################################################################
#### NMI#### 
.NMI_list = function(x) {
    df_cls_label <- x
    unlist(lapply(df_cls_label, function(x) {
        clusters <- x$Clusters
        classes <- x$Classes
        NMI(clusters, classes)
        }
    )) -> return_object
    return(return_object)
}

####  ARI#### 
.ARI_list = function(x) {
	df_cls_label <- x
	unlist(lapply(df_cls_label, function(x) {
    clusters <- x$Clusters
    classes <- x$Classes
    adjustedRandIndex(clusters, classes)
    }
    )) -> return_object
	return(return_object)
}

####  purity#### 
ClusterPurity <- function(clusters, classes) {
  sum(apply(table(classes, clusters), 2, max)) / length(clusters)
}
.purity_list = function(x) {
	df_cls_label <- x
	unlist(lapply(df_cls_label, function(x) {
    clusters <- x$Clusters
    classes <- x$Classes
    ClusterPurity(clusters, classes)
    }
    )) -> return_object
	return(return_object)
}

#### Fmeasure####
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
.Fmeasure_list = function(x) {
	df_cls_label <- x
	unlist(lapply(df_cls_label, function(x) {
    clusters <- x$Clusters
    classes <- x$Classes
    Fmeasure(clusters, classes)
    }
    )) -> return_object
	return(return_object)
}

#### Entropy####
calcEntropy0 <- function(pv){
    p1 <- pv / sum(pv)
    p2 <- p1[p1 !=0]
    - sum(p2 * log2(p2))
}
# 全体のエントロピー（クラスタのデータ数で重み付け平均をとる）
Entropy <- function(cluster, label){
    # クロス集計
    ctbl <- table(cluster, label)
    # 重み
    w <- apply(ctbl, 1, sum) / sum(ctbl)
    # 全体のエントロピー
    sum(w * apply(ctbl, 1, calcEntropy0))
}

.Entropy_list = function(x) {
	df_cls_label <- x
	unlist(lapply(df_cls_label, function(x) {
    clusters <- x$Clusters
    classes <- x$Classes
    Entropy(clusters, classes)
    }
    )) -> return_object
	return(return_object)
}