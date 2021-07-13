#library
##################################################
library(mclust)
library(igraph)
library(tidyverse)
##################################################
.cls_ward = function(x) {
    # clustering
    cls_n <- x
    hclust(d, method = "ward.D2") %>% 
        cutree(., cls_n) %>% 
            as.vector() -> g_cls
    # prepare df
    return_object <- mutate(df_cell_type, cls = c(g_cls))
    return(return_object)
}
#### eval ARI####
.cls_ARI = function(x) {
    cls_n <- x
    first_cls_diff <- set_cutree[1] - 1
    i <- cls_n - first_cls_diff
    df_cls_n <- df_cls[[i]]
    # merge
    df_cls_label <- merge(df_cls_n, 
                          df_label, 
                          by.x = "cell_type", 
                          by.y = "cell_type", 
                          all.x = TRUE)
    clusters <- df_cls_label$cls
    classes <- df_cls_label$label_acf
    return_object <- adjustedRandIndex(clusters, classes)
    return(return_object)
}

.max_eval = function(x) {
  x %>%
    filter(., eval_value == max(eval_value)) %>%
    .$set_cutree  -> return_object
  return(return_object)
}

.min_eval = function(x) {
  x %>%
    filter(., eval_value == min(eval_value)) %>%
    .$set_cutree  -> return_object
  return(return_object)
}


#### eval purity####
ClusterPurity <- function(clusters, classes) {
    sum(apply(table(classes, clusters), 2, max)) / length(clusters)
}
.cls_purity = function(x) {
    cls_n <- x
    first_cls_diff <- set_cutree[1] - 1
    i <- cls_n - first_cls_diff
    df_cls_n <- df_cls[[i]]
    # merge
    df_cls_label <- merge(df_cls_n, 
                          df_label, 
                          by.x = "cell_type", 
                          by.y = "cell_type", 
                          all.x = TRUE)
    clusters <- df_cls_label$cls
    classes <- df_cls_label$label_acf
    return_object <- ClusterPurity(clusters, classes)
    return(return_object)
}

#### eval Fmeasure####
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
.cls_Fmeasure = function(x) {
    cls_n <- x
    first_cls_diff <- set_cutree[1] - 1
    i <- cls_n - first_cls_diff
    df_cls_n <- df_cls[[i]]
    # merge
    df_cls_label <- merge(df_cls_n, 
                          df_label, 
                          by.x = "cell_type", 
                          by.y = "cell_type", 
                          all.x = TRUE)
    clusters <- df_cls_label$cls
    classes <- df_cls_label$label_acf
    return_object <- Fmeasure(clusters, classes)
    return(return_object)
}

#### eval Entropy####
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
.cls_Entropy = function(x) {
    cls_n <- x
    first_cls_diff <- set_cutree[1] - 1
    i <- cls_n - first_cls_diff
    df_cls_n <- df_cls[[i]]
    # merge
    df_cls_label <- merge(df_cls_n, 
                          df_label, 
                          by.x = "cell_type", 
                          by.y = "cell_type", 
                          all.x = TRUE)
    clusters <- df_cls_label$cls
    classes <- df_cls_label$label_acf
    return_object <- Entropy(clusters, classes)
    return(return_object)
}