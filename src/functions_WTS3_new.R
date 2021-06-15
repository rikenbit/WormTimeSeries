#library
##################################################
library(tidyverse)
library(Rtsne)
library(ggrepel)
library(patchwork)
library(igraph)
library(ggpubr)
library(uwot)
library(mclust)
# mclustだけステージ版
##################################################
set.seed(1234)

#### TSNE#### 
wts_tsne = function(x) {
    dist_data <- x
    set.seed(1234)
    tSNE <- Rtsne(dist_data, 
                  is_distance = TRUE, 
                  dims = 2, 
                  perplexity = 15, 
                  verbose = TRUE, 
                  max_iter = 1000)
    df_tSNE <- data.frame(cord_1 = tSNE$Y[,1],
                          cord_2 = tSNE$Y[,2],
                          cell_type = attr(dist_data, "Labels")
                          )
    return(df_tSNE)
}

#### UMAP####
wts_umap = function(x) {
    d <- x
    attr(d, "Labels") %>% length() -> lab_length
    set.seed(1234)
    umap_d <- uwot::umap(d,
                         metric = "precomputed", 
                         nn_method = uwot:::dist_nn(d, k = lab_length),
                         n_neighbors = 15,
                         n_components = 2
                         )
    df_umap <- data.frame(cord_1 = umap_d[,1],
                          cord_2 = umap_d[,2],
                          cell_type = attr(d, "Labels")
                          )
    return(df_umap)
}

#### ggplot Neuron type/group####
gg_n = function(x) {
    # ggplot
    g_col <- x
    gg_n <- eval(parse(text=paste0("ggplot(df_merged, 
                                    aes(x = cord_1,
                                        y = cord_2, 
                                        label = cell_type, 
                                        color = ",factor(g_col),"))"))) +
            geom_point() +
            geom_text_repel(max.overlaps = Inf,
                            min.segment.length = 0) +
            ggtitle(g_col) +
            theme(plot.title = element_text(size = 30, hjust = 0.5))
    return(gg_n)
}

#### ggplot Neuron type/group StimShape####
gg_n_stim = function(x) {
  # ggplot
  g_col <- x
  gg_n <- eval(parse(text=paste0("ggplot(df_merged_stim0, 
                                    aes(x = cord_1,
                                        y = cord_2, 
                                        label = cell_type,
                                        shape = factor(stim),
                                        color = ",factor(g_col),"))"))) +
    geom_point(size = 8.0) +
    scale_shape_manual(values=c(4,15)) +
    geom_text_repel(max.overlaps = Inf,
                    min.segment.length = 0) +
    ggtitle(g_col) +
    theme(plot.title = element_text(size = 30, hjust = 0.5))
  return(gg_n)
}

#### dataframe for evaluation ####
cls_cord = function(x) {
    #### clustering#### 
    cls_n <- x
    hclust(d, method = "ward.D2") %>% 
        cutree(., cls_n) %>% 
        as.vector() -> g_cls
    #### prepare df####
    df_cls <- mutate(df_cord, cls = c(g_cls))
    return(df_cls)
}

#### ggplot ward.D2 clustering####
gg_clusters = function(x) {
    cls_n <- x
    first_cls_diff <- set_cutree[1] - 1
    i <- cls_n - first_cls_diff
    df_cls <- df_cls_cord[[i]]
    #### ggplot#### 
    gg_cls_n <- ggplot(df_cls, 
                       aes(x = cord_1, 
                           y = cord_2, 
                           label = cell_type, 
                           color = factor(cls))) +
                            # color = forcats::fct_explicit_na(factor(cls)))) +
                geom_point() +
                geom_text_repel(max.overlaps = Inf, 
                                min.segment.length = 0)
    eval(parse(text=paste0("title <- ggtitle('cutree_",cls_n,"')")))
    t_1 <- theme(plot.title = element_text(size = 30, hjust = 0.5))
    gg_cls_n <- gg_cls_n +
        title + 
        t_1
    return(gg_cls_n)
}

#### eval purity####
ClusterPurity <- function(clusters, classes) {
    sum(apply(table(classes, clusters), 2, max)) / length(clusters)
}
cls_purity = function(x) {
    cls_n <- x
    first_cls_diff <- set_cutree[1] - 1
    i <- cls_n - first_cls_diff
    df_plot <- df_cls_cord[[i]]
    # merge
    df_plot_stim <- merge(df_plot, 
                       stim_sheet, 
                       by.x = "cell_type", 
                       by.y = "cell_type", 
                       all.x = TRUE)
    # stim列のNAを0に変換する
    df_plot_stim %>% 
        replace_na(., replace = list(stim = 0)) -> df_stim
    classes <- df_stim$stim
    clusters <- df_stim$cls
    ClusterP <- ClusterPurity(clusters, classes)
    return(ClusterP)
}

#### eval ARI####
cls_ARI = function(x) {
    cls_n <- x
    first_cls_diff <- set_cutree[1] - 1
    i <- cls_n - first_cls_diff
    df_plot <- df_cls_cord[[i]]
    # merge
    df_plot_stim <- merge(df_plot, 
                          stim_sheet, 
                          by.x = "cell_type", 
                          by.y = "cell_type", 
                          all.x = TRUE)
    # stim列のNAを0に変換する
    df_plot_stim %>% 
        replace_na(., replace = list(stim = 0)) -> df_stim
    classes <- df_stim$stim
    clusters <- df_stim$cls
    ClusterARI <- adjustedRandIndex(clusters, classes)
    return(ClusterARI)
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
cls_Fmeasure = function(x) {
    cls_n <- x
    first_cls_diff <- set_cutree[1] - 1
    i <- cls_n - first_cls_diff
    df_plot <- df_cls_cord[[i]]
    # merge
    df_plot_stim <- merge(df_plot, 
                          stim_sheet, 
                          by.x = "cell_type", 
                          by.y = "cell_type", 
                          all.x = TRUE)
    # stim列のNAを0に変換する
    df_plot_stim %>% 
        replace_na(., replace = list(stim = 0)) -> df_stim
    classes <- df_stim$stim
    clusters <- df_stim$cls
    ClusterFmeasure <- Fmeasure(clusters, classes)
    return(ClusterFmeasure)
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
cls_Entropy = function(x) {
    cls_n <- x
    first_cls_diff <- set_cutree[1] - 1
    i <- cls_n - first_cls_diff
    df_plot <- df_cls_cord[[i]]
    # merge
    df_plot_stim <- merge(df_plot, 
                          stim_sheet, 
                          by.x = "cell_type", 
                          by.y = "cell_type", 
                          all.x = TRUE)
    # stim列のNAを0に変換する
    df_plot_stim %>% 
        replace_na(., replace = list(stim = 0)) -> df_stim
    classes <- df_stim$stim
    clusters <- df_stim$cls
    ClusterEntropy <- Entropy(clusters, classes)
    return(ClusterEntropy)
}
#### eval Max####
max_eval = function(x) {
    x %>%
        filter(., eval_value == max(eval_value)) %>%
            .$set_cutree %>%
                purrr::map(., gg_clusters) -> gg_cls
    return(gg_cls)
}
#### eval Min####
min_eval = function(x) {
    x %>%
        filter(., eval_value == min(eval_value)) %>%
            .$set_cutree %>%
                purrr::map(., gg_clusters) -> gg_cls
    return(gg_cls)
}