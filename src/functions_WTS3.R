#library
##################################################
options(rgl.useNULL=TRUE)
library(TSclust)
library(tidyverse)
library(Rtsne)
library(ggrepel)
library(patchwork)
library(igraph)
library(dtwclust)
library(ggpubr)
library(umap)
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
# d_nn <- uwot:::dist_nn(d, k = lab_length)
# uwot_umap_d <- uwot::umap(d, 
#                           metric = "precomputed", 
#                           nn_method = list(idx = d_nn$idx, dist = d_nn$dist))
# plot(uwot_umap_d)
attr(d, "Labels") %>% length() -> lab_length
set.seed(1234)
uwot_umap_d <- uwot::umap(d, 
                          metric = "precomputed", 
                          nn_method = uwot:::dist_nn(d, k = lab_length),
                          n_neighbors = 15,
                          n_components = 2)
plot(uwot_umap_d)


wts_umap = function(x) {
    dist_data <- x
    set.seed(1234)

    dist_data %>% 
        as.matrix() %>%
            umap() -> umap_d
    #### umap####
    df_umap <- data.frame(cord_1 = umap_d$layout[,1],
                          cord_2 = umap_d$layout[,2],
                          cell_type = attr(dist_data, "Labels")
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

#### cal purity####
ClusterPurity <- function(clusters, classes) {
    sum(apply(table(classes, clusters), 2, max)) / length(clusters)
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

#### eval purity####
cls_purity = function(x) {
    cls_n <- x
    first_cls_diff <- cls_length[1] - 1
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

#### ggplot ward.D2 clustering####
gg_clusters = function(x) {
    cls_n <- x
    first_cls_diff <- cls_length[1] - 1
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