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
##################################################
set.seed(1234)

gg_clsters = function(x) {
    #### clustering#### 
    cls_n <- x
    hclust(d, method = "ward.D2") %>% 
        cutree(., cls_n) %>% 
        as.vector() -> g_cls
    #### prepare df#### 
    df_tSNE <- data.frame(tsne_1 = tSNE$Y[,1],
                          tsne_2 = tSNE$Y[,2],
                          cell_type = attr(d, "Labels"),
                          cls = c(g_cls)
    )
    #### ggplot#### 
    gg_cls_n <- ggplot(df_tSNE, 
                       aes(x = tsne_1, 
                           y = tsne_2, 
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

gg_n = function(x) {
    # ggplot
    g_col <- x
    # gg_n <- eval(parse(text=paste0("ggplot(df_merged, aes(x = tsne_1, y = tsne_2, label = cell_type, color = forcats::fct_explicit_na(",factor(g_col),")))"))) +
    gg_n <- eval(parse(text=paste0("ggplot(df_merged, 
                                    aes(x = tsne_1,
                                        y = tsne_2, 
                                        label = cell_type, 
                                        color = ",factor(g_col),"))"))) +
            geom_point() +
            geom_text_repel(max.overlaps = Inf,
                            min.segment.length = 0) +
            ggtitle(g_col) +
            theme(plot.title = element_text(size = 30, hjust = 0.5))
    return(gg_n)
}

cls_purity = function(x) {
    #### clustering#### 
    cls_n <- x
    hclust(d, method = "ward.D2") %>% 
        cutree(., cls_n) %>% 
        as.vector() -> g_cls
    #### prepare df#### 
    df_tSNE <- data.frame(tsne_1 = tSNE$Y[,1],
                          tsne_2 = tSNE$Y[,2],
                          cell_type = attr(d, "Labels"),
                          cls = c(g_cls)
                        )
    #### merge df####
    # merge
    df_merged <- merge(df_tSNE, 
                       stim_sheet, 
                       by.x = "cell_type", 
                       by.y = "cell_type", 
                       all.x = TRUE)
    # stim列のNAを0に変換する
    df_merged %>% 
        replace_na(., replace = list(stim = 0)) -> df_tSNE_stim
    #### cal  purity####
    ClusterPurity <- function(clusters, classes) {
        sum(apply(table(classes, clusters), 2, max)) / length(clusters)
    }
    
    classes <- df_tSNE_stim$stim
    clusters <- df_tSNE_stim$cls
    ClusterP <- ClusterPurity(clusters, classes)
    return(ClusterP)
}