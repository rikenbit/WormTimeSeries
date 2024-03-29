#library
##################################################
library(Rtsne)
library(uwot)
library(tidyverse)
library(ggrepel)
library(igraph)
library(patchwork)
library(openxlsx) # read.xlsxを追加
# library(usedist)
library(ggpubr)
set.seed(1234)
##################################################

#### TSNE#### 
.wts_tsne = function(x) {
    dist_data <- x
    set.seed(1234)
    tSNE <- Rtsne(dist_data, 
                  is_distance = TRUE, 
                  dims = 2, 
                  perplexity = 15, 
                  verbose = TRUE, 
                  max_iter = 1000)
    return_object <- data.frame(cord_1 = tSNE$Y[,1],
                          cord_2 = tSNE$Y[,2],
                          cell_type = attr(dist_data, "Labels")
                          )
    return(return_object)
}
#### UMAP####
.wts_umap = function(x) {
    dist_data <- x
    attr(dist_data, "Labels") %>% length() -> lab_length
    set.seed(1234)
    umap_d <- uwot::umap(dist_data,
                         metric = "precomputed", 
                         nn_method = uwot:::dist_nn(dist_data, k = lab_length),
                         n_neighbors = 15,
                         n_components = 2
                         )
    return_object <- data.frame(cord_1 = umap_d[,1],
                          cord_2 = umap_d[,2],
                          cell_type = attr(dist_data, "Labels")
                          )
    return(return_object)
}

#### Dimensionality Reduction####
.df_cords = function(x) {
    d <- D[[x]]
    df_cord <- switch(args_DimReduc,
                      "tsne" = .wts_tsne(d),
                      "umap" = .wts_umap(d),
                      stop("Only can use tsne,")
                      )
    return_object <- df_cord
    return(return_object)
}
#### prepare dataframe & ggplot#####
.df_cord_cls = function(x) {
    df_cord <- DF_cord[[x]]
    cls_n <- C[[x]]
    df_cls <- as.data.frame(cls_n)
    df_cls <- data.frame(
        cell_type = rownames(df_cls),
        Cluster = df_cls$cls_n,
        stringsAsFactors = FALSE
    )
    # アノテーション付きの細胞データにクラスタ情報が付与
    df_cord_cls <- merge(df_cord, 
                         df_cls, 
                         by.x = "cell_type", 
                         by.y = "cell_type", 
                         all.x = TRUE)
    return_object <- df_cord_cls
    return(return_object)
}
.df_cord_cls_nl = function(x) {
    df_cord_cls <- DF_cord_cls[[x]]
    df_cord_cls_nl <- merge(df_cord_cls, 
                            df_nl, 
                            by.x = "cell_type", 
                            by.y = "cell_type", 
                            all.x = TRUE)
    return_object <- df_cord_cls_nl
    return(return_object)
}

.df_cord_cls_nl_be = function(x) {
    df_cord_cls_nl <- DF_cord_cls_nl[[x]]
    df_cord_cls_nl_be <- merge(df_cord_cls_nl, 
                               df_eval_label, 
                               by.x = "cell_type", 
                               by.y = "cell_type", 
                               all.x = TRUE)
    return_object <- df_cord_cls_nl_be
    return(return_object)
}

.df_cord_cls_nl_be_count = function(x) {
    df_cord_cls_nl_be <- DF_cord_cls_nl_be[[x]]
    df_cord_cls_nl_be_count <- merge(df_cord_cls_nl_be, 
                               df_cell_count, 
                               by.x = "cell_type", 
                               by.y = "CellType", 
                               all.x = TRUE)
    # cell countのNAは0に置換
    df_cord_cls_nl_be_count$CellCount[is.na(df_cord_cls_nl_be_count$CellCount)] <- 0
    return_object <- df_cord_cls_nl_be_count
    return(return_object)
}

.plot_dimreduc = function(x) {
    sample_n <- df_weight[x, 2]
    sample_id <- as.numeric(df_weight[x, 1])
    df_merged <- DF_cord_cls_nl_be_count[[sample_id]]
    #### ggplot Cluster####
    if (args_DimReduc == "tsne") {
      cord_x <- c("t-SNE-1")
      cord_y <- c("t-SNE-2")
    } else if (args_DimReduc == "umap") {
      cord_x <- c("UMAP-1")
      cord_y <- c("UMAP-2")
    } else {
      cord_x <-  c("cord-1")
      cord_y <- c("cord-2")
    }
    gg_cls <- ggplot(df_merged, 
                     aes(x = cord_1,
                         y = cord_2, 
                         label = cell_type,
                         color = factor(Cluster)
                         )
                     ) + 
        labs(color = "Cluster") +
        theme(text = element_text(size = 24)) +
        geom_point(size = 6.0, 
                   alpha = 0.6) +
        geom_label_repel(max.overlaps = Inf,
                         min.segment.length = 0,
                         size = 7.0,
                         force = 6.0) + # ラベル間の反発力
        theme(text = element_text(size = 60)) +
        labs(x = cord_x, y = cord_y)
    #### ggplot NeuronType####
    gg_NT <- ggplot(df_merged, 
                    aes(x = cord_1,
                        y = cord_2, 
                        label = cell_type,
                        color = factor(NeuronType)
                    )
    ) + 
        labs(color = "Neuron type") +
        theme(text = element_text(size = 24)) +
        geom_point(size = 6.0, 
                   alpha = 0.6) +
        geom_label_repel(max.overlaps = Inf,
                         min.segment.length = 0,
                         size = 7.0,
                         force = 6.0) + # ラベル間の反発力
        theme(text = element_text(size = 60)) +
        labs(x = cord_x, y = cord_y) 
    #### ggplot eval_label####
    gg_eval_label <- ggplot(df_merged, 
                            aes(x = cord_1,
                                y = cord_2, 
                                label = cell_type,
                                color = factor(Classes)
                            )
    ) + 
        labs(color = "Class") +
        theme(text = element_text(size = 24)) +
        geom_point(size = 6.0, 
                   alpha = 0.6) +
        geom_label_repel(max.overlaps = Inf,
                         min.segment.length = 0,
                         size = 7.0,
                         force = 6.0) + # ラベル間の反発力
        theme(text = element_text(size = 60)) +
        labs(x = cord_x, y = cord_y) 
    #### ggplot cell_count####
    gg_cell_count <- ggplot(df_merged, 
                            aes(x = cord_1,
                                y = cord_2, 
                                label = cell_type,
                                color = CellCount
                                )
                            ) + 
      scale_color_viridis_c(option = "D") +
      labs(color = "No. of cells") +
      geom_point(size = 6.0, 
                 alpha = 0.6) +
      geom_label_repel(max.overlaps = Inf,
                       min.segment.length = 0,
                       size = 9.0,
                       force = 6.0) +# ラベル間の反発力
      theme(text = element_text(size = 60)) +
      labs(x = cord_x,
           y = cord_y) +
      theme(legend.key.height = unit(1.5, "cm")) +
      theme(legend.key.width = unit(1.5, "cm"))
    #### patchwork####
    eval(parse(text=paste0("plot_title <- c('",x,"_SampleNumber_",sample_n,"')")))
    gg <- gg_cls +
        gg_NT +
        gg_eval_label +
        gg_cell_count +
        plot_layout(nrow = 1) +
        plot_annotation(
            title = plot_title,
            caption = 'made with patchwork',
            theme = theme(plot.title = element_text(size = 60, hjust = 0.5))
        )
    return_object <- gg
    return(return_object)
}