#library
##################################################
library(Rtsne)
library(uwot)
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(patchwork)
##################################################
set.seed(1234)

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
    d <- x
    attr(d, "Labels") %>% length() -> lab_length
    set.seed(1234)
    umap_d <- uwot::umap(d,
                         metric = "precomputed", 
                         nn_method = uwot:::dist_nn(d, k = lab_length),
                         n_neighbors = 15,
                         n_components = 2
                         )
    return_object <- data.frame(cord_1 = umap_d[,1],
                          cord_2 = umap_d[,2],
                          cell_type = attr(d, "Labels")
                          )
    return(return_object)
}

.gg_label = function(x,y) {
  title_name <- y
  # ggplot
  return_object <- eval(parse(text=paste0("ggplot(df_label_cord_yshift, 
                            aes(x = cord_1,
                                y = cord_2, 
                                label = cell_type,
                                color = ",factor(x),",
                                shape = ",factor(y),"
                                )
                            )"))) +
    geom_point(size = 6.0, alpha = 0.6) +
    scale_shape_manual(values=c(4,15)) +
    geom_text_repel(max.overlaps = Inf,
                    min.segment.length = 0) +
    ggtitle(title_name) +
    theme(plot.title = element_text(size = 30, hjust = 0.5))
  return(return_object)
}

.gg_label_yshift = function(x,y) {
  df_label_cord_yshift %>% filter(yshift_filter==1) -> df
  title_name <- paste0(y," (yshift_filter)")
  # ggplot
  return_object <- eval(parse(text=paste0("ggplot(df, 
                            aes(x = cord_1,
                                y = cord_2, 
                                label = cell_type,
                                color = ",factor(x),",
                                shape = ",factor(y),"
                                )
                            )"))) +
    geom_point(size = 6.0, alpha = 0.6) +
    scale_shape_manual(values=c(4,15)) +
    geom_text_repel(max.overlaps = Inf,
                    min.segment.length = 0) +
    ggtitle(title_name) +
    theme(plot.title = element_text(size = 30, hjust = 0.5))
  return(return_object)
}