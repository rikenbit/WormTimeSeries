#library
##################################################
library(mclust)
library(igraph)
library(tidyverse)
##################################################
.check_args_shift = function(x) {
  x -> sample_cell_type
  sample_cell_type %>% 
    str_count(., pattern="ASER") %>%
    sum() -> check_ASER
  sample_cell_type %>% 
    str_count(., pattern="BAGR") %>%
    sum() -> check_BAGR
  sample_cell_type %>% 
    str_count(., pattern="BAGL") %>%
    sum() -> check_BAGL
  if (check_ASER >= 1) {
    return_object <- "ASER"
  } else if (check_BAGR >= 1) {
    return_object <- "BAGR"
  } else if (check_BAGL >= 1) {
    return_object <- "BAGL"
  } else {
    return_object <- sample_cell_type[1]
  }
  return(return_object)
}

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

#### eval Max####
.max_eval = function(x) {
  x %>%
    filter(., eval_value == max(eval_value)) %>%
    .$set_cutree  -> return_object
  return(return_object)
}
#### eval Min####
.min_eval = function(x) {
  x %>%
    filter(., eval_value == min(eval_value)) %>%
    .$set_cutree  -> return_object
  return(return_object)
}