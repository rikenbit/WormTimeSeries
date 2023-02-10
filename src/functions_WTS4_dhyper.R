#library
##################################################
library(tidyverse)
##################################################
# # 赤玉:
# # 白玉:
#
# x <- 5 # 取り出された赤玉の数
# m <- 6 # 赤玉の総数
# n <- 194 - 6 # 白玉の総数
# k <- 7 # 取り出された玉の数
#
# p <- dhyper(x, m, n, k)
#
# Log10p <- -log10(p)

dhyper_all_cls = function(x) {
    purrr::map_dbl(sort(unique(table_merged[,"Cluster"])),
                   dhyper_all,
                   y = x)
    }

dhyper_all = function(x,y) {
    table_merged %>%
        dplyr::filter(.,Cluster==x) %>%
        .[,y] %>%
        na.omit() %>%
        length() -> z

    table_merged[,y] %>%
        na.omit() %>%
        length() -> m

    table_merged %>%
        nrow() - m -> n

    table_merged %>%
        dplyr::filter(.,Cluster==x) %>%
        length() -> k

    p <- dhyper(z, m, n, k)
    # Log10p <- -log10(p)
    # return(Log10p)
    return(p)
}