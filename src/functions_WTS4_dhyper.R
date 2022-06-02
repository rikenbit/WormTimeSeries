#library
##################################################
library(tidyverse)
##################################################
# # 赤玉: 注目している細胞のアノテーションがついている細胞（例: PC1_posの細胞）
# # 白玉: 注目している細胞のアノテーションがついていない細胞（例: PC1_pos以外の細胞）
# 
# x <- 5 # 取り出された赤玉の数（例: あるクラスタに属するPC1_pos細胞の数）
# m <- 6 # 赤玉の総数（例: PC1_pos細胞の総数）
# n <- 194 - 6 # 白玉の総数（例: 全細胞数 - PC1_pos細胞の総数）
# k <- 7 # 取り出された玉の数（例: あるクラスタのメンバー数）
# 
# p <- dhyper(x, m, n, k)
# 
# Log10p <- -log10(p)

dhyper_all_PC1_neg = function(x) {
    table_merged %>% 
        dplyr::filter(.,Cluster==x) %>% 
        .[,"PC1_neg"] %>% 
        na.omit() %>% 
        length() -> y
    table_merged[,"PC1_neg"] %>% 
        na.omit() %>% 
        length() -> m
    table_merged %>% 
        nrow() - m -> n
    table_merged %>% 
        dplyr::filter(.,Cluster==x) %>% 
        length() -> k
    p <- dhyper(y, m, n, k)
    Log10p <- -log10(p)
	return(Log10p)
}

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
    Log10p <- -log10(p)
    return(Log10p)
}