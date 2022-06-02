source("src/functions_WTS4_dhyper.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
#### test args####
#
####test####
クラスタリング結果
label_table_cls <- read.csv("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/Merged_tsne/label_table_k6.csv")
# 豊島さんのアノテーション情報
label_table_ann <- read.csv("data/label_table_ann.csv")
table_merged <- merge(label_table_cls, 
                      label_table_ann, 
                      by.x = "Cell_type", 
                      by.y = "Cell_type", 
                      all.x = TRUE
                      )

#### クラスタ5の場合####
table_merged %>% 
    dplyr::filter(.,Cluster==5) %>% 
        .[,"PC1_neg"] %>% 
            na.omit() %>% 
                length() -> x
table_merged[,"PC1_neg"] %>% 
    na.omit() %>% 
        length() -> m
table_merged %>% 
    nrow() - m -> n
table_merged %>% 
    dplyr::filter(.,Cluster==5) %>% 
    length() -> k
p <- dhyper(x, m, n, k)
Log10p <- -log10(p)

#### クラスタnの場合####
table_merged[,"Cluster"] %>% 
    unique() %>% 
        sort() %>% 
    purrr::map_dbl(.,dhyper_all_PC1_neg)
#### クラスタn 項目ぜんぶ####
colnames(label_table_ann) %>% 
    .[3:length(.)] %>% 
        purrr::map(., dhyper_all_cls)-> list_dhyper
bind_cols(list_dhyper) -> table_dhyper
colnames(label_table_ann) %>% 
    .[3:length(.)] -> colnames(table_dhyper)
#### ggsave####