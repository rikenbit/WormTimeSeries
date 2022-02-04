source("src/functions_WTS4_Plot_SampleWeight.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
args_weight_path <- args[1]
args_sample_path <- args[2]
args_output <- args[3]

# #### test args####
# args_weight_path <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/Merged_data")
# args_sample_path <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Distance")
# args_output <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/Plot_SampleWeight.png")

#### sample number####
sample_path_list <- list.files(args_sample_path, pattern="SampleNumber_", full.names=TRUE)
sample_path_list %>% 
    str_remove(., args_sample_path) %>% 
    str_remove(., "/SampleNumber_") %>% 
    str_remove(., ".RData") %>% 
    as.numeric() %>% 
    sort() -> sample_sort_num

#### cluster number####
weight_path_list <- list.files(args_weight_path, pattern="k_Number_", full.names=TRUE)
weight_path_list %>% 
    str_remove(., args_weight_path) %>% 
    str_remove(., "/k_Number_") %>% 
    str_remove(., ".RData") %>% 
    as.numeric() %>% 
    sort() -> cluster_sort_num

#### MCMI weight table all cluster####
# purrr
cluster_sort_num %>% 
    purrr::map_dfr(., .all_k_table) -> df_weight_all

#### df_weight_group####
df_weight_all %>% 
    group_by(SampleNumber) %>% 
        summarise(weight_ave = mean(weight_abs)) -> df_weight_group
#### weight_sort####
df_weight_group %>% 
    dplyr::arrange(desc(weight_ave)) %>% 
        .$SampleNumber %>% 
            c() -> sample_sort_weight

#### list_sample_weight####
# weightの大きい順にリストを作成
# purrr
sample_sort_weight %>% 
    purrr::map(., .list_beeswarm) -> list_sample_weight

#### beeswarm####
# beeswarm(list_sample_weight, labels = sample_sort_weight)

#### ggplot####
gg <- ggplot(df_weight_all, aes(x=factor(SampleNumber, levels = sample_sort_weight), 
                                y=weight_abs, 
                                colour=factor(SampleNumber, levels = sample_sort_weight), 
                                fill = factor(SampleNumber, levels = sample_sort_weight)
                                )
             )
gg <- gg + guides(colour="none")
# https://stats.biopapyrus.jp/r/ggplot/geom-boxplot.html
gg <- gg + geom_boxplot(outlier.shape = NA, alpha =0.8) #外れ値のプロットを省く
gg <- gg + geom_point(position = position_jitter(width=0.05), size = 5.0 ,alpha = 0.7)
gg <- gg + labs(fill="sample") + xlab("sample") + ylab("weight")
gg <- gg + 
    ggtitle("weight beeswarm (sort by weight average)") +
    theme(plot.title = element_text(hjust = 0.5))
# all text size
gg <- gg + theme(text = element_text(size = 60)) 

#### ggsave####
ggsave(filename = args_output, 
       plot = gg,
       dpi = 100, 
       width = 100.0, 
       height = 10.0,
       limitsize = FALSE
       )