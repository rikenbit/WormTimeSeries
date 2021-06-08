source("src/functions_WTS3_plot.R")

#### test args####
args_sample <- c("2")
# input_Neuron Activity ファイル名
args_input_n <- c("data/normalize_1/ReadData_2.RData")
# input_stim ファイル名
args_input_stim <- c("data/stimulation/stim_2.RData")
# input_mCherry ファイル名
args_input_mCherry <- c("data/mCherry/mCherry_2.RData")
# input_Position ファイル名
args_input_Position <- c("data/Position/Position_2.RData")
# input_tempdat ファイル名
args_input_tempdata <- c("output/WTS3/SBD/normalize_1/all/tsne/ARI/cls_tempdata/SampleNumber_2.RData")
# outputファイル名
args_output <- c("output/WTS3/SBD/normalize_1/all/tsne/ARI/plot/SampleNumber_2.png")
# select data データの指定
args_data <- c("normalize_1")
# クラスター評価手法
args_eval <- c("ARI")
# 次元圧縮手法
args_DimRedu <- c("tsne")
# フィルタリング
args_filter <- c("stim_cell")
# y-shift計算対象の細胞
args_shift <- c("ASER")
#######################
#### SBD yshift####
load(args_input_n)
# input_n <- ReadData_2
eval(parse(text=paste0("input_n <- ReadData_",args_sample)))
# 行列っぽいデータを各細胞（列）ごとに分割
input_n.list <- asplit(input_n, 2)
# cal SBD
library(dtwclust)
RNGkind(kind = "Mersenne-Twister")
set.seed(1234)

# dtwclust::SBD()
sbd_y = function(x) {
    shift_2 <- input_n.list[[x]] %>% as.numeric()
    sbd <- dtwclust::SBD(shift_1,
                         shift_2, 
                         znorm = FALSE, 
                         error.check = TRUE, 
                         return.shifted = TRUE)
    return(sbd$yshift)
    # return(sbd)
}
list_cell_type <- colnames(input_n)
eval(parse(text=paste0("shift_1 <- input_n.list$",args_shift," %>% as.numeric()")))
seq(1:length(list_cell_type)) %>% 
    purrr::map(., sbd_y) %>%
        as.data.frame() -> sbd_yshift_df
colnames(sbd_yshift_df) <- list_cell_type
# # dtwclust::tsclust()
# hc <- dtwclust::tsclust(input_n.list,
#                         distance = "sbd",
#                         trace = TRUE)
# hc <- dtwclust::tsclust(input_n.list, 
#                         distance = "sbd", 
#                         trace = TRUE,
#                         return.shifted = TRUE)
# hc@distmat %>% head()

# proxy::dist()
# library(proxy)
# test_list <- input_n.list[1:10]
# sbD <- proxy::dist(test_list, 
#                    test_list, 
#                    method = "SBD", 
#                    znorm = TRUE)