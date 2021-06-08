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
eval(parse(text=paste0("input_n <- ReadData_",args_sample)))
# 行列っぽいデータを細胞ごとにlist化
input_n.list <- asplit(input_n, 2)
# setup SBD
library(dtwclust)
RNGkind(kind = "Mersenne-Twister")
set.seed(1234)

# prepare sbd y-shift
list_cell_type <- colnames(input_n)
eval(parse(text=paste0("shift_1 <- input_n.list$",args_shift," %>% as.numeric()")))
# sbd y-shift
seq(1:length(list_cell_type)) %>% 
    purrr::map(., sbd_y) %>%
        as.data.frame() -> sbd_yshift_df_wide
colnames(sbd_yshift_df_wide) <- list_cell_type
# convert long df
sbd_yshift_df_wide %>% 
    rownames_to_column("time_frame") %>% 
        pivot_longer(-time_frame, 
                     names_to = "cell_type", 
                     values_to = "y_shift") -> sbd_yshift_df