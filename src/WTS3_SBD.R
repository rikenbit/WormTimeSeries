source("src/functions_WTS3_SBD.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
# sample number 
args_sample <- args[1]
# path NeuronActivity Data
args_neuron <- args[2]
# args_time <- c("all")
args_time <- args[3]
# stimtiming
args_stim_xlsx <- args[4]
# y-shift計算対象の細胞
args_shift <- args[5]
# output SBD距離行列
args_SBD <- args[6]
# output SBD yshift
args_yshift <- args[7]

# #### test args####
# # sample number サンプル番号の指定
# args_sample <- c("1")
# # path NeuronActivity Data
# args_neuron <- c("data/normalize_1/ReadData_1.RData")

# # args_time <- c("all")
# args_time <- c("stimAfter")
# # stimtiming
# args_stim_xlsx <- c("data/stimulation/stimulation_timing.xlsx")
# # y-shift計算対象の細胞
# args_shift <- c("ASER")
# # output SBD距離行列
# args_SBD <- c("output/WTS3/normalize_1/stimAfter/SBD/SampleNumber_1/SBD.RData")
# # output SBD yshift
# args_yshift <- c("output/WTS3/normalize_1/stimAfter/SBD/SampleNumber_1/yshift.RData")


#### load NeuronActivity####
load(args_neuron)
# 元データがディレクトリごとではなく，ファイル名で各サンプルがわかれている対応
eval(parse(text=paste0("ReadData <- ReadData_",args_sample)))

#### switch time range & trim ReadData####
ReadData <- switch(args_time,
                   "all" = .ReadData_all(ReadData),
                   "stimAfter" = .ReadData_stimAfter(ReadData, args_stim_xlsx),
                   stop("Only can use all, stimAfter ")
                   )

#### SBD####
# 行列っぽいデータを細胞ごとにlist化
ReadData.list <- asplit(ReadData,2)
# create sbd matrix
hc <- tsclust(ReadData.list, distance = "sbd", trace = TRUE)
# convert dist
d <- stats::as.dist(hc@distmat)
# save SBD dist
save(d, file=args_SBD)

#### yshift####
# prepare shift_1
eval(parse(text=paste0("shift_1 <- ReadData.list$",args_shift," %>% as.numeric()")))
colnames(ReadData) %>% 
    purrr::map(., .sbd_y) %>%
        as.data.frame() -> sbd_yshift_df_wide
colnames(sbd_yshift_df_wide) <- colnames(ReadData)
# convert long df
sbd_yshift_df_wide %>% 
    # rownames_to_column("time_frame") %>% 
    mutate(time_frame = rownames(ReadData)) %>% 
        pivot_longer(-time_frame, 
                     names_to = "cell_type", 
                     values_to = "yshift") -> sbd_yshift_df
# save SBD yshift dataframe
save(sbd_yshift_df, file=args_yshift)