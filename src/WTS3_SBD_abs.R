install.packages("dtwclust_5.5.6.tar.gz", repos = NULL, type = "source", INSTALL_opts = '--no-lock')
# install.packages("dtwclust_abstest_5.5.6.tar.gz", repos = NULL, type = "source", INSTALL_opts = '--no-lock')
# install.packages("dtwclust_abs_5.5.6.tar.gz", repos = NULL, type = "source", INSTALL_opts = '--no-lock')

source("src/functions_WTS3_SBD_abs.R")

# #### args setting####
# args <- commandArgs(trailingOnly = T)
# # sample number 
# args_sample <- args[1]
# # path NeuronActivity Data
# args_neuron <- args[2]
# # args_time <- c("all")
# args_time <- args[3]
# # stimtiming
# args_stim_xlsx <- args[4]
# # y-shift計算対象の細胞
# args_shift <- args[5]
# # output SBD_abs距離行列
# args_SBD_abs <- args[6]
# # output SBD_abs yshift
# args_yshift <- args[7]

#### test args####
# sample number サンプル番号の指定
args_sample <- c("2")
# path NeuronActivity Data
args_neuron <- c("data/normalize_1/ReadData_2.RData")

# args_time <- c("all")
args_time <- c("stimAfter")
# stimtiming
args_stim_xlsx <- c("data/stimulation/stimulation_timing.xlsx")
# y-shift計算対象の細胞
args_shift <- c("ASER")
# output SBD_abs距離行列
args_SBD_abs <- c("output/WTS3/normalize_1/stimAfter/SBD_abs/SampleNumber_2/SBD_abs.RData")
# output SBD_abs yshift
args_yshift <- c("output/WTS3/normalize_1/stimAfter/SBD_abs/SampleNumber_2/yshift.RData")


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

#### SBD_abs####
# # 行列っぽいデータを細胞ごとにlist化
# ReadData.list <- asplit(ReadData,2)
# # create sbd matrix
# # hc <- tsclust(ReadData.list, distance = "sbd", trace = TRUE)

# # sbd(時系列データの1列目, 時系列データの2列目)の値
# sbd_dist_value <- sbd(as.numeric(ReadData.list[[1]]),
#                as.numeric(ReadData.list[[2]]))

# # sbd(時系列データの1列目, すべての時系列データ)のベクトル
# ReadData.list[[1]] %>%
#     as.numeric() -> shift_1
# names(ReadData.list) %>% 
#     purrr::map_dbl(., .sbd_d) -> sbd_dist_vec

# sbd(時系列データの1列目, すべての時系列データ)の列ベクトルをまとめた行列
SBD_zero_mat <- sapply(1:ncol(ReadData), function(x) {
    sapply(1:ncol(ReadData), function(z){
        if(x!=z){
            shift_1 <- ReadData[,x]
            shift_2 <- ReadData[,z]
            return_object <- dtwclust::SBD(shift_1,
                         shift_2, 
                         znorm = FALSE, 
                         error.check = TRUE, 
                         return.shifted = TRUE)
            return_object$dist
        } else{
            return_object <- 0
            return_object
        }
    })
})


# convert dist
d <- stats::as.dist(hc@distmat)
# save SBD_abs dist
save(d, file=args_SBD_abs)

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
# save SBD_abs yshift dataframe
save(sbd_yshift_df, file=args_yshift)