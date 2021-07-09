source("src/functions_WTS3_DTW.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
# sample number 
args_sample <- args[1]
# input NeuronActivity Data
args_neuron <- args[2]
# time range
args_time <- args[3]
# stimtiming
args_stim_xlsx <- args[4]
# output SBD距離行列
args_DTW <- args[5]

#### test args####
# # sample number サンプル番号の指定
# args_sample <- c("1")
# # input NeuronActivity Data
# args_neuron <- c("data/normalize_1/ReadData_1.RData")
# # time range
# # args_time <- c("all")
# args_time <- c("stimAfter")
# # stimtiming
# args_stim_xlsx <- c("data/stimulation/stimulation_timing.xlsx")
# 
# # output SBD距離行列
# # args_DTW <- c("output/WTS3/normalize_1/all/DTW/SampleNumber_1/DTW.RData")
# args_DTW <- c("output/WTS3/normalize_1/stimAfter/DTW/SampleNumber_1/DTW.RData")

#### load NeuronActivity####
load(args_neuron)
# 元データがディレクトリごとではなく，ファイル名で各サンプルがわかれている対応
eval(parse(text=paste0("ReadData <- ReadData_",args_sample)))
#### test####
ReadData <- ReadData_1[,1:3]
#### test####
#### switch time range & DTW calc####
d <- switch(args_time,
            "all" = .dtw_cal_all(ReadData),
            "stimAfter" = .dtw_cal_stimAfter(ReadData, args_stim_xlsx),
            stop("Only can use all, stimAfter ")
            )

#### DTW####
save(d, file=args_DTW)
