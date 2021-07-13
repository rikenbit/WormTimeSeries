source("src/functions_WTS3_EUCL.R")

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
args_EUCL <- args[5]

# #### test args####
# # sample number サンプル番号の指定
# args_sample <- c("1")
# # input NeuronActivity Data
# args_neuron <- c("data/normalize_1/ReadData_1.RData")
# # time range
# # args_time <- c("all")
# args_time <- c("stimAfter")
# # stimtiming
# args_stim_xlsx <- c("data/stimulation/stimulation_timing.xlsx")

# # output EUCL距離行列
# # args_EUCL <- c("output/WTS3/normalize_1/all/EUCL/SampleNumber_1/EUCL.RData")
# args_EUCL <- c("output/WTS3/normalize_1/stimAfter/EUCL/SampleNumber_1/EUCL.RData")

#### load NeuronActivity####
load(args_neuron)
# 元データがディレクトリごとではなく，ファイル名で各サンプルがわかれている対応
eval(parse(text=paste0("ReadData <- ReadData_",args_sample)))

#### switch time range & EUCL calc####
d <- switch(args_time,
            "all" = .eucl_cal_all(ReadData),
            "stimAfter" = .eucl_cal_stimAfter(ReadData, args_stim_xlsx),
            stop("Only can use all, stimAfter ")
            )

#### EUCL####
save(d, file=args_EUCL)