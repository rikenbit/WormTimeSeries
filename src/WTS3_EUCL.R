source("src/functions_WTS3.R")

# #### args####
# args <- commandArgs(trailingOnly = T)
# # select data データの指定
# args_data <- args[1]
# # select timeframe 時系列区間の指定
# args_TF <- args[2]
# # select lag（ラグ） ラグの間隔の指定
# args_lag <- as.numeric(args[3])
# # select animal number 個体番号の指定
# args_sample <- args[4]
# # outputファイル名
# args_output <- args[5]
# #######################
#######################
#### test####
# select data データの指定
args_data <- c("normalize_1")
# select timeframe 時系列区間の指定
args_TF <- c("all")
# select lag（ラグ） ラグの間隔の指定
args_lag <- as.numeric(c("1"))
# select animal number 個体番号の指定
args_sample <- c("1")
# outputファイル名
args_output <- c("output/WTS3/EUCL/normalize_1/all/SampleNumber_1/EUCL.png")
########################

#### load NeuronActivity####
# inputpath <- paste('data', args_data, 'ReadData_1.RData', sep = '/')
eval(parse(text=paste0("inputpath <- paste('data', args_data, 'ReadData_",args_sample,".RData', sep = '/')")))
load(inputpath)
# ReadData <- ReadData_1
eval(parse(text=paste0("ReadData <- ReadData_",args_sample)))

#### test####
# ReadData <- ReadData[,1:50]
########

### EUCL####
# ユークリッド距離で距離行列を作成
# #### test####
# before <- Sys.time()
# tSNE <- Rtsne(d, is_distance = TRUE, dims = 2, perplexity = 5, verbose = TRUE, max_iter = 1000)
# after <- Sys.time()
# Time <- after - before
# save(tSNE, file="output/WTS3/EUCL/normalize_1/all/SampleNumber_1/tSNE.RData")
# save(Time, file="output/WTS3/EUCL/normalize_1/all/SampleNumber_1/Time.RData")
# ########
d <- diss(ReadData, "EUCL")
tSNE <- Rtsne(d, is_distance = TRUE, dims = 2, perplexity = 5, verbose = TRUE, max_iter = 1000)
########
plot(tSNE$Y)