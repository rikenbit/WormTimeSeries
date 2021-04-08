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
args_output <- c("output/WTS2/heatmap/normalize_1/all/SampleNumber_1/τ1.png")
########################

#### load NeuronActivity####
# inputpath <- paste('data', args_data, 'ReadData_1.RData', sep = '/')
eval(parse(text=paste0("inputpath <- paste('data', args_data, 'ReadData_",args_sample,".RData', sep = '/')")))
load(inputpath)
# ReadData <- ReadData_1
eval(parse(text=paste0("ReadData <- ReadData_",args_sample)))

#### test####
ReadData <- ReadData[,1:4]
########

#### DTW####
# DTW 距離で距離行列を作成
d <- diss(ReadData, "DTWARP")
tSNE <- Rtsne(d, is_distance = TRUE, dims = 2, perplexity = 1, verbose = TRUE, max_iter = 500)
########
plot(tSNE$Y)

#### サンプルデータ####
library(Rtsne)
library(vegan)
df = data.frame(A = c(4, 11, 17, 0, 2, 4, 8, 10, 2, 4),
                B = c(6, 10, 7, 2, 21, 3, 3, 0, 2, 17),
                C = c(5, 2, 3, 12, 12, 14, 0, 7, 8, 2),
                D = c(7, 16, 24, 18, 31, 8, 2, 21, 3, 13))
bc <- vegdist(df, method = "bray")
tSNE <- Rtsne(bc, is_distance = TRUE, dims = 2, perplexity = 2, verbose = TRUE, max_iter = 500)
plot(tSNE)
########