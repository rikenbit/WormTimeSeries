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

### EUCL####
# ユークリッド距離で距離行列を作成
### test####
ReadData <- ReadData[1:5,1:5]
########
d <- diss(ReadData, "EUCL")
tSNE <- Rtsne(d, is_distance = TRUE, dims = 2, perplexity = 5, verbose = TRUE, max_iter = 1000)
plot(tSNE$Y)

### test####
df_test <- data.frame(列A = c(0.2,0.3,0.4),
                  列B = c(0.4,0.6,0.8),
                  列C = c(0.19,0.29,0.39),
                  列D = c(0.41,0.61,0.81),
                  列E = c(0.39,0.59,0.79)
                  )
d_test <- diss(df_test, "EUCL")
tSNE_test <- Rtsne(d_test, is_distance = TRUE, dims = 2, perplexity = 1, verbose = TRUE, max_iter = 1000)
df_tSNE_test <- data.frame(tsne_1 = tSNE_test$Y[,1],
                           tsne_2 = tSNE_test$Y[,2],
                           列名 = attr(d_test, "Labels"))
install.packages("ggrepel") 
library("ggrepel")
gg<- ggplot(df_tSNE_test, aes(x = tsne_1, y = tsne_2, label = 列名)) +
    geom_point() +
    geom_text_repel()
########