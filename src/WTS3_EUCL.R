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

#### EUCL####
# ユークリッド距離で距離行列を作成
d <- diss(ReadData, "EUCL")
tSNE <- Rtsne(d, is_distance = TRUE, dims = 2, perplexity = 5, verbose = TRUE, max_iter = 1000)
df_tSNE <- data.frame(tsne_1 = tSNE$Y[,1],
                      tsne_2 = tSNE$Y[,2],
                      celltype = attr(d, "Labels"))

#### ggplot####
# 数字から始まる細胞型名を取り除場合
# df_tSNE$celltype[which(!is.na(str_extract(df_tSNE$celltype, "^\\d")))] <- ""
# 元のデータフレームに戻す
# df_tSNE$celltype <- attr(d, "Labels")

gg <- ggplot(df_tSNE, aes(x = tsne_1, y = tsne_2, label = celltype)) +
    geom_point() +
    geom_text_repel(max.overlaps = Inf,
                   min.segment.length = 0)
    # geom_label_repel(max.overlaps = Inf,
    #                  min.segment.length = 0)
# ggsave
ggsave(filename = args_output, plot = gg, dpi = 100, width = 10.0, height = 10.0)