source("src/functions_WTS2.R")

#### args####
args <- commandArgs(trailingOnly = T)
# select data データの指定
args_data <- args[1]
# select timeframe 時系列区間の指定
args_TF <- args[2]
# select lag（ラグ） ラグの間隔の指定
args_lag <- as.numeric(args[3])
# select animal number 個体番号の指定
args_sample <- args[4]
# outputファイル名
args_output <- args[5]
#######################
# #######################
# ### test####
# # select data データの指定
# args_data <- c("normalize_1")
# # select timeframe 時系列区間の指定
# args_TF <- c("all")
# # select lag（ラグ） ラグの間隔の指定
# args_lag <- as.numeric(c("1"))
# # select animal number 個体番号の指定
# args_sample <- c("1")
# # outputファイル名
# # args_output <- c("output/WTS2/heatmap/normalize_1/all/SampleNumber_1/τ1.png")
# args_output <- c("output/WTS2/heatmap/normalize_1/all/SampleNumber_1/new/τ1.png")
# #######################

#### load NeuronActivity####
# inputpath <- paste('data', args_data, 'ReadData_1.RData', sep = '/')
eval(parse(text=paste0("inputpath <- paste('data', args_data, 'ReadData_",args_sample,".RData', sep = '/')")))
load(inputpath)
# ReadData <- ReadData_1
eval(parse(text=paste0("ReadData <- ReadData_",args_sample)))
mat <- as.matrix(ReadData)
rownames(mat) <- rownames(ReadData)

#### CCF_ACF####
CCF_ACF_mat <- sapply(1:ncol(mat), function(x) {
    sapply(1:ncol(mat), function(z){
        if(x!=z){
            resTmp <- ccf(x = mat[, x],
                          y = mat[, z], plot=F, 
                          na.action = na.contiguous, 
                          lag.max = args_lag)
            resTmp$acf[which.max(resTmp$lag)]
        } else{
            resTmp <- acf(x = mat[, x], 
                          plot=F, 
                          na.action = na.contiguous, 
                          lag.max = args_lag)
            resTmp$acf[which.max(resTmp$lag)]
        }
    })
})

#### reshape CCF_ACF####
# 行をSender 列をReceiverに転値する
CCF_ACF_mat <- t(CCF_ACF_mat)
# rename colnames,rownames
rownames(CCF_ACF_mat) <- colnames(mat)
colnames(CCF_ACF_mat) <- colnames(mat)
# matrix to dataframe
CCF_ACF_df <- data.frame(CCF_ACF_mat)
# rename colnames,rownames
rownames(CCF_ACF_df) <- colnames(mat)
colnames(CCF_ACF_df) <- colnames(mat)
# convert wider to longer
CCF_ACF_df %>%
    rownames_to_column("cell_Sender") %>%
        pivot_longer(-cell_Sender, names_to = "cell_Receiver", values_to = "CCF_ACF") -> df
    
#### ggplot####
titlename <- paste0("TF",args_TF,"_SampleNumber",args_sample,"_τ",args_lag)

ghm <- ggplot_ghm(df) + ggtitle(titlename) + theme(plot.title = element_text(size = 30, hjust = 0.5))
ggsave(filename = args_output, plot = ghm, dpi = 80, width = 24.0, height = 22.0)
