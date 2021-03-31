source("src/functions_WTS2.R")

# #### args####
# args <- commandArgs(trailingOnly = T)
# # select data データの指定
# args_data <- args[1]
# # select timeframe 時系列区間の指定
# args_TF <- args[2]
# # select lag（ラグ） ラグの間隔の指定
# args_lag <- as.numeric(args[3])
# # select Acf 自己相関の指定
# args_Acf <- args[4]
# # select animal number 個体番号の指定
# args_sample <- args[5]
# # select cell number 細胞番号の指定
# args_cell <- args[6]
# # select celltype 細胞型名の指定
# args_celltype <- args[7]
# # outputファイル名
# args_output <- args[8]
# #######################
#######################
#### test####
# select data データの指定
args_data <- c("normalize_1")
# select timeframe 時系列区間の指定
args_TF <- c("all")
# select lag（ラグ） ラグの間隔の指定
# args_lag <- c("300")
args_lag <- as.numeric(c("300"))
# select Acf 自己相関の指定
args_Acf <- c("Acf")
# select animal number 個体番号の指定
args_sample <- c("1")
# select cell number 細胞番号の指定
args_cell <- c("1")
# select celltype 細胞型名の指定
args_celltype <- c("102")
# outputファイル名
args_output <- c("output/WTS2/correlogram/raw_CFP/all/50/Acf/SampleNumber_1/CellNumber_87_CellType_ADAR.png")
########################

#### load NeuronActivity####
# inputpath <- paste('data', args_data, 'ReadData_1.RData', sep = '/')
eval(parse(text=paste0("inputpath <- paste('data', args_data, 'ReadData_",args_sample,".RData', sep = '/')")))
load(inputpath)
# ReadData <- ReadData_1
eval(parse(text=paste0("ReadData <- ReadData_",args_sample)))
mat <- as.matrix(ReadData)
rownames(mat) <- rownames(ReadData)
#### first stim TimeFrame####

#### real data####
#test
args_lag <- as.numeric(c("1"))
mat <-mat[1:6000,1:10]
#
CCF_ACF_mat <- sapply(1:ncol(mat), function(x) {
    sapply(1:ncol(mat), function(z){
        if(x!=z){
            resTmp <- ccf(x = mat[, x],
                          y = mat[, z], plot=F, 
                          na.action = na.contiguous, 
                          lag.max = args_lag)
            resTmp$acf[which.max(resTmp$lag)]
        } else{
            resTmp <- acf(x = mat[, 1], 
                          plot=F, 
                          na.action = na.contiguous, 
                          lag.max = args_lag)
            resTmp$acf[which.max(resTmp$lag)]
        }
    })
})
# 行をSender 列をReceiverに転値する
CCF_ACF_mat <- t(CCF_ACF_mat)
rownames(CCF_ACF_mat) <- colnames(mat)
colnames(CCF_ACF_mat) <- colnames(mat)
# matrix to dataframe
CCF_ACF_df <- data.frame(CCF_ACF_mat)
#colnamesとrownames
rownames(CCF_ACF_df) <- colnames(mat)
colnames(CCF_ACF_df) <- colnames(mat)
# convert wider to longer
CCF_ACF_df %>%
    rownames_to_column("cell_Sender") %>%
        pivot_longer(-cell_Sender, names_to = "cell_Receiver", values_to = "CCF_ACF") -> df
    
#### ggplot####
# ghm <- ggplot(df, aes(x = cell_Receiver, y = cell_Sender, fill = CCF_ACF))
# ghm <- ghm + geom_tile()

ghm <- ggplot_ghm(df) + ggtitle(args_sample) +theme(plot.title = element_text(size = 30, hjust = 0.5))
#### reorder#### 
# CCF_ACF_mat %>%
#     heatmap(., scale = "none", 
#             hclustfun = function(x) {hclust(x, method = "ward.D2")})-> clr


# cell_Receiver.idx <- colnames(CCF_ACF_mat)[clr$colInd]
# df$cell_Receiver <- factor(df$cell_Receiver, levels = cell_Receiver.idx)
# 
# cell_Sender.idx <- colnames(CCF_ACF_mat)[clr$rowInd]
# df$cell_Sender <- factor(df$cell_Sender, levels = cell_Sender.idx)
#### ラグごとにアニメーションにするので並び替えしない####