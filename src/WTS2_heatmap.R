source("src/functions_WTS2.R")

# #### args####
args <- commandArgs(trailingOnly = T)
# select data データの指定
args_data <- args[1]
# select timeframe 時系列区間の指定
args_TF <- args[2]
# select lag（ラグ） ラグの間隔の指定
args_lag <- as.numeric(args[3])
# # select Acf 自己相関の指定
# args_Acf <- args[4]
# select animal number 個体番号の指定
args_sample <- args[5]
# # select cell number 細胞番号の指定
# args_cell <- args[6]
# # select celltype 細胞型名の指定
# args_celltype <- args[7]
# outputファイル名
args_output <- args[8]
#######################
#######################
#### test####
# select data データの指定
args_data <- c("normalize_1")
# select timeframe 時系列区間の指定
args_TF <- c("all")
# select lag（ラグ） ラグの間隔の指定
args_lag <- as.numeric(c("300"))
# # select Acf 自己相関の指定
# args_Acf <- c("Acf")
# select animal number 個体番号の指定
args_sample <- c("1")
# # select cell number 細胞番号の指定
# args_cell <- c("1")
# # select celltype 細胞型名の指定
# args_celltype <- c("102")
# outputファイル名
args_output <- c("output/WTS2/heatmap/Data_normalize_1/TF_all/SampleNumber_1/τ1.png")
########################

#### load NeuronActivity####
# inputpath <- paste('data', args_data, 'ReadData_1.RData', sep = '/')
eval(parse(text=paste0("inputpath <- paste('data', args_data, 'ReadData_",args_sample,".RData', sep = '/')")))
load(inputpath)
# ReadData <- ReadData_1
eval(parse(text=paste0("ReadData <- ReadData_",args_sample)))
mat <- as.matrix(ReadData)
rownames(mat) <- rownames(ReadData)

#### load stimulation timing####
# inputpath <- paste('data','stimulation', 'stim_1.RData', sep = '/')
eval(parse(text=paste0("inputpath <- paste('data','stimulation', 'stim_",args_sample,".RData', sep = '/')")))
load(inputpath)
# stim <- stim_1
eval(parse(text=paste0("stim <- stim_",args_sample)))
stim %>% 
    as.numeric() -> stimtiming
timeframe <- as.numeric(rownames(ReadData))
#dataframe
data.frame(
  TimeFrame = timeframe,
  # Nactivity = nactivity,
  StimTiming = stimtiming,
  stringsAsFactors = FALSE
) -> g_all

#### first stim TimeFrame####
length(stimtiming) -> stim_all
g_all %>%
    filter(StimTiming != 0) %>%
        slice_head() %>%
            .$TimeFrame -> stim_after
stim_after -1 -> stim_before
eval(parse(text=paste0("TF <- stim_",args_TF)))
mat <- mat[1:TF,]
#### CCF_ACF####
#test
args_lag <- as.numeric(c("1"))
mat <- mat[,1:10]
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
tilename <- paste0("TF",args_TF,"_SampleNumber",args_sample,"_τ",args_lag)

ghm <- ggplot_ghm(df) + ggtitle(tilename) +theme(plot.title = element_text(size = 30, hjust = 0.5))
ggsave(filename = args_output, plot = ghm, dpi = 100, width = 10.0, height = 7.0)
