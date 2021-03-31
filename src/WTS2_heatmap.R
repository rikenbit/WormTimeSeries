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
# TimeFrame
timeframe <- as.numeric(rownames(ReadData))
# NeuronActivity 
eval(parse(text=paste0("nactivity <- ReadData[,'",args_celltype,"']")))

#### dataframe####
data.frame(
  TimeFrame = timeframe,
  Nactivity = nactivity,
  # StimTiming = stimtiming,
  stringsAsFactors = FALSE
) -> g_all
ReadData %>%
  rowid_to_column("timeframe") %>%
  pivot_longer(-timeframe, names_to = "celltype", values_to = "nactivity") -> g_all_cell
#### first stim TimeFrame####

#### Acf/Ccf####
# ccf(x, y, plot = FALSE, lag.max = 1, type = c("correlation"))
# g_all$Nactivity %>% acf(plot = FALSE,
#                         lag.max = 1, 
#                         type = c("correlation")) -> acf_lag1
####sample dat CCF value fix####
mat <- matrix(rnorm(100), ncol=10)
mat[sample(1:length(mat), 10)] <- NA 
# args_lag <- c("5")
args_lag <- as.numeric(c("3"))

res_fix_if <- sapply(1:ncol(mat), function(x) {
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
        # print(resTmp)
        # print(resTmp$lag)
        # print(which.max(resTmp$lag))
        # print(resTmp$acf[which.max(resTmp$lag)])
    })
})
########

#### real data####
args_lag <- as.numeric(c("3"))

res_fix_if <- sapply(1:ncol(mat), function(x) {
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
        # print(resTmp)
        # print(resTmp$lag)
        # print(which.max(resTmp$lag))
        # print(resTmp$acf[which.max(resTmp$lag)])
    })
})


#### ggplot####
