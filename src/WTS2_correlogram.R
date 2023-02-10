source("src/functions_WTS2.R")

#### 6data 3TF 28sample 200cell ####

#### args####
args <- commandArgs(trailingOnly = T)
# select data データの指定
args_data <- args[1]
# select timeframe 時系列区間の指定
args_TF <- args[2]
# select lag（ラグ） ラグの間隔の指定
args_lag <- as.numeric(args[3])
# select Acf 自己相関の指定
args_Acf <- args[4]
# select animal number 個体番号の指定
args_sample <- args[5]
# select cell number 細胞番号の指定
args_cell <- args[6]
# select celltype 細胞型名の指定
args_celltype <- args[7]
# outputファイル名
args_output <- args[8]
#######################

#### load NeuronActivity####
# inputpath <- paste('data', args_data, 'ReadData_1.RData', sep = '/')
eval(parse(text=paste0("inputpath <- paste('data', args_data, 'ReadData_",args_sample,".RData', sep = '/')")))
load(inputpath)
# ReadData <- ReadData_1
eval(parse(text=paste0("ReadData <- ReadData_",args_sample)))
# TimeFrame
timeframe <- as.numeric(rownames(ReadData))
# NeuronActivity
# nactivity <- ReadData[,'細胞型名']
eval(parse(text=paste0("nactivity <- ReadData[,'",args_celltype,"']")))

#### load stimulation timing####
# inputpath <- paste('data','stimulation', 'stim_1.RData', sep = '/')
eval(parse(text=paste0("inputpath <- paste('data','stimulation', 'stim_",args_sample,".RData', sep = '/')")))
load(inputpath)
# stim <- stim_1
eval(parse(text=paste0("stim <- stim_",args_sample)))
stim %>%
    .[1:length(nactivity)] %>%
        as.numeric() -> stimtiming

#### dataframe####
data.frame(
  TimeFrame = timeframe,
  Nactivity = nactivity,
  StimTiming = stimtiming,
  stringsAsFactors = FALSE
) -> g_all

#### first stim TimeFrame####
g_all %>%
    filter(StimTiming != 0) %>%
        slice_head() %>%
            .$TimeFrame -> stim_after
stim_before <- stim_after -1
head(g_all,eval(parse(text=paste0(stim_before)))) -> g_before
g_all %>%
    .[eval(parse(text=paste0(stim_after))):nrow(g_all),] -> g_after

#### ggAcf####
# g_TF <- g_all
eval(parse(text=paste0("g_TF <- g_",args_TF)))
p_type <- switch(args_Acf,
              "Acf" = "correlation",
              "pAcf" = "partial",
              stop("Only can use Acf, pAcf")
)
# calculate Acf
g_TF$Nactivity %>%
  ggAcf(lag.max = args_lag, type = p_type, plot = TRUE) -> p
  # ggAcf(lag.max = 300, type = p_type, plot = TRUE) -> p

#### ggplot####
# title name
# titlename <- raw_CFP_SampleNumber1_102_Acf_all_stim124_lag50
eval(parse(text=paste0("titlename <- paste(args_data,
                       'SampleNumber",args_sample,"',
                        args_celltype,
                        args_Acf,
                        args_TF,
                        'stim",stim_after,"',
                        'lag",args_lag,"', sep = '_')")))
title <- ggtitle(titlename)
# title theme
t_1 <- theme(plot.title = element_text(size = 16, hjust = 0.5))
t_2 <- theme(axis.title = element_text(size = 16))
gg <- p +
  title +
  t_1 +
  t_2
# ggsave
ggsave(filename = args_output, plot = gg, dpi = 100, width = 10.0, height = 7.0)