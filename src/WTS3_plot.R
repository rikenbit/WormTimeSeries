source("src/functions_WTS3_plot.R")

# # args setting
# ##################################################
# args <- commandArgs(trailingOnly = T)
# # select animal number 個体番号の指定
# args_sample <- args[1]
# # input_Neuron Activity ファイル名
# args_input_n <- args[2]
# # input_stim ファイル名
# args_input_stim <- args[3]
# # input_mCherry ファイル名
# args_input_mCherry <- args[4]
# # input_Position ファイル名
# args_input_Position <- args[5]
# # input_tempdat ファイル名
# args_input_tempdata <- args[6]
# # outputファイル名
# args_output <- args[7]
# # select data データの指定
# args_data <- c("normalize_1")
# # クラスター評価手法
# args_eval <- args[8]
# # 次元圧縮手法
# args_DimRedu <- args[9]
# ##################################################
#### test args####
args_sample <- c("1")
# input_Neuron Activity ファイル名
args_input_n <- c("data/normalize_1/ReadData_1.RData")
# input_stim ファイル名
args_input_stim <- c("data/stimulation/stim_1.RData")
# input_mCherry ファイル名
args_input_mCherry <- c("data/mCherry/mCherry_1.RData")
# input_Position ファイル名
args_input_Position <- c("data/Position/Position_1.RData")
# input_tempdat ファイル名
args_input_tempdata <- c("output/WTS3/SBD/normalize_1/all/tsne/ARI/cls_tempdata/SampleNumber_1.RData")
# outputファイル名
args_output <- c("output/WTS3/SBD/normalize_1/all/tsne/ARI/plot/SampleNumber_1.png")
# select data データの指定
args_data <- c("normalize_1")
# クラスター評価手法
args_eval <- c("ARI")
# 次元圧縮手法
args_DimRedu <- c("tsne")
#######################

# Neuron Activity Data
##################################################
load(args_input_n)
# input_n <- ReadData_1
eval(parse(text=paste0("input_n <- ReadData_",args_sample)))
# ReadData_1[,'102'] %>% as.numeric() -> nactivity
#### test####
input_n[,1] %>% 
    as.numeric() -> nactivity
########
#TimeFrame
# rownames(ReadData_1) %>% as.numeric() -> timeframe
rownames(input_n) %>% 
    as.numeric() -> timeframe
##################################################

# Stimulation Data
##################################################
load(args_input_stim)
# input_stim <- stim_1
eval(parse(text=paste0("input_stim <- stim_",args_sample)))
# stim_1 %>% as.numeric() -> stimtiming
input_stim %>%
    as.numeric() -> stimtiming
# cut TimeFrame
stimtiming[1:length(timeframe)] -> stimtiming
##################################################

# mCherry
##################################################
load(args_input_mCherry)
# input_mCherry <- mCherry_1
eval(parse(text=paste0("input_mCherry <- mCherry_",args_sample)))
# mCherry_1[,1] %>% as.numeric() -> mcherry
eval(parse(text=paste0("mCherry_",args_sample,"[,'",args_celltype,"'] %>% as.numeric() -> mcherry")))
input_mCherry[,1] %>% 
    as.numeric() -> mcherry
# cut TimeFrame
mcherry[1:length(timeframe)] -> mcherry
##################################################

# Position
##################################################
load(args_input_Position)
# input_Position <- Position_1
eval(parse(text=paste0("input_Position <- Position_",args_sample)))
# Position_1$MoveX %>% as.numeric() -> position
input_Position$MoveX %>% 
    as.numeric() -> position
# cut TimeFrame
position[1:length(timeframe)] -> position
##################################################

# dataframe for ggplot
##################################################
data.frame(
        TimeFrame = timeframe,
        Nactivity = nactivity,
        StimTiming = stimtiming,
        mCherry = mcherry,
        Position = position,
        stringsAsFactors = FALSE
) -> g

# rollmean
g %>% mutate(N_roll = roll_meanr(Nactivity, n=51, align="right", fill=NA)) -> g_roll

# diff
diff_value <- 50
n_diff <- append(rep(NA, diff_value), diff(g$Nactivity, diff_value))
g_roll %>% 
    mutate(N_diff = n_diff) -> g_roll_diff
##################################################

# ggplot
##################################################
p_1 <- ggplot(data = g_roll_diff, aes(TimeFrame))
p_2 <- p_1 +
        geom_line(aes(y = Nactivity, colour = "Nactivity"), size = 0.5) +
        geom_line(aes(y = N_roll, colour = "N_rollmean"), size = 1.5, alpha = 0.7) +
        geom_line(aes(y = N_diff, colour = "N_diff"), size = 1.5, alpha = 0.7, linetype = "dotted")
p_3 <- p_1 +        
         geom_line(aes(y = StimTiming, colour = "StimTiming") , size = 1.5)
p_4 <- p_1 +        
         geom_line(aes(y = mCherry, colour = "mCherry") , size = 1.5)
p_5 <- p_1 +        
         geom_line(aes(y = Position, colour = "Position") , size = 1.5)


sX <- scale_x_continuous(name = "TimeFrame(1frame/0.2sec)",    # 軸の名前を変える
                         breaks = seq(0, length(timeframe), by= 1000),     # 軸の区切りを0,2,4にする
                        )
s_2 <- scale_color_manual(values = c("orange", "green", "black"))
s_3 <- scale_color_manual(values = c("purple"))
s_4 <- scale_color_manual(values = c("red"))
s_5 <- scale_color_manual(values = c("blue"))


# title <- ggtitle('SampleNumber1_CellNumber1_X1_raw_CFP')
eval(parse(text=paste0("title <- ggtitle('SampleNumber",args_sample,"_CellNumber",args_cell,"_",args_celltype,"_",args_datadir,"')")))
t_1 <- theme(plot.title = element_text(size = 30, hjust = 0.5))
t_2 <- theme(axis.title = element_text(size = 20))
t_3 <- theme(legend.title = element_text(size = 28),
             legend.text = element_text(size = 20))


gg2 <- p_2 +
    s_2 +
    sX +
    title +
    t_1 +
    t_2 +
    t_3 +
    labs(colour="each data")
gg3 <- p_3 +
    s_3 +
    sX +
    t_2 +
    t_3 +
    labs(colour="each data")
gg4<- p_4 +
    s_4 +
    sX +
    t_2 +
    t_3 +
    labs(colour="each data")
gg5 <- p_5 +
    s_5 +
    sX +
    t_2 +
    t_3 +
    labs(colour="each data")
gg <- gg2 + gg3 + gg4 + gg5 + plot_layout(ncol = 1, heights = c(2, 1, 1, 1))
##################################################

# ggsave
##################################################
ggsave(filename = args_output, plot = gg, dpi = 100, width = 20.0, height = 15.0)
##################################################