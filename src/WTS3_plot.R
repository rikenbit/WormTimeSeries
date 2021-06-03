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

#### input Neuron Activity Data####
load(args_input_n)
# input_n <- ReadData_1
eval(parse(text=paste0("input_n <- ReadData_",args_sample)))
# rownamesを一行目に追加したデータフレーム作成
input_n %>% 
    rownames_to_column("time_frame") %>% 
    pivot_longer(-time_frame, 
                 names_to = "cell_type", 
                 values_to = "n_activity") -> df_input_n
# convert type
as.numeric(df_input_n$time_frame) -> df_input_n$time_frame

#### input other Data####
# TimeFrame
rownames(input_n) %>% 
    as.numeric() -> timeframe
# Stimulation Data
load(args_input_stim)
eval(parse(text=paste0("input_stim <- stim_",args_sample)))
input_stim %>%
    as.numeric() -> stimtiming
stimtiming[1:length(timeframe)] -> stimtiming

# mCherry
load(args_input_mCherry)
eval(parse(text=paste0("input_mCherry <- mCherry_",args_sample)))
input_mCherry[,1] %>% 
    as.numeric() -> mcherry
mcherry[1:length(timeframe)] -> mcherry

# Position
load(args_input_Position)
eval(parse(text=paste0("input_Position <- Position_",args_sample)))
input_Position$MoveX %>% 
    as.numeric() -> position
position[1:length(timeframe)] -> position

data.frame(
    time_frame = timeframe,
    stim_timing = stimtiming,
    m_cherry = mcherry,
    position = position,
    stringsAsFactors = FALSE
) -> df_input_other

#### merge input####
df_input_n %>% 
    # test 3細胞にフィルタ#
    # filter(., cell_type %in% colnames(input_n)[1:3]) %>% 
    ##
        merge(., 
              df_input_other, 
              by.x = "time_frame", 
              by.y = "time_frame", 
              all.x = TRUE) -> df_merged_other
#### merge tempdata####
load(args_input_tempdata)
merge(df_merged_other,
      df_tempdata,
      by.x = "cell_type", 
      by.y = "cell_type", 
      all.x = TRUE) -> df_merged_temp
#### filter####
df_merged_temp %>% 
    filter(., stim == 1) %>% 
        mutate(.,
           stim_timing = if_else(stim_timing == 1, 
                                 max(.$n_activity), 
                                 min(.$n_activity))) -> df_merged

#### ggplot Neuron Activity Dat####
p_1 <- ggplot(data = df_merged, aes(x = time_frame))
p_2 <- p_1 + 
    geom_line(aes(y = n_activity, 
                  group = cell_type, 
                  colour = cell_type)
    ) +
    scale_color_viridis(option = "D", discrete = T) +
    geom_text_repel(data = subset(df_merged, 
                                  time_frame == max(time_frame)),
                    aes(x = time_frame,
                        y = n_activity,
                        label = cell_type),
                    nudge_x = 50,
                    segment.alpha = 0.5,
                    size = 3,
                    max.overlaps = Inf,
                    min.segment.length = 0
    ) +
    # lims(x = c(min(df_merged$time_frame), max(df_merged$time_frame)*1.1)) +
    geom_line(aes(y = stim_timing) , linetype = "dashed", alpha = 0.3)
sX <- scale_x_continuous(name = "TimeFrame(1frame/0.2sec)",    # 軸の名前を変える
                         breaks = seq(0, length(timeframe), by= 1000),     # 軸の区切りを0,2,4にする
                        )
eval(parse(text=paste0("title <- ggtitle('SampleNumber",args_sample,"_",args_data,"')")))
t_1 <- theme(plot.title = element_text(size = 30, hjust = 0.5))
t_2 <- theme(axis.title = element_text(size = 20))
t_3 <- theme(legend.title = element_text(size = 28),
             legend.text = element_text(size = 20))

gg2 <- p_2 +
    sX +
    title +
    t_1 +
    t_2 +
    t_3 +
    labs(colour="each data")
########

#### ggplot other Data####
p_4 <- p_1 +        
         geom_line(aes(y = m_cherry, colour = "m_cherry") , size = 1.5)
p_5 <- p_1 +        
         geom_line(aes(y = position, colour = "position") , size = 1.5)

s_4 <- scale_color_manual(values = c("red"))
s_5 <- scale_color_manual(values = c("blue"))

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
##################################################

# ggsave
##################################################
gg <- gg2  + gg4 + gg5 + plot_layout(ncol = 1, heights = c(2, 1, 1))
ggsave(filename = args_output, plot = gg, dpi = 100, width = 30.0, height = 15.0)
##################################################