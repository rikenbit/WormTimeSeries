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
# # フィルタリング
# args_filter <- args[10]
# ##################################################
# #### test args####
# args_sample <- c("1")
# # input_Neuron Activity ファイル名
# args_input_n <- c("data/normalize_1/ReadData_1.RData")
# # input_stim ファイル名
# args_input_stim <- c("data/stimulation/stim_1.RData")
# # input_mCherry ファイル名
# args_input_mCherry <- c("data/mCherry/mCherry_1.RData")
# # input_Position ファイル名
# args_input_Position <- c("data/Position/Position_1.RData")
# # input_tempdat ファイル名
# args_input_tempdata <- c("output/WTS3/SBD/normalize_1/all/tsne/ARI/cls_tempdata/SampleNumber_1.RData")
# # outputファイル名
# args_output <- c("output/WTS3/SBD/normalize_1/all/tsne/ARI/plot/SampleNumber_1.png")
# # select data データの指定
# args_data <- c("normalize_1")
# # クラスター評価手法
# args_eval <- c("ARI")
# # 次元圧縮手法
# args_DimRedu <- c("tsne")
# # フィルタリング
# args_filter <- c("stim_cell")
# #######################
#### test args####
args_sample <- c("2")
# input_Neuron Activity ファイル名
args_input_n <- c("data/normalize_1/ReadData_2.RData")
# input_stim ファイル名
args_input_stim <- c("data/stimulation/stim_2.RData")
# input_mCherry ファイル名
args_input_mCherry <- c("data/mCherry/mCherry_2.RData")
# input_Position ファイル名
args_input_Position <- c("data/Position/Position_2.RData")
# input_tempdat ファイル名
args_input_tempdata <- c("output/WTS3/SBD/normalize_1/all/tsne/ARI/cls_tempdata/SampleNumber_2.RData")
# outputファイル名
args_output <- c("output/WTS3/SBD/normalize_1/all/tsne/ARI/plot/SampleNumber_2.png")
# select data データの指定
args_data <- c("normalize_1")
# クラスター評価手法
args_eval <- c("ARI")
# 次元圧縮手法
args_DimRedu <- c("tsne")
# フィルタリング
args_filter <- c("stim_cell")
# y-shift計算対象の細胞
args_shift <- c("ASER")
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
df_merged <- switch(args_filter,
              "stim_cell" = filter_stim(df_merged_temp),
              stop("Only can use stim_cell")
)
# List of filtered cells 
df_merged %>%
    .$cell_type %>%
        unique() -> list_cell_type
#### check ASER or BAGR or BAGL####
list_cell_type %>% 
    str_count(., pattern="ASER") %>%
        sum() -> check_ASER
list_cell_type %>% 
    str_count(., pattern="BAGR") %>%
        sum() -> check_BAGR
list_cell_type %>% 
    str_count(., pattern="BAGL") %>%
        sum() -> check_BAGL
if (check_ASER >= 1) {
    args_shift <- "ASER"
} else if (check_BAGR >= 1) {
    args_shift <- "BAGR"
} else if (check_BAGL >= 1) {
    args_shift <- "BAGL"
} else {
    args_shift <- list_cell_type[1]
}
#### WTS3_yshift#####
# 行列っぽいデータを細胞ごとにlist化
input_n.list <- asplit(input_n, 2)
# prepare shift_1
eval(parse(text=paste0("shift_1 <- input_n.list$",args_shift," %>% as.numeric()")))
# sbd y-shift
list_cell_type %>% 
    purrr::map(., sbd_y) %>%
        as.data.frame() -> sbd_yshift_df_wide
colnames(sbd_yshift_df_wide) <- list_cell_type
# convert long df
sbd_yshift_df_wide %>% 
    rownames_to_column("time_frame") %>% 
    pivot_longer(-time_frame, 
                 names_to = "cell_type", 
                 values_to = "y_shift") -> sbd_yshift_df
#### merge yshift#####
merge(x=df_merged, 
      y=sbd_yshift_df, 
      by.x=c("time_frame", "cell_type"), 
      by.y=c("time_frame", "cell_type"),
      all.x = TRUE) -> df_merged_yshift
df_merged_yshift %>% 
    mutate(., 
           shifted_TF = time_frame + y_shift) -> df_merged_TF
#### ggplot test####
# ggplot theme
sX <- scale_x_continuous(name = "TimeFrame(1frame/0.2sec)",
                         breaks = seq(0, length(timeframe), by= 1000)
                         )
t_1 <- theme(plot.title = element_text(size = 30, hjust = 0.5))
t_2 <- theme(axis.title = element_text(size = 20))
t_3 <- theme(legend.title = element_text(size = 28),
             legend.text = element_text(size = 20))

# plot Neuron Activity Data
seq(1:length(list_cell_type)) %>% 
    purrr::map(., plot_one_cell) -> gg_cells

# plot other data
p_1 <- ggplot(data = df_merged,
              aes(x = time_frame))
gg_m <- p_1 +        
        geom_line(aes(y = m_cherry, colour = "m_cherry")) +
        scale_color_manual(values = c("red")) +
        t_2 +
        t_3 +
        sX
gg_p <- p_1 +
        geom_line(aes(y = position, colour = "position")) +
        scale_color_manual(values = c("blue")) +
        t_2 +
        t_3 +
        sX

#### wrap_plots####
# list n_activity & other data
append(gg_cells, list(gg_m)) %>% 
    append(., list(gg_p)) -> gg_list
# wrap
eval(parse(text=paste0("plot_title <- c('SampleNumber_",args_sample,"_",args_data,"')")))
gg_list %>% 
    wrap_plots(., ncol = 1) +
    plot_annotation(
        title = plot_title,
        caption = 'made with patchwork::wrap_plots',
        theme = theme(plot.title = element_text(size = 48, hjust = 0.5))
    ) -> gg
#### ggsave####
ggsave(filename = args_output, 
       plot = gg, 
       dpi = 100, 
       width = 25.0, 
       height = 50.0, 
       limitsize = FALSE)