source("src/functions_WTS3_tsPlot.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
# select animal number 個体番号の指定
args_sample <- args[1]
# input_Neuron Activity
args_input_n <- args[2]
# input_stim
args_input_stim <- args[3]
# input_mCherry
args_input_mCherry <- args[4]
# input_Position
args_input_Position <- args[5]

# input SBD yshift Neuron
args_yshift <- args[6]
# input yshift filter
args_yshift_value <- args[7]
# input label filter df_label
args_label_table <- args[8]
# input select label
args_label <- args[9]
# input select shift cell_type
args_shift <- args[10]
# output tsPlot
args_output <- args[11]

# # stimtiming
args_stim_xlsx <- args[12]
# time range
args_time <- args[13]



#### load Neuron Activity Data####
load(args_input_n)
eval(parse(text=paste0("input_n <- ReadData_",args_sample)))
#### switch time range & trim ReadData####
input_n <- switch(args_time,
                   "all" = .ReadData_all(input_n),
                   "stimAfter" = .ReadData_stimAfter(input_n, args_stim_xlsx),
                   stop("Only can use all, stimAfter ")
                   )
input_n %>%
    rownames_to_column("time_frame") %>%
        pivot_longer(-time_frame,
                     names_to = "cell_type",
                     values_to = "n_activity") -> df_input_n
as.numeric(df_input_n$time_frame) -> df_input_n$time_frame
#### load other Data####
# TimeFrame
rownames(input_n) %>%
    as.numeric() -> timeframe
# Stimulation Data
load(args_input_stim)
eval(parse(text=paste0("input_stim <- stim_",args_sample)))
input_stim %>%
    as.numeric() -> stimtiming
stimtiming[timeframe] -> stimtiming
# mCherry
load(args_input_mCherry)
eval(parse(text=paste0("input_mCherry <- mCherry_",args_sample)))
input_mCherry[,"ASER"] %>%
    as.numeric() -> mcherry
mcherry[timeframe] -> mcherry
# Position
load(args_input_Position)
eval(parse(text=paste0("input_Position <- Position_",args_sample)))
input_Position$MoveX %>%
    as.numeric() -> position
position[timeframe] -> position
data.frame(
    time_frame = timeframe,
    stim_timing = stimtiming,
    m_cherry = mcherry,
    position = position,
    stringsAsFactors = FALSE
    ) -> df_input_other
#### load yshift Neuron Activity Data#####
load(args_yshift)
#### load yshift filter####
load(args_yshift_value)
#### load label####
load(args_label_table)
#### prepare ggplot table####
#### merge input####
df_input_n %>%
    merge(.,
        df_input_other,
        by.x = "time_frame",
        by.y = "time_frame",
        all.x = TRUE) -> df_input
#### merge input , yshift Neuron Activity Data####
df_input %>%
  merge(.,
        sbd_yshift_df,
        by.x = c("time_frame", "cell_type"),
        by.y = c("time_frame", "cell_type"),
        all.x = TRUE) -> df_tsPlot
#### merge yshift filter, label####
yshift_value_table %>%
  merge(.,
        df_label,
        by.x = "cell_type",
        by.y = "cell_type",
        all.x = TRUE) -> df_label_filter
#### merge all####
df_tsPlot %>%
  merge(.,
        df_label_filter,
        by.x = "cell_type",
        by.y = "cell_type",
        all.x = TRUE) -> df_tsPlot
df_tsPlot$time_frame <- as.numeric(df_tsPlot$time_frame)

#### ggtheme####
sX <- scale_x_continuous(name = "TimeFrame(1frame/0.2sec)",
                         breaks = seq(0, length(timeframe), by= 1000)
                         )
t_1 <- theme(plot.title = element_text(size = 30, hjust = 0.5, family ="HiraKakuPro-W3"))
t_2 <- theme(axis.title = element_text(size = 20))
t_3 <- theme(legend.title = element_text(size = 28),
             legend.text = element_text(size = 20)
             )

#### ggplot labeled cell####
# select label
label_filter_df <- switch(args_label,
          "label_acf" = filter(df_label_filter, label_acf==1, yshift_filter==1),
          "label_cls" = filter(df_label_filter, label_cls==1, yshift_filter==1),
          stop("Only can use label_acf,label_cls")
          )
# sort yshift_abs
label_filter_df %>%
    dplyr::arrange(yshift_abs) %>%
        .$cell_type -> label_filter_list
match(args_shift, label_filter_list) -> yshift_order
c(label_filter_list[yshift_order], label_filter_list[-yshift_order]) -> label_filter_list
# plot cells
seq(1:length(label_filter_list)) %>%
    purrr::map(., .plot_yshift) -> gg_cells

#### plot other data####
p_1 <- ggplot(data = df_input_other,
              aes(x = time_frame))
gg_m <- p_1 +
    geom_line(aes(y = m_cherry, colour = "m_cherry")) +
    scale_color_manual(values = c("red")) +
    t_2 +
    t_3 +
    sX +
    t_1 +
    ggtitle(args_shift)
gg_p <- p_1 +
    geom_line(aes(y = position, colour = "position")) +
    scale_color_manual(values = c("blue")) +
    t_2 +
    t_3 +
    sX
#### wrap_plots####
append(gg_cells, list(gg_m)) %>%
    append(., list(gg_p)) -> gg_list
# wrap
eval(parse(text=paste0("plot_title <- c('SampleNumber_",args_sample,"_",args_label,"')")))
gg_list %>%
    wrap_plots(., ncol = 1) +
    plot_annotation(title = plot_title,
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