source("src/functions_WTS5_tsPlot.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
args_sample <- args[1]
args_celltype <- args[2]
args_s_celltype <- args[3]
args_input_n <- args[4]
args_input_stim <- args[5]
args_yshift <- args[6]
args_yshift_value <- args[7]
args_stim_xlsx <- args[8]
args_output <- args[9]

#### load Neuron Activity Data####
load(args_input_n)
eval(parse(text=paste0("input_n <- ReadData_", args_sample)))
#### switch time range & trim ReadData####
input_n <- .ReadData_stimAfter(input_n, args_stim_xlsx)
input_n %>%
    rownames_to_column("time_frame") %>%
        pivot_longer(-time_frame,
                     names_to = "cell_type",
                     values_to = "n_activity") -> df_input_n
df_input_n$time_frame <- as.numeric(df_input_n$time_frame)
### load other Data####
# TimeFrame
rownames(input_n) %>%
    as.numeric() -> timeframe
# Stimulation Data
load(args_input_stim)
eval(parse(text=paste0("input_stim <- stim_",args_sample)))
input_stim %>%
    as.numeric() -> stimtiming
stimtiming <- stimtiming[timeframe]
df_input_other <- data.frame(time_frame = timeframe,
                             stim_timing = stimtiming,
                             stringsAsFactors = FALSE)
#### merge input####
df_input_n %>%
    merge(.,
          df_input_other,
          by.x = "time_frame",
          by.y = "time_frame",
          all.x = TRUE) -> df_input
#### load yshift Neuron Activity Data#####
# sbd_yshift_df 神経活動値（yshift後）
load(args_yshift)
#### merge input , yshift Neuron Activity Data####
df_input %>%
    merge(.,
          sbd_yshift_df,
          by.x = c("time_frame", "cell_type"),
          by.y = c("time_frame", "cell_type"),
          all.x = TRUE) -> df_tsPlot

#### load yshift filter####
# yshift_value_table y_shift値
load(args_yshift_value)

#### merge all####
df_tsPlot %>%
    merge(.,
          yshift_value_table,
          by.x = "cell_type",
          by.y = "cell_type",
          all.x = TRUE) -> df_tsPlot
df_tsPlot$time_frame <- as.numeric(df_tsPlot$time_frame)
#### ggtheme####
sX <- scale_x_continuous(name = "TimeFrame(1frame/0.2sec)",
                         breaks = seq(0,
                                      # length(timeframe),
                                      6000,
                                      by= 1000)
                         )
t_1 <- theme(plot.title = element_text(size = 30, hjust = 0.5, family ="HiraKakuPro-W3"))
t_2 <- theme(axis.title = element_text(size = 20))
t_3 <- theme(legend.title = element_text(size = 28),
             legend.text = element_text(size = 20)
)

label_filter_list <-c(args_celltype,
                      args_s_celltype)

seq(1:length(label_filter_list)) %>%
    purrr::map(., .plot_yshift) -> gg_cells

#### wrap_plots####
gg_cells %>%
    wrap_plots(., ncol = 1) +
    plot_annotation(title = "",
                    theme = theme(plot.title = element_text(size = 48, hjust = 0.5))
    ) -> gg
#### ggsave####
ggsave(filename = args_output,
       plot =gg,
       dpi = 100,
       width = 50.0,
       height = 20.0,
       limitsize = FALSE)