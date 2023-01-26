source("src/functions_WTS5_tsPlot.R")

# #### args sample21 vsASER####
# args_sample <- c("21")
# args_input_n <- c("data/normalize_1/ReadData_21.RData")
# args_output <- c("output/WTS5/normalize_1/stimAfter/SBD_abs/tsPlot/SampleNumber_21_ASER_1000.eps")
# args_stim_xlsx <- c("data/stimulation/stimulation_timing.xlsx")
# args_input_stim <- c("data/stimulation/stim_21.RData")
# args_yshift <- c("output/WTS3/normalize_1/stimAfter/SBD_abs_manual/SampleNumber_21/yshift.RData")
# args_yshift_value <- c("output/WTS3/normalize_1/stimAfter/SBD_abs_manual/SampleNumber_21/yshift_value.RData")
# #### args sample2 vsASER####
# args_sample <- c("2")
# args_input_n <- c("data/normalize_1/ReadData_2.RData")
# args_stim_xlsx <- c("data/stimulation/stimulation_timing.xlsx")
# args_input_stim <- c("data/stimulation/stim_2.RData")
# args_yshift <- c("output/WTS3/normalize_1/stimAfter/SBD_abs_manual/SampleNumber_2/yshift.RData")
# args_yshift_value <- c("output/WTS3/normalize_1/stimAfter/SBD_abs_manual/SampleNumber_2/yshift_value.RData")
# 
# # args_output <- c("output/WTS5/normalize_1/stimAfter/SBD_abs/tsPlot/SampleNumber_2.eps")
# # args_output <- c("output/WTS5/normalize_1/stimAfter/SBD_abs/tsPlot/SampleNumber_2_BAGR.eps")
# args_output <- c("output/WTS5/normalize_1/stimAfter/SBD_abs/tsPlot/SampleNumber_2_ASEL.eps")

# #### args sample2 vs AVAR####
# args_sample <- c("2")
# args_input_n <- c("data/normalize_1/ReadData_2.RData")
# # args_output <- c("output/WTS5/normalize_1/stimAfter/SBD_abs/tsPlot/SampleNumber_2_AVAR_RIML.eps")
# args_output <- c("output/WTS5/normalize_1/stimAfter/SBD_abs/tsPlot/SampleNumber_2_AVAR_RIMR.eps")
# args_stim_xlsx <- c("data/stimulation/stimulation_timing.xlsx")
# args_input_stim <- c("data/stimulation/stim_2.RData")
# args_yshift <- c("output/WTS3/normalize_1/stimAfter/SBD_abs_manual/SampleNumber_2/yshift_AVAR.RData")
# args_yshift_value <- c("output/WTS3/normalize_1/stimAfter/SBD_abs_manual/SampleNumber_2/yshift_value_AVAR.RData")

#### args sample2 vs RIMR####
args_sample <- c("2")
args_input_n <- c("data/normalize_1/ReadData_2.RData")
# args_output <- c("output/WTS5/normalize_1/stimAfter/SBD_abs/tsPlot/SampleNumber_2_AVAR_RIML.eps")
args_output <- c("output/WTS5/normalize_1/stimAfter/SBD_abs/tsPlot/SampleNumber_2_RIMR_AVAR.eps")
args_stim_xlsx <- c("data/stimulation/stimulation_timing.xlsx")
args_input_stim <- c("data/stimulation/stim_2.RData")
args_yshift <- c("output/WTS3/normalize_1/stimAfter/SBD_abs_manual/SampleNumber_2/yshift_RIMR.RData")
args_yshift_value <- c("output/WTS3/normalize_1/stimAfter/SBD_abs_manual/SampleNumber_2/yshift_value_RIMR.RData")


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

# label_filter_list <- c("ASER","ASEL")
# label_filter_list <- c("ASER","BAGR")
# label_filter_list <- c("AVAR","RIML")
label_filter_list <- c("RIMR","AVAR")

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