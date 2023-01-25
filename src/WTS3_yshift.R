source("src/functions_WTS3_yshift.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
# sample number
args_sample <- args[1]
# input SBD yshift df
args_yshift <- args[2]
# input stim timing extra
args_stim_xlsx <- args[3]
# output SBD yshift df
args_yshift_value <- args[4]

# #### test args####
# # sample number
# args_sample <- c("1")
# # input SBD yshift df
# args_yshift <- c("output/WTS3/normalize_1/stimAfter/SBD/SampleNumber_1/yshift.RData")
# # input stim timing extra
# args_stim_xlsx <- c("data/stimulation/stimulation_timing.xlsx")
# # output SBD yshift df
# args_yshift_value <- c("output/WTS3/normalize_1/stimAfter/SBD/SampleNumber_1/yshift_value.RData")

# ### test args####
# # sample number
# args_sample <- c("2")
# # input SBD yshift df
# args_yshift <- c("output/WTS3/normalize_1/stimAfter/SBD_abs_manual/SampleNumber_2/yshift_AVAR.RData")
# # input stim timing extra
# args_stim_xlsx <- c("data/stimulation/stimulation_timing.xlsx")
# # output SBD yshift df
# args_yshift_value <- c("output/WTS3/normalize_1/stimAfter/SBD_abs_manual/SampleNumber_2/yshift_value_AVAR.RData")

#### load y-shift####
load(args_yshift)
as.numeric(sbd_yshift_df$time_frame) -> sbd_yshift_df$time_frame
#### load stim timing extra####
stimtimng_sheet <- read.xlsx(args_stim_xlsx,
                             sheet = "Sheet1",
                             rowNames = FALSE,
                             colNames =TRUE)
#### get shift max####
stimtimng_sheet %>%
    dplyr::select(sample_number = 1,
                  frame_sec = 6,
                  stim_first = 7,
                  half_period = 8
                  ) %>%
        mutate(period = half_period*2,
               shift_max = half_period*2 + 50,
               shift_min = -half_period*2 - 50
               ) -> stimtimng_sheet
stimtimng_sheet %>%
    filter(sample_number == as.numeric(args_sample)) %>%
        .$shift_max -> yshift_max

#### cell_typeのリストゲット####
sbd_yshift_df %>%
    .$cell_type %>%
        unique() -> list_cell_type
#### cell_typeごとのy-shift値計算####
seq(1:length(list_cell_type)) %>%
    purrr::map(., .df_yshift) -> yshift_list
seq(1:length(list_cell_type)) %>%
    purrr::map_dbl(., .yshift_value) -> yshift_value
#### cell_typeごとのy-shift値の絶対値計算####
yshift_value %>%
    abs() -> yshift_abs
#### cell_typeごとのy-shift値絶対値がyshift_max以下かどうか判定####
seq(1:length(list_cell_type)) %>%
    purrr::map_dbl(., .yshift_filter) -> yshift_filter
#### create y-shift_value table####
data.frame(
    cell_type = list_cell_type,
    yshift_value = yshift_value,
    yshift_abs = yshift_abs,
    yshift_filter = yshift_filter,
    stringsAsFactors = FALSE
    ) -> yshift_value_table
#### save yshift_table####
save(yshift_value_table, file=args_yshift_value)