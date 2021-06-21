source("src/functions_WTS3_plot.R")

### args setting###
args <- commandArgs(trailingOnly = T)
# select animal number 個体番号一覧取得
args_numbers <- c("1","2","4","5","6","7","9","10","11","12","13","14","15","16","17","18","19","21","22","23","24","26","27","28")
# y-shift計算対象の細胞
args_shift <- c("ASER")
# input_dir
args_tempdata_dir <- args[1]
# ラベリング
args_filter <- args[2]
# output table stim_cell
args_output_neuron <- args[3]
# output table stim_cell
args_output_rmSensory <- args[4]
# #### args####
# # サンプル番号一覧取得
# args_numbers <- c("1","2","4","5","6","7","9","10","11","12","13","14","15","16","17","18","19","21","22","23","24","26","27","28")
# # y-shift計算対象の細胞
# args_shift <- c("ASER")
# # input
# args_tempdata_dir <- c("output/WTS3/normalize_1/all/SBD/ARI/tsne/cls_tempdata/SampleNumber_")
# # ラベリング
# args_filter <- c("stim_cell")
# # output table stim_cell
# args_output_neuron <- c("output/WTS3/normalize_1/all/SBD/ARI/tsne/table/stim_cell/table_neuron.csv")
# args_output_rmSensory <- c("output/WTS3/normalize_1/all/SBD/ARI/tsne/table/stim_cell/table_rmSensory.csv")
# ########

#### prepare dataframe####
# df_temp_list 24サンプルのdfのリストを作成
seq(1,length(args_numbers)) %>% 
    purrr::map(.,load_tempdata) -> df_temp_list
# bind 24 sample list
plyr::join_all(df_temp_list, type = 'full') -> df_table
#### create table####
output_table <- switch(args_filter,
                    "stim_cell" = table_cell(df_table),
                    "stim_cluster" = table_cls(df_table),
                    stop("Only can use stim_cell,stim_cluster")
)
# write_excel_csv
write_excel_csv(output_table[[1]], args_output_neuron) 
write_excel_csv(output_table[[2]], args_output_rmSensory)
