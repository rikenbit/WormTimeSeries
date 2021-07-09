source("src/functions_WTS3_label_table.R")

#### args setting####
# select animal number 個体番号一覧取得
args_numbers <- c("1","2","5","6","7","10","11","12","13","14","17","18","19","21","23","24","27")

args <- commandArgs(trailingOnly = T)
# input df_label
args_label_dir <- args[1]
# input df_label
args_yshift_value_dir <- args[2]
# input select label
args_label <- args[3]
# output table stim_cell
args_output_neuron <- args[4]
args_output_rmSensory <- args[5]
# output yshift
args_output_yshift <- args[6]

# #### test args####
# # select animal number 個体番号一覧取得
# args_numbers <- c("1","2","5","6","10","11","12","13","14","17","18","19","21","23","24","27")
# # input df_label
# args_label_dir <- c("output/WTS3/normalize_1/all/SBD/ARI/SampleNumber_")
# # input df_label
# args_yshift_value_dir <- c("output/WTS3/normalize_1/all/SBD/SampleNumber_")
# # input select label
# args_label <- c("label_acf")
# # output table stim_cell
# args_output_neuron <- c("output/WTS3/normalize_1/all/SBD/ARI/table/label_acf/table_neuron.csv")
# args_output_rmSensory <- c("output/WTS3/normalize_1/all/SBD/ARI/table/label_acf/table_rmSensory.csv")
# # output yshift
# args_output_yshift <- c("output/WTS3/normalize_1/all/SBD/ARI/table/label_acf/table_yshift.csv")

#### prepare label_table####
# 16サンプルのdfのリストを作成
seq(1,length(args_numbers)) %>% 
    purrr::map(., .load_label_table) -> df_label_list
# bind 16 sample list
plyr::join_all(df_label_list, type = 'full') -> df_table

#### prepare yshift_value####
# 16サンプルのdfのリストを作成
seq(1,length(args_numbers)) %>% 
  purrr::map(., .load_yshift_value) -> df_yshift_value_list
# bind 16 sample list
plyr::join_all(df_yshift_value_list, type = 'full') -> df_yshift_table

#### create 塩刺激 table####
output_table <- switch(args_label,
                       "label_acf" = .table_acf(df_table, df_yshift_table),
                       "label_cls" = .table_cls(df_table, df_yshift_table),
                       stop("Only can use label_acf,label_cls")
)
#### 加工####
output_table[[1]] %>% 
    .add_sum_sort() -> output_neuron
output_table[[2]] %>% 
    .add_sum_sort() -> output_rmSensory
#### save 塩刺激 write_excel_csv####
write_excel_csv(output_neuron, args_output_neuron) 
write_excel_csv(output_rmSensory, args_output_rmSensory)

#### create yshift table####
output_yshift_table <- switch(args_label,
                       "label_acf" = .table_acf_yshift(df_table, df_yshift_table),
                       "label_cls" = .table_cls_yshift(df_table, df_yshift_table),
                       stop("Only can use label_acf,label_cls")
                       )
output_yshift_table %>% 
    as.data.frame() %>% 
    column_to_rownames(var = "cell_type") -> output_yshift_table
output_yshift_table[,order(as.numeric(colnames(output_yshift_table)))] -> output_yshift_table
#### sort inorm####
seq(1, nrow(output_yshift_table)) %>% 
  purrr::map_dbl(., .yshift_norm) %>% 
      trunc()  -> yshift_i_norms
output_yshift_table %>%
    mutate(i_norm = yshift_i_norms) %>%
        dplyr::arrange(i_norm) -> output_yshift_table_sort
output_yshift_table_sort %>% 
    rownames_to_column("cell_type") -> output_yshift_table_sort
#### save yshift write_excel_csv####
write_excel_csv(output_yshift_table_sort, args_output_yshift) 