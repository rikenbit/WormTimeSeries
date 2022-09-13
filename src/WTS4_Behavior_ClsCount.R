source("src/functions_WTS4_Behavior_ClsCount.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
args_input_cls <- args[1]
args_output_all <- args[2]
args_output_sum <- args[3]
args_k <- args[4]
args_eval_label <- args[5]
args_sample_path <- args[6]
  
# #### test args####
# # No. of Clusters
# args_input_cls <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Cluster_sample/k_Number_6/sample_cls.RData")
# args_output_all <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/ClsCount/k_Number_6/df_count_all.RData")
# args_output_sum <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/ClsCount/k_Number_6/df_count_sum.RData")
# 
# args_k <- c("6")
# 
# args_eval_label <- c("data/WTS4_Eval_behavior_ACF.xlsx")
# args_sample_path <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Distance")

#### No. of Clusters####
k <- as.numeric(args_k)

#### sample number####
sample_path_list <- list.files(args_sample_path, pattern="SampleNumber_", full.names=TRUE)
sample_path_list %>% 
  str_remove(., args_sample_path) %>% 
  str_remove(., "/SampleNumber_") %>% 
  str_remove(., ".RData") %>% 
  as.numeric() %>% 
  sort() -> sample_sort_num

##### load WTS4_Eval_behavior.xlsx####
read.xlsx(args_eval_label,
          sheet = "Sheet1",
          rowNames = FALSE,
          colNames =TRUE) %>% 
  dplyr::rename(CellType = celltype, 
                Classes = class) -> df_eval_label

#### laod sample cls result####
# クラスタリング結果のリスト
load(args_input_cls)
# Cがloadされる

#### 各個体の行動ラベル回数データフレームをlistで保存####
# purrr::map でlistにしてsave
seq(1:length(C)) %>%
  purrr::map(., .df_count) -> df_count_list

#### 全個体の累積行動ラベル回数データフレーム作成####
# purrr::map_dfr でrbind
seq(1:length(df_count_list)) %>% 
    purrr::map_dfr(., .df_count_rbind) -> df_count_all
#### save####
save(df_count_all, file=args_output_all)

#### groupby count_sum####
df_count_all %>% 
  group_by(CellType) %>% 
  summarise(Count_sum = sum(Count)) -> df_count_sum

#### save####
save(df_count_sum, file=args_output_sum)