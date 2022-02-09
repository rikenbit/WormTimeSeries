source("src/functions_WTS4_Eval_behavior.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
# input merged_cls
args_input_cls <- args[1]
# output merged_data
args_output_data <- args[2]
# ReClustering Method
args_method <- args[3]
# Evaluation Method
args_eval_method <- args[4]
# Evaluation label list
args_eval_label <- args[5]

# #### test args NMI####
# # input merged_cls
# args_input_cls <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/Merged_cls/k_Number_5.RData")
# # output merged_data
# args_output_data <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/Eval/NMI_behavior/k_Number_5.RData")
# # ReClustering Method
# args_method <- c("MCMIHOOI")
# # Evaluation Method
# args_eval_method <- c("NMI_behavior")
# # Evaluation label list
# args_eval_label <- c("data/WTS4_Eval_behavior_fix.xlsx")

##### load merged_cls####
load(args_input_cls)
data.frame(CellType = names(merged_cls),
           Clusters = merged_cls,
           stringsAsFactors = FALSE,
           row.names = NULL
           ) -> df_cls
##### load WTS4_Eval_behavior_fix.xlsx####
read.xlsx(args_eval_label,
          sheet = "Sheet1",
          rowNames = FALSE,
          colNames =TRUE) %>% 
    dplyr::rename(CellType = celltype, 
                  Classes = class) -> df_eval_label
#### merge dataframe####
df_merged <- merge(df_cls, 
                      df_eval_label, 
                      by.x = "CellType", 
                      by.y = "CellType", 
                      all.x = TRUE)
#### translation NA####
df_merged %>% 
    mutate_at(c("Classes"), 
              ~replace(., 
                       is.na(.), 
                       "others")
              ) -> df_cls_label
# チルダでpurrrの処理。各行読み込んではreplace関数を実行
# https://www.delftstack.com/ja/howto/r/replace-na-with-0-in-r/#r-データフレームのサブセット内の-na-を-0-に置き換える
# https://qiita.com/five-dots/items/361a42baf1e94edf5846

#### Eval Dataframe####
clusters <- df_cls_label$Clusters
classes <- df_cls_label$Classes
#### Evaluation####
eval_result <- switch(args_eval_method,
                      "ARI_behavior" = adjustedRandIndex(clusters, classes),
                      "purity_behavior" = ClusterPurity(clusters, classes),
                      "Fmeasure_behavior" = Fmeasure(clusters, classes),
                      "Entropy_behavior" = Entropy(clusters, classes), # Entropyは小さいほど、良い。他の評価は値が高いほどお良い。
                      "NMI_behavior" = NMI(clusters, classes),
                      stop("Only can use all, ARI, purity, Fmeasure, Entropy, NMI")
                      )

#### ggsave####
save(eval_result, file=args_output_data)
