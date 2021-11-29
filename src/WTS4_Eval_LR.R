source("src/functions_WTS4_Eval_LR.R")

# #### args setting####
# args <- commandArgs(trailingOnly = T)
# # input merged_cls
# args_input_cls <- args[1]
# # output merged_data
# args_output_data <- args[2]

# # ReClustering Method
# args_method <- args[3]
# # Evaluation Method
# args_eval_method <- args[4]

#### test args####
# input merged_cls
args_input_cls <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/5_Clusters/MCMIHOOI/merged_cls.RData")
# output merged_data
args_output_data <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/5_Clusters/MCMIHOOI/ARI/eval_result.RData")
# ReClustering Method
args_method <- c("MCMIHOOI")
# Evaluation Method
args_eval_method <- c("ARI")

##### load####
load(args_input_cls)

#### create dataframe####
data.frame(CellType = names(merged_cls),
           Clusters = merged_cls,
           stringsAsFactors = FALSE,
           row.names = NULL
           ) %>% 
    dplyr::arrange(CellType) -> df_cls

df_cls %>% 
    mutate(Check_L = str_ends(.$CellType,"L")) %>% # add LRcheck column
    mutate(Check_R = str_ends(.$CellType,"R")) %>% 
    mutate(rm_LR = str_sub(.$CellType, start = 1, end = -2)) %>% 
        dplyr::filter(Check_L==TRUE | Check_R==TRUE) %>% # filter LR pair # ref. https://id.fnshr.info/2018/05/18/duplicated/
            group_by(rm_LR) %>% 
                filter(n()>1) %>% 
                    ungroup() -> df_cls_LR

#### Eval Dataframe####
df_cls_LR %>%
    select(1,5,2) -> df_cls_label

clusters <- df_cls_label$Clusters
classes <- df_cls_label$rm_LR

#### Evaluation####
eval_result <- switch(args_eval_method,
                      "ARI" = adjustedRandIndex(clusters, classes),
                      "purity" = ClusterPurity(clusters, classes),
                      "Fmeasure" = Fmeasure(clusters, classes),
                      stop("Only can use all, ARI, purity, Fmeasure")
                      )
# Entropyは小さいほど、クラスタリングがうまくいく。他の評価は値が高いほどお良い。他の評価と混ざるとわかりにくいので後回し。

#### ggsave####
save(eval_result, file=args_output_data)