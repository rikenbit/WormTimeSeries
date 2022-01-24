source("src/functions_WTS4_Eval_sample_nolabel.R")

#### args setting####
# args <- commandArgs(trailingOnly = T)
args_input_sample <- args[1]
args_input_cls <- args[2]
args_time <- args[3]
args_k <- args[4]
args_eval_method <- args[5]
args_input_path <- args[6]
args_stim_xlsx <- args[7]
args_output <- args[8]

# #### test args####
# # sample matrix data
# args_input_sample <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Distance")
# # sample cluster
# args_input_cls <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Cluster_sample/k_Number_5/sample_cls.RData")

# # args_time
# args_time <- c("stimAfter")
# # k
# args_k <- c("5")
# # eval method
# args_eval_method <- c("kNN")

# # param path
# args_input_path <- c("data/normalize_1")
# # param stimtiming
# args_stim_xlsx <- c("data/stimulation/stimulation_timing.xlsx")

# # args_output
# args_output <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Eval_sample/kNN_nolabel/k_number_5.RData")

#### No. of Clusters####
k <- as.numeric(args_k)
########

#### load sample cls####
# C
load(args_input_cls)
lapply(C, function(x) {
    data.frame(CellType = attr(x, "names"),
               Clusters = as.numeric(x),
               stringsAsFactors = FALSE,
               row.names = NULL
               ) -> df_cls
    cellnames <- df_cls$CellType
    # ここで数字除去するかは微妙なところ
    df_cls[grep("^[0-9]", cellnames, invert=TRUE),]
    }
) -> df_cls_list

#### fix sample number sort####
input_path_list <- list.files(args_input_sample, pattern="SampleNumber_", full.names=TRUE)
# 24サンプルの順にサンプル番号並び替え
input_path_list %>% 
    str_remove(., args_input_sample) %>% 
    str_remove(., "/SampleNumber_") %>% 
    str_remove(., ".RData") %>% 
    as.numeric() %>% 
    sort() -> sample_sort_num

input_path <- c()
for(i in sample_sort_num){
  eval(parse(text=paste0("path <- c('",args_input_path,"/ReadData_",i,".RData')")))
  input_path <- c(input_path, path)
}
input_path_list <- input_path

#### load sample matrix data####
# inputファイル名のリスト
sample_data_list <- list()
for(i in 1:length(input_path_list)){
    # load ReadData_1
    load(input_path_list[i])
    args_sample <- sample_sort_num[i]
    # ReadData <- ReadData_1
    eval(parse(text=paste0("ReadData <- ReadData_", args_sample)))
    #### switch time range & trim ReadData####
    sample_data <- switch(args_time,
                   "all" = .ReadData_all(ReadData),
                   "stimAfter" = .ReadData_stimAfter(ReadData, args_stim_xlsx),
                   stop("Only can use all, stimAfter ")
                   )
    # 数字の細胞列の除去
    sample_data_rm <- sample_data[, grep("^[0-9]", colnames(sample_data), invert=TRUE)]
    # 行と列の入れ替え
    sample_data_rm <- t(as.matrix(sample_data_rm))
    # test code
    # sample_data_rm <- sample_data_rm[,1:100]
    sample_data_list <- c(sample_data_list, list(sample_data_rm))
}
####


#### Evaluation####
# PseudoF(sample_data_list[[1]], df_cls_list[[1]]$Clusters)
# input 「sample_data_list」, 「df_cls_list」
eval_result <- switch(args_eval_method,
                    # eval method
                    "PseudoF" = purrr::map_dbl(seq(1:length(sample_data_list)), .PseudoF_list),
                    "Connectivity" = purrr::map_dbl(seq(1:length(sample_data_list)), .Connectivity_list),
                    "kNN" = purrr::map_dbl(seq(1:length(sample_data_list)), .kNN_list),
                    stop("Only can use all, PseudoF, Connectivity,kNN")
                    )

save(eval_result, file=args_output)