source("src/functions_WTS4_Silhouette_Sample.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
# input merged_distance
args_params_distance <- args[1]
# input merged_cls
args_input_cls <- args[2]
# output eval_result
args_output_value <- args[3]
# output silhouette png
args_params_plot <- args[4]
# output silhouette gg object
args_output_gg <- args[5]

# #### test args####
# args_params_distance <- c("output/WTS4/normalize_1/stimAfter/EUCL/Distance")
# args_input_cls <-  c("output/WTS4/normalize_1/stimAfter/EUCL/Cluster_sample/k_Number_2/sample_cls.RData")
# args_output_value <- c("output/WTS4/normalize_1/stimAfter/EUCL/Eval_sample/Silhouette/k_Number_2.RData")
# # args_params_plot <- c("output/WTS4/normalize_1/stimAfter/EUCL/DimReduc_sample/k_Number_2/Sil_plot")
# # args_output_gg <- c("output/WTS4/normalize_1/stimAfter/EUCL/DimReduc_sample/k_Number_2/Sil_gg.RData")
# args_params_plot <- c("output/WTS4/normalize_1/stimAfter/EUCL/DimReduc_sample/k_Number_2/Sil_plot")
# args_output_gg <- c("output/WTS4/normalize_1/stimAfter/EUCL/DimReduc_sample/k_Number_2/Sil_plot/Sil_gg.RData")

#### fix sample number sort####
input_path_list <- list.files(args_params_distance, pattern="SampleNumber_", full.names=TRUE)
# 24サンプルの順にサンプル番号並び替え
input_path_list %>% 
  str_remove(., args_params_distance) %>% 
  str_remove(., "/SampleNumber_") %>% 
  str_remove(., ".RData") %>% 
  as.numeric() %>% 
  sort() -> sample_sort_num

#### load cls####
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

#### load distance####
sample_sort_num %>% 
  purrr::map(., .list_distance) -> D

#### save eval_result####
seq(1:length(D)) %>% 
  purrr::map_dbl(., .list_eval_result) -> eval_result
save(eval_result, file=args_output_value)

#### save ggplot####
seq(1:length(D)) %>% 
  purrr::map(., .list_gg_sil) -> gg_sil


#### save png####
# SampleNumber_*.pngで保存
for(i in seq(1:length(D))){
  if (is.na(gg_sil[[i]])) {
    # シルエット係数がNAなら何もしない
  } else {
  sample_num <- sample_sort_num[i]
  eval(parse(text=paste0("args_output_png <- c('",args_params_plot,"/SampleNumber_",sample_num,".png')")))
  gg <- gg_sil[[i]]
  #### ggsave fviz_silhouette####
  ggsave(filename = args_output_png, 
         plot = gg,
         dpi = 100, 
         width = 30.0, 
         height = 20.0,
         limitsize = FALSE)
  }
}

#### save ggplot####
save(gg_sil, file=args_output_gg)