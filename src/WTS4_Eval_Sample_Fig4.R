source("src/functions_WTS4_Eval_Sample_Fig4.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
args_weight_path <- args[1]
args_sample_path <- args[2]
args_input_cls <- args[3]
args_output_k <- args[4]
args_k <- args[5]

# # #### test args####
# # args_weight_path <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/Merged_data")
# # args_sample_path <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Distance")
# # args_input_cls <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Cluster_sample/k_Number_9/sample_cls.RData")
# # args_output <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/DimReduc_sample/k_Number_9/Eval_Sample_Fig4.png")
# # args_output_k <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/DimReduc_sample/k_Number_9/Eval_Sample_Fig4_k9.png")
# # args_k <- c("9")

# args_weight_path <- c("output/WTS4/n1_28sample/stimAfter/SBD_abs/MCMIHOOI/Merged_data")
# args_sample_path <- c("output/WTS4/n1_28sample/stimAfter/SBD_abs/Distance")
# args_input_cls <- c("output/WTS4/n1_28sample/stimAfter/SBD_abs/Cluster_sample/k_Number_9/sample_cls.RData")
# args_output_k <- c("output/WTS4/n1_28sample/stimAfter/SBD_abs/DimReduc_sample/k_Number_9/Eval_Sample_Fig4_k.png")
# args_k <- c("9")

#### No. of Clusters####
k <- as.numeric(args_k)
########

#### sample number####
sample_path_list <- list.files(args_sample_path, pattern="SampleNumber_", full.names=TRUE)
sample_path_list %>% 
    str_remove(., args_sample_path) %>% 
    str_remove(., "/SampleNumber_") %>% 
    str_remove(., ".RData") %>% 
    as.numeric() %>% 
    sort() -> sample_sort_num

#### cluster number####
weight_path_list <- list.files(args_weight_path, pattern="k_Number_", full.names=TRUE)
weight_path_list %>% 
    str_remove(., args_weight_path) %>% 
    str_remove(., "/k_Number_") %>% 
    str_remove(., ".RData") %>% 
    as.numeric() %>% 
    sort() -> cluster_sort_num

#### MCMI weight table all cluster####
# purrr
cluster_sort_num %>% 
    purrr::map_dfr(., .all_k_table) -> df_weight_all

#### df_weight_group####
df_weight_all %>% 
    group_by(SampleNumber) %>% 
    summarise(weight_ave = mean(weight_abs)) -> df_weight_group

#### weight_sort####
df_weight_group %>% 
    dplyr::arrange(desc(weight_ave)) %>% 
    .$SampleNumber %>% 
    c() -> sample_sort_weight

##### load merged_cls####
load(args_input_cls)
#### add annotated count####
#数字を除く
# annotated_count <- attr(C[[1]], "names") %>% grep("^[0-9]", ., invert=TRUE)  %>% length()
annotated_count <- unlist(lapply(C, function(x){length(grep("^[0-9]", attr(x, "names"), invert=TRUE))}))
df_eval_wide <- data.frame(SampleNumber = as.character(sample_sort_num),
                           annotated_count = annotated_count,
                           stringsAsFactors = FALSE,
                           row.names = NULL
                           )

# gg <- ggplot(df_eval_wide, aes(x = SampleNumber, y= annotated_count, group=1))
# gg <- gg + scale_x_discrete(limits=sample_sort_weight)
# gg <- gg + geom_line(color = "black", size= 2)
# gg <- gg + geom_smooth(method="lm", size =0.5, se = TRUE, alpha = 0.4, color = "black")
# gg <- gg + labs(x = "Sample No.", y = "No. of annotated cells") + theme(text = element_text(size = 60))

# #### ggsave####
# ggsave(filename = args_output, 
#        plot = gg,
#        dpi = 100, 
#        width = 20.0, 
#        height = 20.0,
#        limitsize = FALSE)

#### sort one k#####
df_weight_all %>% 
  filter(Cluster==k) %>%
  dplyr::arrange(desc(weight_abs)) %>% 
  .$SampleNumber %>% 
  c() -> sample_sort_weight_k

gg_k <- ggplot(df_eval_wide, aes(x = SampleNumber, y= annotated_count, group=1))
gg_k <- gg_k + scale_x_discrete(limits=sample_sort_weight_k)
gg_k <- gg_k + geom_line(color = "black", size= 2)
gg_k <- gg_k + geom_smooth(method="lm", size =0.5, se = TRUE, alpha = 0.4, color = "black")
gg_k <- gg_k + labs(x = "Sample No.", y = "No. of identified cells") + theme(text = element_text(size = 60))

### ggsave####
ggsave(filename = args_output_k, 
       plot = gg_k,
       dpi = 100, 
       # width = 30.0, 
       # height = 30.0,
       width = 20.0, 
       height = 20.0,
       limitsize = FALSE)