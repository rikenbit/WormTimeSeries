source("src/functions_WTS4_membership_df.R")

#### args setting####
#### test args####
args_sample_path <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Distance")
args_input <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Membership_F/k_Number_6.RData")

load(args_input)

# リスト化につける名前
sample_path_list <- list.files(args_sample_path, pattern="SampleNumber_", full.names=TRUE)
sample_path_list %>% 
    str_remove(., args_sample_path) %>% 
    str_remove(., "/SampleNumber_") %>% 
    str_remove(., ".RData") %>% 
    as.numeric() %>% 
    sort() -> sample_sort_num


.animal_df_list= function(x) {
    names(x) <- c("animal_1")
    return(return_object)
}

lapply(newHs, .animal_df_list)

# listだけエクスポート
names(mem_df_list) <- as.character(sample_sort_num)
save(mem_df_list, file=args_output_list)

# 各個体だけエクスポート
save(mem_df_list, file=args_output_df)