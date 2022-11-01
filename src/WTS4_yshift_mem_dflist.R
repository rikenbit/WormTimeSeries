source("src/functions_WTS4_yshift_mem_dflist.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
args_sample_path <- args[1]
args_output <- args[2]
args_k <- args[3]
#### test args####
# args_sample_path <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Distance")
# args_output <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Membership_df/k_Number_6/DFs.RData")
# args_k <- c("6)

#### サンプル####
sample_path_list <- list.files(args_sample_path, pattern="SampleNumber_", full.names=TRUE)
sample_path_list %>% 
    str_remove(., args_sample_path) %>% 
    str_remove(., "/SampleNumber_") %>% 
    str_remove(., ".RData") %>% 
    as.numeric() %>% 
    sort() -> sample_sort_num
#### yshift####
input_path_list <- c()
for(i in sample_sort_num){
    eval(parse(text=paste0("path <- c('output/WTS4/normalize_1/stimAfter/SBD_abs/Membership_df/k_Number_",args_k,"/SampleNumber_",i,".RData')")))
    input_path_list <- c(input_path_list, path)
}
DFs <- list()
for(i in seq(length(sample_sort_num))){
    load(input_path_list[i])
    x <- sample_sort_num[i]
    eval(parse(text=paste0("DFs <- c(DFs, animal_",x,"=list(yshift_mem_df))")))
}

#### save####
save(DFs, file=args_output)