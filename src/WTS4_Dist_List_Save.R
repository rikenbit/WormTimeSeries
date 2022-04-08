source("src/functions_WTS4_Dist_List_Save.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
args_sample_path <- args[1]
args_output <- args[2]

# #### test args####
# args_sample_path <- c("output/WTS4/n1_28sample/stimAfter/SBD_abs/Distance")
# args_output <- c("output/WTS4/n1_28sample/stimAfter/SBD_abs/Distance/Ds.RData")

# リスト化するサンプルの数字を取得
sample_path_list <- list.files(args_sample_path, pattern="SampleNumber_", full.names=TRUE)
sample_path_list %>% 
    str_remove(., args_sample_path) %>% 
    str_remove(., "/SampleNumber_") %>% 
    str_remove(., ".RData") %>% 
    as.numeric() %>% 
    sort() -> sample_sort_num

input_path_list <- c()
for(i in sample_sort_num){
    eval(parse(text=paste0("path <- c('",args_sample_path,"/SampleNumber_",i,".RData')")))
    input_path_list <- c(input_path_list, path)
}

Ds <- list()

for(i in sample_sort_num){
    load(input_path_list[i])
    # Ds <- c(Ds, animal_i=list(d))
    eval(parse(text=paste0("Ds <- c(Ds, animal_",i,"=list(d))")))
}
#### save####
save(Ds, file=args_output)