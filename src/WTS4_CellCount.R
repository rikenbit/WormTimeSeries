source("src/functions_WTS4_CellCount.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
args_input_path <- args[1]
args_output <- args[2]

#### test args####
# args_input_path <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Distance")
# args_output <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Distance/CellCount.RData")

#### fix sample number sort####
# inputファイル名のリスト
input_path_list <- list.files(args_input_path, pattern="SampleNumber_", full.names=TRUE)
input_path_list %>% 
    str_remove(., args_input_path) %>% 
    str_remove(., "/SampleNumber_") %>% 
    str_remove(., ".RData") %>% 
    as.numeric() %>% 
    sort() -> sample_sort_num
input_path <- c()
for(i in sample_sort_num){
    eval(parse(text=paste0("path <- c('",args_input_path,"/SampleNumber_",i,".RData')")))
    input_path <- c(input_path, path)
}

input_path_list <- input_path 

#### load dist data####
# 空の行列を格納するファイルを作成
D <- list()
# ファイルを読み込んで，リストに加える．各リストのattr(*, "Labels")に細胞型名が残っている
for(i in 1:length(input_path_list)){
    load(input_path_list[i])
    D <- c(D, list(d))
}

#### filter annotated####
#文字の細胞名
cn_list <- list()
for(i in 1:length(D)){
    D[[i]] %>% attr(., "Labels") -> tmp
    cn_list[[i]] <- tmp[grep("^[0-9]", tmp, invert=TRUE)]
}

#### count cell####
# https://www.fixes.pub/program/167651.html
cn_list %>% 
    unlist() %>% 
    as.data.frame() %>% 
    dplyr::select(., CellType = 1, dplyr::everything()) %>% 
    dplyr::count(., CellType) %>% 
    dplyr::select(., CellType = 1, CellCount =  2, dplyr::everything()) -> df_cell_count

save(df_cell_count, file=args_output)