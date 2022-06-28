source("src/functions_WTS4_CellCount_F.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
args_input <- args[1]
args_output <- args[2]

#### test args####
# args_input <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Distance/Ds_F.RData")
# args_output <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Distance/CellCount.RData")

# 距離行列のリスト,数字の細胞フィルター済み
load(args_input)

#### filter annotated####
#文字の細胞名
cn_list <- list()
for(i in 1:length(Ds_F)){
    Ds_F[[i]] %>% 
    	attr(., "Labels") -> cn_list[[i]]
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