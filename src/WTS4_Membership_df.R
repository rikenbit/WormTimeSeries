source("src/functions_WTS4_Membership_df.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
args_input <- args[1]
args_output <- args[2]
args_animal <- args[3]
args_sample_path <- args[4]
#### test args####
# args_input <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Membership_F/k_Number_6.RData")
# args_output <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Membership_df/k_Number_6/SampleNumber_28.RData")
# args_animal <-c("28")
# args_sample_path <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Distance")

load(args_input)

#### リストにつける名前####
sample_path_list <- list.files(args_sample_path, pattern="SampleNumber_", full.names=TRUE)
sample_path_list %>% 
    str_remove(., args_sample_path) %>% 
    str_remove(., "/SampleNumber_") %>% 
    str_remove(., ".RData") %>% 
    as.numeric() %>% 
    sort() -> sample_sort_num

#### list####
mem_list <- newHs
names(mem_list) <- as.character(sample_sort_num)
#### 各個体####
mem <- mem_list[[args_animal]]
mem_mat <- mem %*% t(mem)

#### matrix to df####
indices <- t(combn(seq(nrow(mem_mat)), 2))
purrr::map_dfr(1:nrow(indices), .mem_df) -> mem_df

#### merge yshift####
args_yshift <- paste0("output/WTS4/normalize_1/stimAfter/SBD_abs/yShift_df/SampleNumber_",args_animal,".RData")
load(args_yshift)
yshift_mem_df <- merge(merge_df, 
                       mem_df,
                       by.x = "cell_cell", 
                       by.y = "cell_cell",
                       all.x = TRUE)
########
save(yshift_mem_df, file=args_output)