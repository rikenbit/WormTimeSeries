source("src/functions_WTS4_Cluster_sample_F.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
args_input_path <- args[1]
args_output <- args[2]
args_k <- args[3]

#### No. of Clusters####
k <- as.numeric(args_k)

#### load dist data####
# 空の行列を格納するファイルを作成
D <- list()
# inputファイル名のリスト
input_path_list <- list.files(args_input_path, pattern="SampleNumber_", full.names=TRUE)

#### fix sample number sort####
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
####

for(i in 1:length(input_path_list)){
    load(input_path_list[i])
    D <- c(D, list(d))
}

#### Clustering against each distance matrix####
C <- lapply(D, function(d, k) {
    cutree(hclust(d, method="ward.D2"), k)
    }, k=k)

save(C, file=args_output)