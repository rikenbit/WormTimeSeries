source("src/functions_WTS4_Cluster_sample.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
args_input_path <- args[1]
args_output <- args[2]
args_k <- args[3]

# # #### test args####
# args_input_path <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Distance")
# args_output <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Cluster_sample/k_Number_5/sample_cls.RData")
# args_k <- c("5")

#### No. of Clusters####
k <- as.numeric(args_k)

#### load dist data####
# 空の行列を格納するファイルを作成
D <- list()
# inputファイル名のリスト
input_path_list <- list.files(args_input_path, pattern="SampleNumber_", full.names=TRUE)
# ファイルを読み込んで，リストに加える．各リストのattr(*, "Labels")に細胞型名が残っている
for(i in 1:length(input_path_list)){
    load(input_path_list[i])
    D <- c(D, list(d))
}

#### Clustering against each distance matrix####
C <- lapply(D, function(d, k) {
    cutree(hclust(d, method="ward.D2"), k)
    }, k=k)

save(C, file=args_output)