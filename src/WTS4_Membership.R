source("src/functions_WTS4_Membership.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
# No. of Clusters 
args_k <- args[1]
# input distance matrix
args_input_path <- args[2]
# output Membership matrix
args_output_membership <- args[3]
#### test args####
# # No. of Clusters 
# args_k <- c("3")
# args_input_path <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Distance")
# args_output_membership  <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Membership/k_Number_3.RData")

#### No. of Clusters####
k <- as.numeric(args_k)
########

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

# 空の行列を格納するファイルを作成
D <- list()
# ファイルを読み込んで，リストに加える．各リストのattr(*, "Labels")に細胞型名が残っている
for(i in 1:length(input_path_list)){
    load(input_path_list[i])
    D <- c(D, list(d))
}

# Clustering against each distance matrix
C <- lapply(D, function(d, k) {
    cutree(hclust(d, method="ward.D2"), k)
    }, k=k)

# Cluster Labels → Indicator Matrices
Hs <- lapply(C, function(x) {
    out <- matrix(0, nrow=length(x), ncol=length(unique(x)))
  	for(i in seq_along(x)) {
  		  out[i,x[i]] <- 1
	  }
	  rownames(out) <- names(x)
      out
    })

# fix Membership
# 全個体で取りうる細胞名を取得。数字名も細胞名も両方
cellnames <- unique(unlist(lapply(Hs, rownames)))
# 上記細胞名のうち数字でない細胞名のみを残す
cellnames <- cellnames[grep("^[0-9]", cellnames, invert=TRUE)]

newHs <- list()
for(i in seq_along(Hs)){
    H <- Hs[[i]]
    out <- matrix(0, nrow=length(cellnames), ncol=ncol(H))
    rownames(out) <- cellnames
    for(j in seq_along(cellnames)){
        target <- which(cellnames[j] == rownames(H))
        if(length(target) != 0){
            out[j, ] <- H[target, ]
        }
    }
    newHs[[i]] <- out
}

#### ggsave####
save(newHs, file=args_output_membership)