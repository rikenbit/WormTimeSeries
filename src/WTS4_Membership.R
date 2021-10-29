# libraryインストールがないので
# source("src/functions_WTS4_Membership.R")

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
# args_input_path <- c("output/WTS4/normalize_1/stimAfter/SBD_abs")
# args_output_membership  <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/3_Clusters/Membership.RData")

#### No. of Clusters####
k <- as.numeric(args_k)
########

# 空の行列を格納するファイルを作成
D <- list()

# inputファイル名のリスト
input_path_list <- list.files(args_input_path, pattern="SampleNumber_", full.names=TRUE)
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
cellnames <- unique(unlist(lapply(Hs, rownames)))
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