source("src/functions_WTS4_Membership_F.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
# No. of Clusters 
args_k <- args[1]
# input
args_input <- args[2]
# output
args_output <- args[3]

#### test args####
# # No. of Clusters 
# args_k <- c("3")
# args_input <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Distance/Ds_F.RData")
# args_output <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Membership_F/k_Number_3.RData")

#### No. of Clusters####
k <- as.numeric(args_k)
########

# 距離行列のリスト,数字の細胞フィルター済み
load(args_input)

# Clustering against each distance matrix
C <- lapply(Ds_F, function(d, k){
    cutree(hclust(d, method="ward.D2"), k)
}, k=k)

# Cluster Labels → Indicator Matrices
Hs <- lapply(C, function(x){
    out <- matrix(0, nrow=length(x), ncol=length(unique(x)))
  	for(i in seq_along(x)){
  	    out[i,x[i]] <- 1
	}
	rownames(out) <- names(x)
    out
})

# 全個体で取りうる細胞名を取得。
cellnames <- unique(unlist(lapply(Hs, rownames)))

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
save(newHs, file=args_output)