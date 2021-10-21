source("src/functions_WTS4_Membership.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
# No. of Clusters 
args_k <- args[1]
# input distance matrix
args_input_dist <- args[2]
# output Membership matrix
args_output_membership <- args[3]

#### test args####
# No. of Clusters 
args_k <- c("3")
args_input_dist <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/SampleNumber_1/SBD_abs.RData")
args_output_membership  <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/SampleNumber_1/3_Clusters/Membership.RData")


#### No. of Clusters####
k <- args_k
########

##################################################
# Approach 1: Consensus Clustering
##################################################
# Distance matrix: Cell × Cell
# D1 <- dist(X1)
# D2 <- dist(X2)
# D3 <- dist(X3)

# load Distance matrix
load(args_input_dist)
D1 <- d

# D <- list(D1, D2, D3)
D <- list(D1)

# Clustering against each distance matrix
C <- lapply(D, function(d, k){
	cutree(hclust(d, method="ward.D2"), k)
}, k=k)

# Cluster Labels → Indicator Matrices
Hs <- lapply(C, function(x){
	out <- matrix(0, nrow=length(x), ncol=length(unique(x)))
	for(i in seq_along(x)){
		out[i,x[i]] <- 1
	}
	out
})

#### ggsave####
save(Hs, file=args_output_membership)