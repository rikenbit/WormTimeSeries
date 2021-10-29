source("src/functions_WTS4_ReClustering.R")

### args setting####
args <- commandArgs(trailingOnly = T)
# No. of Clusters 
args_k <- args[1]
# Method of ReClustering
args_method <- args[2]
# input
args_input_membership <- args[3]
# output merged_data
args_output_data <- args[4]
# output merged_distance
args_output_distance <- args[5]
# output merged_cls
args_output_cls <- args[6]

# #### test args####
# # No. of Clusters 
# args_k <- c("3")
# # Method of ReClustering
# args_method <- c("MCMIHOOI")
# # input
# args_input_membership <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/3_Clusters/Membership.RData")
# # output merged_data
# args_output_data <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/3_Clusters/MCMIHOOI/merged_data.RData")
# # output merged_distance
# args_output_distance <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/3_Clusters/MCMIHOOI/merged_distance.RData")
# # output merged_cls
# args_output_cls <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/3_Clusters/MCMIHOOI/merged_cls.RData")


#### No. of Clusters####
k <- as.numeric(args_k)
########

# newHs
load(args_input_membership)

# S
S <- lapply(newHs, function(h){h %*% t(h)})

# Perform MC-MI-HOOI
A <- array(0, dim=c(nrow(S[[1]]), ncol(S[[1]]), length(S)))
for(i in seq_len(dim(A)[3])){
    A[,,i] <- S[[i]]
}

#MC-MI-HOOI
if(args_method == "MCMIHOOI") {
    # options(repos="https://cran.ism.ac.jp/")
    # if (!requireNamespace("BiocManager", quietly = TRUE))
    #  install.packages("BiocManager",update = TRUE)
    # BiocManager::install("rTensor",update = FALSE)
    # BiocManager::install("einsum",update = FALSE)
    library("rTensor")
    install.packages("data/einsum_0.1.0.tar", repos = NULL, type = "source")
    library("einsum")
}


# merged_data 
merged_data <- switch(args_method,
                    # Perform Consensus Clustering
                    "CSPA" = CSPA(newHs),
                    "OINDSCAL" = OINDSCAL(S, k),
                    "MCMIHOOI" = MCMIHOOI(A, k),
                    stop("Only can use all, CSPA, OINDSCAL, MCMIHOOI")
                    )

# merged_distance
merged_distance <- switch(args_method,
                        # for t-SNE/UMAP/Clustering
                        "CSPA" = as.dist(1 - merged_data),
                        "OINDSCAL" = dist(merged_data$X),
                        "MCMIHOOI" = dist(merged_data$U),
                        stop("Only can use all, CSPA, OINDSCAL, MCMIHOOI")
                        )

# add CellType
attr(merged_distance, "Labels") <- rownames(newHs[[1]])

# merged_cls
merged_cls <- switch(args_method,
                    # Perform Consensus Clustering
                    "CSPA" = cutree(hclust(merged_distance, method="ward.D2"), k),
                    "OINDSCAL" = cutree(hclust(merged_distance, method="ward.D2"), k),
                    "MCMIHOOI" = cutree(hclust(merged_distance, method="ward.D2"), k),
                    stop("Only can use all, CSPA, OINDSCAL, MCMIHOOI")
                    )

#### ggsave####
save(merged_data, file=args_output_data)
save(merged_distance, file=args_output_distance)
save(merged_cls, file=args_output_cls)