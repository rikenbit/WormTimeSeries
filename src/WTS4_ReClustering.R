source("src/functions_WTS4_ReClustering.R")

#### args setting####
# args <- commandArgs(trailingOnly = T)
# # No. of Clusters 
# args_k <- args[1]
# # Method of ReClustering
# args_method <- args[2]
# # input
# args_input_membership <- args[3]
# # output merged_data
# args_output_data <- args[4]
# # output merged_distance
# args_output_distance <- args[5]
# # output merged_cls
# args_output_cls <- args[6]

#### test args####
# No. of Clusters 
args_k <- c("3")
# Method of ReClustering
args_method <- c("CSPA")
# input
args_input_membership <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/3_Clusters/Membership.RData")
# output merged_data
args_output_data <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/3_Clusters/CSPA/merged_data.RData")
# output merged_distance
args_output_distance <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/3_Clusters/CSPA/merged_distance.RData")
# output merged_cls
args_output_cls <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/3_Clusters/CSPA/merged_cls.RData")


#### No. of Clusters####
k <- args_k
########

# Hs
load(args_input_membership)

# S
S <- lapply(Hs, function(h){h %*% t(h)})

# merged_data 
merged_data <- switch(args_method,
                    # Perform Consensus Clustering
                    "CSPA" = CSPA(Hs),
                    "OINDSCAL" = OINDSCAL(input_n),
                    "MCMIHOOI" = MCMIHOOI(input_n, args_stim_xlsx),
                    stop("Only can use all, CSPA, OINDSCAL, MCMIHOOI")
                    )

# merged_distance
merged_distance <- switch(args_method,
                        # for t-SNE/UMAP/Clustering
                        "CSPA" = as.dist(1 - merged_data),
                        "OINDSCAL" = OINDSCAL(input_n),
                        "MCMIHOOI" = MCMIHOOI(input_n, args_stim_xlsx),
                        stop("Only can use all, CSPA, OINDSCAL, MCMIHOOI")
                        )
# merged_cls
merged_cls <- switch(args_method,
                    # Perform Consensus Clustering
                    "CSPA" = cutree(hclust(merged_distance, method="ward.D2"), k),
                    "OINDSCAL" = OINDSCAL(input_n),
                    "MCMIHOOI" = MCMIHOOI(input_n, args_stim_xlsx),
                    stop("Only can use all, CSPA, OINDSCAL, MCMIHOOI")
                    )

#### ggsave####
save(merged_data, file=args_output_data)
save(merged_distance, file=args_output_distance)
save(merged_cls, file=args_output_cls)