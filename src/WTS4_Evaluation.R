source("src/functions_WTS4_Evaluation.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
# input merged_data
args_input_data <- args[1]
# input merged_cls
args_input_cls <- args[2]
# output merged_data
args_output_data <- args[3]

# ReClustering Method
args_method <- args[4]
# Evaluation Method
args_eval_method <- args[5]


# #### test args####
# # input merged_data
# args_input_data <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/3_Clusters/MCMIHOOI/merged_data.RData")
# # input merged_cls
# args_input_cls <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/3_Clusters/MCMIHOOI/merged_cls.RData")
# # output merged_data
# args_output_data <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/3_Clusters/MCMIHOOI/PseudoF/eval_result.RData")
# # ReClustering Method
# args_method <- c("MCMIHOOI")
# # Evaluation Method
# args_eval_method <- c("PseudoF")

# load
load(args_input_data)
load(args_input_cls)

# data
eval_data <- switch(args_method,
                        # for t-SNE/UMAP/Clustering
                        "CSPA" = 1 - merged_data,
                        "OINDSCAL" = merged_data$X,
                        "MCMIHOOI" = merged_data$U,
                        stop("Only can use all, CSPA, OINDSCAL, MCMIHOOI")
                        )
# evaluation
eval_result <- switch(args_eval_method,
                        # for t-SNE/UMAP/Clustering
                        "PseudoF" = .PseudoF(eval_data, merged_cls),
                        "Connectivity" = .Connectivity(eval_data, merged_cls),
                        stop("Only can use all, PseudoF, Connectivity")
                        )
#### ggsave####
save(eval_result, file=args_output_data)