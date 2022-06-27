source("src/functions_WTS4_Dist_Filter.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
args_input <- args[1]
args_output <- args[2]
#### test args####
# args_input <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Distance/Ds.RData")
# args_output <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Distance/Ds_F.RData")

# load Ds
load(args_input)
# filter Ds
Ds_F <- lapply(Ds, filter_cellnames)
#### save####
save(Ds_F, file=args_output)