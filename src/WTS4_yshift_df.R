source("src/functions_WTS4_yshift_df.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
args_input_yshift <- args[1]
#### test args####
# INPUT yshift
args_input_yshift <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Shift_F/SampleNumber_1.RData")

#### load####
load(args_input_yshift)
#### ?####
indices <- t(combn(seq(nrow(shift_matrix)), 2))
purrr::map_dfr(1:nrow(indices), .yshift_dfr) -> yshift_df

#### ?####
#### ?####
#### save####
