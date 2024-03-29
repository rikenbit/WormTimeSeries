source("src/functions_WTS4_yshift_df.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
args_input_yshift <- args[1]
args_input_mSBD <- args[2]
args_output <- args[3]

#### test args####
# # INPUT yshift
# args_input_yshift <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Shift_F/SampleNumber_1.RData")
# # INPUT mSBD
# args_input_mSBD <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Distance/SampleNumber_1.RData")
# # OUTPUT merge dataframe
# args_output <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/yShift_df/SampleNumber_1.RData")

#### load####
load(args_input_yshift)

#### df yshift####
indices <- t(combn(seq(nrow(shift_matrix)), 2))
purrr::map_dfr(1:nrow(indices), .yshift_dfr) -> yshift_df
yshift_df[order(yshift_df$cell_cell),] -> yshift_df_sort

#### df mSBD####
load(args_input_mSBD)
d_f <- .filter_cellnames(d)
mat_d_f <- as.matrix(d_f)
purrr::map_dfr(1:nrow(indices), .mSBD_dfr) -> mSBD_df
mSBD_df[order(mSBD_df$cell_cell),] -> mSBD_df_sort

#### merge df####
merge_df <- merge(yshift_df_sort, 
                  mSBD_df_sort,
                  by.x = "cell_cell", 
                  by.y = "cell_cell"
                  )

#### save####
save(merge_df, file=args_output)