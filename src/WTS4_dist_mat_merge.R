source("src/functions_WTS4_dist_mat_merge.R")

#### args setting####
# args <- commandArgs(trailingOnly = T)
# # input matrix
# args_input <- args[1]
# # labels
# args_output <- args[2]
#### test args####
args_input <- ("output/WTS4/normalize_1/stimAfter/SBD_abs/Distance/Ds_F.RData")
args_output <- ("output/WTS4/normalize_1/stimAfter/SBD_abs/Ds_F_mat_mean/SampleNumber_ALL.RData")
    
load(args_input)
#### 各個体のinput行列 簡易版####
for(i in 1:length(Ds_F)){
    shift_matrix <- as.matrix(Ds_F[[i]])
    ani_n <- str_remove(names(Ds_F)[i],"animal_")
    output <- paste0("output/WTS4/normalize_1/stimAfter/SBD_abs/Ds_F_mat_mean/SampleNumber_",ani_n,".RData")
    save(shift_matrix, file=output)
}

#### 全個体のinput行列####
Ds_F |>
    lapply(function(x) {
        as.matrix(x)
    }) -> Ds

#### cellname####
uni_names <- .union_cellnames(Ds)

#### arr####
I <- length(uni_names)
M <- length(Ds)
# arr template
arr <- array(NA, dim = c(I, I, M))
dimnames(arr) <- list(
    uni_names,
    uni_names,
    names(Ds)
)
# arr value
for(i in 1:M) {
    for(j in 1:I) {
        for(k in 1:I) {
            if(uni_names[j] %in%  colnames(Ds[[i]])) {
                if(uni_names[k] %in%  rownames(Ds[[i]])) {
                    arr[j,k,i] <- Ds[[i]][uni_names[j], uni_names[k]]
                }
            }
        }
    }
}

#### arr to matrix####
shift_matrix <- arr2mat_stat[["mean"]](arr)

#### save####
save(shift_matrix, file=args_output)