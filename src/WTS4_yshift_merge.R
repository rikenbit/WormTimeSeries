source("src/functions_WTS4_yshift_merge.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
args_input <- args[1]
args_output <- args[2]
args_stat <- args[3]

# #### test args####
# args_input <- ("output/WTS4/normalize_1/stimAfter/SBD_abs/Shift_F")
# args_output <- ("output/WTS4/normalize_1/stimAfter/SBD_abs/Shift_FM_mean/SampleNumber_ALL.RData")
# args_stat <- ("mean")

#### load&list####
sample_path_list <- list.files(args_input, pattern="SampleNumber_", full.names=TRUE)
sample_path_list |> 
    str_remove(args_input) |> 
    str_remove("/SampleNumber_") |> 
    str_remove(".RData") |> 
    as.numeric() |> 
    sort() -> sample_sort_num
input_path_list <- c()
for(i in sample_sort_num){
    eval(parse(text=paste0("path <- c('",args_input,"/SampleNumber_",i,".RData')")))
    input_path_list <- c(input_path_list, path)
}

Ds <- list()
for(i in seq(length(sample_sort_num))){
    load(input_path_list[i])
    x <- sample_sort_num[i]
    eval(parse(text=paste0("Ds <- c(Ds, animal_",x,"=list(shift_matrix))")))
}

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
shift_matrix <- arr2mat_stat[[args_stat]](arr)

#### save####
save(shift_matrix, file=args_output)