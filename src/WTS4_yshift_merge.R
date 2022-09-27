source("src/functions_WTS4_yshift_merge.R")

#### args setting####
#### test args####
args_input <- ("output/WTS4/normalize_1/stimAfter/SBD_abs/Shift_F")
args_output <- ("output/WTS4/normalize_1/stimAfter/SBD_abs/Shift_FM_ave/SampleNumber_ALL.RData")

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

#### template matrix####
# # 細胞名の和集合の行列、値はNA
# shift_matrix<- matrix(NA,
#                       nrow = length(uni_names),
#                       ncol = length(uni_names)
#                       )
# colnames(shift_matrix) <- uni_names
# rownames(shift_matrix) <- uni_names

#### matrix purrr###
lapply(Ds, function(z) {
    sapply(1:length(uni_names), function(x) {
        sapply(1:length(uni_names), function(y) {
            if(TRUE){
                return_object <- z[uni_names[x],uni_names[y]]
            }else{
                # return_object <- NA
                return_object <- 0
            }
            return(return_object)
        })
    })
}) -> Ds_list

#### purrr####
# colnames(shift_matrix)[1]
# rownames(shift_matrix)[1]
# 
# sapply(uni_names, function(x){
#     sapply(uni_names, function(y){
#         .value_list = function(x) {
#         	return(return_object)
#         }
#         purrr::map_dbl(Ds,
#                        .[x,y])
#     })
# })

#### arr####
I <- length(uni_names)
M <- length(Ds)
arr <- array(NA, dim = c(I, I, M))
dimnames(arr) <- list(
    uni_names,
    uni_names,
    names(Ds)
)
# if(uni_names[24] %in%  colnames(Ds[[1]])) {
#     if(uni_names[24] %in%  rownames(Ds[[1]])) {
#         #tensor に値を代入
#         arr[24,24,1] <- Ds[[1]][uni_names[24],uni_names[24]]
#     }
# }

# for(i in 1:M) {
#     D <- Ds[[i]]
#     for(j in 1:I) {
#         for(k in 1:I) {
#             if(uni_names[j] %in%  colnames(Ds[[i]])) {
#                 if(uni_names[k] %in%  rownames(Ds[[i]])) {
#                     arr[j,k,i] <- D[uni_names[j], uni_names[k]]
#                 }
#             }
#         }
#     }
# }
# arr_Ds <- arr





# .for_jk <- function(x){
#     i <- x
#     D <- Ds[[i]]
#     for(j in 1:I) {
#         for(k in 1:I) {
#             if(uni_names[j] %in%  colnames(Ds[[i]])) {
#                 if(uni_names[k] %in%  rownames(Ds[[i]])) {
#                     arr[j,k,i] <- D[uni_names[j], uni_names[k]]
#                 }
#             }
#         }
#     }
# }
# purrr::map(1:M,.for_jk) -> res_for_jk

#### purrr to matrix###
#### t(matrix)###
#### save####