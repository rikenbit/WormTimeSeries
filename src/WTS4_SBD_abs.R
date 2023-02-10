source("src/functions_WTS4_SBD_abs.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
# sample number
args_sample <- args[1]
# path NeuronActivity Data
args_neuron <- args[2]
# args_time <- c("all")
args_time <- args[3]
# stimtiming
args_stim_xlsx <- args[4]
# output SBD_abs距離行列
args_SBD <- args[5]

#### load NeuronActivity####
load(args_neuron)
# 元データがディレクトリごとではなく，ファイル名で各サンプルがわかれている対応
# eval(parse(text=paste0("ReadData <- ReadData_",args_sample)))

#### switch time range & trim ReadData####
ReadData <- switch(args_time,
                   "all" = .ReadData_all(ReadData),
                   "stimAfter" = .ReadData_stimAfter(ReadData, args_stim_xlsx),
                   stop("Only can use all, stimAfter ")
                   )

#### SBD_abs####
# create sbd matrix
SBD_zero_mat <- sapply(1:ncol(ReadData), function(x) {
    sapply(1:ncol(ReadData), function(z){
        if(x!=z){
            shift_1 <- ReadData[,x]
            shift_2 <- ReadData[,z]
            return_object <- dtwclust::SBD(shift_1,
                         shift_2,
                         znorm = FALSE,
                         error.check = TRUE,
                         return.shifted = TRUE)
            return_object$dist
        } else{
            return_object <- 0
            return_object
        }
    })
})
# add col/row names
colnames(SBD_zero_mat) <- colnames(ReadData)
rownames(SBD_zero_mat) <- colnames(ReadData)

# convert matrix(symmetrix ) to dist
d <- stats::as.dist(SBD_zero_mat)

# save SBD_abs dist
save(d, file=args_SBD)