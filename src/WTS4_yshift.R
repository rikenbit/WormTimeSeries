source("src/functions_WTS4_yshift.R")

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
# output y-shiftの総当たり結果の行列
args_output<- args[5]
# 値の形式　数字の細胞除去の有無
args_in_mat <- args[6]

#### load NeuronActivity####
load(args_neuron)
eval(parse(text=paste0("ReadData <- ReadData_",args_sample)))
#### switch time range & trim ReadData####
ReadData <- switch(args_time,
                   "all" = .ReadData_all(ReadData),
                   "stimAfter" = .ReadData_stimAfter(ReadData, args_stim_xlsx),
                   stop("Only can use all, stimAfter ")
)
# row is cell, column is time
t(as.matrix(ReadData)) -> ReadData_t

#### yshift value matrix####
# sheet of combination
indices <- t(combn(seq(nrow(ReadData_t)), 2))
indices <- indices[, 2:1]

# get yshift value
shift_value <- apply(indices, 1, function(xx) {
    .mSBD(ReadData_t[xx[1], ], ReadData_t[xx[2], ])$shift_value
})
out <- matrix(0, nrow = nrow(ReadData_t), ncol = nrow(ReadData_t))
out[indices] <- shift_value

# convert to symmetric matrix
out_t <- -t(out)
out_all <- out + out_t

# set cellname
rownames(out_all) <- rownames(ReadData_t)
colnames(out_all) <- rownames(ReadData_t)

#### save####
shift_matrix <- switch (args_in_mat,
                        "Shift" = out_all,
                        "Shift_F" = .filter_cellnames(out_all)
                        )
save(shift_matrix, file=args_output)