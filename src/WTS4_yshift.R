source("src/functions_WTS4_yshift.R")

#### args setting####
#### test args####
# sample number サンプル番号の指定
args_sample <- c("2")
# path NeuronActivity Data
args_neuron <- c("data/normalize_1/ReadData_2.RData")
# time range
args_time <- c("stimAfter")
# stimtiming
args_stim_xlsx <- c("data/stimulation/stimulation_timing.xlsx")
# output y-shiftの総当たり結果
args_output <-("output/WTS4/normalize_1/stimAfter/SBD_abs/Shift/SampleNumber_3.RData")

#### load NeuronActivity####
load(args_neuron)
eval(parse(text=paste0("ReadData <- ReadData_",args_sample)))
#### switch time range & trim ReadData####
ReadData <- switch(args_time,
                   "all" = .ReadData_all(ReadData),
                   "stimAfter" = .ReadData_stimAfter(ReadData, args_stim_xlsx),
                   stop("Only can use all, stimAfter ")
)
ReadData %>% as.matrix() %>% t() -> ReadData_t

####test###
ReadData_TEST <- ReadData_t[101:105,1:30]
indices <- t(combn(seq(nrow(ReadData_TEST)), 2))
indices <- indices[, 2:1]
# yshift <- apply(indices, 1, function(xx) {
#     .mSBD(ReadData_TEST[xx[1], ], ReadData_TEST[xx[2], ])$yshift
# })
shift_value <- apply(indices, 1, function(xx) {
    .mSBD(ReadData_TEST[xx[1], ], ReadData_TEST[xx[2], ])$shift_value
})
out <- matrix(0, nrow = nrow(ReadData_TEST), ncol = nrow(ReadData_TEST))
out[indices] <- shift_value
out_t <- -t(out)
out_all <- out + out_t
rownames(out_all) <- rownames(ReadData_TEST)
colnames(out_all) <- rownames(ReadData_TEST)
# out_allは行方向の細胞を基準（ゼロ）としたさいの列方向の細胞の平行移動
# 下記コマンドでどれだけ平行移動したか取得できる
# out_all["HYPL10VL",]やout_all[2,]やout_all[c("GLRR","HYPL10VL"),]
####test###



indices <- t(combn(seq(nrow(ReadData_t)), 2))
indices <- indices[, 2:1]
shift_value <- apply(indices, 1, function(xx) {
    .mSBD(ReadData_t[xx[1], ], ReadData_t[xx[2], ])$shift_value
})
out <- matrix(0, nrow = nrow(ReadData_t), ncol = nrow(ReadData_t))
out[indices] <- shift_value
out_t <- -t(out)
out_all <- out + out_t
rownames(out_all) <- rownames(ReadData_t)
colnames(out_all) <- rownames(ReadData_t)
