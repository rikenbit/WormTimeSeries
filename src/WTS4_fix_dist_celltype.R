source("src/functions_WTS4_fix_dist_celltype.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
args_dist_sample <- args[1]

#### test args####
# input output SBD_abs距離行列
args_dist_sample <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Distance/SampleNumber_14.RData")

# d 距離行列データ
load(args_dist_sample)

# if文で判定 FLP
d_celltype <- attr(d, "Labels")
if_FLP <- length(which(!is.na(match(d_celltype,"FLP"))))
if(if_FLP==1){
    attr(d, "Labels")[which(!is.na(match(d_celltype,"FLP")))] <- c("FLPR")
}
# if文で判定 HYPL9VR4
d_celltype <- attr(d, "Labels")
if_FLP <- length(which(!is.na(match(d_celltype,"HYPL9VR4"))))
if(if_FLP==1){
    attr(d, "Labels")[which(!is.na(match(d_celltype,"HYPL9VR4")))] <- c("HYPL9VR")
}

#### save####
save(d, file=args_dist_sample)
