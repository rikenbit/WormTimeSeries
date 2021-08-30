#library
##################################################
library(dtwclust)
RNGkind(kind = "Mersenne-Twister")
library(tidyverse)
library(openxlsx)
##################################################
.sbd_y = function(x) {
    shift_2 <- ReadData.list[[x]] %>% as.numeric()
    sbd <- dtwclust::SBD(shift_1,
                         shift_2, 
                         znorm = FALSE, 
                         error.check = TRUE, 
                         return.shifted = TRUE)
    return(sbd$yshift)
}
.ReadData_stimAfter = function(x,y) {
    #### load stim timing extra####
    stimtimng_sheet <- read.xlsx(y,
                                 sheet = "Sheet1",
                                 rowNames = FALSE,
                                 colNames =TRUE)
    stimtimng_sheet %>% 
        dplyr::select(sample_number = 1, 
                      stim_first = 7,
                      ) %>% 
            filter(sample_number == args_sample) %>% 
                .$stim_first %>% 
                    trunc() -> stimtimng #切り捨て
                    # ceiling() -> stimtimng #切り上げ
    return_object <- x[stimtimng:nrow(x),]
    return(return_object)
}
.ReadData_all = function(x) {
    return_object <- x
    return(return_object)
}

########
SBD_abs <- function(x, y, znorm = FALSE, error.check = TRUE, return.shifted = TRUE) {
    if (is_multivariate(list(x, y))) stop("SBD does not support multivariate series.")
    if (error.check) {
        check_consistency(x, "ts")
        check_consistency(y, "ts")
        x <- as.numeric(x)
        y <- as.numeric(y)
    }
    nx <- length(x)
    ny <- length(y)

    if (nx > ny) {
        # The order in which I provide the arguments to NCCc affects 'shift'
        swap <- x
        x <- y
        y <- swap
    }
    else {
        swap <- NULL
    }

    if (znorm)
        CCseq <- NCCc(zscore(x, error.check = FALSE),
                      zscore(y, error.check = FALSE),
                      error.check = FALSE)
    else
        CCseq <- NCCc(x, y, error.check = FALSE)

    # m <- max(CCseq)
    m <- max(abs(CCseq))
    if (!return.shifted) return(1 - m) # nocov
    # shift <- which.max(CCseq) - max(nx, ny)
    shift <- which.max(abs(CCseq))- max(nx, ny)

    if (is.null(swap)) {
        if (shift < 0L)
            yshift <- y[(-shift + 1L):ny]
        else
            yshift <- c( rep(0, shift), y )
    }
    else {
        # Remember, if I swapped them, then I have to shift what is now saved in 'x'
        if (shift < 0L)
            yshift <- c( rep(0, -shift), x )
        else
            yshift <- x[(shift + 1L):ny]
    }

    nys <- length(yshift)
    if (nys < nx)
        yshift <- c( yshift, rep(0, nx-nys) )
    else
        yshift <- yshift[1L:nx]
    # return
    list(dist = 1 - m, yshift = yshift)
}

#' @rdname SBD
#' @export
#'
sbd_abs <- SBD_abs

########

# ########
env <- getNamespace("dtwclust")
assignInNamespace("SBD", SBD_abs, ns="dtwclust")
assignInNamespace("sbd", sbd_abs, ns="dtwclust")

# assignInNamespace("SBD", SBD, ns="dtwclust", envir=env)
# assignInNamespace("sbd", sbd, ns="dtwclust", envir=env)

# assignInNamespace("SBD", SBD_abs, ns="dtwclust")
# assignInNamespace("SBD", SBD_abs, ns="dtwclust", envir=env)
# assignInNamespace("sbd", sbd_abs, ns="dtwclust", envir=env)

# library(proxy)
# env <- getNamespace("proxy")
# unlockBinding("SBD", env = env)
# assignInNamespace("SBD", SBD_abs, ns="proxy", envir=env)
# assignInNamespace("sbd", sbd_abs, ns="proxy", envir=env)
# ########

source("src/functions_WTS3_SBD_abs.R")

#### test args####
# sample number サンプル番号の指定
args_sample <- c("2")
# path NeuronActivity Data
args_neuron <- c("data/normalize_1/ReadData_2.RData")

# args_time <- c("all")
args_time <- c("stimAfter")
# stimtiming
args_stim_xlsx <- c("data/stimulation/stimulation_timing.xlsx")
# y-shift計算対象の細胞
args_shift <- c("ASER")
# output SBD_abs距離行列
args_SBD <- c("output/WTS3/normalize_1/stimAfter/SBD_abs/SampleNumber_2/SBD_abs.RData")
# output SBD_abs yshift
args_yshift <- c("output/WTS3/normalize_1/stimAfter/SBD_abs/SampleNumber_2/yshift.RData")

#### load NeuronActivity####
load(args_neuron)
# 元データがディレクトリごとではなく，ファイル名で各サンプルがわかれている対応
eval(parse(text=paste0("ReadData <- ReadData_",args_sample)))

#### switch time range & trim ReadData####
ReadData <- switch(args_time,
                   "all" = .ReadData_all(ReadData),
                   "stimAfter" = .ReadData_stimAfter(ReadData, args_stim_xlsx),
                   stop("Only can use all, stimAfter ")
                   )

#### SBD####
# 行列っぽいデータを細胞ごとにlist化
ReadData.list <- asplit(ReadData,2)
# create sbd matrix
hc <- dtwclust::tsclust(ReadData.list, distance = "sbd", trace = TRUE)
# hc <- dtwclust::tsclust(ReadData.list, distance = "SBD", trace = TRUE)
# convert dist
d <- stats::as.dist(hc@distmat)

#### SBD_abs####
# create sbd matrix
SBD_zero_mat <- sapply(1:ncol(ReadData), function(x) {
    sapply(1:ncol(ReadData), function(z){
        if(x!=z){
            shift_1 <- ReadData[,x]
            shift_2 <- ReadData[,z]
            return_object <- dtwclust::SBD(shift_1,
            # return_object <- SBD(shift_1,
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

# SBD関数にabsを追加し，getNamespaceで置き換え SBD <- function()のみ変更
SBD <- function(x, y, znorm = FALSE, error.check = TRUE, return.shifted = TRUE) {
    if (is_multivariate(list(x, y))) stop("SBD does not support multivariate series.")
    if (error.check) {
        check_consistency(x, "ts")
        check_consistency(y, "ts")
        x <- as.numeric(x)
        y <- as.numeric(y)
    }
    nx <- length(x)
    ny <- length(y)

    if (nx > ny) {
        # The order in which I provide the arguments to NCCc affects 'shift'
        swap <- x
        x <- y
        y <- swap
    }
    else {
        swap <- NULL
    }

    if (znorm)
        CCseq <- NCCc(zscore(x, error.check = FALSE),
                      zscore(y, error.check = FALSE),
                      error.check = FALSE)
    else
        CCseq <- NCCc(x, y, error.check = FALSE)

    # m <- max(CCseq)
    m <- max(abs(CCseq))
    if (!return.shifted) return(1 - m) # nocov
    # shift <- which.max(CCseq) - max(nx, ny)
    shift <- which.max(abs(CCseq))- max(nx, ny)

    if (is.null(swap)) {
        if (shift < 0L)
            yshift <- y[(-shift + 1L):ny]
        else
            yshift <- c( rep(0, shift), y )
    }
    else {
        # Remember, if I swapped them, then I have to shift what is now saved in 'x'
        if (shift < 0L)
            yshift <- c( rep(0, -shift), x )
        else
            yshift <- x[(shift + 1L):ny]
    }

    nys <- length(yshift)
    if (nys < nx)
        yshift <- c( yshift, rep(0, nx-nys) )
    else
        yshift <- yshift[1L:nx]
    # return
    list(dist = 1 - m, yshift = yshift)
}
sbd <- SBD
# 入れ替え
env <- getNamespace("dtwclust")
assignInNamespace("SBD", SBD, ns = "dtwclust", envir = env)

#### SBD_abs####
# create sbd matrix
SBD_zero_mat <- sapply(1:ncol(ReadData), function(x) {
    sapply(1:ncol(ReadData), function(z){
        if(x!=z){
            shift_1 <- ReadData[,x]
            shift_2 <- ReadData[,z]
            return_object <- dtwclust::SBD(shift_1,
            # return_object <- SBD(shift_1,
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