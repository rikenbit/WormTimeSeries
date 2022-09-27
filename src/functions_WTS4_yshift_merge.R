#library
##################################################
library(tidyverse)
##################################################
.union_cellnames <- function(Ds) {
    Ds |>
        lapply(function(x) {
            colnames(x)
        }) |>
        unlist() |>
        unique() |>
        sort() -> res
    res
}

#### mean matrix#####
.mean_na <- function(x) {
    mean(x, na.rm=TRUE)
}
stat_mean <- function(arr) {
    apply(arr, c(1,2), .mean_na)
}
#### sd matrix#####
.sd_na = function(x) {
    sd(x, na.rm=TRUE)
}
stat_sd <- function(arr) {
    apply(arr, c(1,2), .sd_na) 
}
#### No. of cell matrix#####
.length_na = function(x) {
    # NAとNaNを除去
    length(x[!is.na(x)])
}
stat_count <- function(arr) {
    apply(arr, c(1,2), .length_na)
}

arr2mat_stat <- list(
    "mean" = stat_mean,
    "sd" = stat_sd,
    "count" = stat_count
)