library(usedist)
library(tidyverse)
.filter_cellnames <- function(X){
    D <- X
    D_cell <- attr(D, "Labels")
    D_cell_f <- D_cell[grep("^[0-9]", D_cell, invert=TRUE)]
    D_f <- dist_subset(D, D_cell_f)
    D_f
}

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

.mean_na <- function(x) {
    mean(x, na.rm=TRUE)
}
stat_mean <- function(arr) {
    apply(arr, c(1,2), .mean_na)
}
arr2mat_stat <- list(
    "mean" = stat_mean
)