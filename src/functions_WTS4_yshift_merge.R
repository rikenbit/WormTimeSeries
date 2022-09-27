#library
##################################################
library(tidyverse)
##################################################
# .func = function(x) {
# 	return(return_object)
# }
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