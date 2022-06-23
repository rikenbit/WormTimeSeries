#library
##################################################
library(tidyverse)
library(usedist)
##################################################
.func = function(x) {
	return(return_object)
}

.filter_cellnames <- function(X){
    D <- X
    D_cell <- attr(D, "Labels")
    D_cell_f <- D_cell[grep("^[0-9]", D_cell, invert=TRUE)]
    D_f <- dist_subset(D, D_cell_f)
    D_f
}