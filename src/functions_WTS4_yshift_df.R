#library
##################################################
library(tidyverse)
library(usedist)
##################################################
.func = function(x) {
	return(return_object)
}

.yshift_dfr = function(x) {
    cell_2 <-colnames(shift_matrix)[indices[x,][1]]
    cell_1 <- rownames(shift_matrix)[indices[x,][2]]
    shift_matrix[indices[x,][2],indices[x,][1]]
    
    data.frame(cell_cell = paste0(cell_1, "_", cell_2),
               yshift = shift_matrix[indices[x,][2],indices[x,][1]],
               stringsAsFactors = FALSE
    ) -> return_object
    return(return_object)
}
.mSBD_dfr = function(x) {
    cell_2 <-colnames(mat_d_f)[indices[x,][1]]
    cell_1 <- rownames(mat_d_f)[indices[x,][2]]
    mat_d_f[indices[x,][2],indices[x,][1]]
    
    data.frame(cell_cell = paste0(cell_1, "_", cell_2),
               mSBD = mat_d_f[indices[x,][2],indices[x,][1]],
               stringsAsFactors = FALSE
    ) -> return_object
    return(return_object)
}

.filter_cellnames <- function(X){
    D <- X
    D_cell <- attr(D, "Labels")
    D_cell_f <- D_cell[grep("^[0-9]", D_cell, invert=TRUE)]
    D_f <- dist_subset(D, D_cell_f)
    D_f
}
