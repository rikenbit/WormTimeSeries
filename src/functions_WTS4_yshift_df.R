#library
##################################################
library(tidyverse)
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
