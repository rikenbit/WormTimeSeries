#library
##################################################
library(tidyverse)
##################################################
.mem_df = function(x) {
    cell_2 <-colnames(mem_mat)[indices[x,][1]]
    cell_1 <- rownames(mem_mat)[indices[x,][2]]
    mem_mat[indices[x,][2],indices[x,][1]]
    
    data.frame(cell_cell = paste0(cell_1, "_", cell_2),
               member = mem_mat[indices[x,][2],indices[x,][1]],
               stringsAsFactors = FALSE
    ) -> return_object
    return(return_object)
}
.mem_vec = function(x) {
    mem_mat[merge_df_sep[x,1],merge_df_sep[x,2]]
}