#library
##################################################
library(dtwclust)
library(tidyverse)
##################################################
sbd_y = function(x) {
    shift_2 <- ReadData.list[[x]] %>% as.numeric()
    sbd <- dtwclust::SBD(shift_1,
                         shift_2, 
                         znorm = FALSE, 
                         error.check = TRUE, 
                         return.shifted = TRUE)
    return(sbd$yshift)
}

check_args_shift = function(x) {
    x -> sample_cell_type
    sample_cell_type %>% 
        str_count(., pattern="ASER") %>%
            sum() -> check_ASER
    sample_cell_type %>% 
        str_count(., pattern="BAGR") %>%
            sum() -> check_BAGR
    sample_cell_type %>% 
        str_count(., pattern="BAGL") %>%
            sum() -> check_BAGL
    if (check_ASER >= 1) {
        args_shift <- "ASER"
    } else if (check_BAGR >= 1) {
        args_shift <- "BAGR"
    } else if (check_BAGL >= 1) {
        args_shift <- "BAGL"
    } else {
        args_shift <- sample_cell_type[1]
    }
    return(args_shift)
}