#library
##################################################
library(tidyverse)
library(openxlsx)
library(dtwclust)
# library(Matrix)
##################################################
# setting Random number generator
 RNGkind(kind = "Mersenne-Twister")
# dtwclust_rngkind <- "Mersenne-Twister"

#
# すべてsourceで読み込み
########
# f <- function(dat){
#     str_c("R/", dat) 
# }
# list.files("R") %>% 
#     str_subset("\\.R") %>% 
#     map_chr(., f) %>% 
#     map(~source(.x))
########

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