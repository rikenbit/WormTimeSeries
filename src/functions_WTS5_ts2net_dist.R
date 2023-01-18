#library
##################################################
library(tidyverse)
library(openxlsx)
# install.packages("ts2net")
remotes::install_github("lnferreira/ts2net")
# 1でALLを選択
library(ts2net)
##################################################
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


