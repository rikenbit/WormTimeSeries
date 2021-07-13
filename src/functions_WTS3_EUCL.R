#library
##################################################
options(rgl.useNULL=TRUE)
library(TSclust)
library(openxlsx)
library(tidyverse)
##################################################
set.seed(1234)

.eucl_cal_stimAfter = function(x,y) {
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
  	tsReadData <- x[stimtimng:nrow(x),]
  	return_object <- diss(tsReadData, "EUCL")
  	return(return_object)
}
.eucl_cal_all = function(x) {
  	tsReadData <- x
	return_object <- diss(tsReadData, "EUCL")
	return(return_object)
}