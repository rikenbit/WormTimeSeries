# library
##################################################
library(openxlsx)
library(tidyverse)
library(RcppRoll)
##################################################

# func
##################################################
# WTS1 load rawData
xlsx_to_ReadData <- function(x, FileList = n_filename, Excelsheet = excelsheet) {
  fullpath <- paste(path, FileList[x], sep = '/')
  ReadData <- read.xlsx(fullpath, sheet = Excelsheet, rowNames = TRUE, colNames =TRUE)
  return(ReadData)
}
##################################################