# library
##################################################
library(openxlsx)
library(tidyverse)
library(RcppRoll)
library(patchwork)
##################################################

# func
##################################################
# WTS1 load rawData
xlsx_to_ReadData <- function(x, Path = inputdir, FileList = n_filename, Excelsheet = excelsheet) {
  fullpath <- paste(Path, FileList[x], sep = '/')
  ReadData <- read.xlsx(fullpath, sheet = Excelsheet, rowNames = TRUE, colNames =TRUE)
  return(ReadData)
}
##################################################