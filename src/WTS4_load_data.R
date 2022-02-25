source("src/functions_WTS4_load_data.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
args_sample <- args[1]
args_AnimalName <- args[2]
args_output <- args[3]

#### test args####
# args_sample <- c("3")
# args_AnimalName <- c("data/normalize_1/AnimalName.csv")
# args_output <- c("data/n1_28sample/ReadData_3.RData")
  
#### load AnimalName.csv####
read_csv("data/normalize_1/AnimalName.csv") %>% 
    as.data.frame() -> AnimalName

#### prepare uniqueNames.csv####
dir_csv <-c("data/normalize1_withoutremove/")
name_ID <- AnimalName[args_sample,2]
eval(parse(text=paste0("args_name <- c('",dir_csv,name_ID,"_uniqueNames.csv')")))
#### load uniqueNames.csv####
read_csv(args_name, 
         col_names = FALSE) %>% 
    as.data.frame() %>% .[,1] -> uniqueNames

#### prepare ratio.csv####
ratio_csv <- AnimalName[args_sample,3]
eval(parse(text=paste0("args_ratio <- c('",dir_csv,ratio_csv,"')")))
#### load ratio.csv & make df####
read_csv(args_ratio, 
         col_names = uniqueNames) %>% 
    as.data.frame() -> df_ratio_csv

#### save data matrix####
ReadData <- df_ratio_csv
save(ReadData, file=args_output)