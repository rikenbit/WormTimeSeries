source("src/functions_WTS4_sample_sheet.R")

#### args setting####
# elwood上のRで１回dだけ実行
#### test args####
args_output <- c("data/n1_28sample/WTS4_sample_sheet.csv")
args_input_path <- c("data/n1_28sample")

#### inputファイル名のリスト####
input_path_list <- list.files(args_input_path, pattern="ReadData", full.names=TRUE)
#### fix sample number sort####
input_path_list %>% 
    str_remove(., args_input_path) %>% 
    str_remove(., "/ReadData_") %>% 
    str_remove(., ".RData") %>% 
    as.numeric() %>% 
    sort() -> sample_sort_num
input_path <- c()
for(i in sample_sort_num){
    eval(parse(text=paste0("path <- c('",args_input_path,"/ReadData_",i,".RData')")))
    input_path <- c(input_path, path)
}
##### 格納ベクトル####
n_sample <- c()
n_cell <- c()
celltype <- c()

#### ロードしてリスト####
for (i in 1:length(sample_sort_num)) {
  # load ReadData
    load(input_path[i])
  # create Sample.number
    rep(sample_sort_num[i], ncol(ReadData)) %>% 
        as.character() %>% 
            append(n_sample, .) -> n_sample
    # create Cell.number
    seq(1:length(ReadData)) %>% 
        as.character() %>% 
            append(n_cell, .) -> n_cell
    # create Cell.type
    ReadData %>% 
        colnames() %>% 
            append(celltype, .) -> celltype
}
#### １つのデータフレームにまとめる####
sample_sheet <- data.frame(
    SampleNumber = n_sample,
    CellNumber = n_cell,
    CellType = celltype,
    stringsAsFactors = FALSE
    )

#### write.csv#####
write.csv(sample_sheet, args_output, row.names=FALSE)
