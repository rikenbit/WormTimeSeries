source("src/functions_WTS1.R")

# # load and export neuron activity data
# ##################################################
# # number of c.elegans neuron activity file
# list.files("data/cleandata/", all.files=TRUE, recursive=TRUE) %>%
#     length() %>%
#         seq() -> celegans
# # list of load files
# loadfiles <- paste0("data/cleandata/", celegans, "_ratio.csv")
# # TimePoints 0.2sec/Point
# Time <- seq(1, 6000, 1)

# # load csv 
# # Rの仕様で数字の細胞型名は先頭に「X」がついた名前に変換される
# for(i in celegans){
# 	eval(parse(text=paste0("data_",i," <- read.csv(loadfiles[",i,"], header=TRUE)")))
# }
# # Manual Fix sample2 Cellname　FLP→FLPR
# colnames(data_2)[match("FLP",colnames(data_2))] <- c("FLPR")

# for(i in celegans){
#     eval(parse(text=paste0("data_",i," %>% as.matrix() -> mat_data_",i)))
#     eval(parse(text=paste0("Time -> rownames(mat_data_",i,")")))
#     eval(parse(text=paste0("t(mat_data_",i,") -> matrix_",i)))
#     eval(parse(text=paste0("save(matrix_",i,", file ='data/cleandata_mat/matrix_",i,".RData')")))
# }
# ##################################################

# # export sample sheet
# ##################################################
# # 行列データをリストにまとめた
# matrix_list <- list()
# for (i in 1:15) {
#     eval(parse(text=paste0("matrix_list <- c(matrix_list, list(matrix_",i,"))")))
# }
# # 個体No.，細胞No.，細胞型名のデータフレーム作成
# n_sample <- c()
# n_cell <- c()
# celltype <- c()
# for (i in 1:15) {
#     # create Sample.number
#     # rep(celegans[1],nrow(matrix_list[[1]])) %>% as.character()
#     rep(celegans[i],nrow(matrix_list[[i]])) %>%
#         as.character() %>%
#             append(n_sample, .) -> n_sample
#     # create Cell.number
#     # seq(nrow(matrix_1))
#     seq(nrow(matrix_list[[i]])) %>%
#         as.character() %>%
#             append(n_cell, .) -> n_cell
#     # create Cell.type
#     rownames(matrix_list[[i]]) %>%
#             append(celltype, .) -> celltype
# }
# # create dataframe
# data.frame(
#         SampleNumber = n_sample,
#         CellNumber = n_cell,
#         CellType = celltype,
#         stringsAsFactors = FALSE
# ) -> sample_sheet

# # export
# write.csv(sample_sheet, "data/WTS1_sample_sheet.csv", row.names=FALSE)
# ##################################################


# # load and export stimulation timing data
# ##################################################
# # load csv 
# stimdata <- read.csv("data/stimulation_timing.csv")

# # stim celegans1~15
# # １列目のタイムフレームを削除したデータフレーム作成
# stimdata %>% select(., -1) -> stim15

# # timeframe
# timeframe <- stimdata[,1]

# # save each stimulation timing
# for(i in celegans){
#     # stimtiming <- stim15[,1]
#     eval(parse(text=paste0("stimtiming <- stim15[,",i,"]")))
#     # dataframe
#     data.frame(
#         TimeFrame = timeframe,
#         StimTiming = stimtiming,
#         stringsAsFactors = FALSE
#     ) -> stimdf
#     # stimdf -> stim_1
#     eval(parse(text=paste0("stimdf -> stim_",i)))
#     # save(stim_1, file ='data/stimulation/stim_1.RData')
#     eval(parse(text=paste0("save(stim_",i,", file ='data/stimulation/stim_",i,".RData')")))
# }
# ##################################################

# raw CFP
####################################################################################
# Neuron Activity Data
##################################################
# config
path <- "data"
excelsheet <- "pi_k_Ch2"
outputdir <- "raw_CFP"

# input file list
n_filename <- list.files(path, pattern=".xlsx") 
seq(length(n_filename)) %>%
    map(., xlsx_to_ReadData)  -> ReadDataList

# export
for(i in 1:length(n_filename)){
    #outputpath <- paste(path, outputdir, 'ReadData_1.RData', sep = '/')
    eval(parse(text=paste0("outputpath <- paste(path, outputdir, 'ReadData_",i,".RData', sep = '/')")))
    #ReadData_1 <- ReadDataList[[1]]
    eval(parse(text=paste0("ReadData_",i," <- ReadDataList[[",i,"]]")))
    #save(ReadData_1, file = outputpath)
    eval(parse(text=paste0("save(ReadData_",i,", file = outputpath)")))
}
##################################################

# Sample Sheet
##################################################
# config
path <- "data"
outputdir <- "raw_CFP"

n_sample <- c()
n_cell <- c()
celltype <- c()

for (i in 1:15) {
  # create Sample.number
  rep(seq(length(n_filename))[i], ncol(ReadDataList[[i]])) %>%
    as.character() %>%
    append(n_sample, .) -> n_sample
  # create Cell.number
  seq(ncol(ReadDataList[[i]])) %>%
    as.character() %>%
    append(n_cell, .) -> n_cell
  # create Cell.type
  colnames(ReadDataList[[i]]) %>%
    append(celltype, .) -> celltype
}
data.frame(
  SampleNumber = n_sample,
  CellNumber = n_cell,
  CellType = celltype,
  stringsAsFactors = FALSE
) -> sample_sheet

# export
# write.csv(sample_sheet, "data/WTS1_sample_sheet.csv", row.names=FALSE)
sample_sheet_path <- paste(path, outputdir, 'WTS1_sample_sheet.csv', sep = '/')
write.csv(sample_sheet, sample_sheet_path, row.names=FALSE)
##################################################

# Number_to_Animalname
##################################################
# config
path <- "data"
outputdir <- "raw_CFP"

samplenumber <- seq(length(n_filename))
animalname <- c()
for (i in 1:length(n_filename)) {
  animalname[i] <- str_sub(n_filename[i], start = 1, end = 7)
}
# Number_to_Animalname df
data.frame(
  SampleNumber = samplenumber,
  AnimaleName = animalname,
  N_FileName = n_filename,
  stringsAsFactors = FALSE
) -> AnimalName

# export
AnimalName_path <- paste(path, outputdir, 'AnimalName.csv', sep = '/')
write.csv(AnimalName, AnimalName_path, row.names=FALSE)
##################################################

# mCherry
##################################################
# config
path <- "data"
excelsheet <- "pi_k_Ch1"
outputdir <- "mCherry"

# input file list
n_filename <- list.files(path, pattern=".xlsx") 

seq(length(n_filename)) %>%
    map(., xlsx_to_ReadData)  -> mCherryList

# export
for(i in 1:length(n_filename)){
    #outputpath <- paste(path, outputdir, 'mCherry_1.RData', sep = '/')
    eval(parse(text=paste0("outputpath <- paste(path, outputdir, 'mCherry_",i,".RData', sep = '/')")))
    #mCherry_1 <- mCherryList[[1]]
    eval(parse(text=paste0("mCherry_",i," <- mCherryList[[",i,"]]")))
    #save(mCherry_1, file = outputpath)
    eval(parse(text=paste0("save(mCherry_",i,", file = outputpath)")))
}
##################################################
# Position
##################################################
# config
path <- "data"
excelsheet <- "Sheet1"
outputdir <- "Position"
# outputpath <- paste(path, outputdir, sep = '/')

# input file list
n_filename <- list.files(path, pattern=".xlsx") 

seq(length(n_filename)) %>%
    map(., xlsx_to_ReadData)  -> PositionList

# export
for(i in 1:length(n_filename)){
    #outputpath <- paste(path, outputdir, 'Position_1.RData', sep = '/')
    eval(parse(text=paste0("outputpath <- paste(path, outputdir, 'Position_",i,".RData', sep = '/')")))
    #Position_1 <- PositionList[[1]]
    eval(parse(text=paste0("Position_",i," <- PositionList[[",i,"]]")))
    #save(Position_1, file = outputpath)
    eval(parse(text=paste0("save(Position_",i,", file = outputpath)")))
}
##################################################
# Stimulation Data
##################################################
# config
path <- "data/stimulation"

StimData <- read.xlsx("data/stimulation/stimulation_timing.xlsx",
                      sheet = "Sheet2",
                      rowNames = TRUE,
                      colNames =TRUE)
StimData_Edit <- StimData[7:6006,]

# export
for(i in 1:ncol(StimData_Edit)){
  # stim_1 <- StimData_Edit[,1]
  eval(parse(text=paste0("Stim_",i," <- StimData_Edit[,",i,"]")))
  eval(parse(text=paste0("save(Stim_",i,", file = 'data/stimulation/Stim_",i,".RData')")))
}

# Number_to_Animalname
samplenumber_stim <- colnames(StimData_Edit)
animalname_stim <- unname(unlist(StimData[2,]))
data.frame(
  SampleNumber = samplenumber_stim,
  AnimaleName = animalname_stim,
  stringsAsFactors = FALSE
) -> AnimalName_Stim

# export
AnimalName_path <- paste(path, 'AnimalName_Stim.csv', sep = '/')
write.csv(AnimalName_Stim, AnimalName_path, row.names=FALSE)
##################################################
####################################################################################