source("src/functions_WTS1.R")

# raw CFP
####################################################################################
args <- commandArgs(trailingOnly = T)
args_data <- args[1]
args_sheet <- args[2]
# Neuron Activity Data
##################################################
# config
path <- "data"
# excelsheet <- "pi_k_Ch2"
excelsheet <- args_sheet
outputdir <- args_data
inputdir <- "data/raw"

# export
n_filename <- list.files(inputdir, pattern=".xlsx") 
for(i in 1:length(n_filename)){
    fullpath <- paste(inputdir, n_filename[i], sep = '/')
    ReadData <- read.xlsx(fullpath, sheet = excelsheet, rowNames = TRUE, colNames =TRUE)
    #ReadData_1 <- ReadDataList[[1]]
    eval(parse(text=paste0("ReadData_",i," <- ReadData")))
    #outputpath <- paste(path, outputdir, 'ReadData_1.RData', sep = '/')
    eval(parse(text=paste0("outputpath <- paste(path, outputdir, 'ReadData_",i,".RData', sep = '/')")))
    #save(ReadData_1, file = outputpath)
    eval(parse(text=paste0("save(ReadData_",i,", file = outputpath)")))
}

##################################################

# Sample Sheet
##################################################
# config
path <- "data"
outputdir <- args_data

n_sample <- c()
n_cell <- c()
celltype <- c()
ReadDataList <- list()
for(i in 1:length(n_filename)) {
    eval(parse(text=paste0("ReadDataList <- c(ReadDataList,list(ReadData_",i,"))")))
}

for (i in 1:length(n_filename)) {
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
outputdir <- args_data

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
inputdir <- "data/raw"

# export
n_filename <- list.files(inputdir, pattern=".xlsx") 
for(i in 1:length(n_filename)){
    fullpath <- paste(inputdir, n_filename[i], sep = '/')
    mCherry <- read.xlsx(fullpath, sheet = excelsheet, rowNames = TRUE, colNames =TRUE)
    #mCherry_1 <- mCherryList[[1]]
    eval(parse(text=paste0("mCherry_",i," <- mCherry")))
    #outputpath <- paste(path, outputdir, 'mCherry_1.RData', sep = '/')
    eval(parse(text=paste0("outputpath <- paste(path, outputdir, 'mCherry_",i,".RData', sep = '/')")))
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
inputdir <- "data/raw"

# export
n_filename <- list.files(inputdir, pattern=".xlsx") 
for(i in 1:length(n_filename)){
    fullpath <- paste(inputdir, n_filename[i], sep = '/')
    Position <- read.xlsx(fullpath, sheet = excelsheet, rowNames = TRUE, colNames =TRUE)
    #Position_1 <- PositionList[[1]]
    eval(parse(text=paste0("Position_",i," <- Position")))
    #outputpath <- paste(path, outputdir, 'Position_1.RData', sep = '/')
    eval(parse(text=paste0("outputpath <- paste(path, outputdir, 'Position_",i,".RData', sep = '/')")))
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