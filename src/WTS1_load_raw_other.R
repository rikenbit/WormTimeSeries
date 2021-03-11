source("src/functions_WTS1.R")

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