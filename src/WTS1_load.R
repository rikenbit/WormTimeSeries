source("src/functions_WTS1.R")

# load and export neuron activity data
##################################################
# number of c.elegans neuron activity file
list.files("data/cleandata/", all.files=TRUE, recursive=TRUE) %>%
    length() %>%
        seq() -> celegans
# list of load files
loadfiles <- paste0("data/cleandata/", celegans, "_ratio.csv")
# TimePoints 0.2sec/Point
Time <- seq(1, 6000, 1)

# load csv 
# Rの仕様で数字の細胞型名は先頭に「X」がついた名前に変換される
for(i in celegans){
	eval(parse(text=paste0("data_",i," <- read.csv(loadfiles[",i,"], header=TRUE)")))
}
# Manual Fix sample2 Cellname　FLP→FLPR
colnames(data_2)[match("FLP",colnames(data_2))] <- c("FLPR")

for(i in celegans){
    eval(parse(text=paste0("data_",i," %>% as.matrix() -> mat_data_",i)))
    eval(parse(text=paste0("Time -> rownames(mat_data_",i,")")))
    eval(parse(text=paste0("t(mat_data_",i,") -> matrix_",i)))
    eval(parse(text=paste0("save(matrix_",i,", file ='data/cleandata_mat/matrix_",i,".RData')")))
}
##################################################

# export sample sheet
##################################################
# 行列データをリストにまとめた
matrix_list <- list()
for (i in 1:15) {
    eval(parse(text=paste0("matrix_list <- c(matrix_list, list(matrix_",i,"))")))
}
# 個体No.，細胞No.，細胞型名のデータフレーム作成
n_sample <- c()
n_cell <- c()
celltype <- c()
for (i in 1:15) {
    # create Sample.number
    # rep(celegans[1],nrow(matrix_list[[1]])) %>% as.character()
    rep(celegans[i],nrow(matrix_list[[i]])) %>%
        as.character() %>%
            append(n_sample, .) -> n_sample
    # create Cell.number
    # seq(nrow(matrix_1))
    seq(nrow(matrix_list[[i]])) %>%
        as.character() %>%
            append(n_cell, .) -> n_cell
    # create Cell.type
    rownames(matrix_list[[i]]) %>%
            append(celltype, .) -> celltype
}
# create dataframe
data.frame(
        Sample.number = n_sample,
        Cell.number = n_cell,
        Cell.type = celltype,
        stringsAsFactors = FALSE
) -> sample_sheet

# export
write.csv(sample_sheet, "data/WTS1_sample_sheet.csv", row.names=FALSE)
##################################################


# load and export stimulation timing data
##################################################
# load csv 
stimdata <- read.csv("data/stimulation_timing.csv")

# stim celegans1~15
# １列目のタイムフレームを削除したデータフレーム作成
stimdata %>% select(., -1) -> stim15

# timeframe
timeframe <- stimdata[,1]

# save each stimulation timing
for(i in celegans){
    # stimtiming <- stim15[,1]
    eval(parse(text=paste0("stimtiming <- stim15[,",i,"]")))
    # dataframe
    data.frame(
        TimeFrame = timeframe,
        StimTiming = stimtiming,
        stringsAsFactors = FALSE
    ) -> stimdf
    # stimdf -> stim_1
    eval(parse(text=paste0("stimdf -> stim_",i)))
    # save(stim_1, file ='data/stimulation/stim_1.RData')
    eval(parse(text=paste0("save(stim_",i,", file ='data/stimulation/stim_",i,".RData')")))
}
##################################################