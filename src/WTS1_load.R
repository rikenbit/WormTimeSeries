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
        SampleNumber = n_sample,
        CellNumber = n_cell,
        CellType = celltype,
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

# raw CFP
####################################################################################
# Neuron Activity Data
##################################################
#List_Number_to_Samplename.csvを使って，読み込むエクセルを指定
#ファイル名一覧取得
path <- 'data/raw/'
NeuronList <- list.files(path, pattern='.xlsx') 

# animalname 7char
str_sub(NeuronList[1], start = 1, end = 7) <- read.xlsx(NeuronList[1], sheet = 'pi_k_Ch2', rowNames = TRUE, colNames =TRUE)

df_test <- read.xlsx("data/raw/180712e_4d_catv10_190824045149_edited_retrack=4847.xlsx",
                 sheet = "pi_k_Ch2", #読み込むシート名を指定できる。sheet = 1，のように番号でシートを指定することもできる。
                 rowNames = TRUE, #一行目を行名として扱う
                 colNames =TRUE) #一列目を列名として扱う

str_sub(NeuronList[i], start = 1, end = 7) <- read.xlsx(NeuronList[i], sheet = 'pi_k_Ch2', rowNames = TRUE, colNames =TRUE)


str_sub(NeuronList[1], start = 1, end = 7) %>% assign(sprintf(.), 1) 

eval(parse(text=paste0("assign('",str_sub(NeuronList[1], start = 1, end = 7),"',1)")))

# ファイル名一覧取得
path <- 'data/raw/'
NeuronList <- list.files(path, pattern='.xlsx')
# 読み込みオブジェクト用意（動的）
for (i in 1:3) {
  # v_name <- str_sub(NeuronList[i], start = 1, end = 7)
  eval(parse(text=paste0("v_name <- str_sub(NeuronList[",i,"], start = 1, end = 7)")))
  assign(v_name, i) #作成した変数名(v1,v2,v3)にnを格納
}

for (n in 1:3) {
  v_name <- c(str_sub(NeuronList[n], start = 1, end = 7)) #変数名作成(v1,v2,v3)
  assign(v_name, n) #作成した変数名(v1,v2,v3)にnを格納
}
# 読み込み


##################################################
# Sample Sheet
##################################################
##################################################
# Stimulation Data
##################################################
##################################################

####################################################################################