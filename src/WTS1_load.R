source("src/functions_WTS1.R")

# load & save cleandata
##################################################
# number of c.elegans file
list.files("data/cleandata/", all.files = TRUE, recursive = TRUE) %>%
    length() %>%
        seq() -> celegans
# list of load files
loadfiles <- paste0("data/cleandata/", celegans, "_ratio.csv")
# TimePoints 0.2sec/Point
Time <- seq(1, 6000, 1)

# load csv 
# Rの仕様で数字の細胞型名は先頭に「X」がついた名前に変換される
for(i in celegans){
	eval(parse(text = paste0("data_",i," <- read.csv(loadfiles[",i,"], header=TRUE)")))
}
# Manual Fix sample2 Cellname　FLP→FLPR
colnames(data_2)[match("FLP",colnames(data_2))] <- c("FLPR")

for(i in celegans){
    eval(parse(text = paste0("data_",i," %>% as.matrix() -> mat_data_",i)))
    eval(parse(text = paste0("Time -> rownames(mat_data_",i,")")))
    eval(parse(text = paste0("t(mat_data_",i,") -> matrix_",i)))
    eval(parse(text = paste0("save(matrix_",i,", file ='data/cleandata_mat/matrix_",i,".RData')")))
}
##################################################



# load stimulation_timing.xlsx
##################################################
##################################################