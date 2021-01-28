source("src/functions_WTS1.R")

# args setting
## args
##################################################
# # args = test
# args <- commandArgs(trailingOnly = T)
# test <- args
##################################################
## test 個体No.1の1つ目の細胞のグラフを描く
##################################################
# testargs 
test <- c("1","1")
# select celegans 個体の指定
args_celegans <- test[1]
# select cell 細胞の指定
args_cell <- test[2]
##################################################

# neuron activity data(matrix)
##################################################
# eval load('data/cleandata_mat/matrix_1.RData')
eval(parse(text = paste0("load('data/cleandata_mat/matrix_",args_celegans,".RData')")))
# eval n_matrix <- matrix_1
eval(parse(text = paste0("n_matrix <- matrix_",args_celegans)))

# Time
timeframe <- colnames(n_matrix)

# CellType
# eval cellType <- rownames(n_matrix)[1]
eval(parse(text = paste0("celltype <- rownames(n_matrix)[",args_cell,"]")))

# Neuron activity ：matrix  to 1cell dataframe
# nactivity <- n_matrix[1,]
eval(parse(text = paste0("nactivity <- n_matrix[",args_cell,",]")))
##################################################

# stimulation timing data(dataframe)
##################################################
# eval load('data/stimulation/stim_1.RData')
eval(parse(text = paste0("load('data/stimulation/stim_",args_celegans,".RData')")))

# StimTiming
# eval stim <- stim_1
eval(parse(text = paste0("stim <- stim_",args_celegans)))
stimtiming <- stim[,2]
##################################################

# dataframe for ggplot
data.frame(
        TimeFrame = timeframe,
        Nactivity= nactivity,
        StimTiming = stimtiming,
        stringsAsFactors = FALSE
) -> g

# ggplot
##################################################
# plot
p_0 <- ggplot(data = g, mapping = aes(x = timeframe, y = nactivity, group = 1))
p_1 <- p_0 +
    layer(geom = "point", stat = "identity", position = "identity")
p_1_2 <- p_1 +
  layer(geom = "line", stat = "identity", position = "identity")





# ggsave(filename = 'output/WTS1/1_1_X1.png', plot = gg, dpi = 100, width = 7.0, height = 7.0)
eval(parse(text = paste0("ggsave(filename = 'output/WTS1/",args_celegans,"_",args_cell,"_",celltype,".png', plot = p_1_2, dpi = 100, width = 21.0, height = 7.0)")))
##################################################

