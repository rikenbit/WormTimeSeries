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
timeframe <- as.numeric(colnames(n_matrix))

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
stimtiming <- as.numeric(stim[,2])
##################################################

# dataframe for ggplot
data.frame(
        TimeFrame = timeframe,
        Nactivity = nactivity,
        StimTiming = stimtiming,
        stringsAsFactors = FALSE
) -> g

# ggplot
##################################################
# p_1 <- ggplot(data = g, mapping = aes(x = TimeFrame, y = Nactivity, group = 1)) +
#     geom_point() +
#     geom_line()
# p_2 <- ggplot(data = g, mapping = aes(x = TimeFrame, y = StimTiming, group = 1)) +
#     geom_point() +
#     geom_line()

# sX <- scale_x_continuous(name = "TimeFrame(1frame/0.2sec)",    # 軸の名前を変える
#                          breaks = seq(0, 6000, by=1000),     # 軸の区切りを0,2,4にする
#                          #labels = c("zero","two","four"), # 区切りを名付ける
#                          # limits = c(0,4)       # 0から4までしか表示しない
#                         )
# gg <- p_1 + 
#     p_2 +
#     sX


# p_1_2 <- ggplot(data = g, aes(TimeFrame)) +
#     geom_line(aes(y = Nactivity, colour = "Nactivity")) +
#     geom_line(aes(y = StimTiming, colour = "StimTiming")) +
#     scale_color_manual(values = c("black", "red"))
# gg <- p_1_2 +
#     sX

p_1 <- ggplot(data = g, aes(TimeFrame)) +
    geom_line(aes(y = Nactivity, colour = "Nactivity"))
p_1_2 <- p_1 + geom_line(aes(y = StimTiming, colour = "StimTiming") ,size = 1.5)
s_1 <- scale_color_manual(values = c("black", "red"))
sX <- scale_x_continuous(name = "TimeFrame(1frame/0.2sec)",    # 軸の名前を変える
                         breaks = seq(0, 6000, by= 1000),     # 軸の区切りを0,2,4にする
                         #labels = c("zero","two","four"), # 区切りを名付ける
                         # limits = c(0,4)       # 0から4までしか表示しない
                        )
# title <- ggtitle('celegans1_cell1_X1')
eval(parse(text = paste0("title <- ggtitle('celegans",args_celegans,"_cell",args_cell,"_",celltype,"')")))
t_1 <- theme(plot.title = element_text(size = 30, hjust = 0.5))
t_2 <- theme(axis.title = element_text(size = 20))
t_3 <- theme(legend.title = element_text(size = 28),
			 legend.text = element_text(size = 20))

gg <- p_1_2 +
    s_1 +
    sX +
    title +
    t_1 +
    t_2 +
    t_3 +
    labs(colour="each data")

# ggsave(filename = 'output/WTS1/celegans1/1_1_X1.png', plot = gg, dpi = 100, width = 7.0, height = 7.0)
eval(parse(text = paste0("ggsave(filename = 'output/WTS1/celegans",args_celegans,"/",args_celegans,"_",args_cell,"_",celltype,".png', plot = gg, dpi = 100, width = 21.0, height = 7.0)")))
##################################################