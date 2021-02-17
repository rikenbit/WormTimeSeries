source("src/functions_WTS2.R")

# # args setting
# args <- commandArgs(trailingOnly = T)
# # select animal number 個体番号の指定
# args_sample <- args[1]
# # select cell number 細胞番号の指定
# args_cell <- args[2]
# # select celltype 細胞型名の指定
# args_celltype <- args[3]

# test args
##################################################
# select animal number 個体番号の指定
args_sample <- c("1")
# select cell number 細胞番号の指定
args_cell <- c("1")
# select celltype 細胞型名の指定
args_celltype <-c("X1")
##################################################

# load
##################################################
# load('data/cleandata_mat/matrix_1.RData')
eval(parse(text=paste0("load('data/cleandata_mat/matrix_",args_sample,".RData')")))
# n_matrix <- matrix_1
eval(parse(text=paste0("n_matrix <- matrix_",args_sample)))
# TimeFrame
timeframe <- as.numeric(colnames(n_matrix))
# NeuronActivity nactivity <- n_matrix[1,]
eval(parse(text=paste0("nactivity <- n_matrix[",args_cell,",]")))
##################################################


# dataframe
##################################################
data.frame(
        TimeFrame = timeframe,
        Nactivity = nactivity,
        stringsAsFactors = FALSE
) -> g
##################################################

# acf
##################################################
# acf(g$Nactivity, plot=T)
g$Nactivity %>% acf(plot=T, lag.max=50)
# library(magrittr)
# g$Nactivity %>% acf(plot = T) %T>% print()
##################################################
# plot
##################################################
# png('output/WTS2/SampleNumber_1/CellNumber_1_CellType_X1.png', width = 100, height = 100)
eval(parse(text=paste0("png('output/WTS2/SampleNumber_",args_sample,"/CellNumber_",args_cell,"_CellType_",args_celltype,".png', width = 500, height = 500)")))
g$Nactivity %>% acf(plot=T, lag.max=50) %T>% print()
dev.off() 
##################################################
# ggplot acf
##################################################
correlo <- acf(g$Nactivity, plot = FALSE)
correlo_df <- with(correlo, data.frame(lag, acf))
# gg <- ggplot(data=correlo_df, mapping=aes(x=lag, y=acf)) +
#        geom_bar(stat = "identity", position = "identity")
gg <- ggplot(data = correlo_df, mapping = aes(x = lag, y = acf)) +
       geom_hline(aes(yintercept = 0)) +
       geom_segment(mapping = aes(xend = lag, yend = 0))

gg <- ggplot(data = correlo_df, mapping = aes(x = lag, y = acf)) +
       geom_hline(aes(yintercept = 0)) +
       geom_segment(mapping = aes(xend = lag, yend = 0))
##################################################

# ggacf
##################################################
# ggacf
g$Nactivity %>%
    ggAcf(lag.max = 50, type = c("correlation"), plot = TRUE) -> p_acf
# title name
eval(parse(text=paste0("title <- ggtitle('SampleNumber",args_sample,"_CellNumber",args_cell,"_",args_celltype,"')")))
# title theme
t_1 <- theme(plot.title = element_text(size = 24, hjust = 0.5))

gg <- p_acf +
    title +
    t_1
##################################################