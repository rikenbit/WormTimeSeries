source("src/functions_WTS2.R")

# args setting
##################################################
args <- commandArgs(trailingOnly = T)
# select animal number 個体番号の指定
args_sample <- args[1]
# select cell number 細胞番号の指定
args_cell <- args[2]
# select celltype 細胞型名の指定
args_celltype <- args[3]
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

# ggAcf
##################################################
g$Nactivity %>%
    ggAcf(lag.max = 50, type = c("correlation"), plot = TRUE) -> p_Acf
# title name
eval(parse(text=paste0("title <- ggtitle('SampleNumber",args_sample,"_CellNumber",args_cell,"_",args_celltype,"_Acf')")))
# title theme
t_1 <- theme(plot.title = element_text(size = 24, hjust = 0.5))
t_2 <- theme(axis.title = element_text(size = 16))

gg <- p_Acf +
    title +
    t_1+
    t_2

eval(parse(text=paste0("ggsave(filename = 'output/WTS2/SampleNumber_",args_sample,"/CellNumber_",args_cell,"_CellType_",args_celltype,"_Acf.png', plot = gg, dpi = 100, width = 7.0, height = 7.0)")))
##################################################

# ggAcf partial
##################################################
g$Nactivity %>%
    ggAcf(lag.max = 50, type = c("partial"), plot = TRUE) -> p_pAcf
# title name
eval(parse(text=paste0("title <- ggtitle('SampleNumber",args_sample,"_CellNumber",args_cell,"_",args_celltype,"_pAcf')")))
# title theme
t_1 <- theme(plot.title = element_text(size = 24, hjust = 0.5))
t_2 <- theme(axis.title = element_text(size = 16))

gg <- p_pAcf +
    title +
    t_1+
    t_2

eval(parse(text=paste0("ggsave(filename = 'output/WTS2/SampleNumber_",args_sample,"/CellNumber_",args_cell,"_CellType_",args_celltype,"_pAcf.png', plot = gg, dpi = 100, width = 7.0, height = 7.0)")))
##################################################