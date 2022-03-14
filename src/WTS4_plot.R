source("src/functions_WTS4_plot.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
# select animal number 個体番号の指定
args_sample <- args[1]
# select cell number 細胞番号の指定
args_cell <- args[2]
# select celltype 細胞型名の指定
args_celltype <- args[3]
# select datadir ディレクトリ 名の指定
args_datadir <- args[4]
# outputファイル名
args_output <- args[5]
#### test args####
# select animal number 個体番号の指定
args_sample <- c("")
# select cell number 細胞番号の指定
args_cell <- c("")
# select celltype 細胞型名の指定
args_celltype <- c("")
# select datadir ディレクトリ 名の指定
args_datadir <- c("")
# outputファイル名
args_output <- c("")
#### ggsave####