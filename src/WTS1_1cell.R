source("src/functions_WTS1.R")

# args setting
args <- commandArgs(trailingOnly = T)

# select animal number 個体番号の指定
args_sample <- args[1]
# select cell number 細胞番号の指定
args_cell <- args[2]
# select celltype 細胞型名の指定
args_celltype <- args[3]

# neuron activity data(matrix)
##################################################
# load('data/cleandata_mat/matrix_1.RData')
eval(parse(text=paste0("load('data/cleandata_mat/matrix_",args_sample,".RData')")))
# n_matrix <- matrix_1
eval(parse(text=paste0("n_matrix <- matrix_",args_sample)))
# Time
timeframe <- as.numeric(colnames(n_matrix))
# Neuron activity ：matrix  to 1cell dataframe
# nactivity <- n_matrix[1,]
eval(parse(text=paste0("nactivity <- n_matrix[",args_cell,",]")))
##################################################

# stimulation timing data(dataframe)
##################################################
# load('data/stimulation/stim_1.RData')
eval(parse(text=paste0("load('data/stimulation/stim_",args_sample,".RData')")))
# StimTiming
# stim <- stim_1
eval(parse(text=paste0("stim <- stim_",args_sample)))
stimtiming <- as.numeric(stim[,2])
##################################################

# dataframe for ggplot
##################################################
data.frame(
        TimeFrame = timeframe,
        Nactivity = nactivity,
        StimTiming = stimtiming,
        stringsAsFactors = FALSE
) -> g

# rollmean
g %>% mutate(N_roll = roll_meanr(Nactivity, n=51, align="right", fill=NA)) -> g_roll

# diff
diff_value <- 50
n_diff <- append(rep(NA, diff_value), diff(g$Nactivity, diff_value))
g_roll %>% 
    mutate(N_diff = n_diff) -> g_roll_diff
##################################################

## ggplot
##################################################
p_1 <- ggplot(data = g_roll_diff, aes(TimeFrame)) +
    geom_line(aes(y = Nactivity, colour = "Nactivity"))
p_2 <- p_1 +
         geom_line(aes(y = StimTiming, colour = "StimTiming") , size = 1.5, alpha = 0.7) +
         geom_line(aes(y = N_roll, colour = "N_rollmean"), size = 1.5, alpha = 0.7) +
         geom_line(aes(y = N_diff, colour = "N_diff"), size = 1.5, alpha = 0.4, linetype = "dotted")

s_1 <- scale_color_manual(values = c("blue", "red" , "black", "green"))
sX <- scale_x_continuous(name = "TimeFrame(1frame/0.2sec)",    # 軸の名前を変える
                         breaks = seq(0, 6000, by= 1000),     # 軸の区切りを0,2,4にする
                        )

# title <- ggtitle('SampleNumber1_CellNumber1_X1')
eval(parse(text=paste0("title <- ggtitle('SampleNumber",args_sample,"_CellNumber",args_cell,"_",args_celltype,"')")))
t_1 <- theme(plot.title = element_text(size = 30, hjust = 0.5))
t_2 <- theme(axis.title = element_text(size = 20))
t_3 <- theme(legend.title = element_text(size = 28),
             legend.text = element_text(size = 20))

gg <- p_2 +
    s_1 +
    sX +
    title +
    t_1 +
    t_2 +
    t_3 +
    labs(colour="each data")

eval(parse(text=paste0("ggsave(filename = 'output/WTS1/Samplenumber",args_sample,"/Cellnumber",args_cell,"_Celltype",args_celltype,".png', plot = gg, dpi = 100, width = 21.0, height = 7.0)")))
##################################################