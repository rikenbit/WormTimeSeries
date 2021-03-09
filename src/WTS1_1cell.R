source("src/functions_WTS1.R")

# raw CFP
####################################################################################
# args setting
##################################################
# args <- commandArgs(trailingOnly = T)
# # select animal number 個体番号の指定
# args_sample <- args[1]
# # select cell number 細胞番号の指定
# args_cell <- args[2]
# # select celltype 細胞型名の指定
# args_celltype <- args[3]

# select animal number 個体番号の指定
args_sample <- c("1")
# select cell number 細胞番号の指定
args_cell <- c("1")
# select celltype 細胞型名の指定
args_celltype <-c("X1")
##################################################

# Neuron Activity Data
##################################################
path <- "data"
inputdir <- "raw_CFP"

# inputpath <- paste(path, inputdir, 'ReadData_1.RData', sep = '/')
eval(parse(text=paste0("inputpath <- paste(path, inputdir, 'ReadData_",args_sample,".RData', sep = '/')")))
load(inputpath)
# ReadData_1[,1] %>% as.numeric() -> nactivity
eval(parse(text=paste0("ReadData_",args_sample,"[,",args_cell,"] %>% as.numeric() -> nactivity")))

#TimeFrame
# rownames(ReadData_1) %>% as.numeric() -> timeframe
eval(parse(text=paste0("rownames(ReadData_",args_sample,") %>% as.numeric() -> timeframe")))
##################################################

# Stimulation Data
##################################################
path <- "data"
inputdir <- "stimulation"

# inputpath <- paste(datapath, 'Stim_1.RData', sep = '/')
eval(parse(text=paste0("inputpath <- paste(path, inputdir, 'Stim_",args_sample,".RData', sep = '/')")))
load(inputpath)
# Stim_1 %>% as.numeric() -> stimtiming
eval(parse(text=paste0("Stim_",args_sample," %>% as.numeric() -> stimtiming")))
##################################################
# mCherry
##################################################
path <- "data"
inputdir <- "mCherry"

# inputpath <- paste(path, inputdir, 'mCherry_1.RData', sep = '/')
eval(parse(text=paste0("inputpath <- paste(path, inputdir, 'mCherry_",args_sample,".RData', sep = '/')")))
load(inputpath)
# mCherry_1[,1] %>% as.numeric() -> mcherry
eval(parse(text=paste0("mCherry_",args_sample,"[,",args_cell,"] %>% as.numeric() -> mcherry")))
##################################################
# Position
##################################################
path <- "data"
inputdir <- "Position"

# inputpath <- paste(path, inputdir, 'Position_1.RData', sep = '/')
eval(parse(text=paste0("inputpath <- paste(path, inputdir, 'Position_",args_sample,".RData', sep = '/')")))
load(inputpath)
# Position_1$MoveX %>% as.numeric() -> position
eval(parse(text=paste0("Position_",args_sample,"$MoveX %>% as.numeric() -> position")))
##################################################

# dataframe for ggplot
##################################################
data.frame(
        TimeFrame = timeframe,
        Nactivity = nactivity,
        StimTiming = stimtiming,
        mCherry = mcherry,
        Position = position,
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
p_1 <- ggplot(data = g_roll_diff, aes(TimeFrame))
p_2 <- p_1 +
        geom_line(aes(y = Nactivity, colour = "Nactivity")) +
         # geom_line(aes(y = StimTiming, colour = "StimTiming") , size = 1.5, alpha = 0.7) +
        geom_line(aes(y = N_roll, colour = "N_rollmean"), size = 1.5, alpha = 0.7) +
        geom_line(aes(y = N_diff, colour = "N_diff"), size = 1.5, alpha = 0.4, linetype = "dotted")

s_1 <- scale_color_manual(values = c("orange", "green", "black"))
sX <- scale_x_continuous(name = "TimeFrame(1frame/0.2sec)",    # 軸の名前を変える
                         breaks = seq(0, 6000, by= 1000),     # 軸の区切りを0,2,4にする
                        )


outputdir <- "raw_CFP"
# title <- ggtitle('SampleNumber1_CellNumber1_X1_raw_CFP')
eval(parse(text=paste0("title <- ggtitle('SampleNumber",args_sample,"_CellNumber",args_cell,"_",args_celltype,"_",outputdir,"')")))
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
##################################################

## ggsave
##################################################
path <- "output/WTS1/plot"
outputdir <- "raw_CFP"
# outputpath <- paste(path, outputdir, 'SampleNumber_",args_sample,"', 'CellNumber_",args_cell,"_CellType_",args_celltype,".png', sep = '/')
eval(parse(text=paste0("outputpath <- paste(path, outputdir, 'SampleNumber_",args_sample,"', 'CellNumber_",args_cell,"_CellType_",args_celltype,".png', sep = '/')")))
ggsave(filename = outputpath, plot = gg, dpi = 100, width = 21.0, height = 7.0)
##################################################