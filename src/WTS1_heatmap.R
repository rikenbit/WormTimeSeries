source("src/functions_WTS1.R")

#### test 1data 1sample ####
#### test args####
# args <- commandArgs(trailingOnly = T)

# select animal number 個体番号の指定
args_sample <- c("1")
# select cell number 細胞番号の指定
args_cell <- c("1")
# select celltype 細胞型名の指定
args_celltype <- c("102")
# select datadir ディレクトリ 名の指定
args_datadir <- c("raw_CFP")
# outputファイル名
args_output <- c("output/WTS1/heatmap/raw_CFP/SampleNumber_1.png")

#### Neuron Activity####
#### load matrix####
path <- "data"
inputdir <- args_datadir
eval(parse(text=paste0("inputpath <- paste(path, inputdir, 'ReadData_",args_sample,".RData', sep = '/')")))
load(inputpath)

#### convert wider to longer####
# ReadData <- ReadData_1
eval(parse(text=paste0("ReadData <- ReadData_",args_sample)))
ReadData %>%
    rowid_to_column("timeframe") %>%
        pivot_longer(-timeframe, names_to = "celltype", values_to = "nactivity") -> df

#### ggplot####
ghm <- ggplot(df, aes(x = timeframe, y = celltype, fill = nactivity))
ghm <- ghm + geom_tile()
ghm <- ghm + theme_bw()
ghm <- ghm + theme(plot.background = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.background = element_blank(),
                   axis.line = element_blank(),
                   axis.ticks = element_blank(),
                   strip.background = element_rect(fill = "white", colour = "white"),
                   axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ghm <- ghm + scale_fill_gradientn("value", colours = rev(brewer.pal(9, "Spectral")), na.value = "white")
ghm <- ghm + xlab("timeframe") + ylab("celltype")

title <- ggtitle("NeuronActivity")
t_1 <- theme(plot.title = element_text(size = 30, hjust = 0.5))
sX <- scale_x_continuous(name = "TimeFrame(1frame/0.2sec)",    # 軸の名前を変える
                         breaks = seq(0, nrow(ReadData), by= 1000),     # 軸の区切りを0,2,4にする
)
ghm <- ghm +
    title +
    t_1 +
    sX
########################

#### Stimulation & Position####
#### load Stimulation####
eval(parse(text=paste0("rownames(ReadData_",args_sample,") %>% as.numeric() -> timeframe")))
path <- "data"
inputdir <- "stimulation"
# inputpath <- paste(path, inputdir, 'ReadData_1.RData', sep = '/')
eval(parse(text=paste0("inputpath <- paste(path, inputdir, 'stim_",args_sample,".RData', sep = '/')")))
load(inputpath)
eval(parse(text=paste0("stim_",args_sample," %>% as.numeric() -> stimtiming")))
# cut TimeFrame
stimtiming[1:length(timeframe)] -> stimtiming

#### load Position####
path <- "data"
inputdir <- "Position"
# inputpath <- paste(path, inputdir, 'Position_1.RData', sep = '/')
eval(parse(text=paste0("inputpath <- paste(path, inputdir, 'Position_",args_sample,".RData', sep = '/')")))
load(inputpath)
# Position_1$MoveX %>% as.numeric() -> position
eval(parse(text=paste0("Position_",args_sample,"$MoveX %>% as.numeric() -> position")))
# cut TimeFrame
position[1:length(timeframe)] -> position

#### dataframe for ggplot####
data.frame(
  TimeFrame = timeframe,
  StimTiming = stimtiming,
  Position = position,
  stringsAsFactors = FALSE
) -> g

#### ggplot####
p_1 <- ggplot(data = g, aes(TimeFrame))
p_3 <- p_1 +        
  geom_line(aes(y = StimTiming, colour = "StimTiming") , size = 1.5)
p_5 <- p_1 +        
  geom_line(aes(y = Position, colour = "Position") , size = 1.5)
sX <- scale_x_continuous(name = "TimeFrame(1frame/0.2sec)",    # 軸の名前を変える
                         breaks = seq(0, length(timeframe), by= 1000),     # 軸の区切りを0,2,4にする
)
s_3 <- scale_color_manual(values = c("purple"))
s_5 <- scale_color_manual(values = c("blue"))
# title <- ggtitle('SampleNumber1_raw_CFP')
eval(parse(text=paste0("title <- ggtitle('SampleNumber",args_sample,"_",args_datadir,"')")))
t_1 <- theme(plot.title = element_text(size = 30, hjust = 0.5))
t_2 <- theme(axis.title = element_text(size = 20))
t_3 <- theme(legend.title = element_text(size = 28),
             legend.text = element_text(size = 20))
gg3 <- p_3 +
  s_3 +
  sX +
  title +
  t_1 +
  t_2 +
  t_3 +
  labs(colour="each data")
gg5 <- p_5 +
  s_5 +
  sX +
  t_2 +
  t_3 +
  labs(colour="each data")
########################

#### mCherry####
#### load####
eval(parse(text=paste0("rownames(ReadData_",args_sample,") %>% as.numeric() -> timeframe")))
path <- "data"
inputdir <- "mCherry"

# inputpath <- paste(path, inputdir, 'mCherry_1.RData', sep = '/')
eval(parse(text=paste0("inputpath <- paste(path, inputdir, 'mCherry_",args_sample,".RData', sep = '/')")))
load(inputpath)
# cut TimeFrame
eval(parse(text=paste0("mCherry_1[",args_sample,":length(timeframe),] -> mCherry")))

mCherry %>%
  rowid_to_column("timeframe") %>%
      pivot_longer(-timeframe, names_to = "celltype", values_to = "mcherry") -> df

#### ggplot####
ghm_m <- ggplot(df, aes(x = timeframe, y = celltype, fill = mcherry))
ghm_m <- ghm_m + geom_tile()
ghm_m <- ghm_m + theme_bw()
ghm_m <- ghm_m + theme(plot.background = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.background = element_blank(),
                   axis.line = element_blank(),
                   axis.ticks = element_blank(),
                   strip.background = element_rect(fill = "white", colour = "white"),
                   axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ghm_m <- ghm_m + scale_fill_gradientn("value", colours = rev(brewer.pal(9, "Spectral")), na.value = "white")
ghm_m <- ghm_m + xlab("timeframe") + ylab("celltype")
title <- ggtitle("mCherry")
t_1 <- theme(plot.title = element_text(size = 30, hjust = 0.5))
sX <- scale_x_continuous(name = "TimeFrame(1frame/0.2sec)",    # 軸の名前を変える
                         breaks = seq(0, nrow(ReadData), by= 1000),     # 軸の区切りを0,2,4にする
)
ghm_m <- ghm_m +
  title +
  t_1 +
  sX
########################

#### patchwork####
ggp <- (gg3 + plot_spacer()) / (gg5 + plot_spacer()) / (ghm + ghm_m) + plot_layout(heights = c(1, 1, 5))
ggsave(filename = args_output, plot = ggp, dpi = 100, width = 30.0, height = 30.0)
########################