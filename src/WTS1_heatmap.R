source("src/functions_WTS1.R")

#### 6data All 28sample ####
#### args####
args <- commandArgs(trailingOnly = T)
# select animal number 個体番号の指定
args_sample <- args[1]
# select datadir ディレクトリ 名の指定
args_datadir <- args[2]
# outputファイル名
args_output <- args[3]
########################


#### Neuron Activity####
#### load matrix####
path <- "data"
inputdir <- args_datadir
eval(parse(text=paste0("inputpath <- paste(path, inputdir, 'ReadData_",args_sample,".RData', sep = '/')")))
load(inputpath)
# ReadData <- ReadData_1
eval(parse(text=paste0("ReadData <- ReadData_",args_sample)))

#### convert wider to longer####
ReadData %>%
    rowid_to_column("timeframe") %>%
        pivot_longer(-timeframe, names_to = "celltype", values_to = "nactivity") -> df

#### merge igraph####
# load igraphdata
load("data/igraph/Fig1_HNS.RData")
# node information convert dataframe
ig_Fig1_HNS %>% 
    igraph::as_data_frame(., what="vertices") -> df_node
# remove space from colnames
df_node %>% 
    names() %>% 
        str_replace_all(., c(" " = "")) -> names(df_node)
# add node information
merge(df, df_node, by.x = "celltype", by.y = "name",all.x = TRUE) -> df_merged
df_merged %>% .$NeuronType %>% unique() %>% na.omit() -> NeuronType_Name

#### hit NA####
# create group df
df_merged %>% 
    filter(.,is.na(NeuronType)) -> df_onegroup
# filter matrix
df_onegroup %>% 
    .$celltype %>% 
        unique()-> v_celltype
ReadData[,v_celltype] %>% 
    as.matrix() -> ReadData_NA
# clustering
ReadData_NA %>%
    heatmap(., Rowv = NA, scale = "none", 
            hclustfun = function(x) { hclust(x, method = "ward.D2")}) -> clr
# reorder
celltype.idx <- colnames(ReadData_NA)[clr$colInd]
df_onegroup$celltype <- factor(df_onegroup$celltype, levels = celltype.idx)
# ggplot
groupname <- "other"
ghm <- ggplot_ghm(df_onegroup) + ggtitle("NeuronActivity") +theme(plot.title = element_text(size = 30, hjust = 0.5))
# NeuronType length NA
df_merged %>%
    filter(.,is.na(NeuronType)) %>% 
      .$celltype %>%
          unique() %>% length() -> NeuronType_length

#### hit Neuron Type#### 
ReadData_flist <- list()
for (i in 1:length(NeuronType_Name)) {
    # create group df
    df_merged %>% 
        filter(.,NeuronType==NeuronType_Name[i]) -> df_onegroup
    # filter matrix
    df_onegroup %>% 
        .$celltype %>% 
            unique() -> v_celltype
    ReadData[,v_celltype] %>% 
        as.matrix() -> ReadData_f
    append(ReadData_flist, list(ReadData_f)) -> ReadData_flist
    # clustering
    ReadData_flist[[i]] %>%
        heatmap(., Rowv = NA, scale = "none", 
                hclustfun = function(x) { hclust(x, method = "ward.D2")}) -> clr
    # reorder
    celltype.idx <- colnames(ReadData_flist[[i]])[clr$colInd]
    df_onegroup$celltype <- factor(df_onegroup$celltype, levels = celltype.idx)
    # ggplot
    groupname <- NeuronType_Name[i]
    ghm <- ghm / ggplot_ghm(df_onegroup)
    # NeuronType length
    df_onegroup$celltype %>%
        unique() %>% length() %>%
            c(NeuronType_length,.) -> NeuronType_length
}
########################


#### mCherry####
#### load####
# Neuron Activityデータのtimeframeの範囲で，mCherryのtimeframeを取得
eval(parse(text=paste0("rownames(ReadData_",args_sample,") %>% as.numeric() -> timeframe")))
path <- "data"
inputdir <- "mCherry"
# inputpath <- paste(path, inputdir, 'mCherry_1.RData', sep = '/')
eval(parse(text=paste0("inputpath <- paste(path, inputdir, 'mCherry_",args_sample,".RData', sep = '/')")))
load(inputpath)
# cut TimeFrame
eval(parse(text=paste0("mCherry_",args_sample,"[1:length(timeframe),] -> mCherry")))

#### convert wider to longer####
mCherry %>%
  rowid_to_column("timeframe") %>%
      pivot_longer(-timeframe, names_to = "celltype", values_to = "mCherry") -> df

#### merge igraph####
# add node information
merge(df, df_node, by.x = "celltype", by.y = "name",all.x = TRUE) -> df_merged
df_merged %>% .$NeuronType %>% unique() %>% na.omit() -> NeuronType_Name

#### hit NA####
# create group df
df_merged %>% 
  filter(.,is.na(NeuronType)) -> df_onegroup
# filter matrix
# clustering
ReadData_NA %>%
  heatmap(., Rowv = NA, scale = "none", 
          hclustfun = function(x) { hclust(x, method = "ward.D2")}) -> clr
# reorder
celltype.idx <- colnames(ReadData_NA)[clr$colInd]
df_onegroup$celltype <- factor(df_onegroup$celltype, levels = celltype.idx)
# ggplot
groupname <- "other"
ghm_m <- ggplot_ghm_m(df_onegroup) + ggtitle("mCherry") +theme(plot.title = element_text(size = 30, hjust = 0.5))

#### hit Neuron Type#### 
for (i in 1:length(NeuronType_Name)) {
  # create group df
  df_merged %>% 
    filter(.,NeuronType==NeuronType_Name[i]) -> df_onegroup
  # clustering
  ReadData_flist[[i]] %>%
    heatmap(., Rowv = NA, scale = "none", 
            hclustfun = function(x) { hclust(x, method = "ward.D2")}) -> clr
  # reorder
  celltype.idx <- colnames(ReadData_flist[[i]])[clr$colInd]
  df_onegroup$celltype <- factor(df_onegroup$celltype, levels = celltype.idx)
  # ggplot
  groupname <- NeuronType_Name[i]
  ghm_m <- ghm_m / ggplot_ghm_m(df_onegroup)
}
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


#### patchwork####
ggp <- (gg3 + plot_spacer()) /
    (gg5 + plot_spacer()) / 
    (ghm + plot_layout(heights = NeuronType_length) | 
       ghm_m + plot_layout(heights = NeuronType_length)) + 
    plot_layout(heights = c(1, 1, 5))
ggsave(filename = args_output, plot = ggp, dpi = 100, width = 30.0, height = 40.0)
########################