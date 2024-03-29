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

# # WTS2 correlogram τ50
# #####################################################################################
# # ggAcf autocorrelation
# ##################################################
# g$Nactivity %>%
#     ggAcf(lag.max = 50, type = c("correlation"), plot = TRUE) -> p_Acf
# # title name
# eval(parse(text=paste0("title <- ggtitle('SampleNumber",args_sample,"_CellNumber",args_cell,"_",args_celltype,"_Acf')")))
# # title theme
# t_1 <- theme(plot.title = element_text(size = 20, hjust = 0.5))
# t_2 <- theme(axis.title = element_text(size = 16))

# gg <- p_Acf +
#     title +
#     t_1 +
#     t_2

# eval(parse(text=paste0("ggsave(filename = 'output/WTS2/ACFτ50/SampleNumber_",args_sample,"/CellNumber_",args_cell,"_CellType_",args_celltype,"_Acf.png', plot = gg, dpi = 100, width = 7.0, height = 7.0)")))
# ##################################################

# # ggAcf Partial autocorrelation
# ##################################################
# g$Nactivity %>%
#     ggAcf(lag.max = 50, type = c("partial"), plot = TRUE) -> p_pAcf
# # title name
# eval(parse(text=paste0("title <- ggtitle('SampleNumber",args_sample,"_CellNumber",args_cell,"_",args_celltype,"_pAcf')")))
# # title theme
# t_1 <- theme(plot.title = element_text(size = 20, hjust = 0.5))
# t_2 <- theme(axis.title = element_text(size = 16))

# gg <- p_pAcf +
#     title +
#     t_1 +
#     t_2

# eval(parse(text=paste0("ggsave(filename = 'output/WTS2/ACFτ50/SampleNumber_",args_sample,"/CellNumber_",args_cell,"_CellType_",args_celltype,"_pAcf.png', plot = gg, dpi = 100, width = 7.0, height = 7.0)")))
# ##################################################
# #####################################################################################

# # WTS2 correlogram τ500
# #####################################################################################
# # ggAcf autocorrelation
# ##################################################
# g$Nactivity %>%
#     ggAcf(lag.max = 500, type = c("correlation"), plot = TRUE) -> p_Acf
# # title name
# eval(parse(text=paste0("title <- ggtitle('SampleNumber",args_sample,"_CellNumber",args_cell,"_",args_celltype,"_Acf')")))
# # title theme
# t_1 <- theme(plot.title = element_text(size = 20, hjust = 0.5))
# t_2 <- theme(axis.title = element_text(size = 16))

# gg <- p_Acf +
#     title +
#     t_1 +
#     t_2 +
#     scale_x_continuous(breaks=seq(0,500,length=26),limits=c(0,500))

# eval(parse(text=paste0("ggsave(filename = 'output/WTS2/ACFτ500/SampleNumber_",args_sample,"/CellNumber_",args_cell,"_CellType_",args_celltype,"_Acf.png', plot = gg, dpi = 300, width = 14.0, height = 7.0)")))
# ##################################################

# # ggAcf Partial autocorrelation
# ##################################################
# g$Nactivity %>%
#     ggAcf(lag.max = 500, type = c("partial"), plot = TRUE) -> p_pAcf
# # title name
# eval(parse(text=paste0("title <- ggtitle('SampleNumber",args_sample,"_CellNumber",args_cell,"_",args_celltype,"_pAcf')")))
# # title theme
# t_1 <- theme(plot.title = element_text(size = 20, hjust = 0.5))
# t_2 <- theme(axis.title = element_text(size = 16))

# gg <- p_pAcf +
#     title +
#     t_1 +
#     t_2 +
#     scale_x_continuous(breaks=seq(0,500,length=26),limits=c(0,500))

# eval(parse(text=paste0("ggsave(filename = 'output/WTS2/ACFτ500/SampleNumber_",args_sample,"/CellNumber_",args_cell,"_CellType_",args_celltype,"_pAcf.png', plot = gg, dpi = 300, width = 14.0, height = 7.0)")))
# ##################################################
# #####################################################################################

# WTS2 correlogram stim τ50
#####################################################################################
# stimulation timing data(dataframe)
##################################################
eval(parse(text=paste0("load('data/stimulation/stim_",args_sample,".RData')")))
eval(parse(text=paste0("stim <- stim_",args_sample)))
stimtiming <- as.numeric(stim[,2])
##################################################
# dataframe
##################################################
data.frame(
        TimeFrame = timeframe,
        Nactivity = nactivity,
        StimTiming = stimtiming,
        stringsAsFactors = FALSE
) -> g
##################################################
# first stim TimeFrame
##################################################
g %>%
    filter(StimTiming != 0) %>%
        slice_head() %>%
            .$TimeFrame -> stimtiming_AF
stimtiming_BF <- stimtiming_AF -1
head(g,eval(parse(text=paste0(stimtiming_BF)))) -> g_BF
g %>%
    .[eval(parse(text=paste0(stimtiming_AF))):nrow(g),]  -> g_AF
##################################################
# ggAcf autocorrelation BF
##################################################
g_BF$Nactivity %>%
    ggAcf(lag.max = 50, type = c("correlation"), plot = TRUE) -> p_Acf
# title name
eval(parse(text=paste0("title <- ggtitle('SampleNumber",args_sample,"_CellNumber",args_cell,"_",args_celltype,"_Acf_BF')")))
# title theme
t_1 <- theme(plot.title = element_text(size = 20, hjust = 0.5))
t_2 <- theme(axis.title = element_text(size = 16))
gg <- p_Acf +
    title +
    t_1 +
    t_2
# ggsave
eval(parse(text=paste0("ggsave(filename = 'output/WTS2/ACFτ50_stim/ACF/SampleNumber_",args_sample,"/CellNumber_",args_cell,"_CellType_",args_celltype,"_BF.png', plot = gg, dpi = 100, width = 7.0, height = 7.0)")))
##################################################
# ggAcf autocorrelation AF
##################################################
g_AF$Nactivity %>%
    ggAcf(lag.max = 50, type = c("correlation"), plot = TRUE) -> p_Acf
# title name
eval(parse(text=paste0("title <- ggtitle('SampleNumber",args_sample,"_CellNumber",args_cell,"_",args_celltype,"_Acf_AF')")))
# title theme
t_1 <- theme(plot.title = element_text(size = 20, hjust = 0.5))
t_2 <- theme(axis.title = element_text(size = 16))
gg <- p_Acf +
    title +
    t_1 +
    t_2
# ggsave
eval(parse(text=paste0("ggsave(filename = 'output/WTS2/ACFτ50_stim/ACF/SampleNumber_",args_sample,"/CellNumber_",args_cell,"_CellType_",args_celltype,"_AF.png', plot = gg, dpi = 100, width = 7.0, height = 7.0)")))
##################################################
# ggAcf Partial autocorrelation BF
##################################################
g_BF$Nactivity %>%
    ggAcf(lag.max = 50, type = c("partial"), plot = TRUE) -> p_pAcf
# title name
eval(parse(text=paste0("title <- ggtitle('SampleNumber",args_sample,"_CellNumber",args_cell,"_",args_celltype,"_pAcf_BF')")))
# title theme
t_1 <- theme(plot.title = element_text(size = 20, hjust = 0.5))
t_2 <- theme(axis.title = element_text(size = 16))

gg <- p_pAcf +
    title +
    t_1 +
    t_2

eval(parse(text=paste0("ggsave(filename = 'output/WTS2/ACFτ50_stim/pACF/SampleNumber_",args_sample,"/CellNumber_",args_cell,"_CellType_",args_celltype,"_BF.png', plot = gg, dpi = 100, width = 7.0, height = 7.0)")))
##################################################
# ggAcf Partial autocorrelation AF
##################################################
g_AF$Nactivity %>%
    ggAcf(lag.max = 50, type = c("partial"), plot = TRUE) -> p_pAcf
# title name
eval(parse(text=paste0("title <- ggtitle('SampleNumber",args_sample,"_CellNumber",args_cell,"_",args_celltype,"_pAcf_AF')")))
# title theme
t_1 <- theme(plot.title = element_text(size = 20, hjust = 0.5))
t_2 <- theme(axis.title = element_text(size = 16))

gg <- p_pAcf +
    title +
    t_1 +
    t_2

eval(parse(text=paste0("ggsave(filename = 'output/WTS2/ACFτ50_stim/pACF/SampleNumber_",args_sample,"/CellNumber_",args_cell,"_CellType_",args_celltype,"_AF.png', plot = gg, dpi = 100, width = 7.0, height = 7.0)")))
##################################################
#####################################################################################