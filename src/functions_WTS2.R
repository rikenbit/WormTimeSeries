#library
##################################################
library(tidyverse)
library(forecast)
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(viridis)
library(igraph)
##################################################
ggplot_ghm = function(x) {
# ggplot_ghm = function(x, N_ReadData=ReadData, GroupName=groupname) {
	ghm <- ggplot(df, aes(x = cell_Receiver, y = cell_Sender, fill = CCF_ACF))
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
	ghm <- ghm + scale_fill_viridis(na.value = "white") #heatmap color is viridis
	# # ghm <- ghm + xlab("timeframe") + ylab(GroupName)

	# sX <- scale_x_continuous(name = "TimeFrame(1frame/0.2sec)",    # 軸の名前を変える
	#                          breaks = seq(0, nrow(N_ReadData), by= 1000),     # 軸の区切りを0,2,4にする
	# )
	# ghm <- ghm +
	#     sX
	return(ghm)
}