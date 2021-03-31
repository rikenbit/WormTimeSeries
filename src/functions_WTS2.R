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
	return(ghm)
}