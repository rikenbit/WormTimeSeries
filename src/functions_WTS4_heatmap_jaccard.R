#library
##################################################
library(tidyverse)
library(RColorBrewer)
library(viridis)
##################################################
#### jaccard####
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

ggplot_ghm = function(x) {
  ghm <- ggplot(x, aes(x = col_celltype, y = row_celltype, fill = dist_value))
  ghm <- ghm + geom_tile()
  ghm <- ghm + theme_bw()
  ghm <- ghm + theme(
    # plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    strip.background = element_rect(fill = "white", colour = "white"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  ghm <- ghm + scale_fill_viridis(na.value = "grey", direction = -1) # heatmap color is viridis http://www.okadajp.org/RWiki/?色見本
  ghm <- ghm + scale_y_discrete(limits = rev(levels(x$row_celltype))) # reverse y axis
  return(ghm)
}