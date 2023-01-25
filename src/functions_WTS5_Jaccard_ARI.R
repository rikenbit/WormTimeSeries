#library
##################################################
library(tidyverse)
library(RColorBrewer)
library(viridis)

library(openxlsx) # read.xlsxを追加
library(svglite) # svgでの保存できるようにする
library(mclust) #ARI
##################################################
#### jaccard####
jaccard <- function(a, b) {
    intersection = length(intersect(a, b))
    union = length(a) + length(b) - intersection
    return(intersection/union)
}

ggplot_ghm = function(x) {
  ghm <- ggplot(x, aes(x = label, y = cls, fill = jaccard))
    # ghm <- ggplot(x, aes(x = label, y = cls, fill = ari))
    ghm <- ghm + geom_tile()
    ghm <- ghm + theme_bw()
    ghm <- ghm + theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    strip.background = element_rect(fill = "white", colour = "white"),
    plot.title = element_text(size = 60, hjust = 0.5),
    axis.title = element_text(size = 60),
    axis.text = element_text(size = 60, angle = 90),
    legend.key.height = unit(2.5, "cm"),
    legend.key.width = unit(1.5, "cm"),
    legend.text = element_text(size = 60),
    legend.title = element_text(size = 60)
    )
    ghm <- ghm + scale_fill_viridis(na.value = "grey", direction = 1) # heatmap color is viridis http://www.okadajp.org/RWiki/?色見本
    ghm <- ghm + labs(x = "",
                    y = "Cluster Number")
    ghm <- ghm + labs(fill = "Jaccard")
    # ghm <- ghm + labs(fill = "ARI")
    return(ghm)
}