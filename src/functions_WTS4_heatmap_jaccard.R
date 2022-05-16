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
  ghm <- ggplot(x, aes(x = CSPA_k, y = MCMI_k, fill = jaccard))
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
    axis.text = element_text(size = 60),
    legend.key.height = unit(2.5, "cm"),
    legend.key.width = unit(1.5, "cm"),
    legend.text = element_text(size = 60),
    legend.title = element_text(size = 60)
    )
  ghm <- ghm + scale_fill_viridis(na.value = "grey", direction = 1) # heatmap color is viridis http://www.okadajp.org/RWiki/?色見本
  ghm <- ghm + labs(x = "CSPA",
                    y = "MCMI")
  ghm <- ghm + labs(fill = "Jaccard")
  return(ghm)
}