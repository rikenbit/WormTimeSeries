#library
##################################################
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(openxlsx)
##################################################
ggplot_ghm <- function(x) {
    ghm <- ggplot(x, aes(x = col_celltype, y = row_celltype, fill = value))
    ghm <- ghm + geom_tile()
    ghm <- ghm + theme_bw()
    ghm <- ghm + theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_rect(fill = "white", colour = "white"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    ghm <- ghm + scale_y_discrete(limits = colnames(input_mat_F_S))
    ghm <- ghm + scale_x_discrete(limits = colnames(input_mat_F_S))
    return(ghm)
}

con_heatmap_df <- function(x) {
    x |> 
        as.data.frame() |> 
            rownames_to_column("row_celltype") |> 
                pivot_longer(-row_celltype, 
                             names_to = "col_celltype", 
                             values_to = "value") -> long_celltype
    return(long_celltype)
}

vis_z <- function(x) {
    x |> 
        con_heatmap_df() |> 
            ggplot_ghm() -> ghm
    ghm <- ghm + scale_fill_viridis(na.value = "grey", direction = 1) 
    return(ghm)
}

vis_abs <- function(x) {
    x |> 
        abs() |> 
            con_heatmap_df() |> 
                ggplot_ghm() -> ghm
    ghm <- ghm + scale_fill_viridis(na.value = "grey", direction = 1, limits=c(0, args_thrr))
    # if (args_stat=="count") {
    #     ghm <- ghm + scale_fill_viridis(na.value = "grey", direction = 1, limits=c(0, args_stat_limit)) 
    # } else {
    #     ghm <- ghm + scale_fill_viridis(na.value = "grey", direction = -1, limits=c(0, args_stat_limit))
    # }
    # return(ghm)
}

vis_ghm <- list(
    "zahlen" = vis_z,
    "abs" = vis_abs
    )

vis_labels <- list(
    "ALL" = c("NaCl", "PC1_neg", "PC1_pos", "PC2", "PC3"),
    "NaCl" = c("NaCl"),
    "1n" = c("PC1_neg"),
    "1p" = c("PC1_pos"),
    "1np" = c("PC1_neg", "PC1_pos")
    )