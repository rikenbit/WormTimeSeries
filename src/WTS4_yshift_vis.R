source("src/functions_WTS4_yshift_vis.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
args_input <- args[1]
args_output <- args[2]

#### test args####
# OUTPUT merge dataframe
# args_input <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/yShift_df/SampleNumber_1.RData")
# args_output <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/yShift_vis/SampleNumber_1.png")

args_stat_limit <- c("6000")
args_stat_limit <- as.numeric(args_stat_limit)

#1個体のvis
#1個体のデータフレームをload merge_df
load(args_input)
#ggplot
#### ggplot####
cord_x <- c("yshift")
cord_y <- c("1 - mSBD")
args_output |> 
    str_remove("output/WTS4/normalize_1/stimAfter/SBD_abs/yShift_vis/") |> 
    str_remove(".png") -> ggplot_title
gg <- ggplot(merge_df, 
                aes(x = yshift,
                    y = 1-mSBD, 
                    label = cell_cell
                )
                ) + 
    geom_point(size = 6.0, 
               alpha = 0.6) +
    # geom_label_repel(max.overlaps = Inf,
    #                  min.segment.length = 0,
    #                  size = 9.0,
    #                  force = 6.0) +# ラベル間の反発力
    theme(text = element_text(size = 60)) +
    labs(x = cord_x,
         y = cord_y,
         title = ggplot_title) +
    ylim(c(0, 1)) +
    xlim(c(-6000,6000)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text = element_text(size = 90))
#### ggsave####
ggsave(filename = args_output, 
       plot = gg,
       dpi = 50, 
       width = 50.0, 
       height = 50.0,
       limitsize = FALSE)