source("src/functions_WTS4_yshift_vis_all.R")
#### args setting####
args <- commandArgs(trailingOnly = T)
args_input <- args[1]
args_output <- args[2]
#### test args####
# args_input <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Membership_df/k_Number_6/DFs.RData")
# args_output <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/yShift_vis_all/k_Number_6.png")

#### load####
load(args_input)
#### yshift list to df####
# DFs_yshift <- purrr::map(names(DFs), .DFs_yshift)
dfr_yshift <- purrr::map_dfr(names(DFs), .DFs_yshift)

#### ggplot####
cord_x <- c("yshift")
cord_y <- c("1 - mSBD")
ggplot_title <- c("yShift_vis_all")
gg <- ggplot(dfr_yshift, 
             aes(x = yshift,
                 y = 1-mSBD, 
                 label = cell_cell,
                 color = factor(animal, levels = names(DFs))
                 )
             ) + 
    labs(color = "Animal No.") +
    geom_point(size = 6.0, 
               alpha = 0.6) +
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
       width = 100.0, 
       height = 50.0,
       limitsize = FALSE)