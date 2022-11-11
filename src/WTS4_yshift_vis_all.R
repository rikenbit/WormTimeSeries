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

# #### ggplot####
# cord_x <- c("yshift")
# cord_y <- c("1 - mSBD")
# ggplot_title <- c("Animal_all")
# gg <- ggplot(dfr_yshift, 
#              aes(x = yshift,
#                  y = 1-mSBD,
#                  # label = cell_cell,
#                  # color = factor(animal, levels = names(DFs))
#                  )
#              ) + 
#     # labs(color = "Animal No.") +
#     guides(colour="none") +
#     geom_point(size = 3.0, 
#                alpha = 0.6) +
#     # theme(text = element_text(size = 60)) +
#     labs(x = cord_x,
#          y = cord_y,
#          title = ggplot_title) +
#     ylim(c(0, 1)) +
#     xlim(c(-6000,6000)) +
#     theme(plot.title = element_text(hjust = 0.5)) +
#     theme(text = element_text(size = 90))

#### ggplot density group####
cord_x <- c("yshift")
cord_y <- c("density")
ggplot_title <- c("Animal_all")
gg <- ggplot(dfr_yshift,
             aes(x = yshift,
                 fill=factor(animal, levels = names(DFs)))
             ) +
    geom_density(alpha=0.2) +
    labs(x = cord_x,
         y = cord_y,
         title = ggplot_title,
         fill = "Animal_No.") +
    ylim(c(0, 0.05)) +
    xlim(c(-1000, 1000)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text = element_text(size = 90))
# #### ggplot density####
# cord_x <- c("yshift")
# cord_y <- c("density")
# ggplot_title <- c("Animal_all")
# gg <- ggplot(dfr_yshift,
#              aes(x = yshift)
#              ) +
#     guides(colour="none") +
#     geom_density(fill="blue", alpha=0.2) +
#     labs(x = cord_x,
#          y = cord_y,
#          title = ggplot_title) +
#     ylim(c(0, 0.015)) +
#     xlim(c(-1000, 1000)) +
#     theme(plot.title = element_text(hjust = 0.5)) +
#     theme(text = element_text(size = 90))

# #### ggplot density 2d####
# # https://mukkujohn.hatenablog.com/entry/2016/09/25/164940
# cord_x <- c("yshift")
# cord_y <- c("1-mSBD")
# ggplot_title <- c("Animal_all")
# gg <- ggplot(dfr_yshift, aes(x = yshift,y = 1-mSBD)) +
#     stat_density2d(aes(fill=..density..),
#                    geom = "raster", 
#                    contour = FALSE) +
#     scale_fill_continuous(type = "viridis", 
#                           limits=c(0, 0.002)) +
#     labs(x = cord_x,
#          y = cord_y,
#          title = ggplot_title) +
#     ylim(c(0, 1)) +
#     xlim(c(-6000, 6000)) +
#     theme(plot.title = element_text(hjust = 0.5),
#           text = element_text(size = 90),
#           legend.key.height = unit(2.0, "cm"),
#           legend.key.width = unit(1.0, "cm")
#           )

#### ggsave####
ggsave(filename = args_output, 
       plot = gg,
       dpi = 50, 
       width = 50.0, 
       height = 40.0,
       limitsize = FALSE)