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

# #### ggplot density group####
# cord_x <- c("yshift")
# cord_y <- c("density")
# ggplot_title <- c("Animal_all")
# gg <- ggplot(dfr_yshift,
#              aes(x = yshift,
#                  fill=factor(animal, levels = names(DFs)))
#              ) +
#     geom_density(alpha=0.2) +
#     labs(x = cord_x,
#          y = cord_y,
#          title = ggplot_title,
#          fill = "Animal_No.") +
#     ylim(c(0, 0.05)) +
#     xlim(c(-1000, 1000)) +
#     theme(plot.title = element_text(hjust = 0.5)) +
#     theme(text = element_text(size = 90))
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

#### ggplot density group####
cord_x <- c("1-mSBD")
cord_y <- c("density")
ggplot_title <- c("Animal_all")

# # case_when https://www.jaysong.net/RBook/visualization4.html
# dfr_yshift %>% mutate(yshift_abs = case_when(abs(yshift) < 100 ~ "100",
#                                              abs(yshift) < 200 ~ "101~200",
#                                              abs(yshift) < 300 ~ "201~300",
#                                              TRUE              ~ "300以上"),
#                       yshift_abs = factor(yshift_abs, 
#                                           levels = c("0~100",
#                                                      "101~200",
#                                                      "201~300",
#                                                      "300以上"))
#                       ) -> dfr_yshift_group

# case_when https://www.jaysong.net/RBook/visualization4.html
dfr_yshift %>% mutate(yshift_abs = case_when(abs(yshift) <= 100 ~ "0~100",
                                             abs(yshift) <= 200 ~ "101~200",
                                             abs(yshift) <= 300 ~ "201~300",
                                             abs(yshift) <= 400 ~ "301~400",
                                             abs(yshift) <= 500 ~ "401~500",
                                             abs(yshift) <= 600 ~ "501~600",
                                             abs(yshift) <= 700 ~ "601~700",
                                             abs(yshift) <= 800 ~ "701~800",
                                             abs(yshift) <= 900 ~ "801~900",
                                             abs(yshift) <= 1000 ~ "901~1000",
                                             TRUE              ~ "1001~"),
                      # yshift_abs = factor(yshift_abs,
                      #                     levels = c("0~100",
                      #                                "101~200",
                      #                                "201~300",
                      #                                "301~400",
                      #                                "401~500",
                      #                                "501~600",
                      #                                "601~700",
                      #                                "701~800",
                      #                                "801~900",
                      #                                "901~1000",
                      #                                "1001以上"))
                      yshift_abs = factor(yshift_abs,
                                          levels = c("1001~",
                                                     "901~1000",
                                                     "801~900",
                                                     "701~800",
                                                     "601~700",
                                                     "501~600",
                                                     "401~500",
                                                     "301~400",
                                                     "201~300",
                                                     "101~200",
                                                     "0~100"))
                      ) -> dfr_yshift_group

gg <- ggplot(dfr_yshift_group,
             # dfr_yshift_group[dfr_yshift_group$yshift_abs!="1000以上",],
             aes(x = 1-mSBD,
                 fill=yshift_abs)
             ) +
    geom_density(alpha=0.8) +
    labs(x = cord_x,
         y = cord_y,
         title = ggplot_title,
         fill = "yshift_abs") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(text = element_text(size = 90))

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