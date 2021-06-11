# library
##################################################
library(tidyverse)
library(patchwork)
library(dtwclust)
# library(dtwclust)でデフォルトの変数生成器じゃなくなるので元に戻す
RNGkind(kind = "Mersenne-Twister")
set.seed(1234)
##################################################
# check args_shift
check_args_shift = function(x) {
    x -> sample_cell_type
    sample_cell_type %>% 
        str_count(., pattern="ASER") %>%
        sum() -> check_ASER
    sample_cell_type %>% 
        str_count(., pattern="BAGR") %>%
        sum() -> check_BAGR
    sample_cell_type %>% 
        str_count(., pattern="BAGL") %>%
        sum() -> check_BAGL
    if (check_ASER >= 1) {
        args_shift <- "ASER"
    } else if (check_BAGR >= 1) {
        args_shift <- "BAGR"
    } else if (check_BAGL >= 1) {
        args_shift <- "BAGL"
    } else {
        args_shift <- sample_cell_type[1]
    }
    return(args_shift)
}

# filter stim
filter_stim = function(x) {
    x %>% 
        filter(., stim == 1) -> df_filter
    return(df_filter)
}

# filter same_clusters
filter_clusters = function(x,y) {
    x %>% 
        filter(., cell_type == y) %>% 
            .$cls %>% 
                unique() -> same_cls
    x %>% 
        filter(., cls == same_cls) -> df_filter
    return(df_filter)
}

# dtwclust::SBD()
sbd_y = function(x) {
    shift_2 <- input_n.list[[x]] %>% as.numeric()
    sbd <- dtwclust::SBD(shift_1,
                         shift_2, 
                         znorm = FALSE, 
                         error.check = TRUE, 
                         return.shifted = TRUE)
    return(sbd$yshift)
}

# plot one cell 3coloar
plot_yshift = function(x) {
    df_merged_yshift %>% 
        filter(., cell_type == list_cell_type[x]) %>% 
        mutate(.,
               stim_timing = if_else(stim_timing == 1, 
                                     max(.$n_activity), 
                                     min(.$n_activity))
        ) -> data_shifted
    p_1 <- ggplot(data = data_shifted)
    p_2 <- p_1 + 
        geom_line(aes(x = time_frame, 
                      y = n_activity, 
                      colour = "n_activity")
        ) +
        geom_line(aes(x = time_frame, 
                      y = y_shift, 
                      colour = "n_yshift")
        ) +
        geom_line(aes(x = time_frame, 
                      y = stim_timing, 
                      colour = "stim_timing"),
                  linetype = "dotted", 
                  alpha = 0.5
        ) +
        scale_colour_manual(values = c("black", "red", "purple"),
                            breaks = c("n_activity", "n_yshift", "stim_timing")) +
        eval(parse(text=paste0("ggtitle('celltype_",list_cell_type[x],"')"))) +
        t_1 +
        t_2 +
        t_3 +
        sX
    return(p_2)
}

# # plot one cell 4coloar
# plot_yshift = function(x) {
#     df_merged_yshift %>% 
#         filter(., cell_type == list_cell_type[x]) %>% 
#         mutate(.,
#                stim_timing = if_else(stim_timing == 1, 
#                                      max(.$n_activity), 
#                                      min(.$n_activity))
#         ) -> data_shifted
#     df_merged_yshift %>% 
#         filter(., cell_type == list_cell_type[match(args_shift,list_cell_type)]) %>% 
#         mutate(.,
#                stim_timing = if_else(stim_timing == 1, 
#                                      max(.$n_activity), 
#                                      min(.$n_activity))
#         ) -> data_shifted_ASER
#     p_1 <- ggplot(data = data_shifted)
#     p_2 <- p_1 + 
#         geom_line(aes(x = time_frame, 
#                       y = n_activity, 
#                       colour = "n_activity")
#         ) +
#         geom_line(aes(x = time_frame, 
#                       y = y_shift, 
#                       colour = "n_yshift")
#         ) +
#         geom_line(data = data_shifted_ASER, 
#                   aes(x = time_frame, 
#                       y = n_activity, 
#                       colour = args_shift),
#                   linetype = "dashed"
#         ) +
#         geom_line(aes(x = time_frame, 
#                       y = stim_timing, 
#                       colour = "stim_timing"),
#                   linetype = "dotted", 
#                   alpha = 0.5
#         ) +
#         scale_colour_manual(values = c("black", "red", "green", "purple"),
#                             breaks = c("n_activity", "n_yshift", args_shift,"stim_timing")) +
#         eval(parse(text=paste0("ggtitle('celltype_",list_cell_type[x],"')"))) +
#         t_1 +
#         t_2 +
#         t_3 +
#         sX
#     return(p_2)
# }