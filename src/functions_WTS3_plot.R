# library
##################################################
library(tidyverse)
library(patchwork)
library(dtwclust)
# library(dtwclust)でデフォルトの変数生成器じゃなくなるので元に戻す
RNGkind(kind = "Mersenne-Twister")
set.seed(1234)
##################################################
#filter stim
filter_stim = function(x) {
    x %>% 
        filter(., stim == 1) -> df_filter
    return(df_filter)
}

# plot one cell 4coloar
plot_yshift = function(x) {
    df_merged_yshift %>% 
        filter(., cell_type == list_cell_type[x]) %>% 
        mutate(.,
               stim_timing = if_else(stim_timing == 1, 
                                     max(.$n_activity), 
                                     min(.$n_activity))
        ) -> data_shifted
    df_merged_yshift %>% 
        filter(., cell_type == list_cell_type[match(args_shift,list_cell_type)]) %>% 
        mutate(.,
               stim_timing = if_else(stim_timing == 1, 
                                     max(.$n_activity), 
                                     min(.$n_activity))
        ) -> data_shifted_ASER
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
        geom_line(data = data_shifted_ASER, 
                  aes(x = time_frame, 
                      y = n_activity, 
                      colour = args_shift),
                  linetype = "dashed"
        ) +
        geom_line(aes(x = time_frame, 
                      y = stim_timing, 
                      colour = "stim_timing"),
                  linetype = "dotted", 
                  alpha = 0.5
        ) +
        scale_colour_manual(values = c("black", "red", "green", "purple"),
                            breaks = c("n_activity", "n_yshift", args_shift,"stim_timing")) +
        eval(parse(text=paste0("ggtitle('celltype_",list_cell_type[x],"')"))) +
        t_1 +
        t_2 +
        t_3 +
        sX
    return(p_2)
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