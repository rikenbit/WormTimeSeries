# library
##################################################
library(tidyverse)
library(patchwork)
library(ggrepel)
library(dtwclust)
RNGkind(kind = "Mersenne-Twister")
set.seed(1234)
##################################################
#filter stim
filter_stim = function(x) {
    x %>% 
        filter(., stim == 1) -> df_filter
    return(df_filter)
}
# plot one cell
plot_one_cell = function(x) {
    df_merged %>% 
        filter(., cell_type == list_cell_type[x]) %>% 
            mutate(.,
                   stim_timing = if_else(stim_timing == 1, 
                                         max(.$n_activity), 
                                         min(.$n_activity))) %>% 
                ggplot(., aes(x = time_frame)) -> p_1
  
    p_2 <- p_1 + 
        geom_line(aes(y = n_activity, colour = "n_activity")) +
        geom_line(aes(y = stim_timing, colour = "stim_timing"),
                  linetype = "dashed", 
                  alpha = 0.5) +
        scale_colour_manual(values = c("black", "purple")) +
        # ggtitle(list_cell_type[x]) +
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