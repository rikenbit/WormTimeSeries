# library
##################################################
library(tidyverse)
library(patchwork)
library(ggrepel)
##################################################
filter_stim = function(x) {
    x %>% 
    filter(., stim == 1) %>% 
        mutate(.,
           stim_timing = if_else(stim_timing == 1, 
                                 max(.$n_activity), 
                                 min(.$n_activity))) -> df_filter
    return(df_filter)
}