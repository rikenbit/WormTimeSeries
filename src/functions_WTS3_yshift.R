#library
##################################################
library(openxlsx)
library(tidyverse)
##################################################
.df_yshift = function(x) {
    sbd_yshift_df %>% 
        filter(., cell_type == list_cell_type[x]) -> return_object
    return(return_object)
}

.yshift_value = function(x) {
    data_shifted <- yshift_list[[x]]
    if (data_shifted$yshift[1]==0) {
      # はじめて0でない要素番号の取得
        data_shifted %>% 
            filter(yshift != 0) %>% 
                slice_head() %>% 
                    .$time_frame -> first_not_0
        # yshift_valueの計算
        return_object <- first_not_0 -1
    } else {
      # はじめて0になる要素番号の取得
        data_shifted %>% 
            filter(yshift == 0) %>% 
                slice_head() %>% 
                    .$time_frame -> first_0
        # yshift_valueの計算
        return_object <- first_0 -1 -length(data_shifted$time_frame)
    }
    # yshiftが0の場合，numeric(0)になるので,0を代入
    if (length(return_object) == 0) {
        return_object <- 0
    }
    return(return_object)
}

.yshift_filter = function(x) {
    if (yshift_abs[x] <= yshift_max) {
        return_object <- 1
    } else {
        return_object <- 0
    }
    return(return_object)
}