#library
##################################################
library(tidyverse)
library(patchwork)
library(openxlsx)
##################################################
.plot_yshift= function(x) {
    df_tsPlot %>%
        filter(., cell_type == label_filter_list[x]) %>%
            mutate(.,
                   stim_timing = if_else(stim_timing == 1,
                                         max(.$n_activity),
                                         min(.$n_activity))
                   ) -> data_shifted
    data_shifted %>%
        dplyr::arrange(time_frame) -> data_shifted
    p_1 <- ggplot(data = data_shifted)
    return_object <- p_1 +
        geom_line(aes(x = time_frame,
                      y = n_activity,
                      colour = "n_activity")
        ) +
        geom_line(aes(x = time_frame,
                      y = yshift,
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
        eval(parse(text=paste0("ggtitle('celltype_",label_filter_list[x],"_平行移動 ",data_shifted$yshift_value[1],"')"))) +
        t_1 +
        t_2 +
        t_3 +
        sX
    return(return_object)
}

.ReadData_stimAfter = function(x,y) {
    #### load stim timing extra####
    stimtimng_sheet <- read.xlsx(y,
                                 sheet = "Sheet1",
                                 rowNames = FALSE,
                                 colNames =TRUE)
    stimtimng_sheet %>% 
        dplyr::select(sample_number = 1, 
                      stim_first = 7,
                      ) %>% 
            filter(sample_number == args_sample) %>% 
                .$stim_first %>% 
                    trunc() -> stimtimng #切り捨て
    return_object <- x[stimtimng:nrow(x),]
    return(return_object)
}
.ReadData_all = function(x) {
    return_object <- x
    return(return_object)
}