# library
##################################################
library(openxlsx)
library(tidyverse)
library(patchwork)
library(dtwclust)
# library(dtwclust)でデフォルトの変数生成器じゃなくなるので元に戻す
RNGkind(kind = "Mersenne-Twister")
set.seed(1234)
##################################################
#### check args_shift####
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

#### filter stim####
filter_stim = function(x) {
    x %>% 
        filter(., stim == 1) -> df_filter
    return(df_filter)
}

#### filter same_clusters####
filter_clusters = function(x,y) {
    x %>% 
        filter(., cell_type == y) %>% 
            .$cls %>% 
                unique() -> same_cls
    x %>% 
        filter(., cls == same_cls) -> df_filter
    return(df_filter)
}

#### dtwclust::SBD()####
sbd_y = function(x) {
    shift_2 <- input_n.list[[x]] %>% as.numeric()
    sbd <- dtwclust::SBD(shift_1,
                         shift_2, 
                         znorm = FALSE, 
                         error.check = TRUE, 
                         return.shifted = TRUE)
    return(sbd$yshift)
}
#### func df_yshift####
df_yshift = function(x) {
    df_merged_yshift %>% 
        filter(., cell_type == list_cell_type[x]) %>% 
            mutate(.,
                   stim_timing = if_else(stim_timing == 1, 
                                         max(.$n_activity), 
                                         min(.$n_activity))
        ) -> data_shifted
    data_shifted %>% 
        dplyr::arrange(time_frame) -> data_shifted
    return(data_shifted)
}
#### func y_shift_dbl####
y_shift_dbl = function(x) {
    data_shifted <- data_shifted_list[[x]]
    if (data_shifted$y_shift[1]==0) {
        # はじめて0でない要素番号の取得
        data_shifted %>% 
            filter(y_shift != 0) %>% 
            slice_head() %>% 
            .$time_frame -> first_not_0
        # y_shift_valueの計算
        y_shift_value <- first_not_0 -1
    } else {
        # はじめて0になる要素番号の取得
        data_shifted %>% 
            filter(y_shift == 0) %>% 
            slice_head() %>% 
            .$time_frame -> first_0
        # y_shift_valueの計算
        y_shift_value <- first_0 -1 -length(data_shifted$time_frame)
    }
    # yshiftが0の場合，numeric(0)になるので,0を代入
    if (length(y_shift_value) == 0) {
        y_shift_value <- 0
    }
    return(y_shift_value)
}
# #### plot one cell 3coloar y-shift value####
plot_yshift = function(x) {
    df_merged_yshift %>%
        filter(., cell_type == list_cell_type[x]) %>%
        mutate(.,
               stim_timing = if_else(stim_timing == 1,
                                     max(.$n_activity),
                                     min(.$n_activity))
        ) -> data_shifted
    # get y-shift value
    data_shifted %>%
        dplyr::arrange(time_frame) -> data_shifted
    if (data_shifted$y_shift[1]==0) {
        # はじめて0でない要素番号の取得
        data_shifted %>% 
            filter(y_shift != 0) %>% 
            slice_head() %>% 
            .$time_frame -> first_not_0
        # y_shift_valueの計算
        y_shift_value <- first_not_0 -1
    } else {
        # はじめて0になる要素番号の取得
        data_shifted %>% 
            filter(y_shift == 0) %>% 
            slice_head() %>% 
            .$time_frame -> first_0
        # y_shift_valueの計算
        y_shift_value <- first_0 -1 -length(data_shifted$time_frame)
    }
    # yshiftが0の場合，numeric(0)になるので,0を代入
    if (length(y_shift_value) == 0) {
        y_shift_value <- 0
    }

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
        eval(parse(text=paste0("ggtitle('celltype_",list_cell_type[x],"_平行移動 ",y_shift_value,"')"))) +
        t_1 +
        t_2 +
        t_3 +
        sX
    return(p_2)
}


##### サンプル番号を読み込んで，loadする関数#####
load_tempdata = function(x) {
    args_number <- args_numbers[x]
    args_tempdata <- paste(args_tempdata_dir,args_number, sep = "")
    eval(parse(text=paste0("tempdata <- c('",args_tempdata,".RData')")))
    load(tempdata)
    df_tempdata$sample_number <- rep(args_number,nrow(df_tempdata))
    df_tempdata$cell_type %>% 
        check_args_shift() -> args_shift
    filter_clusters(df_tempdata,args_shift) -> cls_n
    cls_n$cls %>% unique() -> cls_num
    df_tempdata$cls_number <- rep(cls_num,nrow(df_tempdata))
    df <- df_tempdata
    return(df)
}

#### table_cell####
table_cell = function(x) {
    df_table <- x
    # cell_typeのリストを取得(stim = 1になっている行のみ)
    df_table %>% 
        filter(stim==1) %>%
            select(cell_type) %>% 
                .$cell_type %>% 
                    unique() %>% 
                        sort() -> stim_cell_type
    # フィルター 取得したcell_type列のみの行にする，sample_number，細胞名，stimのみ残す
    df_table %>% 
        filter(cell_type %in% stim_cell_type) %>%
            select(sample_number, 
                   cell_type, 
                   stim, 
                   NeuronType) -> df_stim_table
    # create table neuron
    df_stim_table %>%
        # 数字の細胞除去
        filter(NeuronType %in% c("Sensory","Interneuron","Endorgan","Motorneuron")) %>%
            pivot_wider(id_cols = sample_number,
                        names_from = cell_type,
                        values_from = stim) -> table_neuron
    # create table rmSensory
    df_stim_table %>%
        # 数字の細胞除去
        filter(NeuronType %in% c("Interneuron","Endorgan","Motorneuron")) %>%
            pivot_wider(id_cols = sample_number,
                        names_from = cell_type,
                        values_from = stim) -> table_rmSensory
    output_table <- list(table_neuron,table_rmSensory)
    return(output_table)
}

#### table_cls####
table_cls = function(x) {
    df_table <- x
    df_table %>% 
        filter(cls==cls_number) %>%
            select(sample_number,
                   cell_type,
                   cls,
                   NeuronType) -> df_cls_table
    df_cls_table$cls <- 1
    # create table neuron
    df_cls_table %>% 
        # 数字の細胞除去
        filter(NeuronType %in% c("Sensory","Interneuron","Endorgan","Motorneuron")) %>%
            pivot_wider(id_cols = sample_number,
                        names_from = cell_type,
                        values_from = cls) -> table_neuron
    table_neuron[is.na(table_neuron)] <- 0
    # create table rmSensory
    df_cls_table %>% 
        # 数字の細胞除去
        filter(NeuronType %in% c("Interneuron","Endorgan","Motorneuron")) %>%
            pivot_wider(id_cols = sample_number,
                        names_from = cell_type,
                        values_from = cls) -> table_rmSensory
    table_rmSensory[is.na(table_rmSensory)] <- 0
    output_table <- list(table_neuron,table_rmSensory)
    return(output_table)
}

#### func load_shift_table####
load_shift_table = function(x) {
    args_number <- args_numbers[x]
    args_shift_table <- paste(args_shift_dir,args_number, sep = "")
    eval(parse(text=paste0("tempdata <- c('",args_shift_table,".RData')")))
    load(tempdata)
    y_shift_table %>% 
        mutate(sample_number = x) %>% 
        select(sample_number, cell_type, y_shift)-> y_shift_table
    return(y_shift_table)
}

# ##### First Stim TimeFrame####
# range_all = function(x) {
#     first_stim <- 1
#     return(first_stim)
# }
# range_after = function(x) {
#     data.frame(
#         TimeFrame = timeframe,
#         StimTiming = x,
#         stringsAsFactors = FALSE
#     ) -> input_stim_df
#     input_stim_df %>% 
#         filter(StimTiming != 0) %>%
#             slice_head() %>%
#                 .$TimeFrame -> first_stim
#     return(first_stim)
# }
