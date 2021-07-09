#library
##################################################
library(tidyverse)
##################################################
.load_label_table = function(x) {
    args_number <- args_numbers[x]
    label_dir <- paste(args_label_dir,args_number, sep = "")
    eval(parse(text=paste0("label_path <- c('",label_dir,"/label_table.RData')")))
    load(label_path)
    df_label$sample_number <- rep(args_number,nrow(df_label))
    return_object <- df_label
    return(return_object)
}
.load_yshift_value= function(x) {
    args_number <- args_numbers[x]
    label_dir <- paste(args_yshift_value_dir,args_number, sep = "")
    eval(parse(text=paste0("label_path <- c('",label_dir,"/yshift_value.RData')")))
    load(label_path)
    yshift_value_table$sample_number <- rep(args_number,nrow(yshift_value_table))
    return_object <- yshift_value_table
    return(return_object)
}

.table_acf = function(x,y) {
    df_table <- x
    df_yshift <- y
    merge(df_table, 
          df_yshift, 
          by.x = c("sample_number", "cell_type"), 
          by.y = c("sample_number", "cell_type"), 
          all.x = TRUE) -> df_input
    # create table neuron
    df_input %>%
        filter(label_acf==1, yshift_filter==1) %>%
        # 数字の細胞除去
        filter(NeuronType %in% c("Sensory","Interneuron","Endorgan","Motorneuron")) %>%
        pivot_wider(id_cols = sample_number,
                    names_from = cell_type,
                    values_from = label_acf) -> table_neuron
        table_neuron$sample_number <- as.numeric(table_neuron$sample_number)
        table_neuron %>% dplyr::arrange(sample_number)-> table_neuron
    # create table rmSensory
    df_input %>%
        filter(label_acf==1, yshift_filter==1) %>%
        # 数字の細胞除去
        filter(NeuronType %in% c("Interneuron","Endorgan","Motorneuron")) %>%
        pivot_wider(id_cols = sample_number,
                    names_from = cell_type,
                    values_from = label_acf) -> table_rmSensory
    as.numeric(table_rmSensory$sample_number) -> table_rmSensory$sample_number
    table_rmSensory %>% 
        dplyr::arrange(sample_number) -> table_rmSensory
    
    return_object <- list(table_neuron,table_rmSensory)
    return(return_object)
}
.table_cls = function(x,y) {
    df_table <- x
    df_yshift <- y
    merge(df_table, 
          df_yshift, 
          by.x = c("sample_number", "cell_type"), 
          by.y = c("sample_number", "cell_type"), 
          all.x = TRUE) -> df_input
    # create table neuron
    df_input %>%
        filter(label_cls==1, yshift_filter==1) %>%
        # 数字の細胞除去
        filter(NeuronType %in% c("Sensory","Interneuron","Endorgan","Motorneuron")) %>%
        pivot_wider(id_cols = sample_number,
                    names_from = cell_type,
                    values_from = label_cls) -> table_neuron
    table_neuron$sample_number <- as.numeric(table_neuron$sample_number)
    table_neuron %>% dplyr::arrange(sample_number)-> table_neuron
    # create table rmSensory
    df_input %>%
        filter(label_cls==1, yshift_filter==1) %>%
        # 数字の細胞除去
        filter(NeuronType %in% c("Interneuron","Endorgan","Motorneuron")) %>%
        pivot_wider(id_cols = sample_number,
                    names_from = cell_type,
                    values_from = label_cls) -> table_rmSensory
    as.numeric(table_rmSensory$sample_number) -> table_rmSensory$sample_number
    table_rmSensory %>% 
        dplyr::arrange(sample_number) -> table_rmSensory
    
    return_object <- list(table_neuron,table_rmSensory)
    return(return_object)
}

.table_acf_yshift = function(x,y) {
    df_table <- x
    df_yshift <- y
    merge(df_table, 
          df_yshift, 
          by.x = c("sample_number", "cell_type"), 
          by.y = c("sample_number", "cell_type"), 
          all.x = TRUE) -> df_input
    # create table neuron
    df_input %>%
        filter(label_acf==1, yshift_filter==1) %>%
        # 数字の細胞除去
        filter(NeuronType %in% c("Sensory","Interneuron","Endorgan","Motorneuron")) %>%
        pivot_wider(id_cols = cell_type,
                    names_from = sample_number,
                    values_from = yshift_value
                    ) -> table_neuron
    return_object <- table_neuron
    return(return_object)
}
.table_cls_yshift = function(x,y) {
    df_table <- x
    df_yshift <- y
    merge(df_table, 
          df_yshift, 
          by.x = c("sample_number", "cell_type"), 
          by.y = c("sample_number", "cell_type"), 
          all.x = TRUE) -> df_input
    # create table neuron
    df_input %>%
        filter(label_cls==1, yshift_filter==1) %>%
        # 数字の細胞除去
        filter(NeuronType %in% c("Sensory","Interneuron","Endorgan","Motorneuron")) %>%
        pivot_wider(id_cols = cell_type,
                    names_from = sample_number,
                    values_from = yshift_value
        ) -> table_neuron
    return_object <- table_neuron
    return(return_object)
}

.yshift_norm = function(x) {
    output_yshift_table[x,] %>% 
        .[is.na(.)==F] %>% 
        as.matrix() %>% 
        # norm(.,"F") -> norm_value #フロベニウスノルム
        norm(.,"I") -> norm_value 
    return(norm_value)
}
.add_sum_sort = function(x) {
    x %>% 
        column_to_rownames(var = "sample_number") %>% 
            t() %>% 
                as.data.frame() -> output_neuron
    output_neuron %>% 
        rowSums(.,na.rm = TRUE,dims = 1) -> sort_vec
    output_neuron %>% 
        mutate(sum_sample = sort_vec) %>% 
            dplyr::arrange(desc(sum_sample)) %>% 
                rownames_to_column(var ="cell_type") -> return_object
    return(return_object)
}