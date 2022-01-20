source("src/functions_WTS4_DimReduc_sample.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
# input merged_distance
args_input_MCMIHOOI <- args[1]
# DimReduc
args_output <- args[2]
args_input_path <- args[3]
args_DimReduc <- args[4]
# No. of Clusters
args_k <- args[5]
# add anotation data
args_NL <- args[6]
args_eval_label <- args[7]

args_input_cls  <- args[8]

# #### test args####
# # input merged_distance
# args_input_MCMIHOOI <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/Merged_data/k_Number_5.RData")
# # DimReduc
# args_output <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/DimReduc_sample/k_Number_5/tsne/table.png")
# args_input_path <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Distance")
# args_DimReduc <- c("tsne")
# args_k <- c("5")
# # add anotation data
# args_NL <- c("data/igraph/Fig1_HNS.RData")
# args_eval_label <- c("data/WTS4_Eval_behavior_fix.xlsx")

# args_input_cls <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Cluster_sample/k_Number_5/sample_cls.RData")

#### No. of Clusters####
k <- as.numeric(args_k)

#### load dist data####
# inputファイル名のリスト
input_path_list <- list.files(args_input_path, pattern="SampleNumber_", full.names=TRUE)

#### fix sample number sort####
input_path_list %>% 
    str_remove(., args_input_path) %>% 
    str_remove(., "/SampleNumber_") %>% 
    str_remove(., ".RData") %>% 
    as.numeric() %>% 
    sort() -> sample_sort_num
input_path <- c()
for(i in sample_sort_num){
    eval(parse(text=paste0("path <- c('",args_input_path,"/SampleNumber_",i,".RData')")))
    input_path <- c(input_path, path)
}

input_path_list <- input_path 
########

# 空の行列を格納するファイルを作成
D <- list()
# ファイルを読み込んで，リストに加える．各リストのattr(*, "Labels")に細胞型名が残っている
for(i in 1:length(input_path_list)){
    load(input_path_list[i])
    D <- c(D, list(d))
}

#### Clustering against each distance matrix####
load(args_input_cls)

#### Dimensionality Reduction####
seq(1:length(D)) %>%
    purrr::map(., .df_cords) -> DF_cord

#### merge cord and cls####
seq(1:length(DF_cord)) %>%
    purrr::map(., .df_cord_cls) -> DF_cord_cls

#### merge igraph NeuronType####
load(args_NL)
# node information convert dataframe
ig_Fig1_HNS %>% 
    igraph::as_data_frame(., what="vertices") -> df_node
# remove space from colnames
df_node %>% 
    names() %>% 
        str_replace_all(., c(" " = "")) -> names(df_node)
df_nl <- data.frame(cell_type = df_node$name,
                    NeuronType = df_node$NeuronType,
                    NeuronGroup = df_node$NeuronGroup,
                    stringsAsFactors = FALSE
                    )
# purrr
seq(1:length(DF_cord_cls)) %>%
    purrr::map(., .df_cord_cls_nl) -> DF_cord_cls_nl

#### merge xlsx behavoir label####
# load eval_label
read.xlsx(args_eval_label,
          sheet = "Sheet1",
          rowNames = FALSE,
          colNames =TRUE) %>% 
    dplyr::rename(cell_type = celltype, 
                  Classes = class) -> df_eval_label
seq(1:length(DF_cord_cls_nl)) %>%
    purrr::map(., .df_cord_cls_nl_be) -> DF_cord_cls_nl_be

### MCMI weight table####
# merged_data
load(args_input_MCMIHOOI)
data.frame(weight = merged_data$W,
		   SampleNumber = sample_sort_num,
           stringsAsFactors = FALSE) %>% 
    rownames_to_column("ID") %>%
        mutate(weight_abs =abs(weight)) %>% 
            dplyr::arrange(desc(weight_abs)) %>% 
            	dplyr::select(1,3,2,4) -> df_weight

#### gg_weight####
# ggtexttable
df_weight %>% 
dplyr::select(2,3,4) %>% 
    rownames_to_column("Ranking") %>% 
        mutate_if(is.numeric, round, digits = 3) %>% 
            ggtexttable(rows = NULL, theme = ttheme(base_size = 48)) -> gg_weight

#### gg_list purrr####
seq(1:length(DF_cord_cls_nl_be)) %>%
    purrr::map(., .plot_dimreduc) -> gg_list

#### ggplot for####
for(x in 1:length(gg_list)){
    # weight table
    gg_weight_bg <- table_cell_bg(gg_weight, 
                               row = x + 1,
                               column = 1:4, 
                               linewidth = 5,
                               fill="darkolivegreen1", 
                               color = "darkolivegreen4")

    # patchwork add gg_weight
    gg <- gg_list[[x]] +
        gg_weight_bg +
        plot_layout(nrow = 1)

    # filename
    args_output %>% 
        str_remove(., "/table.png") -> args_output_name
    sample_n <- df_weight[x, 2]
    eval(parse(text=paste0("plot_title <- c('",x,"_SampleNumber_",sample_n,".png')")))
    eval(parse(text=paste0("args_output_sample <- c('",args_output_name,"/",plot_title,"')")))
    
    # ggsave
    ggsave(filename = args_output_sample, 
           plot = gg,
           dpi = 100, 
           width = 80.0, 
           height = 20.0,
           limitsize = FALSE)
}
#### ggplot output file####
ggsave(filename = args_output, 
       plot = gg_weight,
       dpi = 100, 
       width = 80.0, 
       height = 20.0,
       limitsize = FALSE)