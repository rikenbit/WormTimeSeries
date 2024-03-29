source("src/functions_WTS4_DimReduc_MCMI.R")

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


#### No. of Clusters####
k <- as.numeric(args_k)

### MCMI weight table####
# merged_data
load(args_input_MCMIHOOI)
data.frame(weight = merged_data$W,
           stringsAsFactors = FALSE) %>%
    rownames_to_column("SampleNumber") %>%
        mutate(weight_abs =abs(weight)) %>%
            dplyr::arrange(desc(weight_abs)) -> df_weight

#### load dist data####
# 空の行列を格納するファイルを作成
D <- list()
# inputファイル名のリスト
input_path_list <- list.files(args_input_path, pattern="SampleNumber_", full.names=TRUE)
# ファイルを読み込んで，リストに加える．各リストのattr(*, "Labels")に細胞型名が残っている
for(i in 1:length(input_path_list)){
    load(input_path_list[i])
    D <- c(D, list(d))
}

#### Clustering against each distance matrix####
C <- lapply(D, function(d, k) {
    cutree(hclust(d, method="ward.D2"), k)
    }, k=k)

#### Filter dist matrix####
# # purrr
# seq(1:length(D)) %>%
#     purrr::map(., .dist_f) -> D_f
D -> D_f

#### Dimensionality Reduction####
# purrr
seq(1:length(D_f)) %>%
    purrr::map(., .df_cords) -> DF_cord

#### merge cord and cls####
# purrr
seq(1:length(DF_cord)) %>%
    purrr::map(., .df_cord_cls) -> DF_cord_cls
# save sample_cls
sample_cls <- lapply(DF_cord_cls, function(x){dplyr::select(x, cell_type, Cluster)})
args_output %>%
    str_remove(., "/table.png") %>%
        str_remove(., args_DimReduc) -> args_output_sample_cls
eval(parse(text=paste0("args_output_sample_cls <- c('",args_output_sample_cls,"sample_cls.RData')")))
save(sample_cls, file=args_output_sample_cls)
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

#### gg_weight####
# ggtexttable
df_weight %>%
    rownames_to_column("Ranking") %>%
        mutate_if(is.numeric, round, digits = 3) %>%
            ggtexttable(rows = NULL, theme = ttheme(base_size = 48)) -> gg_weight

#### gg_list purrr####
seq(1:length(DF_cord_cls_nl_be)) %>%
    purrr::map(., .plot_dimreduc) -> gg_list

#### ggplot for####
for(x in 1:length(gg_list)){
# for(x in 1:2 ){
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
    sample_n <- df_weight[x,1]
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