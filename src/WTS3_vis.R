source("src/functions_WTS3_new.R")

#### args####
args <- commandArgs(trailingOnly = T)
# select animal number 個体番号の指定
args_sample <- args[1]
# output
args_output <- args[2]
# output 中間データ
args_tempdata <- args[3]
# input 距離データ
args_dist <- args[4]
# input igraph
args_igraph <- c("data/igraph/Fig1_HNS.RData")
# input PeriodicACF
args_periodic <- c("output/WTS2/WTS2_PeriodicACF.csv")
# select data データの指定
args_data <- c("normalize_1")
# クラスター評価手法
args_eval <- args[5]
# 次元圧縮手法
args_DimRedu <- args[6]
# time range
args_range <- args[7]
# #### test args####
# # select animal number 個体番号の指定
# args_sample <- c("1")
# # output
# args_output <- c("output/WTS3/normalize_1/all/SBD/ARI/tsne/cls_plot/SampleNumber_1.png")
# # output 中間データ
# args_tempdata <- c("output/WTS3/normalize_1/all/SBD/ARI/tsne/cls_tempdata/SampleNumber_1.RData")
# # input 距離データ
# args_dist <- c("output/WTS3/normalize_1/all/SBD/SampleNumber_1/SBD.RData")
# # input igraph
# args_igraph <- c("data/igraph/Fig1_HNS.RData")
# # input PeriodicACF
# args_periodic <- c("output/WTS2/WTS2_PeriodicACF.csv")
# # select data データの指定
# args_data <- c("normalize_1")
# # クラスター評価手法
# args_eval <- c("ARI")
# # 次元圧縮手法
# args_DimRedu <- c("tsne")
# # time range
# args_range <- c("all")

### load SBD####
load(args_dist)

#### Dimensionality Reduction####
df_cord <- switch(args_DimRedu,
          "tsne" = wts_tsne(d),
          "umap" = wts_umap(d),
          stop("Only can use tsne,")
          )

#### merge igraph####
load(args_igraph)
# node information convert dataframe
ig_Fig1_HNS %>% 
    igraph::as_data_frame(., what="vertices") -> df_node
# remove space from colnames
df_node %>% 
    names() %>% 
        str_replace_all(., c(" " = "")) -> names(df_node)
# add node information
df_merged <- merge(df_cord, 
                   df_node, 
                   by.x = "cell_type", 
                   by.y = "name", 
                   all.x = TRUE)

#### import stim sheet####
periodic_sheet <- read.csv(args_periodic, 
                           colClasses=c("numeric", 
                                        "character", 
                                        rep("numeric",2)))
periodic_sheet %>% 
    filter(.,
           stim == 1,
           sample_number == args_sample) -> stim_sheet
#### Stim Cell Checker####
df_merged_stim <- merge(df_merged, 
                        stim_sheet, 
                        by.x = "cell_type", 
                        by.y = "cell_type", 
                        all.x = TRUE)
df_merged_stim %>% 
    replace_na(., replace = list(stim = 0)) -> df_merged_stim0
#### GGplot neuron group####
g_col <- c('NeuronType')
gg_nt <- gg_n_stim(g_col)
g_col <- c('NeuronGroup')
gg_ng <- gg_n_stim(g_col)
#### clustering evaluation####
set_cutree <- seq(3,10)
set_cutree %>% 
    purrr::map(., cls_cord) -> df_cls_cord
# select method
eval_type <- switch(args_eval,
              "purity" = cls_purity,
              "ARI" = cls_ARI,
              "Fmeasure" = cls_Fmeasure,
              "Entropy" = cls_Entropy,
              stop("Only can use cls_purity,ARI,Fmeasure,Entropy")
)
# evaluation
set_cutree %>% 
    purrr::map_dbl(., eval_type) -> eval_value
ClusterP_df <- data.frame(set_cutree = set_cutree, 
                          eval_value = eval_value
                          )

#### GGplot clustering group####
gg_cls <- switch(args_eval,
                 "purity" = max_eval(ClusterP_df),
                 "ARI" = max_eval(ClusterP_df),
                 "Fmeasure" = max_eval(ClusterP_df),
                 "Entropy" = min_eval(ClusterP_df),
                 stop("Only can use cls_purity,ARI,Fmeasure,Entropy")
)
            
#### GGpubr table of eval_value####
ClusterP_df %>% 
    dplyr::summarise_all(list(round), digits=5) %>% 
        # ggpubr
        ggtexttable(rows = NULL, theme = ttheme(base_size = 50)) -> gg_cls_table

#### output temp data####
# merge cell group vs clustering group
df_cell_cls <- merge(df_merged, 
                      gg_cls[[1]]$data, 
                      by.x = "cell_type", 
                      by.y = "cell_type", 
                      all.x = TRUE)
# merge stim sheet
df_cell_cls_stim <- merge(df_cell_cls, 
                     stim_sheet, 
                     by.x = "cell_type", 
                     by.y = "cell_type", 
                     all.x = TRUE)
# stim列のNAを0に変換する
df_cell_cls_stim %>% 
    replace_na(., replace = list(stim = 0)) -> df_cell_cls_stim0
# remove column
df_cell_cls_stim0 %>% 
    dplyr::select(., 
                  cell_type, 
                  NeuronType, 
                  cls, 
                  stim, 
                  NeuronGroup, 
                  Color, 
                  acf) -> df_tempdata
# save tempdata
save(df_tempdata, file=args_tempdata)
#### patchwork & save####
# append(list(gg_nt), list(gg_ng)) -> gg_cls
append(list(gg_nt), list(gg_ng)) %>% 
    append(., gg_cls) -> gg_cls
eval(parse(text=paste0("plot_title <- c('",args_eval,"_SampleNumber_",args_sample,"')")))
gg <- wrap_plots(gg_cls) +
    plot_annotation(
        title = plot_title,
        caption = 'made with patchwork',
        theme = theme(plot.title = element_text(size = 48, hjust = 0.5))
    )
gg <- gg + gg_cls_table
ggsave(filename = args_output, 
       plot = gg, 
       dpi = 100, 
       width = 40.0, 
       height = 30.0)