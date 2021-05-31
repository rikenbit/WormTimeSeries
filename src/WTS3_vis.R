source("src/functions_WTS3_new.R")

#### args####
args <- commandArgs(trailingOnly = T)
# select animal number 個体番号の指定
args_sample <- args[1]
# outputファイル名
args_output <- args[2]
# 中間データファイル名
args_dist <- args[3]
# select data データの指定
args_data <- c("normalize_1")
# クラスター評価手法
args_eval <- args[4]
# 次元圧縮手法
args_DimRedu <- args[5]
# #######################
# #### test args####
# args_sample <- c("1")
# # outputファイル名
# # args_output <- c("output/WTS3/SBD/normalize_1/all/umap/purity/SampleNumber_1.png")
# args_output <- c("output/WTS3/SBD/normalize_1/all/umap/ARI/SampleNumber_1.png")
# # args_output <- c("output/WTS3/SBD/normalize_1/all/tsne/Fmeasure/SampleNumber_1.png")
# # args_output <- c("output/WTS3/SBD/normalize_1/all/tsne/Entropy/SampleNumber_1.png")

# # 中間データファイル名
# args_dist <- c("output/WTS3/SBD/normalize_1/all/SampleNumber_1/SBD.RData")

# # select data データの指定
# args_data <- c("normalize_1")

# # クラスター評価手法
# # args_eval <- c("purity")
# args_eval <- c("ARI")
# # args_eval <- c("Fmeasure")
# # args_eval <- c("Entropy")

# # 次元圧縮手法
# args_DimRedu <- c("umap")
# #######################

### SBD####
load(args_dist)

#### Dimensionality Reduction####
df_cord <- switch(args_DimRedu,
          "tsne" = wts_tsne(d),
          "umap" = wts_umap(d),
          stop("Only can use tsne,")
          )

#### merge igraph####
load("data/igraph/Fig1_HNS.RData")
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
periodic_sheet <- read.csv("output/WTS2/WTS2_PeriodicACF.csv", 
                           colClasses=c("numeric", 
                                        "character", 
                                        rep("numeric",2)))
periodic_sheet %>% 
    filter(.,
           stim == 1,
           sample_number == args_sample) -> stim_sheet

#### ggplot neuron group####
g_col <- c('NeuronType')
gg_nt <- gg_n(g_col)
g_col <- c('NeuronGroup')
gg_ng <- gg_n(g_col)

#### clustering evaluation####
cls_length <- seq(3,10)
cls_length %>% 
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
cls_length %>% 
    purrr::map_dbl(., eval_type) -> ClusterP_n
ClusterP_df <- data.frame(cls_length = cls_length, 
                          cls_eval = ClusterP_n
                          )

#### ggplot clustering group####
gg_cls <- switch(args_eval,
                 "purity" = max_eval(ClusterP_df),
                 "ARI" = max_eval(ClusterP_df),
                 "Fmeasure" = max_eval(ClusterP_df),
                 "Entropy" = min_eval(ClusterP_df),
                 stop("Only can use cls_purity,ARI,Fmeasure,Entropy")
)
            
#### table of cls_eval####
ClusterP_df %>% 
    dplyr::summarise_all(list(round), digits=5) %>% 
        # ggpubr
        ggtexttable(rows = NULL, theme = ttheme(base_size = 50)) -> gg_cls_table

#### patchwork####
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