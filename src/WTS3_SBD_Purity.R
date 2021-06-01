source("src/functions_WTS3.R")

#### args####
args <- commandArgs(trailingOnly = T)
# select animal number 個体番号の指定
args_sample <- args[1]
# outputファイル名
args_output <- args[2]
# 中間データファイル名
args_SBD <- args[3]
# option para
args_op1 <- as.numeric(args[4])
# # select data データの指定
args_data <- c("normalize_1")
#######################
# #### test args####
# args_sample <- c("1")
# # outputファイル名
# args_output <- c("output/WTS3/SBD/normalize_1/all/SampleNumber_1/SBD_Purity_15.png")
# # 中間データファイル名
# args_SBD <- c("output/WTS3/SBD/normalize_1/all/SampleNumber_1/SBD.RData")
# # option para
# args_op1 <- as.numeric(15)
# # select data データの指定
# args_data <- c("normalize_1")
# #######################

#### SBD####
load(args_SBD)

#### Rtsne####
set.seed(1234)
tSNE <- Rtsne(d, 
              is_distance = TRUE, 
              dims = 2, 
              perplexity = args_op1, 
              verbose = TRUE, 
              max_iter = 1000)
df_tSNE <- data.frame(tsne_1 = tSNE$Y[,1],
                      tsne_2 = tSNE$Y[,2],
                      cell_type = attr(d, "Labels")
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
df_merged <- merge(df_tSNE, 
                   df_node, 
                   by.x = "cell_type", 
                   by.y = "name", 
                   all.x = TRUE)

#### import stim sheet####
# load WTS2_PeriodicACF.csv
# periodic_sheet <- read.csv("output/WTS2/WTS2_PeriodicACF.csv", 
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

#### ggplot clustering group####
cls_length <- seq(3,10)
cls_length %>% 
    purrr::map_dbl(., cls_purity) -> ClusterP_n
ClusterP_df <- data.frame(cls_length = cls_length, 
                          purity = ClusterP_n
                          )
ClusterP_df %>% 
    filter(., purity == max(purity)) %>%
        .$cls_length %>% 
            purrr::map(., gg_clusters) -> gg_cls
# table of purity
ClusterP_df %>% 
    # dplyr::summarise_all(list(round), digits=3) %>% 
        ggtexttable(rows = NULL, theme = ttheme(base_size = 50)) -> gg_cls_table

#### patchwork####
append(list(gg_nt), list(gg_ng)) %>% 
    append(., gg_cls) -> gg_cls
eval(parse(text=paste0("plot_title <- c('SBD_SampleNumber_",args_sample,"')")))
# eval(parse(text=paste0("plot_title <- c('SBD_SampleNumber_",args_sample,"_",args_op1,"')")))
gg <- wrap_plots(gg_cls) +
    plot_annotation(
        title = plot_title,
        caption = 'made with patchwork',
        theme = theme(plot.title = element_text(size = 48, hjust = 0.5))
    )
gg <- gg + gg_cls_table
ggsave(filename = args_output, plot = gg, dpi = 100, width = 40.0, height = 30.0)