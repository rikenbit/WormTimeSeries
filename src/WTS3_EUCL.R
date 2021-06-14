source("src/functions_WTS3.R")

#### args####
args <- commandArgs(trailingOnly = T)
# select animal number 個体番号の指定
args_sample <- args[1]
# outputファイル名
args_output <- args[2]
# 中間データファイル名
args_DTW <- args[3]
# # select data データの指定
args_data <- c("normalize_1")
#######################

#### load NeuronActivity####
# inputpath <- paste('data', args_data, 'ReadData_1.RData', sep = '/')
eval(parse(text=paste0("inputpath <- paste('data', args_data, 'ReadData_",args_sample,".RData', sep = '/')")))
load(inputpath)
eval(parse(text=paste0("ReadData <- ReadData_",args_sample)))

#### EUCL####
d <- diss(ReadData, "EUCL")
save(d, file=args_DTW)

# #### Rtsne####
# tSNE <- Rtsne(d, is_distance = TRUE, dims = 2, perplexity = 15, verbose = TRUE, max_iter = 1000)
# # tSNE <- Rtsne(d, is_distance = TRUE, dims = 2, perplexity = 5, verbose = TRUE, max_iter = 1000)
# df_tSNE <- data.frame(tsne_1 = tSNE$Y[,1],
#                       tsne_2 = tSNE$Y[,2],
#                       celltype = attr(d, "Labels")
#                       )

# #### merge igraph####
# load("data/igraph/Fig1_HNS.RData")
# # node information convert dataframe
# ig_Fig1_HNS %>% 
#     igraph::as_data_frame(., what="vertices") -> df_node
# # remove space from colnames
# df_node %>% 
#     names() %>% 
#         str_replace_all(., c(" " = "")) -> names(df_node)
# # add node information
# df_merged <- merge(df_tSNE, 
#                    df_node, 
#                    by.x = "celltype", 
#                    by.y = "name", 
#                    all.x = TRUE)

# #### ggplot neuron group####
# g_col <- c('NeuronType')
# gg_nt <- gg_n(g_col)

# g_col <- c('NeuronGroup')
# gg_ng <- gg_n(g_col)

# #### ggplot clustering group####
# seq(3,10) %>% 
#     purrr::map(., gg_clsters) -> gg_cls

# #### patchwork####
# append(gg_cls, list(gg_nt)) %>% 
#     append(., list(gg_ng)) -> gg_cls
# eval(parse(text=paste0("plot_title <- c('EUCL_SampleNumber_",args_sample,"')")))

# gg <- wrap_plots(gg_cls) +
#     plot_annotation(
#         title = plot_title,
#         caption = 'made with patchwork',
#         theme = theme(plot.title = element_text(size = 48, hjust = 0.5))
#     )
# ggsave(filename = args_output, plot = gg, dpi = 100, width = 40.0, height = 30.0)