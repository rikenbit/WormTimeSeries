source("src/functions_WTS4_DimReduc.R")

# #### args setting####
# args <- commandArgs(trailingOnly = T)
# # input merged_distance
# args_input_distance <- args[1]
# # input merged_cls
# args_input_cls <- args[2]
# # dimentionaly reduction
# args_DimReduc <- args[3]
# # output plot
# args_output <- args[4]

#### test args####
args_input_distance <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/10_Clusters/CSPA/merged_distance.RData")
args_input_cls <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/10_Clusters/CSPA/merged_cls.RData")
args_DimReduc <- c("tsne")
args_output <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/10_Clusters/CSPA/tsne_plot.png")
########

args_igraph <- c("data/igraph/Fig1_HNS.RData")

### load dist object####
load(args_input_distance)
d <- merged_distance

#### Dimensionality Reduction####
df_cord <- switch(args_DimReduc,
          "tsne" = .wts_tsne(d),
          "umap" = .wts_umap(d),
          stop("Only can use tsne,")
          )

### load cls number object####
load(args_input_cls)
df_cls <- as.data.frame(merged_cls)
df_cls <- data.frame(
  cell_type = rownames(df_cls),
  cls = df_cls$merged_cls,
  stringsAsFactors = FALSE
)

### merge cord and cls####
df_cord_cls <- merge(df_cord, 
                     df_cls, 
                     by.x = "cell_type", 
                     by.y = "cell_type", 
                     all.x = TRUE)
# 連続値→離散値
# df_cord_cls$cls <- as.character(df_cord_cls$cls)

# #### load Neuron Label####
# load(args_igraph)
# # node information convert dataframe
# ig_Fig1_HNS %>% 
#     igraph::as_data_frame(., what="vertices") -> df_node
# # remove space from colnames
# df_node %>% 
#     names() %>% 
#         str_replace_all(., c(" " = "")) -> names(df_node)
# df_NeuronType <- data.frame(cell_type = df_node$name,
#                             NeuronType = df_node$NeuronType,
#                             NeuronGroup = df_node$NeuronGroup,
#                             stringsAsFactors = FALSE
#                             )

# ### merge Neuron Label####
# df_cord_cls_label <- merge(df_cord_cls, 
#                            df_NeuronType, 
#                            by.x = "cell_type", 
#                            by.y = "cell_type", 
#                            all.x = TRUE)
# # unknown label is NA

#### ggplot cls####
gg <- ggplot(df_cord_cls, 
       aes(x = cord_1,
           y = cord_2, 
           label = cell_type,
           color = factor(cls)
           )
       ) + 
      geom_point(size = 6.0, 
                 alpha = 0.6) +
      geom_text_repel(max.overlaps = Inf,
                      min.segment.length = 0) 

#### ggplot NeuronType####
# coming soon

#### patchwork 2plot####
# coming soon

#### ggsave####
ggsave(filename = args_output, 
       plot = gg,
       dpi = 100, 
       width = 25.0, 
       height = 20.0,
       limitsize = FALSE)