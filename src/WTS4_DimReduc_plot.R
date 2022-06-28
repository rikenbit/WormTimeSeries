source("src/functions_WTS4_DimReduc_plot.R")

#### args setting####
args <- commandArgs(trailingOnly = T)

# input Dimensionality Reduction cord
args_input_cord <- args[1]
# input merged_cls
args_input_cls <- args[2]
# dimentionaly reduction
args_DimReduc <- args[3]

# Neuron Label Path
args_NL <- args[4]
# Evaluation label list
args_eval_label <- args[5]
# df cell count
args_cell_count <- args[6]
# df_count_sum
args_count_sum <- args[7]

args_output_gg_cls <- args[8]
args_output_gg_NT <- args[9]
args_output_gg_eval_label <- args[10]
args_output_gg_count_sum <- args[11]
args_output_gg_cell_count <- args[12]
# Additional File 1: Cellular labels to interpret the clustering results
args_output_csv <- args[13]

# #### test args####
# args_input_cord <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/CSPA/Merged_tsne_cord/k_Number_5.RData")
# args_input_cls <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/CSPA/Merged_cls/k_Number_5.RData")
# args_DimReduc <- c("tsne")

# args_NL <- c("data/igraph/Fig1_HNS.RData")
# args_eval_label <- c("data/WTS4_Eval_behavior_ACF.xlsx")
# args_cell_count <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Distance/CellCount.RData")
# args_count_sum <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/ClsCount/k_Number_5/df_count_sum.RData")

# args_output_gg_cls <- c('output/WTS4/normalize_1/stimAfter/SBD_abs/CSPA/Merged_tsne_plot/k_Number_5_gg_cls.png')
# args_output_gg_NT <- c('output/WTS4/normalize_1/stimAfter/SBD_abs/CSPA/Merged_tsne_plot/k_Number_5_gg_NT.png')
# args_output_gg_eval_label <- c('output/WTS4/normalize_1/stimAfter/SBD_abs/CSPA/Merged_tsne_plot/k_Number_5_gg_eval_label.png')
# args_output_gg_count_sum <- c('output/WTS4/normalize_1/stimAfter/SBD_abs/CSPA/Merged_tsne_plot/k_Number_5_gg_count_sum.png')
# args_output_gg_cell_count <- c('output/WTS4/normalize_1/stimAfter/SBD_abs/CSPA/Merged_tsne_plot/k_Number_5_gg_cell_count.png')
# args_output_csv <- c('output/WTS4/normalize_1/stimAfter/SBD_abs/CSPA/Merged_tsne_plot/label_table_k9.csv')

#### load Dimensionality Reduction cord####
load(args_input_cord)

#### load cls number object####
load(args_input_cls)
df_cls <- as.data.frame(merged_cls)
df_cls <- data.frame(
    cell_type = rownames(df_cls),
    cls = df_cls$merged_cls,
    stringsAsFactors = FALSE
    )

#### merge cord and cls####
df_cord_cls <- merge(df_cord, 
                     df_cls, 
                     by.x = "cell_type", 
                     by.y = "cell_type", 
                     all.x = TRUE)

#### load Neuron Label####
load(args_NL)
# node information convert dataframe
ig_Fig1_HNS %>% 
    igraph::as_data_frame(., what="vertices") -> df_node
# remove space from colnames
df_node %>% 
    names() %>% 
        str_replace_all(., c(" " = "")) -> names(df_node)
df_NL <- data.frame(cell_type = df_node$name,
                    NeuronType = df_node$NeuronType,
                    NeuronGroup = df_node$NeuronGroup,
                    stringsAsFactors = FALSE
                    )

#### merge Neuron Label####
df_cord_cls_NL <- merge(df_cord_cls, 
                           df_NL, 
                           by.x = "cell_type", 
                           by.y = "cell_type", 
                           all.x = TRUE
                        )
# unknown label is NA
######## load&merge cell_count########
# df_cell_count
load(args_cell_count)

df_cord_cls_NL_count <- merge(df_cord_cls_NL, 
                              df_cell_count, 
                              by.x = "cell_type", 
                              by.y = "CellType",
                              all.x = TRUE
                              )
######## load&merge df_count_sum########
load(args_count_sum)

df_cord_cls_NL_count_sum <- merge(df_cord_cls_NL_count, 
                                  df_count_sum, 
                                  by.x = "cell_type", 
                                  by.y = "CellType",
                                  all.x = TRUE
                                  )
merge(df_cord_cls_NL_count, 
      df_count_sum, 
      by.x = "cell_type", 
      by.y = "CellType",
      all.x = TRUE
      ) %>% 
    # NAを0に変換
    replace_na(., replace = list(Count_sum = 0)) -> df_cord_cls_NL_count_sum
#### load eval_label####
read.xlsx(args_eval_label,
          sheet = "Sheet1",
          rowNames = FALSE,
          colNames =TRUE
          ) %>% 
    dplyr::rename(CellType = celltype, 
                  Classes = class) -> df_eval_label

#### merge eval_label####
df_merged <- merge(df_cord_cls_NL_count_sum, 
                   df_eval_label, 
                   by.x = "cell_type", 
                   by.y = "CellType", 
                   all.x = TRUE
                   )

#### ggplot axis label name####
if (args_DimReduc == "tsne") {
  cord_x <- c("t-SNE-1")
  cord_y <- c("t-SNE-2")
} else if (args_DimReduc == "umap") {
  cord_x <- c("UMAP-1")
  cord_y <- c("UMAP-2")
} else {
  cord_x <- c("cord-1")
  cord_y <- c("cord-2")
}


#### ggplot cls####
gg_cls <- ggplot(df_merged, 
                 aes(x = cord_1,
                     y = cord_2, 
                     label = cell_type,
                     color = factor(cls)
                    )
                ) + 
    labs(color = "Cluster") +
    geom_point(size = 6.0, 
               alpha = 0.6) +
    geom_label_repel(max.overlaps = Inf,
                     min.segment.length = 0,
                     size = 9.0,
                     force = 6.0) +# ラベル間の反発力
    theme(text = element_text(size = 60)) +
    labs(x = cord_x,
         y = cord_y)

#### ggplot NeuronType####
gg_NT <- ggplot(df_merged, 
                 aes(x = cord_1,
                     y = cord_2, 
                     label = cell_type,
                     color = factor(NeuronType)
                     )
                ) + 
    labs(color = " Neuron type") +
    geom_point(size = 6.0, 
               alpha = 0.6) +
    geom_label_repel(max.overlaps = Inf,
                     min.segment.length = 0,
                     size = 9.0,
                     force = 6.0) +# ラベル間の反発力
    theme(text = element_text(size = 60)) +
    labs(x = cord_x,
         y = cord_y) 

#### ggplot eval_label####
gg_eval_label <- ggplot(df_merged, 
                 aes(x = cord_1,
                     y = cord_2, 
                     label = cell_type,
                     color = factor(Classes)
                     )
                ) + 
    labs(color = "Class") +
    geom_point(size = 6.0, 
               alpha = 0.6) +
    geom_label_repel(max.overlaps = Inf,
                     min.segment.length = 0,
                     size = 9.0,
                     force = 6.0) +# ラベル間の反発力
    theme(text = element_text(size = 60)) +
    labs(x = cord_x,
         y = cord_y) 
#### ggplot cell_count####
gg_cell_count <- ggplot(df_merged, 
                        aes(x = cord_1,
                            y = cord_2, 
                            label = cell_type,
                            color = CellCount
                            )
                        ) + 
  scale_color_viridis_c(option = "D") +
  labs(color = "No. of cells") +
  geom_point(size = 6.0, 
             alpha = 0.6) +
  geom_label_repel(max.overlaps = Inf,
                   min.segment.length = 0,
                   size = 9.0,
                   force = 6.0) +# ラベル間の反発力
    theme(text = element_text(size = 60)) +
    labs(x = cord_x,
         y = cord_y) +
  theme(legend.key.height = unit(1.5, "cm")) +
  theme(legend.key.width = unit(1.5, "cm"))
#### ggplot count_sum####
gg_count_sum <- ggplot(df_merged, 
                        aes(x = cord_1,
                            y = cord_2, 
                            label = cell_type,
                            color = Count_sum
                            )
                        ) + 
  scale_color_viridis_c(option = "D") +
  labs(color = "Consistency") +
  geom_point(size = 6.0, 
             alpha = 0.6) +
  geom_label_repel(max.overlaps = Inf,
                   min.segment.length = 0,
                   size = 9.0,
                   force = 6.0) +# ラベル間の反発力
    theme(text = element_text(size = 60)) +
    labs(x = cord_x,
         y = cord_y) +
  theme(legend.key.height = unit(1.5, "cm")) +
  theme(legend.key.width = unit(1.5, "cm"))

#### ggsave####
ggsave(filename = args_output_gg_cls, 
       plot = gg_cls,
       dpi = 100, 
       width = 25.0, 
       height = 20.0,
       limitsize = FALSE)
ggsave(filename = args_output_gg_NT, 
       plot = gg_NT,
       dpi = 100, 
       width = 25.0, 
       height = 20.0,
       limitsize = FALSE)
ggsave(filename = args_output_gg_eval_label, 
       plot = gg_eval_label,
       dpi = 100, 
       width = 25.0, 
       height = 20.0,
       limitsize = FALSE)
ggsave(filename = args_output_gg_count_sum, 
       plot = gg_count_sum,
       dpi = 100, 
       width = 25.0, 
       height = 20.0,
       limitsize = FALSE)
ggsave(filename = args_output_gg_cell_count, 
       plot = gg_cell_count,
       dpi = 100, 
       width = 25.0, 
       height = 20.0,
       limitsize = FALSE)

#### write.csv#####
label_table <- data.frame(
  Cell_type = df_merged$cell_type,
  Cluster = df_merged$cls,
  Neuron_type = df_merged$NeuronType,
  Class = df_merged$Classes,
  Consistency = df_merged$Count_sum,
  No_cells = df_merged$CellCount,
  stringsAsFactors = FALSE
  )

write.csv(label_table, 
          args_output_csv, 
          row.names=FALSE)