source("src/functions_WTS4_MCMI_pairs.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
args_MCMI_matrix <- args[1]
args_count_sum <- args[2]
args_cell_count <- args[3]
args_merged_cls <- args[4]
args_NL <- args[5]
args_eval_label <- args[6]
args_output_Cluster <- args[7]
args_output_Neuron_type <- args[8]
args_output_Class <- args[9]

# #### test args####
# # matrix U
# args_MCMI_matrix <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/Merged_data/k_Number_9.RData")
# # Consistency
# args_count_sum <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/ClsCount/k_Number_9/df_count_sum.RData")
# # No_of_cells
# args_cell_count <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Distance/CellCount.RData")
# # Cluster
# args_merged_cls <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/Merged_cls/k_Number_9.RData")
# # Neuron type
# args_NL <- c("data/igraph/Fig1_HNS.RData")
# # Class
# args_eval_label <- c("data/WTS4_Eval_behavior_fix.xlsx")
# # Output Cluster
# args_output_Cluster <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/PairPlot/k_Number_9/Cluster.png")
# # Output Neuron_type
# args_output_Neuron_type <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/PairPlot/k_Number_9/Neuron_type.png")
# # Output Class
# args_output_Class <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/PairPlot/k_Number_9/Class.png")

#### load merged_data####
load(args_MCMI_matrix)
# matirx U to dataframe
merged_data$U %>% 
    as.data.frame() -> U_df
#### load Consistency####
# df_count_sum
load(args_count_sum)
df_count_sum %>% 
  dplyr::rename(., Consistency=Count_sum) -> df_Consistency
#### load No_of_cells####
# df_cell_count
load(args_cell_count)
df_cell_count %>% 
     dplyr::rename(., No_of_cells = CellCount) -> df_No_of_cells

#### load Cluster####
# merged_cls
load(args_merged_cls)
as.data.frame(merged_cls) %>% 
    .$merged_cls -> vec_Cluster

#### load Neuron type#####
# ig_Fig1_HNS
load(args_NL)
# node information convert dataframe
ig_Fig1_HNS %>% 
    igraph::as_data_frame(., what="vertices") -> df_node
# remove space from colnames
df_node %>% 
    names() %>% 
        str_replace_all(., c(" " = "")) -> names(df_node)
data.frame(cell_type = df_node$name,
           Neuron_type = df_node$NeuronType,
           # Neuron_group = df_node$NeuronGroup,
           stringsAsFactors = FALSE
           ) -> df_Neuron_type

#### load Class####
# WTS4_Eval_behavior_fix.xlsx
read.xlsx(args_eval_label,
          sheet = "Sheet1",
          rowNames = FALSE,
          colNames =TRUE) %>% 
  dplyr::rename(cell_type = celltype, 
                Class = class) -> df_Class

#### merge cluster####
U_df %>%
  dplyr::mutate(cell_type = rownames(as.data.frame(merged_cls))) %>% 
  dplyr::mutate(Cluster = as.data.frame(merged_cls)$merged_cls) -> U_df_CellType #細胞名 あとで除去
# クラスターを数字から文字列(離散値)
U_df_CellType$Cluster <- as.character(U_df_CellType$Cluster)
#### merge Consistency####
merge(U_df_CellType, 
      df_Consistency, 
      by.x = "cell_type", 
      by.y = "CellType", 
      all.x = TRUE) -> U_df_Consistency
#### merge No_of_cells####
merge(U_df_Consistency, 
      df_No_of_cells, 
      by.x = "cell_type", 
      by.y = "CellType", 
      all.x = TRUE) -> U_df_CellType

#### merge Neuron_type####
merge(U_df_CellType, 
      df_Neuron_type, #離散
      by.x = "cell_type", 
      by.y = "cell_type", 
      all.x = TRUE) -> U_df_NT
#### merge Class####
merge(U_df_NT, 
      df_Class, #離散
      by.x = "cell_type", 
      by.y = "cell_type", 
      all.x = TRUE) -> U_df_Class
U_df_Class$Class <- as.character(U_df_Class$Class)

#### U_input####
# 細胞名はggpairsで使わないので削除
U_input <- dplyr::select(U_df_Class, -cell_type)
U_col <- ncol(U_input) - 5
U_input %>% 
  dplyr::select(1:U_col, 
                Consistency,
                No_of_cells, 
                Cluster, 
                Class, 
                Neuron_type) -> U_input
# PC3 to NA
U_input <- mutate(U_input, Class = na_if(Class, "PC3"))

# set data col
ncol_gg <- ncol(U_input) - 3
#### ggpairs Cluster####
# ggpairs
U_input %>% 
  ggpairs(columns = 1:ncol_gg,
          aes_string(colour="Cluster", alpha=0.5),
          upper=list(continuous=wrap("cor", size=5))
          ) -> gg_Cluster
# ggsave
ggsave(filename = args_output_Cluster, 
       plot = gg_Cluster,
       dpi = 100, 
       width = 35.0, 
       height = 35.0,
       limitsize = FALSE)
#### ggpairs Neuron_type####
# ggpairs
U_input %>% 
  ggpairs(columns = 1:ncol_gg,
          aes_string(colour="Neuron_type", alpha=0.5),
          upper=list(continuous=wrap("cor", size=5))
          ) -> gg_Neuron_type
# ggsave
ggsave(filename = args_output_Neuron_type, 
       plot = gg_Neuron_type,
       dpi = 100, 
       width = 35.0, 
       height = 35.0,
       limitsize = FALSE)
#### ggpairs Class####
# ggpairs
U_input %>% 
  ggpairs(columns = 1:ncol_gg,
          aes_string(colour="Class", alpha=0.5),
          # upper = list(continuous = "blank")
          upper=list(continuous=wrap("cor", size=5))
          ) -> gg_Class
# ggsave
ggsave(filename = args_output_Class, 
       plot = gg_Class,
       dpi = 100, 
       width = 35.0, 
       height = 35.0,
       limitsize = FALSE)

# #### test ggplot####
# library(ggrepel)
# gg <- ggplot(U_df_Class, 
#        aes(x = V1,
#            y = V2, 
#            label = cell_type,
#            color = Class
#            ) 
#        ) +
#   geom_point(alpha = 0.6) +
#   geom_label_repel(max.overlaps = Inf,
#                    min.segment.length = 0,
#                    force = 6.0) 
# 
# gg <- ggplot(U_df_Class, 
#              aes(x = V1,
#                  y = V2, 
#                  label = cell_type,
#                  color = Cluster
#                  ) 
#              ) +
#   geom_point(alpha = 0.6) +
#   geom_label_repel(max.overlaps = Inf,
#                    min.segment.length = 0,
#                    force = 6.0) 
