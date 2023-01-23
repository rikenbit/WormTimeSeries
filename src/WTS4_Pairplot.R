source("src/functions_WTS4_Pairplot.R")

#### args setting####
args <- commandArgs(trailingOnly = T)

#### test args####
# Cluster
args_merged_cls <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/Merged_cls/k_Number_6.RData")
# Neuron type
args_NL <- c("data/igraph/Fig1_HNS.RData")
# Class
args_eval_label <- c("data/WTS4_Eval_behavior_fix.xlsx")
# Consistency
args_count_sum <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/ClsCount/k_Number_6/df_count_sum.RData")
# No_of_cells
args_cell_count <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Distance/CellCount.RData")
# Output
args_output <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/PairPlot/k_Number_6/Pairplot.svg")
args_output_Cluster <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/PairPlot/k_Number_6/Pairplot_Cluster.svg")
args_output_NT <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/PairPlot/k_Number_6/Pairplot_Neuron_type.svg")
args_output_Class <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/PairPlot/k_Number_6/Pairplot_Class.svg")
#### load Cluster####
load(args_merged_cls)
as.data.frame(merged_cls) %>%
  rownames_to_column("cell_type") -> U_df_Cluster
data.frame(cell_type = U_df_Cluster$cell_type,
           Cluster = as.character(U_df_Cluster$merged_cls),
           stringsAsFactors = FALSE
           ) -> U_df_CellType
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
#### merge Consistency####
merge(U_df_Class,
      df_Consistency,
      by.x = "cell_type",
      by.y = "CellType",
      all.x = TRUE) -> U_df_Consistency
#### merge No_of_cells####
merge(U_df_Consistency,
      df_No_of_cells,
      by.x = "cell_type",
      by.y = "CellType",
      all.x = TRUE) -> U_df_cells

#### U_input####
# PC3 to NA
U_input <- mutate(U_df_cells, Class = na_if(Class, "PC3"))
#### ggpairs####
# ggpairs
U_input %>%
  ggpairs(columns = 2:6,
          upper=list(continuous=wrap("cor", size=6),
                     combo="facethist",
                     discrete = "facetbar"),
          lower=list(continuous="smooth_loess",
                     combo ="box",
                     discrete ="facetbar")
          ) -> gg
ggsave(filename = args_output,
       plot = gg,
       dpi = 100,
       width = 19.0,
       height = 19.0,
       limitsize = FALSE)

#### ggpairs Cluster####
# ggpairs
U_input %>%
  ggpairs(columns = 2:6,
          aes_string(colour="Cluster", alpha=0.5),
          upper=list(continuous=wrap("cor", size=6),
                     combo="facethist",
                     discrete = "facetbar"),
          lower=list(continuous="points",
                     combo ="box",
                     discrete ="facetbar")
          ) -> gg_Cluster
ggsave(filename = args_output_Cluster,
       plot = gg_Cluster,
       dpi = 100,
       width = 19.0,
       height = 19.0,
       limitsize = FALSE)
#### ggpairs Neuron_type####
U_input %>%
  ggpairs(columns = 2:6,
          aes_string(colour="Neuron_type", alpha=0.5),
          upper=list(continuous=wrap("cor", size=6),
                     combo="facethist",
                     discrete = "facetbar"),
          lower=list(continuous="points",
                     combo ="box",
                     discrete ="facetbar")
          ) -> gg_NT
ggsave(filename = args_output_NT,
       plot = gg_NT,
       dpi = 100,
       width = 19.0,
       height = 19.0,
       limitsize = FALSE)
#### ggpairs Class####
U_input %>%
  ggpairs(columns = 2:6,
          aes_string(colour="Class", alpha=0.5),
          upper=list(continuous=wrap("cor", size=6),
                     combo="facethist",
                     discrete = "facetbar"),
          lower=list(continuous="points",
                     combo ="box",
                     discrete ="facetbar")
          ) -> gg_Class
ggsave(filename = args_output_Class,
       plot = gg_Class,
       dpi = 100,
       width = 19.0,
       height = 19.0,
       limitsize = FALSE)