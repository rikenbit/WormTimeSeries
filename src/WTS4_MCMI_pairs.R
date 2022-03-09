source("src/functions_WTS4_MCMI_pairs.R")

#### args setting####
#### test args####
args_MCMI_matrix <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/Merged_data/k_Number_9.RData")
args_MCMI_cluster <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/Merged_cls/k_Number_9.RData")
args_output <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/PairPlot/k_Number_9/cluster.png")
args_cell_count <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Distance/CellCount.RData")

#### test iris####
# https://qiita.com/h_kobayashi1125/items/46bc28a26f888d03cee3
pairs(iris, panel = panel.smooth)
install.packages("psych")
library(psych)
psych::pairs.panels(iris)
install.packages("corrplot")
library(corrplot)
corrplot::corrplot(cor(iris[,-5]))
install.packages("qgraph")
library(qgraph)
qgraph(cor(iris[,-5]),edge.labels=T )


install.packages("GGally")
library(ggplot2)
library(GGally)
gg <- ggpairs(iris, aes_string(colour="Species", alpha=0.5))

iris_test <- iris
iris_test$Species <- as.character(U_Cluster[1:150])
ggpairs(iris_test, aes_string(colour="Species", alpha=0.5))
#### test real data####
# load merged_data
load(args_MCMI_matrix)
# matirx U
merged_data$U %>% 
  as.data.frame() -> U_df

#### load merged_cls####
load(args_MCMI_cluster)
# Cluster
as.data.frame(merged_cls) %>% 
  .$merged_cls -> U_Cluster

#### test ggpair####
U_df %>% 
  dplyr::mutate(Cluster=U_Cluster) -> U_input
U_input$Cluster <- as.character(U_input$Cluster)
gg <- ggpairs(U_input, aes_string(colour="Cluster", alpha=0.5))
#### load df_cell_count####
load(args_cell_count)
#### test ggpair####
U_df %>% 
  dplyr::mutate(CellCount=df_cell_count$CellCount) %>% 
  dplyr::mutate(Cluster=U_Cluster) -> U_input
U_input$Cluster <- as.character(U_input$Cluster)
gg <- ggpairs(U_input, aes_string(colour="Cluster", alpha=0.5))
# ggpairs(U_input, aes(colour=CellCount, alpha=0.5))


#### test diamonds####
data(diamonds, package="ggplot2")
diamonds.samp <- diamonds[sample(1:dim(diamonds)[1], 1000), ]

ggpairs(
  diamonds.samp[, 1:5],
  mapping = ggplot2::aes(color = cut),
  upper = list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"),
  lower = list(continuous = wrap("points", alpha = 0.3), combo = wrap("dot_no_facet", alpha = 0.4)),
  title = "Diamonds"
)
#### ggsave####
ggsave(filename = args_output, 
       plot = gg,
       dpi = 100, 
       width = 20.0, 
       height = 20.0,
       limitsize = FALSE)