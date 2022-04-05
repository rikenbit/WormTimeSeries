source("src/functions_WTS4_heatmap_dist.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
# 距離行列
args_input_distance <- args[1]
# 細胞名の和集合取得
args_params_celltype <- args[2]
# ヒートマップ
args_output_heatmap <- args[3]
  
#### test args####
# # 距離行列
# args_input_distance <- c("output/WTS4/n1_28sample/stimAfter/SBD_abs/Distance/SampleNumber_1.RData")
# # 細胞名の和集合取得
# args_params_celltype <- c("output/WTS4/n1_28sample/stimAfter/SBD_abs/MCMIHOOI/Merged_cls/k_Number_9.RData")
# # ヒートマップ
# args_output_heatmap <- c("output/WTS4/n1_28sample/stimAfter/SBD_abs/Dist_heatmap/SampleNumber_1.png")

#### load distance####
load(args_input_distance)

#### filter distance####
d_cell <- attr(d, "Labels")
# 数字の細胞削除
d_cell_f <- d_cell[grep("^[0-9]", d_cell, invert=TRUE)]
d_f <- dist_subset(d, d_cell_f)

#### transform dist -> matrix -> data.frame####
d_f %>% 
  as.matrix() %>% 
  data.frame() %>%
  rownames_to_column("row_celltype") %>% 
  pivot_longer(-row_celltype, 
               names_to = "col_celltype", 
               values_to = "dist_value") -> long_d_f

#### load celltype list####
load(args_params_celltype)

#### get all celltype####
as.data.frame(merged_cls) %>% 
  rownames_to_column("cell_type") %>% 
  .$cell_type -> all_celltype

#### create null matirx####
mat_celltype <- matrix(nrow=length(all_celltype),
                       ncol=length(all_celltype))
colnames(mat_celltype) <- all_celltype
rownames(mat_celltype) <- all_celltype

#### transform null matrix to data.frame####
mat_celltype %>% 
  data.frame() %>%
  rownames_to_column("row_celltype") %>% 
  pivot_longer(-row_celltype, 
               names_to = "col_celltype", 
               values_to = "dist_value") -> long_celltype

#### merge df####
left_join(long_celltype,
          long_d_f,
          by = c("row_celltype","col_celltype")) %>% 
  dplyr::select(row_celltype=1,
                col_celltype=2,
                dist_value=4) -> df_ghm

#### levels cell type####
as.data.frame(merged_cls) %>%
  rownames_to_column("cell_type") -> df_cls
df_cls %>% 
  dplyr::arrange(merged_cls) %>% 
  .$cell_type -> lbs
df_ghm_lbs <- df_ghm
df_ghm_lbs$row_celltype <- factor(df_ghm_lbs$row_celltype, levels = lbs)
df_ghm_lbs$col_celltype <- factor(df_ghm_lbs$col_celltype, levels = lbs)

#### plot_title####
str_remove(args_output_heatmap, 
           "output/WTS4/n1_28sample/stimAfter/") %>% 
  str_remove(., 
             ".png") -> plot_title

#### geom_tile####
ghm <- ggplot_ghm(df_ghm_lbs) + 
      ggtitle(plot_title) + 
      theme(plot.title = element_text(size = 30, hjust = 0.5)) +
      theme(axis.title = element_text(size = 40)) + 
      theme(legend.key.height = unit(1.5, "cm")) +
      theme(legend.key.width = unit(1.5, "cm")) +
      theme(legend.text = element_text(size = 30)) + 
      theme(legend.title = element_text(size = 30))

#### ggsave####
ggsave(filename = args_output_heatmap, 
       plot = ghm, 
       dpi = 80, 
       width = 23.5, 
       height = 22.0
       )
