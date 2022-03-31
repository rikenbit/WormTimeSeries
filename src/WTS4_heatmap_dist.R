source("src/functions_WTS4_heatmap_dist.R")

#### args setting####
#### test args####
# 距離行列
args_input_distance <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Distance/SampleNumber_1.RData")
# 細胞名の和集合取得
args_params_celltype <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/Merged_cls/k_Number_9.RData")
# ヒートマップ
args_output_heatmap <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Dist_heatmap/SampleNumber_1.RData")

#### load distance####
load(args_input_distance)
#### filter distance####
d_cell <- attr(d, "Labels")
# 数字の細胞削除
d_cell_f <- d_cell[grep("^[0-9]", d_cell, invert=TRUE)]
d_f <- dist_subset(d, d_cell_f)

#### transform dist to matrix####
#### transform matrix to data.frame####

#### load celltype list####
load(args_params_celltype)
#### 和集合の細胞名を取得####
as.data.frame(merged_cls) %>% 
  rownames_to_column("cell_type") %>% 
  .$cell_type -> all_celltype
#### transform null matrix to data.frame####


#### test data.frame####
# 空の行列作成
mat_celltype <- matrix(nrow=3,ncol=3)
# 行・列名を挿入
colnames(mat_celltype) <- all_celltype[1:3]
rownames(mat_celltype) <- all_celltype[1:3]

# データフレーム変換
mat_celltype %>% 
  data.frame() %>%
  rownames_to_column("row_celltype") %>% 
  pivot_longer(-row_celltype, 
               names_to = "col_celltype", 
               values_to = "dist_value") -> long_celltype

#### test dist####
as.matrix(d_f) %>% 
  .[2:4,2:4] -> test_d_f
# データフレーム変換
test_d_f %>% 
  data.frame() %>%
  rownames_to_column("row_celltype") %>% 
  pivot_longer(-row_celltype, 
               names_to = "col_celltype", 
               values_to = "dist_value") -> long_d_f

#### test merge####
left_join(long_celltype,
          long_d_f,
          by = c("row_celltype","col_celltype")) %>% 
  dplyr::select(row_celltype=1,
                col_celltype=2,
                dist_value=4) -> df_ghm

#### test geom_tile####
ghm <- ggplot_ghm(df_ghm)