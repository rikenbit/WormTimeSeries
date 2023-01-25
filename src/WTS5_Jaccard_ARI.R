source("src/functions_WTS5_Jaccard_ARI.R")

#### test args####
args_input_label <- c("data/WTS4_Eval_behavior_ACF.xlsx")

#### Jaccard args####
# mSBD MCMI
args_input_cls <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/Merged_cls/k_Number_6.RData")
args_output_heatmap <- c("output/WTS5/normalize_1/stimAfter/Jaccard/Heatmap_mSBD_MCMI_Classes.eps")
# # mSBD CSPA
# args_input_cls <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/CSPA/Merged_cls/k_Number_5.RData")
# args_output_heatmap <- c("output/WTS5/normalize_1/stimAfter/Jaccard/Heatmap_mSBD_CSPA_Classes.eps")
# EUCL MCMI
# args_input_cls <- c("output/WTS4/normalize_1/stimAfter/EUCL/MCMIHOOI/Merged_cls/k_Number_10.RData")
# args_output_heatmap <- c("output/WTS5/normalize_1/stimAfter/Jaccard/Heatmap_EUCL_MCMI_Classes.eps")
# # EUCL CSPA
# args_input_cls <- c("output/WTS4/normalize_1/stimAfter/EUCL/CSPA/Merged_cls/k_Number_4.RData")
# args_output_heatmap <- c("output/WTS5/normalize_1/stimAfter/Jaccard/Heatmap_EUCL_CSPA_Classes.eps")


# #### ARI args####
# # mSBD MCMI
# args_input_cls <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/Merged_cls/k_Number_6.RData")
# args_output_heatmap <- c("output/WTS5/normalize_1/stimAfter/ARI/Heatmap_mSBD_MCMI_Classes.eps")
# mSBD CSPA
# args_input_cls <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/CSPA/Merged_cls/k_Number_5.RData")
# args_output_heatmap <- c("output/WTS5/normalize_1/stimAfter/ARI/Heatmap_mSBD_CSPA_Classes.eps")
# # EUCL MCMI
# args_input_cls <- c("output/WTS4/normalize_1/stimAfter/EUCL/MCMIHOOI/Merged_cls/k_Number_10.RData")
# args_output_heatmap <- c("output/WTS5/normalize_1/stimAfter/ARI/Heatmap_EUCL_MCMI_Classes.eps")
# EUCL CSPA
# args_input_cls <- c("output/WTS4/normalize_1/stimAfter/EUCL/CSPA/Merged_cls/k_Number_4.RData")
# args_output_heatmap <- c("output/WTS5/normalize_1/stimAfter/ARI/Heatmap_EUCL_CSPA_Classes.eps")

#### load clusterin result####
load(args_input_cls)
merged_cls %>%
    as.data.frame() %>%
        rownames_to_column("cell_type") %>%
            dplyr::select(cell_type=1, cluster=2) -> input_cls
#### load eval_label####
read.xlsx(args_input_label,
          sheet = "Sheet1",
          rowNames = FALSE,
          colNames =TRUE
          ) %>%
    dplyr::rename(cell_type = celltype,
                  class = class) -> input_label

merge(input_cls,
      input_label,
      by.x = "cell_type",
      by.y = "cell_type",
      all.x = TRUE) -> merged_input
merged_input$class <- replace_na(merged_input$class,"NA")

merged_input <- dplyr::arrange(merged_input, cluster)
label_order <- c("NaCl", "PC1_pos", "PC1_neg", "PC2", "PC3", "NA")

#### ARI wider####
merged_input %>% 
    dplyr::select(cell_type, cluster) %>% 
    mutate(value = 1) %>%
    pivot_wider(names_from = "cluster", 
                values_from = value, 
                values_fill = 0) %>% 
    column_to_rownames("cell_type") -> input_cls_wider
merged_input %>% 
    dplyr::select(cell_type, class) %>% 
    mutate(value = 1) %>%
    pivot_wider(names_from = "class", 
                values_from = value, 
                values_fill = 0) %>% 
    column_to_rownames("cell_type") -> input_label_wider

# #### jaccard sapply#####
# sapply(unique(merged_input$cluster),
#        function(x) {
#            sapply(label_order,
#                   function(y){
#                       # select MCMI cluster
#                       merged_input %>%
#                           dplyr::filter(cluster==x) -> cls_cluster_x
#                       # select CSPA cluster
#                       merged_input %>%
#                           dplyr::filter(class==y) -> label_cluster_y
#                       # jaccard
#                       jaccard(cls_cluster_x$cell_type,
#                               label_cluster_y$cell_type) -> return_object

#                       return_object
#                   }
#            )
#        }
# ) -> sapply_jaccard

#### ARI sapply####
sapply(colnames(input_cls_wider),
       function(x) {
           sapply(colnames(input_label_wider),
                  function(y){
                      # cluster
                      input_cls_wider[,x] -> cls_x
                      # label
                      input_label_wider[,y]-> label_y
                      # jaccard
                      adjustedRandIndex(cls_x, label_y) -> return_object
                      return_object
                  }
           )
       }
) -> sapply_ari

# #### transform data.frame####
# sapply_jaccard %>%
#     as.data.frame() -> df_jaccard
# colnames(df_jaccard) <- seq(1:ncol(df_jaccard))
# df_jaccard %>%
#     rownames_to_column("label") %>%
#         pivot_longer(-label,
#                      names_to = "cls",
#                      values_to = "jaccard") -> df_ghm
# df_ghm <-as.data.frame(df_ghm)
# # sort label
# df_ghm$label <- factor(df_ghm$label, levels = c("NaCl", "PC1_pos", "PC1_neg", "PC2", "PC3", "NA"))
# df_ghm$cls <- factor(df_ghm$cls, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))

#### ARI transform data.frame####
sapply_ari %>%
    as.data.frame() -> df_ari
colnames(df_ari) <- seq(1:ncol(df_ari))

df_ari %>% 
    rownames_to_column("label") %>% 
    pivot_longer(-label,
                 names_to = "cls", 
                 values_to = "ari") -> df_ghm_ari
df_ghm_ari <-as.data.frame(df_ghm_ari)
# sort label
df_ghm_ari$label <- factor(df_ghm_ari$label, levels = c("NaCl", "PC1_pos", "PC1_neg", "PC2", "PC3", "NA"))
df_ghm_ari$cls <- factor(df_ghm_ari$cls, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))

#### geom_tile####
df_ghm |> 
    dplyr::filter(label!="NA") |> 
        ggplot_ghm() -> ghm
#### ggsave####
# mSBD MCMI
ggsave(filename = args_output_heatmap,
       plot = ghm,
       dpi = 80,
       width = 18.0,
       height = 20.0,
       device="eps")
# mSBD CSPA
# ggsave(filename = args_output_heatmap,
#        plot = ghm,
#        dpi = 80,
#        width = 20.0,
#        height = 20.0,
#        device="eps")
# # EUCL MCMI
# ggsave(filename = args_output_heatmap,
#        plot = ghm,
#        dpi = 80,
#        width = 20.0,
#        height = 30.0,
#        device="eps")
# # EUCL CSPA
# ggsave(filename = args_output_heatmap,
#        plot = ghm,
#        dpi = 80,
#        width = 20.0,
#        height = 18.0,
#        device="eps")

#### ARI ggsave####
# df_ghm_ari |> 
#     dplyr::filter(label!="NA") |> 
#     ggplot_ghm() -> ghm_ari

# mSBD MCMI
# ggsave(filename = args_output_heatmap,
#        plot = ghm_ari,
#        dpi = 80,
#        width = 18.0,
#        height = 20.0,
#        device="eps")
# # mSBD CSPA
# ggsave(filename = args_output_heatmap, 
#        plot = ghm_ari, 
#        dpi = 80, 
#        width = 20.0, 
#        height = 20.0,
#        device="eps")
# # EUCL MCMI
# ggsave(filename = args_output_heatmap, 
#        plot = ghm_ari, 
#        dpi = 80, 
#        width = 20.0, 
#        height = 30.0,
#        device="eps")
# EUCL CSPA
# ggsave(filename = args_output_heatmap, 
#        plot = ghm_ari, 
#        dpi = 80, 
#        width = 20.0, 
#        height = 18.0,
#        device="eps")