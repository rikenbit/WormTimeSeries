source("src/functions_WTS5_Jaccard_ARI.R")

#### test args####
args_input_label <- c("data/WTS4_Eval_behavior_ACF.xlsx")

#### test args####
# # mSBD MCMI
# args_input_cls <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/Merged_cls/k_Number_6.RData")
# args_output_heatmap <- c("output/WTS5/normalize_1/stimAfter/Jaccard/Heatmap_mSBD_MCMI_Classes.eps")
# # mSBD CSPA
# args_input_cls <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/CSPA/Merged_cls/k_Number_5.RData")
# args_output_heatmap <- c("output/WTS5/normalize_1/stimAfter/Jaccard/Heatmap_mSBD_CSPA_Classes.eps")
# # EUCL MCMI
# args_input_cls <- c("output/WTS4/normalize_1/stimAfter/EUCL/MCMIHOOI/Merged_cls/k_Number_10.RData")
# args_output_heatmap <- c("output/WTS5/normalize_1/stimAfter/Jaccard/Heatmap_EUCL_MCMI_Classes.eps")
# EUCL CSPA
args_input_cls <- c("output/WTS4/normalize_1/stimAfter/EUCL/CSPA/Merged_cls/k_Number_4.RData")
args_output_heatmap <- c("output/WTS5/normalize_1/stimAfter/Jaccard/Heatmap_EUCL_CSPA_Classes.eps")

#### MCMI####
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

#### jaccard sapply#####
sapply(unique(merged_input$cluster),
       function(x) {
           sapply(label_order,
                  function(y){
                      # select MCMI cluster
                      merged_input %>%
                          dplyr::filter(cluster==x) -> cls_cluster_x
                      # select CSPA cluster
                      merged_input %>%
                          dplyr::filter(class==y) -> label_cluster_y
                      # jaccard
                      jaccard(cls_cluster_x$cell_type,
                              label_cluster_y$cell_type) -> return_object

                      return_object
                  }
           )
       }
) -> sapply_jaccard

#### transform data.frame####
sapply_jaccard %>%
    as.data.frame() -> df_jaccard
colnames(df_jaccard) <- seq(1:ncol(df_jaccard))

df_jaccard %>% 
    rownames_to_column("label") %>% 
        pivot_longer(-label,
                     names_to = "cls", 
                     values_to = "jaccard") -> df_ghm
df_ghm <-as.data.frame(df_ghm)
# sort label
df_ghm$label <- factor(df_ghm$label, levels = c("NaCl", "PC1_pos", "PC1_neg", "PC2", "PC3", "NA"))
df_ghm$cls <- factor(df_ghm$cls, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))
#### geom_tile####
df_ghm |> 
    dplyr::filter(label!="NA") |> 
        ggplot_ghm() -> ghm
#### ggsave####
# # mSBD MCMI
# ggsave(filename = args_output_heatmap, 
#        plot = ghm, 
#        dpi = 80, 
#        width = 18.0, 
#        height = 20.0,
#        device="eps")
# # mSBD CSPA
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
# EUCL CSPA
ggsave(filename = args_output_heatmap, 
       plot = ghm, 
       dpi = 80, 
       width = 20.0, 
       height = 18.0,
       device="eps")