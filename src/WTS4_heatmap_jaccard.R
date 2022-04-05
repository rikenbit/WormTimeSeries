source("src/functions_WTS4_heatmap_jaccard.R")

#### args setting####
#### test args####

#### MCMI####
args_input_MCMI <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/Merged_cls/k_Number_9.RData")
load(args_input_MCMI)
merged_cls %>% 
    as.data.frame() %>% 
        rownames_to_column("cell_type") %>% 
            dplyr::select(cell_type=1, cluster=2) -> MCMI_merged_cls

#### CSPA####
args_input_CSPA <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/CSPA/Merged_cls/k_Number_5.RData")
load(args_input_CSPA)
merged_cls %>% 
    as.data.frame() %>% 
        rownames_to_column("cell_type") %>% 
            dplyr::select(cell_type=1, cluster=2) -> CSPA_merged_cls

#### jaccard sapply#####
# 行がCSPAのクラスタ、列がMCMIのクラスタ
sapply(unique(MCMI_merged_cls$cluster), 
       function(x) {
         sapply(unique(CSPA_merged_cls$cluster), 
                function(y){
                  # select MCMI cluster
                  MCMI_merged_cls %>% 
                    dplyr::filter(cluster==x) -> MCMI_cluster_x
                  # select CSPA cluster
                  CSPA_merged_cls %>% 
                    dplyr::filter(cluster==y) -> CSPA_cluster_y
                  # jaccard
                  jaccard(MCMI_cluster_x$cell_type, 
                          CSPA_cluster_y$cell_type) -> return_object
                  
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
  rownames_to_column("CSPA_k") %>% 
  pivot_longer(-CSPA_k,
               names_to = "MCMI_k", 
               values_to = "jaccard") -> df_ghm

#### plot_title####
str_remove(args_output_heatmap, 
           "output/WTS4/n1_28sample/stimAfter/") %>% 
  str_remove(., 
             ".png") -> plot_title

#### geom_tile####
ghm <- ggplot_ghm(df_ghm) + 
  # ggtitle(plot_title) + 
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

