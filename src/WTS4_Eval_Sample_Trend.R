source("src/functions_WTS4_Eval_Sample_Trend.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
# input merged_cls
args_input_cls <- args[1]
# output ggplot
args_output <- args[2]
# params merged_cls
args_input_merged_cls <- args[3]
# input merged_distance
args_input_MCMIHOOI <- args[4]
# params dist path
args_input_path<- args[5]

# #### test args####
# # input sample_cls
# args_input_cls <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Cluster_sample/k_Number_9/sample_cls.RData")
# # output ggplot
# args_output <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/DimReduc_sample/k_Number_9/Eval_Sample_trend.png")
# # params merged_cls
# args_input_merged_cls <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/Merged_cls/k_Number_9.RData")
# # input merged_distance
# args_input_MCMIHOOI <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/Merged_data/k_Number_9.RData")
# # params dist path
# args_input_path<- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Distance")

# #### test args n1_28sample####
# # input sample_cls
# args_input_cls <- c("output/WTS4/n1_28sample/stimAfter/SBD_abs/Cluster_sample/k_Number_9/sample_cls.RData")
# # output ggplot
# args_output <- c("output/WTS4/n1_28sample/stimAfter/SBD_abs/DimReduc_sample/k_Number_9/Eval_Sample_trend.png")
# # params merged_cls
# args_input_merged_cls <- c("output/WTS4/n1_28sample/stimAfter/SBD_abs/MCMIHOOI/Merged_cls/k_Number_9.RData")
# # input merged_distance
# args_input_MCMIHOOI <- c("output/WTS4/n1_28sample/stimAfter/SBD_abs/MCMIHOOI/Merged_data/k_Number_9.RData")
# # params dist path
# args_input_path<- c("output/WTS4/n1_28sample/stimAfter/SBD_abs/Distance")

##### load sample_cls list####
load(args_input_cls)
lapply(C, function(x) {
  data.frame(CellType = attr(x, "names"),
             Clusters = as.numeric(x),
             stringsAsFactors = FALSE,
             row.names = NULL
  )
}
) -> df_cls_list

##### load merged_cls####
load(args_input_merged_cls)
data.frame(CellType = names(merged_cls),
           Classes = merged_cls,
           stringsAsFactors = FALSE,
           row.names = NULL
) -> df_merged_cls
lapply(df_cls_list, function(x) {
  merge(x, 
        df_merged_cls, 
        by.x = "CellType", 
        by.y = "CellType", 
        all.y = TRUE
  ) -> df_cls_label_NA
  na.omit(df_cls_label_NA)
}
) -> df_cls_label

#### Eval vector####
unlist(lapply(df_cls_label, function(x) {
  clusters <- x$Clusters
  classes <- x$Classes
  adjustedRandIndex(clusters, classes)
}
)
) -> ARI_value

#### fix sample number sort####
input_path_list <- list.files(args_input_path, pattern="SampleNumber_", full.names=TRUE)
input_path_list %>% 
    str_remove(., args_input_path) %>% 
    str_remove(., "/SampleNumber_") %>% 
    str_remove(., ".RData") %>% 
    as.numeric() %>% 
    sort() -> sample_sort_num

#### create dataframe####
df_eval <- data.frame(SampleNumber = sample_sort_num,
                      ARI = ARI_value,
                      stringsAsFactors = FALSE
                      )
#### add annotated count####
#数字を除く
annotated_count <- unlist(lapply(df_cls_label, function(x){nrow(na.omit(x))}))
df_eval$annotated_count <- annotated_count

df_eval_wide <- df_eval
df_eval_wide$SampleNumber <- as.character(df_eval_wide$SampleNumber)

#### transform long format####
df_eval %>% 
  pivot_longer(col= -SampleNumber, 
               names_to = "Eval", 
               values_to ="Eval_Value") -> df_eval_long

#### MCMI weight table####
# merged_data
load(args_input_MCMIHOOI)
data.frame(SampleNumber = as.character(sample_sort_num),
           weight = merged_data$W,
           stringsAsFactors = FALSE) %>% 
    mutate(weight_abs =abs(weight)) %>% 
    dplyr::arrange(desc(weight_abs)) -> df_weight

#### ggplot ARI####
g1 <- ggplot(df_eval_wide, aes(x = SampleNumber, y= ARI , group=1))
g1 <- g1 + geom_line(color = "red", size= 2)
g1 <- g1 + scale_x_discrete(limits=df_weight$SampleNumber)
g1 <- g1 + theme_half_open()
g1 <- g1 + theme(text = element_text(size = 36))
g1 <- g1 + theme(axis.title.y=element_text(colour = "red",size = 36))
g1 <- g1 + geom_smooth(method="lm", size =0.5, se = TRUE, alpha = 0.4, color = "red")
g1 <- g1 + theme(legend.position = 'none')

g2 <- ggplot(df_eval_wide, aes(x = SampleNumber, y= annotated_count, group=1))
g2 <- g2 + geom_line(color = "black", size= 2)
g2 <- g2 + scale_x_discrete(limits=df_weight$SampleNumber)
g2 <- g2 + scale_y_continuous(position = "right")
g2 <- g2 + geom_smooth(method="lm", size =0.5, se = TRUE, alpha = 0.4, color = "black")
g2 <- g2 + theme_half_open()
# theme_half_openで軸の書式(色・文字サイズ)がリセットされる
g2 <- g2 + theme(text = element_text(size = 36))
g2 <- g2 + theme(axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank())

aligned_plots_result <- cowplot::align_plots(g1, g2, align="hv", axis="tblr")
gg_label_ARI <- cowplot::ggdraw(aligned_plots_result[[1]]) + cowplot::draw_plot(aligned_plots_result[[2]])

#### ggtexttable weight####
df_weight %>% 
  rownames_to_column("Ranking") %>% 
  mutate_if(is.numeric, round, digits = 3) %>% 
  ggtexttable(rows = NULL, theme = ttheme(base_size = 48)) -> gg_weight_table

#### prop.trend.test####
merge(df_weight, 
      df_eval_wide, 
      by.x = "SampleNumber", 
      by.y = "SampleNumber", 
      all.x = TRUE) %>% 
    dplyr::arrange(desc(weight_abs)) -> df_eval_weight
# アノテーション数の場合 CA検定
ann <- df_eval_weight$annotated_count
n_ann <- rep(nrow(df_merged_cls), length=length(ann))
trendtest_anotation <- prop.trend.test(ann, n_ann)$p.value
# ARIの場合 JT検定(ヨンキー検定)
ari <- df_eval_weight$ARI
#ヨンキー検定はケンドールの順位相関と同じ考え方
trendtest_ARI <- cor.test(ari, seq(length(ari)), method="kendall")$ p.value
#### ggtexttable trend####
data.frame(ann_CA_pvalue = signif(trendtest_anotation, digits = 3),
           ARI_JT_pvalue = signif(trendtest_ARI, digits = 3),
           stringsAsFactors = FALSE,
           row.names = NULL) -> trendtest_table
trendtest_table %>% 
    ggtexttable(rows = NULL, 
                theme = ttheme(base_size = 60, base_style ="mBlue")) -> gg_trendtest_table

#### graph title####
str_remove(args_output, 
           "output/WTS4/normalize_1/stimAfter/") %>% 
  str_remove(., 
             ".png") %>% 
  str_remove(., 
             "DimReduc_sample/") -> plot_title

#### patchwork####
gg <- gg_weight_table +
  gg_label_ARI +
  gg_trendtest_table +
  plot_layout(nrow = 1) +
  plot_annotation(title = plot_title,
                  caption = 'made with patchwork',
                  theme = theme(plot.title = element_text(size = 60, hjust = 0.5))
  )

#### ggsave####
ggsave(filename = args_output, 
       plot = gg,
       dpi = 100, 
       width = 60.0,
       height = 20.0,
       limitsize = FALSE)