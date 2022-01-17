source("src/functions_WTS4_EvalPlot.R")

args_output <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Eval_CSPA_MCMI_behavior.png")

#### CSPA####
args_input_path_CSPA <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/CSPA/Eval")
args_label <- c("behavior")

list.files(args_input_path_CSPA, 
           full.names=TRUE, 
           recursive=T, 
           pattern="k_Number_"
           ) -> list_eval_data_CSPA

grep(args_label, list_eval_data_CSPA) %>% 
    list_eval_data_CSPA[.] -> list_cls_method_CSPA
 
grep("ARI", list_cls_method_CSPA) %>% 
  list_cls_method_CSPA[.] -> input_path_ARI_CSPA
ARI_value_CSPA <- numeric()
for(i in 1:length(input_path_ARI_CSPA)){
  load(input_path_ARI_CSPA[i])
  ARI_value_CSPA <- c(ARI_value_CSPA, eval_result)
}

#### MCMIHOOI####
args_input_path_MCMIHOOI <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/Eval")
args_label <- c("behavior")

list.files(args_input_path_MCMIHOOI, 
           full.names=TRUE, 
           recursive=T, 
           pattern="k_Number_"
           ) -> list_eval_data_MCMIHOOI

grep(args_label, list_eval_data_MCMIHOOI) %>% 
    list_eval_data_MCMIHOOI[.] -> list_cls_method_MCMIHOOI
 
grep("ARI", list_cls_method_MCMIHOOI) %>% 
  list_cls_method_MCMIHOOI[.] -> input_path_ARI_MCMIHOOI
ARI_value_MCMIHOOI <- numeric()
for(i in 1:length(input_path_ARI_MCMIHOOI)){
  load(input_path_ARI_MCMIHOOI[i])
  ARI_value_MCMIHOOI <- c(ARI_value_MCMIHOOI, eval_result)
}

#### ggplot####
df_eval <- data.frame(CSPA = ARI_value_CSPA,
                      MCMIHOOI = ARI_value_MCMIHOOI,
                      stringsAsFactors = FALSE
                      )
# クラスタ数は2から始まるため
df_eval$Cluster <- as.numeric(rownames(df_eval)) + 1
# transform long format
df_eval_long <- df_eval %>% 
    pivot_longer(col= -Cluster, 
                 names_to = "method", 
                 values_to ="Eval_Value")
                 
df_eval_long %>% 
        ggplot(., 
               aes(x = Cluster, 
                   y = Eval_Value, 
                   colour = method)
        ) +
        geom_line(size = 3) +
        theme(text = element_text(size = 60)) +
        scale_x_continuous(breaks=seq(2,20,1)) -> gg

df_eval %>% 
  dplyr::select(Cluster,CSPA,MCMIHOOI) %>%
  mutate_if(is.numeric, round, digits = 3) %>% 
  ggtexttable(rows = NULL, theme = ttheme(base_size = 60)) -> gg_eval_table

str_remove(args_output, 
           "output/WTS4/normalize_1/stimAfter/") %>% 
    str_remove(., 
               ".png") -> plot_title

gg <- gg +
    gg_eval_table +
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
       height = 30.0,
       # height = 60.0,
       limitsize = FALSE
       )