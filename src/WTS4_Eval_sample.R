source("src/functions_WTS4_Eval_sample.R")

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
# args_input_cls <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Cluster_sample/k_Number_5/sample_cls.RData")
# # output ggplot
# args_output <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/DimReduc_sample/k_Number_5/Eval_sample.png")
# # params merged_cls
# args_input_merged_cls <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/Merged_cls/k_Number_5.RData")
# # input merged_distance
# args_input_MCMIHOOI <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/Merged_data/k_Number_5.RData")
# # params dist path
# args_input_path<- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Distance")

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

unlist(lapply(df_cls_label, function(x) {
    clusters <- x$Clusters
    classes <- x$Classes
    ClusterPurity(clusters, classes)
    }
    )
) -> purity_value

unlist(lapply(df_cls_label, function(x) {
    clusters <- x$Clusters
    classes <- x$Classes
    Fmeasure(clusters, classes)
    }
    )
) -> Fmeasure_value

unlist(lapply(df_cls_label, function(x) {
    clusters <- x$Clusters
    classes <- x$Classes
    Entropy(clusters, classes)
    }
    )
) -> Entropy_value


#### fix sample number sort####
input_path_list <- list.files(args_input_path, pattern="SampleNumber_", full.names=TRUE)
input_path_list %>% 
    str_remove(., args_input_path) %>% 
    str_remove(., "/SampleNumber_") %>% 
    str_remove(., ".RData") %>% 
    as.numeric() %>% 
    sort() -> sample_sort_num
########

#### create dataframe####
df_eval <- data.frame(SampleNumber = sample_sort_num,
                      ARI = ARI_value,
                      purity = purity_value,
                      Fmeasure = Fmeasure_value,
                      Entropy = Entropy_value,
                      # Entropyは小さいほど、良い。他の評価は値が高いほどお良い。
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
    # rownames_to_column("SampleNumber") %>%
        mutate(weight_abs =abs(weight)) %>% 
            dplyr::arrange(desc(weight_abs)) -> df_weight

#### ggplot ARI####
g1 <- ggplot(df_eval_wide, aes(x = SampleNumber, y= ARI , group=1))
g1 <- g1 + geom_line(color = "red", size= 2)
g1 <- g1 + scale_x_discrete(limits=df_weight$SampleNumber)
g1 <- g1 + theme_half_open()
g1 <- g1 + theme(text = element_text(size = 24))
g1 <- g1 + theme(axis.title.y=element_text(colour = "red",size = 24))
g1 <- g1 + geom_smooth(method="lm", size =0.5, se = TRUE, alpha = 0.4, color = "red")
g1 <- g1 + theme(legend.position = 'none')


g2 <- ggplot(df_eval_wide, aes(x = SampleNumber, y= annotated_count, group=1))
g2 <- g2 + geom_line(color = "black", size= 2)
g2 <- g2 + scale_x_discrete(limits=df_weight$SampleNumber)
g2 <- g2 + scale_y_continuous(position = "right")
g2 <- g2 + geom_smooth(method="lm", size =0.5, se = TRUE, alpha = 0.4, color = "black")
g2 <- g2 + theme_half_open()
# theme_half_openで軸の書式(色・文字サイズ)がリセットされる
g2 <- g2 + theme(text = element_text(size = 24))
g2 <- g2 + theme(axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank())

aligned_plots_result <- cowplot::align_plots(g1, g2, align="hv", axis="tblr")
gg_label_ARI <- cowplot::ggdraw(aligned_plots_result[[1]]) + cowplot::draw_plot(aligned_plots_result[[2]])

#### ggplot purity####
g1 <- ggplot(df_eval_wide, aes(x = SampleNumber, y= purity , group=1))
g1 <- g1 + geom_line(color = "red", size= 2)
g1 <- g1 + scale_x_discrete(limits=df_weight$SampleNumber)
g1 <- g1 + theme_half_open()
g1 <- g1 + theme(text = element_text(size = 24))
g1 <- g1 + theme(axis.title.y=element_text(colour = "red",size = 24))
g1 <- g1 + geom_smooth(method="lm", size =0.5, se = TRUE, alpha = 0.4, color = "red")
g1 <- g1 + theme(legend.position = 'none')


g2 <- ggplot(df_eval_wide, aes(x = SampleNumber, y= annotated_count, group=1))
g2 <- g2 + geom_line(color = "black", size= 2)
g2 <- g2 + scale_x_discrete(limits=df_weight$SampleNumber)
g2 <- g2 + scale_y_continuous(position = "right")
g2 <- g2 + geom_smooth(method="lm", size =0.5, se = TRUE, alpha = 0.4, color = "black")
g2 <- g2 + theme_half_open()
# theme_half_openで軸の書式(色・文字サイズ)がリセットされる
g2 <- g2 + theme(text = element_text(size = 24))
g2 <- g2 + theme(axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank())

aligned_plots_result <- cowplot::align_plots(g1, g2, align="hv", axis="tblr")
gg_label_purity <- cowplot::ggdraw(aligned_plots_result[[1]]) + cowplot::draw_plot(aligned_plots_result[[2]])

#### ggplot Fmeasure####
g1 <- ggplot(df_eval_wide, aes(x = SampleNumber, y= Fmeasure , group=1))
g1 <- g1 + geom_line(color = "red", size= 2)
g1 <- g1 + scale_x_discrete(limits=df_weight$SampleNumber)
g1 <- g1 + theme_half_open()
g1 <- g1 + theme(text = element_text(size = 24))
g1 <- g1 + theme(axis.title.y=element_text(colour = "red",size = 24))
g1 <- g1 + geom_smooth(method="lm", size =0.5, se = TRUE, alpha = 0.4, color = "red")
g1 <- g1 + theme(legend.position = 'none')


g2 <- ggplot(df_eval_wide, aes(x = SampleNumber, y= annotated_count, group=1))
g2 <- g2 + geom_line(color = "black", size= 2)
g2 <- g2 + scale_x_discrete(limits=df_weight$SampleNumber)
g2 <- g2 + scale_y_continuous(position = "right")
g2 <- g2 + geom_smooth(method="lm", size =0.5, se = TRUE, alpha = 0.4, color = "black")
g2 <- g2 + theme_half_open()
# theme_half_openで軸の書式(色・文字サイズ)がリセットされる
g2 <- g2 + theme(text = element_text(size = 24))
g2 <- g2 + theme(axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank())

aligned_plots_result <- cowplot::align_plots(g1, g2, align="hv", axis="tblr")
gg_label_Fmeasure <- cowplot::ggdraw(aligned_plots_result[[1]]) + cowplot::draw_plot(aligned_plots_result[[2]])

#### ggplot Entropy####
g1 <- ggplot(df_eval_wide, aes(x = SampleNumber, y= Entropy , group=1))
g1 <- g1 + geom_line(color = "red", size= 2)
g1 <- g1 + scale_x_discrete(limits=df_weight$SampleNumber)
g1 <- g1 + theme_half_open()
g1 <- g1 + theme(text = element_text(size = 24))
g1 <- g1 + theme(axis.title.y=element_text(colour = "red",size = 24))
g1 <- g1 + geom_smooth(method="lm", size =0.5, se = TRUE, alpha = 0.4, color = "red")
g1 <- g1 + theme(legend.position = 'none')


g2 <- ggplot(df_eval_wide, aes(x = SampleNumber, y= annotated_count, group=1))
g2 <- g2 + geom_line(color = "black", size= 2)
g2 <- g2 + scale_x_discrete(limits=df_weight$SampleNumber)
g2 <- g2 + scale_y_continuous(position = "right")
g2 <- g2 + geom_smooth(method="lm", size =0.5, se = TRUE, alpha = 0.4, color = "black")
g2 <- g2 + theme_half_open()
# theme_half_openで軸の書式(色・文字サイズ)がリセットされる
g2 <- g2 + theme(text = element_text(size = 24))
g2 <- g2 + theme(axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank())

aligned_plots_result <- cowplot::align_plots(g1, g2, align="hv", axis="tblr")
gg_label_Entropy <- cowplot::ggdraw(aligned_plots_result[[1]]) + cowplot::draw_plot(aligned_plots_result[[2]])

#### ggtexttable####
# 各Evalのmaxかminの行番号
df_eval_long %>% 
    mutate(num = row_number()) -> df_eval_long_ID
c("Entropy") %>% 
    purrr::map_int(., eval_min) -> eval_id_min
c("ARI","purity","Fmeasure") %>% 
    purrr::map_int(., eval_max) -> eval_id_max
c(eval_id_max, eval_id_min) -> eval_id

df_eval_long_ID[eval_id,] %>% 
    dplyr::select(Eval,SampleNumber,Eval_Value) %>%
        mutate_if(is.numeric, round, digits = 3) %>% 
            ggtexttable(rows = NULL, theme = ttheme(base_size = 60)) -> gg_eval_table

#### graph title####
str_remove(args_output, 
           "output/WTS4/normalize_1/stimAfter/") %>% 
    str_remove(., 
               ".png") -> plot_title
#### patchwork####
gg <- gg_label_ARI +
    gg_label_purity +
    gg_label_Fmeasure +
    gg_label_Entropy +
    gg_eval_table +
    plot_layout(ncol = 1) +
    plot_annotation(title = plot_title,
                    caption = 'made with patchwork',
                    theme = theme(plot.title = element_text(size = 40, hjust = 0.5))
                    )

#### ggsave####
ggsave(filename = args_output, 
       plot = gg,
       dpi = 100, 
       width = 30.0,
       height = 60.0,
       limitsize = FALSE
       )