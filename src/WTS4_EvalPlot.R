source("src/functions_WTS4_EvalPlot.R")

#### args setting####
# input directory path
args_input_path <- args[1]
# output
args_output <- args[2]

# ReClustering Method
args_method <- args[3]
# Evaluation Method
args_eval_method <- args[4]

# #### test args####
# args_input_path <- c("output/WTS4/normalize_1/stimAfter/SBD_abs")
# # output
# args_output <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/CSPA/PseudoF/EvalPlot.png")
# # ReClustering Method
# args_method <- c("CSPA")
# # Evaluation Method
# args_eval_method <- c("PseudoF")


#### prepare input filepath####
list_eval_data <- list.files(args_input_path, 
                             full.names=TRUE, 
                             recursive=T, 
                             pattern="eval_result")
grep(args_method, list_eval_data) %>% 
    list_eval_data[.] -> list_cls_method

grep(args_eval_method, list_cls_method) %>% 
  	list_cls_method[.] -> input_path_list

#### load####
# 評価値を格納する空ファイルを作成
eval_value <- numeric()

for(i in 1:length(input_path_list)){
    load(input_path_list[i])
    eval_value <- c(eval_value, eval_result)
    }

list.files(args_input_path, 
           pattern="_Clusters") %>% 
    str_remove(., "_Clusters") %>% 
        as.numeric() %>% 
            sort() -> eval_cls
#### create dataframe####
df_eval <- data.frame(Cluster = eval_cls,
                      # Evaluation_Value = trunc(eval_value),
                      Evaluation_Value = eval_value,
                      stringsAsFactors = FALSE
                      )

#### ggplot####
str_remove(args_output, 
           "output/WTS4/normalize_1/stimAfter/") %>% 
  	str_remove(., 
  	           "/EvalPlot.png") -> plot_title
# 参考 https://stats.biopapyrus.jp/r/ggplot/geom_bar.html
gg <- ggplot(df_eval, 
            aes(x = Cluster, 
                y = Evaluation_Value, 
                fill = factor(Cluster)
                # fill = Cluster
                )
            ) + 
    geom_bar(stat = "identity") +
    theme(text = element_text(size = 30))+
    plot_annotation(title = plot_title,
                    caption = 'made with patchwork',
                    theme = theme(plot.title = element_text(size = 40, 
                                                            hjust = 0.5)
                                  )
                    )

#### ggsave####
ggsave(filename = args_output, 
       plot = gg,
       dpi = 100, 
       width = 20.0, 
       height = 20.0,
       limitsize = FALSE
       )