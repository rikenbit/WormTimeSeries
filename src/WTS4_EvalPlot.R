source("src/functions_WTS4_EvalPlot.R")

# #### args setting####
# args <- commandArgs(trailingOnly = T)
# # input directory path
# args_input_path <- args[1]
# # output
# args_output <- args[2]
# # ReClustering Method
# args_method <- args[3]


#### test args####
args_input_path <- c("output/WTS4/normalize_1/stimAfter/SBD_abs")
# output
args_output <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/CSPA/EvalPlot.png")
# ReClustering Method
args_method <- c("CSPA")

#### prepare input filepath####
list.files(args_input_path, 
           full.names=TRUE, 
           recursive=T, 
           pattern="eval_result"
           ) -> list_eval_data

grep(args_method, list_eval_data) %>% 
    list_eval_data[.] -> list_cls_method

#### load  value PseudoF####
grep("PseudoF", list_cls_method) %>% 
      list_cls_method[.] -> input_path_PseudoF

PseudoF_value <- numeric()

for(i in 1:length(input_path_PseudoF)){
    load(input_path_PseudoF[i])
    PseudoF_value <- c(PseudoF_value, eval_result)
}
#### load  value Connectivity####
grep("Connectivity", list_cls_method) %>% 
    list_cls_method[.] -> input_path_Connectivity

Connectivity_value <- numeric()

for(i in 1:length(input_path_Connectivity)){
    load(input_path_Connectivity[i])
    Connectivity_value <- c(Connectivity_value, eval_result)
}

#### load  value ARI####
grep("ARI", list_cls_method) %>% 
  list_cls_method[.] -> input_path_ARI

ARI_value <- numeric()

for(i in 1:length(input_path_ARI)){
  load(input_path_ARI[i])
  ARI_value <- c(ARI_value, eval_result)
}
#### load  value purity####
grep("purity", list_cls_method) %>% 
  list_cls_method[.] -> input_path_purity

purity_value <- numeric()

for(i in 1:length(input_path_purity)){
  load(input_path_purity[i])
  purity_value <- c(purity_value, eval_result)
}
#### load  value Fmeasure####
grep("Fmeasure", list_cls_method) %>% 
  list_cls_method[.] -> input_path_Fmeasure

Fmeasure_value <- numeric()

for(i in 1:length(input_path_Fmeasure)){
  load(input_path_Fmeasure[i])
  Fmeasure_value <- c(Fmeasure_value, eval_result)
}
#### load  value Entropy####
grep("Entropy", list_cls_method) %>% 
  list_cls_method[.] -> input_path_Entropy

Entropy_value <- numeric()

for(i in 1:length(input_path_Entropy)){
  load(input_path_Entropy[i])
  Entropy_value <- c(Entropy_value, eval_result)
}

#### eval_cls####
list.files(args_input_path, 
           pattern="_Clusters") %>% 
    str_remove(., "_Clusters") %>% 
        as.numeric() %>% 
            sort() -> eval_cls
#### create dataframe####
df_eval <- data.frame(Cluster = eval_cls,
                      PseudoF = PseudoF_value,
                      Connectivity = Connectivity_value,
                      ARI = ARI_value,
                      purity = purity_value,
                      Fmeasure = Fmeasure_value,
                      Entropy = Entropy_value,
                      stringsAsFactors = FALSE
                      )
# transform long format
df_eval_long <- df_eval %>% 
    pivot_longer(col= -Cluster, 
                 names_to = "Eval", 
                 values_to ="Eval_Value")

#### graph title####
str_remove(args_output, 
           "output/WTS4/normalize_1/stimAfter/") %>% 
    str_remove(., 
               "/EvalPlot.png") -> plot_title
# 参考 https://stats.biopapyrus.jp/r/ggplot/geom_bar.html
#### ggplot PseudoF Connectivity####
df_eval_long %>% 
    dplyr::filter(., Eval=="PseudoF" | Eval=="Connectivity") %>%
        ggplot(., 
               aes(x = Cluster, 
                   y = Eval_Value, 
                   colour = Eval)
        ) +
        geom_line(size = 3) +
        theme(text = element_text(size = 30)) +
        scale_x_continuous(breaks=seq(0,20,2)) -> gg_no_label
#### ggplot label####
df_eval_long %>% 
    dplyr::filter(., Eval=="ARI" | Eval=="purity" | Eval=="Fmeasure") %>%
        ggplot(., 
               aes(x = Cluster, 
                   y = Eval_Value, 
                   colour = Eval)
               ) +
        geom_line(size = 3) +
        theme(text = element_text(size = 30)) +
        scale_x_continuous(breaks=seq(0,20,2)) -> gg_label
#### ggplot Entropy####
df_eval_long %>% 
    dplyr::filter(., Eval=="Entropy") %>%
        ggplot(., 
               aes(x = Cluster, 
                   y = Eval_Value, 
                   colour = Eval)
        ) +
        geom_line(size = 3) +
        theme(text = element_text(size = 30)) +
        scale_x_continuous(breaks=seq(0,20,2)) -> gg_label_Entropy

#### patchwork####
gg <- gg_no_label +
    gg_label +
    gg_label_Entropy +
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
       height = 30.0,
       limitsize = FALSE
       )