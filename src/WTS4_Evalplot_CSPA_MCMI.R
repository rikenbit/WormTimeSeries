source("src/functions_WTS4_Evalplot_CSPA_MCMI.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
args_input <-  args[1]
args_eval_method <- args[2]
args_output <- args[3]
#### test args####
# # args input
# args_input <-c("output/WTS4/normalize_1/stimAfter/SBD_abs")
# # eval method
# args_eval_method <- c("ARI_behavior")
# # output
# args_output <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Evalplot_CSPA_MCMI/ARI_behavior.png")

#### args_k_number####
args_k_number <- as.numeric(2:20)

#### load CSPA####
eval(parse(text=paste0("input_path_CSPA <- c('",args_input,"/CSPA/Eval/",args_eval_method,"')")))
input_path <- c()
for(i in args_k_number){
    eval(parse(text=paste0("path <- c('",input_path_CSPA,"/k_Number_",i,".RData')")))
    input_path <- c(input_path, path)
}
value_CSPA <- c()
for(i in 1:length(input_path)){
    load(input_path[i])
    value_CSPA <- c(value_CSPA, eval_result)
}

#### load MCMI####
eval(parse(text=paste0("input_path_MCMIHOOI <- c('",args_input,"/MCMIHOOI/Eval/",args_eval_method,"')")))
input_path <- c()
for(i in args_k_number){
    eval(parse(text=paste0("path <- c('",input_path_MCMIHOOI,"/k_Number_",i,".RData')")))
    input_path <- c(input_path, path)
}
value_MCMIHOOI <- c()
for(i in 1:length(input_path)){
    load(input_path[i])
    value_MCMIHOOI <- c(value_MCMIHOOI, eval_result)
}

#### df CSPA MCMI####
df_CSPA_MCMI <- data.frame(k_number = as.numeric(args_k_number),
                           CSPA = value_CSPA,
                           MCMIHOOI = value_MCMIHOOI,
                           stringsAsFactors = FALSE
                           )

#### sample number####
eval(parse(text=paste0("args_input_dist <- c('",args_input,"/Distance')")))
input_dist_list <- list.files(args_input_dist, pattern="SampleNumber_", full.names=TRUE)
input_dist_list %>% 
    str_remove(., args_input_dist) %>% 
    str_remove(., "/SampleNumber_") %>% 
    str_remove(., ".RData") %>% 
    as.numeric() %>% 
    sort() -> sample_sort_num

#### load sample####
eval(parse(text=paste0("input_path_sample <- c('",args_input,"/Eval_sample/",args_eval_method,"')")))
input_path <- c()
for(i in args_k_number){
    eval(parse(text=paste0("path <- c('",input_path_sample,"/k_Number_",i,".RData')")))
    input_path <- c(input_path, path)
}

df_sample <- c()
for(i in 1:length(input_path)){
    load(input_path[i])
    df_sample <- rbind(df_sample, eval_result)
}
rownames(df_sample) <- args_k_number
colnames(df_sample) <- sample_sort_num
df_sample %>% 
  as.data.frame() %>% 
    mutate(k_number=args_k_number) -> df_sample

#### df CSPA MCMI sample####
df_all <- merge(df_CSPA_MCMI,
                df_sample,
                by.x = "k_number", 
                by.y = "k_number", 
                all.x = TRUE)

#### transform long format####
df_all %>% 
    pivot_longer(col= -k_number, 
                 names_to = "DataName", 
                 values_to ="Eval_Value") -> df_all_long

df_all_long %>% 
    dplyr::filter(., DataName=="CSPA"| DataName=="MCMIHOOI") -> df_CSPA_MCMI
df_all_long %>%
    dplyr::filter(., DataName != "CSPA") %>%
        dplyr::filter(., DataName != "MCMIHOOI") -> df_sample

#### ggplot####
gg <- ggplot()
# line plot
gg <- gg + geom_line(data = df_CSPA_MCMI, 
                     aes(x = k_number, y = Eval_Value, colour = DataName),
                     size = 3)
# point plot
gg <- gg + geom_point(data = df_all_long, 
                      aes(x = k_number, y = Eval_Value, colour = DataName),
                      size = 6.0)
# all text size
gg <- gg + theme(text = element_text(size = 40)) 

#### ggplot scale_color_manual####
#参考 色見本 http://www.okadajp.org/RWiki/?色見本
gg <- gg + scale_color_manual(values = c("CSPA"="black",
                                         "MCMIHOOI"="red",
                                         "1"="orange",
                                         "2"="slategray",
                                         "4"="lightblue",
                                         "5"="slateblue",
                                         "6"="olivedrab",
                                         "7"="navyblue",
                                         "9"="darkorchid",
                                         "10"="sienna",
                                         "11"="green",
                                         "12"="darkgoldenrod",
                                         "13"="cyan",
                                         "14"="hotpink",
                                         "15"="royalblue",
                                         "16"="lightcoral",
                                         "17"="yellow",
                                         "18"="magenta",
                                         "19"="goldenrod",
                                         "21"="pink",
                                         "22"="brown",
                                         "23"="darkslateblue",
                                         "24"="blue",
                                         "26"="deepskyblue",
                                         "27"="aquamarine",
                                         "28"="lightgray")
                              )

#### graph title####
# str_remove(args_output, "output/WTS4/normalize_1/stimAfter/") %>% 
#   str_remove(., ".png") -> plot_title
args_eval_method -> plot_title
# 参考 https://stats.biopapyrus.jp/r/ggplot/geom_bar.html
#### patchwork####
gg <- gg +
    plot_annotation(title = plot_title,
                    caption = 'made with patchwork',
                    theme = theme(plot.title = element_text(size = 60, hjust = 0.5))
                    )

#### ggsave####
ggsave(filename = args_output, 
       plot = gg,
       dpi = 100, 
       width = 40.0, 
       height = 30.0,
       limitsize = FALSE
)