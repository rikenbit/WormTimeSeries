source("src/functions_WTS4_Evalplot_Dist_Merge.R")

### test args####
# args input EUCL
args_input_EUCL <-c("output/WTS4/normalize_1/stimAfter/EUCL")
# args input SBD_abs
args_input_SBD_abs <-c("output/WTS4/normalize_1/stimAfter/SBD_abs")
# eval method
args_eval_method <- c("Silhouette")
# output
args_output <- c("output/WTS4/normalize_1/stimAfter/Evalplot_Dist_Merge/Silhouette.png")

#### args_k_number####
args_k_number <- as.numeric(2:20)

#### load CSPA EUCL####
eval(parse(text=paste0("input_path_CSPA <- c('",args_input_EUCL,"/CSPA/Eval/",args_eval_method,"')")))

input_path <- c()
for(i in args_k_number){
    eval(parse(text=paste0("path <- c('",input_path_CSPA,"/k_Number_",i,".RData')")))
    input_path <- c(input_path, path)
    }

value_CSPA_EUCL <- c()
for(i in 1:length(input_path)){
    load(input_path[i])
    value_CSPA_EUCL <- c(value_CSPA_EUCL, eval_result)
    }

#### load CSPA SBD_abs####
eval(parse(text=paste0("input_path_CSPA <- c('",args_input_SBD_abs,"/CSPA/Eval/",args_eval_method,"')")))

input_path <- c()
for(i in args_k_number){
    eval(parse(text=paste0("path <- c('",input_path_CSPA,"/k_Number_",i,".RData')")))
    input_path <- c(input_path, path)
    }

value_CSPA_SBD_abs <- c()
for(i in 1:length(input_path)){
    load(input_path[i])
    value_CSPA_SBD_abs <- c(value_CSPA_SBD_abs, eval_result)
    }
#### load MCMIHOOI EUCL####
eval(parse(text=paste0("input_path_MCMIHOOI <- c('",args_input_EUCL,"/MCMIHOOI/Eval/",args_eval_method,"')")))

input_path <- c()
for(i in args_k_number){
    eval(parse(text=paste0("path <- c('",input_path_MCMIHOOI,"/k_Number_",i,".RData')")))
    input_path <- c(input_path, path)
    }

value_MCMIHOOI_EUCL <- c()
for(i in 1:length(input_path)){
    load(input_path[i])
    value_MCMIHOOI_EUCL <- c(value_MCMIHOOI_EUCL, eval_result)
    }

#### load MCMIHOOI SBD_abs####
eval(parse(text=paste0("input_path_MCMIHOOI <- c('",args_input_SBD_abs,"/MCMIHOOI/Eval/",args_eval_method,"')")))

input_path <- c()
for(i in args_k_number){
    eval(parse(text=paste0("path <- c('",input_path_MCMIHOOI,"/k_Number_",i,".RData')")))
    input_path <- c(input_path, path)
    }

value_MCMIHOOI_SBD_abs <- c()
for(i in 1:length(input_path)){
    load(input_path[i])
    value_MCMIHOOI_SBD_abs <- c(value_MCMIHOOI_SBD_abs, eval_result)
}

#### dataframe####
df_CSPA_MCMI <- data.frame(k_number = as.numeric(args_k_number),
                           EUCL_CSPA = value_CSPA_EUCL,
                           EUCL_MCMIHOOI = value_MCMIHOOI_EUCL,
                           mSBD_CSPA = value_CSPA_SBD_abs,
                           mSBD_MCMIHOOI = value_MCMIHOOI_SBD_abs,
                           stringsAsFactors = FALSE
                           )

#### transform long format####
df_CSPA_MCMI %>% 
  pivot_longer(col= -k_number, 
               names_to = "DataName", 
               values_to ="Eval_Value") -> df_long

#### ggplot####
if (args_eval_method == "Silhouette") {
    label_y <- c("Average silhouette coefficient")
} else {
    label_y <-  c("Eval_Value")
}

gg <- ggplot(df_long, 
             aes(x=k_number,
                 y=Eval_Value,
                 colour = DataName)
             ) +
  geom_line(size = 3) +
  geom_point(size = 6) +
  theme(text = element_text(size = 120)) +
  labs(x = "Number of clusters",
       y = label_y,
       color = "Clustering methods") +
  scale_color_hue(labels = c(EUCL_CSPA = "CSPA (Euclid)", 
                             EUCL_MCMIHOOI = "MC-MI-HOOI (Euclid)", 
                             mSBD_CSPA = "CSPA (mSBD)", 
                             mSBD_MCMIHOOI = "WormTensor"))

#### ggsave####
ggsave(filename = args_output, 
       plot = gg,
       dpi = 100, 
       width = 50.0, 
       height = 30.0,
       limitsize = FALSE)