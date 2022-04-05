#library
##################################################
# library(beeswarm)
library(tidyverse)
##################################################
.all_k_table = function(x) {
    eval(parse(text=paste0("args_input_weight <- c('",args_weight_path,"/k_Number_",x,".RData')")))
    load(args_input_weight)
    # weigh table
    data.frame(SampleNumber = as.character(sample_sort_num),
               Cluster = as.character(c(rep(x))),
               weight = merged_data$W,
               stringsAsFactors = FALSE) %>% 
        mutate(weight_abs =abs(weight)) -> return_object
    return(return_object)
}

# .list_beeswarm = function(x) {
#     df_weight_all %>% 
#         dplyr::filter(SampleNumber == x) %>% 
#         dplyr::arrange(as.numeric(Cluster)) %>% 
#         .$weight_abs -> return_object
#     return(return_object)
# }