#library
##################################################
library(tidyverse)
library(patchwork)
library(ggpubr)
##################################################
# 適切なクラスタ数の関数
eval_min = function(x) {
  df_eval_long_ID %>% 
    filter(Eval==x) -> test_df
  return_object <- test_df$num[which.min(test_df$Eval_Value)]
  return(return_object)
}
eval_max = function(x) {
  df_eval_long_ID %>% 
    filter(Eval==x) -> test_df
  return_object <- test_df$num[which.max(test_df$Eval_Value)]
  return(return_object)
}