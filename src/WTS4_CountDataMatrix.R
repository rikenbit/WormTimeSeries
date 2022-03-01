library(tidyverse)

sample_path <- list.files("output/WTS4/n1_28sample/stimAfter/SBD_abs/Distance", pattern="SampleNumber_", full.names=TRUE)
sample_path %>% 
    str_remove(., "output/WTS4/n1_28sample/stimAfter/SBD_abs/Distance/SampleNumber_") %>% 
    str_remove(., ".RData") %>% 
    as.numeric() %>% 
    sort() -> sample_sort_num

count_cell_n <- c()
for(i in sample_sort_num){
    eval(parse(text=paste0("load('output/WTS4/n1_28sample/stimAfter/SBD_abs/Distance/SampleNumber_",i,".RData')")))
    count_cell_n <- c(count_cell_n, length(attr(d, "Labels")))
}

df_CountDataMatrix <- data.frame(SampleNumber = sample_sort_num,
								 CountCell = count_cell_n,
                           		 stringsAsFactors = FALSE
                           		)

write.csv(df_CountDataMatrix, "data/CountDataMatrix.csv")