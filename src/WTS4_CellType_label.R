library(tidyverse)

# #### test args####
args_input_MCMIHOOI <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/Merged_data/k_Number_9.RData")
args_input_path <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Distance")

#### load dist data####
# 空の行列を格納するファイルを作成
D <- list()
# inputファイル名のリスト
input_path_list <- list.files(args_input_path, pattern="SampleNumber_", full.names=TRUE)
#### fix sample number sort####
input_path_list %>% 
  str_remove(., args_input_path) %>% 
  str_remove(., "/SampleNumber_") %>% 
  str_remove(., ".RData") %>% 
  as.numeric() %>% 
  sort() -> sample_sort_num
input_path <- c()
for(i in sample_sort_num){
  eval(parse(text=paste0("path <- c('",args_input_path,"/SampleNumber_",i,".RData')")))
  input_path <- c(input_path, path)
}
input_path_list <- input_path

# ファイルを読み込んで，リストに加える．各リストのattr(*, "Labels")に細胞型名が残っている
for(i in 1:length(input_path_list)){
    load(input_path_list[i])
    D <- c(D, list(d))
}

#### name count####
#文字の細胞名
cn_list <- list()
#数字の列名
digit_list <- list()
for(i in 1:length(D)){
    D[[i]] %>% attr(., "Labels") -> tmp
    cn_list[[i]] <- tmp[grep("^[0-9]", tmp, invert=TRUE)]
    digit_list[[i]] <- tmp[grep("^[0-9]", tmp, invert=FALSE)]
}
#文字列の列数のカウント
cn_num <- unlist(lapply(cn_list, function(x){length(x)}))
#数字の列数のカウント
digit_num <- unlist(lapply(digit_list, function(x){length(x)}))

#### not Annotated####
# 和集合ベクトル
load("output/WTS4/normalize_1/stimAfter/SBD_abs/Membership/k_Number_9.RData")
newHs[[1]] %>%
    rownames() -> all_annotated_name
# not_annotated数
not_annotated_count <- unlist(lapply(cn_list, function(x){length(setdiff(all_annotated_name, x))}))

#### name count dataframe####
df_count <- data.frame(
    SampleNumber = as.character(sample_sort_num),
    annotated = cn_num,
    not_annotated = not_annotated_count,
    digit = digit_num
    )
df_count %>%
    pivot_longer(cols = c(annotated, not_annotated, digit),
                 names_to = "CellType",
                 values_to = "count") -> data
data$CellType <- factor(data$CellType,
                        levels =c("digit",
                                  "not_annotated",
                                  "annotated"))
### MCMI weight table####
# merged_data
load(args_input_MCMIHOOI)
data.frame(weight = merged_data$W,
           stringsAsFactors = FALSE) %>%
  mutate(SampleNumber = as.character(sample_sort_num)) %>%
  mutate(weight_abs = abs(weight)) %>%
  dplyr::arrange(desc(weight_abs)) -> df_weight

#### bar plot####
# sort MCMI$Weight
p <- ggplot(data, aes(x=SampleNumber, y=count)) +
    geom_bar(stat="identity", aes(fill=CellType)) +
    scale_x_discrete(limits=df_weight$SampleNumber)
ggsave(filename = "output/WTS4/normalize_1/stimAfter/SBD_abs/CellType_label/MCMI_weightsort_k9.png",
       plot = p)