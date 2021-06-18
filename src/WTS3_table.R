source("src/functions_WTS3_new.R")
#サンプル番号一覧取得
args_numbers <- c("1","2","4","5","6","7","9","10","11","12","13","14","15","16","17","18","19","21","22","23","24","26","27","28")
args_tempdata_dir <- c("output/WTS3/normalize_1/SBD/ARI/tsne/cls_tempdata/SampleNumber_")
check_args_shift = function(x) {
  x -> sample_cell_type
  sample_cell_type %>% 
    str_count(., pattern="ASER") %>%
    sum() -> check_ASER
  sample_cell_type %>% 
    str_count(., pattern="BAGR") %>%
    sum() -> check_BAGR
  sample_cell_type %>% 
    str_count(., pattern="BAGL") %>%
    sum() -> check_BAGL
  if (check_ASER >= 1) {
    args_shift <- "ASER"
  } else if (check_BAGR >= 1) {
    args_shift <- "BAGR"
  } else if (check_BAGL >= 1) {
    args_shift <- "BAGL"
  } else {
    args_shift <- sample_cell_type[1]
  }
  return(args_shift)
}
# filter same_clusters
filter_clusters = function(x,y) {
  x %>% 
    filter(., cell_type == y) %>% 
    .$cls %>% 
    unique() -> same_cls
  x %>% 
    filter(., cls == same_cls) -> df_filter
  return(df_filter)
}

# df_temp_list <- list()

# サンプル番号を読み込んで，loadする関数
load_tempdata = function(x) {
    args_number <- args_numbers[x]
    args_tempdata <- paste(args_tempdata_dir,args_number, sep = "")
    eval(parse(text=paste0("tempdata <- c('",args_tempdata,".RData')")))
    load(tempdata)
    df_tempdata$sample_number <- rep(args_number,nrow(df_tempdata))
    df_tempdata$cell_type %>% 
      check_args_shift() -> args_shift
    filter_clusters(df_tempdata,args_shift) -> cls_n
    cls_n$cls %>% unique() -> cls_num
    df_tempdata$cls_number <- rep(cls_num,nrow(df_tempdata))
    df <- df_tempdata
    return(df)
}
# 24サンプルのdfのリストを作成
seq(1,length(args_numbers)) %>% 
    purrr::map(.,load_tempdata) -> df_temp_list
# 24のデータフレームをくっつける sample_numberの列作成
plyr::join_all(df_temp_list, type = 'full') -> df_table
#### df_stim_table#####
# フィルター stim = 1になっている行のみにし，cell_type列を取得&重複削除
df_table %>% 
    filter(stim==1) %>%
        select(cell_type) %>% 
            .$cell_type %>% 
                unique() %>% 
                    sort() -> stim_cell_type
# フィルター 取得したcell_type列のみの行にする，sample_number，細胞名，stimのみ残す

df_table %>% 
    filter(cell_type %in% stim_cell_type) %>%
        select(sample_number, cell_type, stim, NeuronType) -> df_stim_table

df_stim_table %>% 
  filter(NeuronType != "Sensory" ) %>%
  group_by(cell_type) %>%
  summarise(sample_sum = sum(stim)) -> stim_tibble_notSensory

# pivot_widerにする
df_stim_table %>%
  # 数字の細胞除去
  filter(NeuronType %in% c("Sensory","Interneuron","Endorgan","Motorneuron")) %>%
    pivot_wider(id_cols = sample_number,
                names_from = cell_type,
                values_from = stim) -> stim_table
df_stim_table %>%
  # 数字の細胞除去
  filter(NeuronType %in% c("Interneuron","Endorgan","Motorneuron")) %>%
  pivot_wider(id_cols = sample_number,
              names_from = cell_type,
              values_from = stim) -> stim_tibble_notSensory
# エクスポート
write_excel_csv(stim_table, "output/WTS3/normalize_1/SBD/ARI/tsne/cls_tempdata/stim_table.csv")
write_excel_csv(stim_tibble_notSensory, "output/WTS3/normalize_1/SBD/ARI/tsne/cls_tempdata/stim_tibble_notSensory.csv")
#########

#### df_cls_table#####
df_table %>% 
  filter(cls==cls_number) %>%
  select(sample_number, cell_type, cls, NeuronType) -> df_cls_table
df_cls_table$cls <- 1

df_cls_table %>% 
  # 数字の細胞除去
  filter(NeuronType %in% c("Sensory","Interneuron","Endorgan","Motorneuron")) %>%
  pivot_wider(id_cols = sample_number,
              names_from = cell_type,
              values_from = cls) -> cls_table
cls_table[is.na(cls_table)] <- 0

df_cls_table %>% 
  # 数字の細胞除去
  filter(NeuronType %in% c("Interneuron","Endorgan","Motorneuron")) %>%
  pivot_wider(id_cols = sample_number,
              names_from = cell_type,
              values_from = cls) -> cls_table_notSensory
cls_table_notSensory[is.na(cls_table_notSensory)] <- 0
# エクスポート
write_excel_csv(cls_table, "output/WTS3/normalize_1/SBD/ARI/tsne/cls_tempdata/cls_table.csv") 

write_excel_csv(cls_table_notSensory, "output/WTS3/normalize_1/SBD/ARI/tsne/cls_tempdata/cls_table_notSensory.csv") 