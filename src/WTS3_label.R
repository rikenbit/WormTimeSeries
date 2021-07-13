source("src/functions_WTS3_label.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
# select animal number 個体番号の指定
args_sample <- args[1]
# input 距離データ
args_dist <- args[2]
# input igraph
args_igraph <- args[3]
# input PeriodicACF
args_periodic <- args[4]
# クラスター評価手法
args_eval <- args[5]
# ラベリングするクラスタに含まれる細胞
args_shift <- args[6]
# output
args_output_label <- args[7]
# output
args_output_cutree <- args[8]

# #### test args####
# # select animal number 個体番号の指定
# args_sample <- c("1")
# # input 距離データ
# args_dist <- c("output/WTS3/normalize_1/stimAfter/EUCL/SampleNumber_1/EUCL.RData")
# # input igraph
# args_igraph <- c("data/igraph/Fig1_HNS.RData")
# # input PeriodicACF
# args_periodic <- c("output/WTS2/WTS2_PeriodicACF.csv")
# # クラスター評価手法
# args_eval <- c("ARI")
# # ラベリングするクラスタに含まれる細胞
# args_shift <- c("ASER")
# # output
# args_output_label <- c("output/WTS3/normalize_1/stimAfter/EUCL/ARI/SampleNumber_1/label_table.RData")
# # output
# args_output_cutree <- c("output/WTS3/normalize_1/stimAfter/EUCL/ARI/SampleNumber_1/cutree_table.RData")

#### load SBD####
load(args_dist)
data.frame(
    cell_type = attr(d, "Labels"),
    stringsAsFactors = FALSE
    ) -> df_cell_type
#### load ACF label####
periodic_sheet <- read.csv(args_periodic, 
                           colClasses=c("numeric", 
                                        "character", 
                                        rep("numeric",2)))
periodic_sheet %>% 
    filter(.,
           stim == 1,
           sample_number == args_sample
           ) %>% 
        select(cell_type, acf) -> stim_sheet
#### add ACF####
df_label <- merge(df_cell_type, 
                  stim_sheet, 
                  by.x = "cell_type", 
                  by.y = "cell_type", 
                  all.x = TRUE)
#### create ACF label####
df_label %>% 
    mutate(., label_acf = if_else(is.na(acf),
                                 true = 0, 
                                 false = 1)
           ) -> df_label
#### clustering####
set_cutree <- seq(3,10)
set_cutree %>% 
  purrr::map(., .cls_ward) -> df_cls
#### eval####
# select method
eval_type <- switch(args_eval,
                    "purity" = .cls_purity,
                    "ARI" = .cls_ARI,
                    "Fmeasure" = .cls_Fmeasure,
                    "Entropy" = .cls_Entropy,
                    stop("Only can use cls_purity, ARI, Fmeasure, Entropy")
                    )
# evaluation
set_cutree %>% 
    purrr::map_dbl(., eval_type) -> eval_value
df_cutree <- data.frame(set_cutree = set_cutree, 
                        eval_value = eval_value
                        )

#### select cutree####
result_cutree <- switch(args_eval,
                        "purity" = .max_eval(df_cutree),
                        "ARI" = .max_eval(df_cutree),
                        "Fmeasure" = .max_eval(df_cutree),
                        "Entropy" = .min_eval(df_cutree),
                        stop("Only can use cls_purity,ARI,Fmeasure,Entropy")
                        )
#### create cutree table####
df_cutree %>% 
    mutate(., select_cutree = if_else(set_cutree == result_cutree,
                                            true = 1, 
                                            false = 0)
           ) -> df_cutree
#### save cutree table####
save(df_cutree, file=args_output_cutree)
#### create cls label####
# merge table
df_label <- merge(df_label, 
                  df_cls[[result_cutree - 2]], 
                  by.x = "cell_type", 
                  by.y = "cell_type", 
                  all.x = TRUE)
# select ASER cls number
df_label %>% 
    filter(., cell_type==args_shift) %>% 
        .$cls -> same_cls
# add label acf
df_label %>% 
    mutate(., label_cls = if_else(cls == same_cls,
                                          true = 1, 
                                          false = 0)
           ) -> df_label
#### load NeuronType &Group####
load(args_igraph)
# node information convert dataframe
ig_Fig1_HNS %>% 
    igraph::as_data_frame(., what="vertices") -> df_node
# remove space from colnames
df_node %>% 
    names() %>% 
        str_replace_all(., c(" " = "")) -> names(df_node)
data.frame(
    cell_type = df_node$name,
    NeuronType = df_node$NeuronType,
    NeuronGroup = df_node$NeuronGroup,
    stringsAsFactors = FALSE
    ) -> df_NeuronType
#### add NeuronType &Group####
df_label <- merge(df_label, 
                  df_NeuronType, 
                  by.x = "cell_type", 
                  by.y = "cell_type", 
                  all.x = TRUE)
#### save label table####
df_label %>% dplyr::select(cell_type,
                           label_acf,
                           label_cls,
                           acf,
                           cls,
                           NeuronType,
                           NeuronGroup
                           ) -> df_label
save(df_label, file=args_output_label)