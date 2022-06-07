source("src/functions_WTS4_dhyper.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
args_input_csv <- args[1]
args_params_csv <- args[2]
args_output_csv <- args[3]
args_output_q_csv <- args[4]

# #### test args####
# args_input_csv <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/Merged_tsne/label_table_k6.csv")
# args_params_csv <- c("data/label_table_ann.csv")
# args_output_csv <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/Merged_tsne/dhyper_table_k6.csv")
# args_output_q_csv <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/Merged_tsne/dhyper_q_table_k6.csv")

#### load####
# clustering results
label_table_cls <- read.csv(args_input_csv)
# toyoshima's annotation
label_table_ann <- read.csv(args_params_csv)
# merge table
table_merged <- merge(label_table_cls, 
                      label_table_ann, 
                      by.x = "Cell_type", 
                      by.y = "Cell_type", 
                      all.x = TRUE
                      )

#### dhyper####
colnames(label_table_ann) %>% 
    .[3:length(.)] %>% 
        purrr::map(., dhyper_all_cls) %>% 
            bind_cols() -> table_dhyper
colnames(label_table_ann) %>% 
    .[3:length(.)] -> colnames(table_dhyper)
-log10(table_dhyper) -> table_dhyper_LOG

#### ggsave####
write.csv(table_dhyper_LOG, 
          args_output_csv, 
          row.names=FALSE)

#### dhyper q value####
1:length(table_dhyper) %>% 
    purrr::map(., 
               function(x) {p.adjust(table_dhyper[[x]], "BH")}
               ) -> list_dhyper_q
names(list_dhyper_q) <- colnames(table_dhyper)
list_dhyper_q %>% 
    as.data.frame() -> table_dhyper_q

#### ggsave q value####
write.csv(table_dhyper_q, 
          args_output_q_csv, 
          row.names=FALSE)