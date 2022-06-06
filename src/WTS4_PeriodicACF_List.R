source("src/functions_WTS4_PeriodicACF_List.R")
#### args setting####
args <- commandArgs(trailingOnly = T)
args_periodic <- args[1]
args_N_observation <- args[2]
args_output <- args[3]
args_output_sample <- args[4]

#### test args####
# args_periodic <- c("output/WTS2/WTS2_PeriodicACF.csv")
# args_N_observation <- c("5")
# args_output <- c("output/WTS4/normalize_1/PeriodicACF/periodic_acf_over5.csv")
# args_output_sample <- c("output/WTS4/normalize_1/PeriodicACF/sample_over5.csv")

#### No. of Clusters####
args_N_obs <- as.numeric(args_N_observation)

#### load WTS2_PeriodicACF.csv####
read.csv(args_periodic, 
         colClasses=c("numeric", 
                      "character", 
                      rep("numeric",2)
                      )
         ) %>% 
    dplyr::filter(stim == 1) -> ACF_table

#### Counting duplicates####
sapply(unique(ACF_table$cell_type), 
       function(x){sum(ACF_table$cell_type==x)}
       ) -> ACF_list

#### selection over 5sample####
data.frame(
    cell_type = attr(ACF_list, "names"),
    count_NaCl_sample = as.numeric(ACF_list),
    stringsAsFactors = FALSE
    ) -> ACF_df

ACF_df %>% 
    dplyr::filter(count_NaCl_sample >= args_N_obs) %>% 
        dplyr::arrange(desc(count_NaCl_sample)) -> ACF_selection

#### save ACF_selection list####
write.csv(ACF_selection, 
          args_output, 
          row.names=TRUE)

#### save sample table####
ACF_table %>% 
    dplyr::filter(cell_type %in% ACF_selection$cell_type) %>% 
        dplyr::select(sample_number, cell_type) -> ACF_sample
write.csv(ACF_sample, 
          args_output_sample, 
          row.names=FALSE)
