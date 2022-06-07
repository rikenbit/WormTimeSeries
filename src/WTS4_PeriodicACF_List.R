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

# for (i in 1:nrow(ACF_sample)) {
#     sample_number <-ACF_sample$sample_number[i]
#     sample_celltype<- ACF_sample$cell_type[i]
#     
#     eval(parse(text=paste0("sample_fullpath <- list.files(path='output/WTS2/correlogram/Data_normalize_1/TF_after/LAG_600/ACF_Acf/SampleNumber_",sample_number,"/',
#                               pattern = sample_celltype,
#                               full.names = TRUE)")))
#     eval(parse(text=paste0("dir.create('output/WTS2/correlogram/Data_normalize_1/TF_after/LAG_600/",sample_celltype,"', showWarnings = TRUE)")))
#     eval(parse(text=paste0("cp_path <- c('output/WTS2/correlogram/Data_normalize_1/TF_after/LAG_600/",sample_celltype,"/SampleNumber_",sample_number,".png')")))
#     file.copy(sample_fullpath, cp_path)
# }

# # #sample_number <- 1
# # sample_number <-ACF_sample$sample_number[1]
# # # sample_celltype <- c("ASER")
# # sample_celltype<- ACF_sample$cell_type[1]
# # 
# # sample_fullpath <- c("output/WTS2/correlogram/Data_normalize_1/TF_after/LAG_600/ACF_Acf/SampleNumber_1/CellNumber_89_CellType_ASKR.png")
# eval(parse(text=paste0("sample_fullpath <- list.files(path='output/WTS2/correlogram/Data_normalize_1/TF_after/LAG_600/ACF_Acf/SampleNumber_",sample_number,"/',
#                               pattern = sample_celltype,
#                               full.names = TRUE)")))
# # dir.create('output/WTS2/correlogram/Data_normalize_1/TF_after/LAG_600/ASKR', showWarnings = TRUE)
# eval(parse(text=paste0("dir.create('output/WTS2/correlogram/Data_normalize_1/TF_after/LAG_600/",sample_celltype,"', showWarnings = TRUE)")))
# # cp_path <- c('output/WTS2/correlogram/Data_normalize_1/TF_after/LAG_600/ASKR/SampleNumber_1.png')
# eval(parse(text=paste0("cp_path <- c('output/WTS2/correlogram/Data_normalize_1/TF_after/LAG_600/",sample_celltype,"/SampleNumber_",sample_number,".png')")))
# file.copy(sample_fullpath,cp_path)

# for (i in 1:nrow(ACF_sample)) {
#     sample_number <-ACF_sample$sample_number[i]
#     sample_celltype<- ACF_sample$cell_type[i]
#     
#     eval(parse(text=paste0("sample_fullpath <- list.files(path='output/WTS1/plot/normalize_1/SampleNumber_",sample_number,"/',
#                               pattern = sample_celltype,
#                               full.names = TRUE)")))
#     eval(parse(text=paste0("dir.create('output/WTS1/plot/normalize_1/",sample_celltype,"', showWarnings = TRUE)")))
#     eval(parse(text=paste0("cp_path <- c('output/WTS1/plot/normalize_1/",sample_celltype,"/SampleNumber_",sample_number,".png')")))
#     file.copy(sample_fullpath, cp_path)
# }
