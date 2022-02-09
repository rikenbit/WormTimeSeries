source("src/functions_WTS4_Eval_sample_behavior.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
# input merged_cls
args_input_cls <- args[1]
args_output <- args[2]
args_eval_method <- args[3]
args_eval_label <- args[4]

# # #### test args####
# # input sample_cls
# args_input_cls <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Cluster_sample/k_Number_5/sample_cls.RData")
# # output eval_result
# args_output <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/Eval_sample/NMI_behavior/k_Number_5.RData")
# # Evaluation Method
# args_eval_method <- c("NMI")
# # Evaluation label list
# args_eval_label <- c("data/WTS4_Eval_behavior_fix.xlsx")

##### load WTS4_Eval_behavior.xlsx####
read.xlsx(args_eval_label,
          sheet = "Sheet1",
          rowNames = FALSE,
          colNames =TRUE) %>% 
    dplyr::rename(CellType = celltype, 
                  Classes = class) -> df_eval_label

##### load sample_cls list####
load(args_input_cls)
lapply(C, function(x) {
    data.frame(CellType = attr(x, "names"),
               Clusters = as.numeric(x),
               stringsAsFactors = FALSE,
               row.names = NULL
               ) -> df_cls
    cellnames <- df_cls$CellType
    df_cls[grep("^[0-9]", cellnames, invert=TRUE),]
    }
) -> df_cls_list

#### merge dataframe####
lapply(df_cls_list, function(x) {
	df_cls_label <- merge(x, 
                      df_eval_label, 
                      by.x = "CellType", 
                      by.y = "CellType", 
                      all.x = TRUE)
	df_cls_label %>% 
		mutate_at(c("Classes"), 
	              ~replace(., 
	                       is.na(.), 
	                       "others")
	    )
	}
) -> df_cls_label

#### Evaluation####
eval_result <- switch(args_eval_method,
                    # eval method
                    "ARI" = .ARI_list(df_cls_label),
                    "purity" = .purity_list(df_cls_label),
                    "Fmeasure" = .Fmeasure_list(df_cls_label),
                    "Entropy" = .Entropy_list(df_cls_label),
                    "NMI" = .NMI_list(df_cls_label),
                    stop("Only can use all, ARI, purity, Fmeasure, Entropy")
                    )
save(eval_result, file=args_output)