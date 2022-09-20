source("src/functions_WTS4_yshift_visualize.R")

#### args setting####

#### test args####
# input matrix
args_input <- ("output/WTS4/normalize_1/stimAfter/SBD_abs/Shift_F/SampleNumber_1.RData")
# labels
args_label <-("data/WTS4_Eval_behavior_ACF.xlsx")
# value type
args_type <- ("abs")
# filter (label combination)
args_label_comb <- c("ALL")
# output 
args_output <- ("output/WTS4/normalize_1/stimAfter/SBD_abs/Shift_F/vis_abs/ALL.RData")

#### load####
load(args_input)
input_mat <- shift_matrix_F

read.xlsx(args_label,
          sheet = "Sheet1",
          rowNames = FALSE,
          colNames =TRUE) |> 
    dplyr::rename(CellType = celltype, 
                  Classes = class) |> 
    dplyr::arrange(Classes,CellType) -> label_df

#### filter cell type####
label_df %>%
    dplyr::filter(Classes %in% vis_labels[[args_label_comb]]) %>%
        .$CellType -> label_F
which_F <- which(!is.na(match(colnames(input_mat),label_F)))
input_mat_F <- input_mat[which_F, which_F]

#### sort cell type####
label_df %>%
    dplyr::filter(CellType %in% colnames(input_mat_F)) %>%
    .$CellType -> input_mat_order
input_mat_F_S <- input_mat_F[input_mat_order, input_mat_order]
#### visualize ####
ghm <- vis_ghm[[args_type]](input_mat_F_S)
#### save####
ggsave(filename = args_output, 
       plot = ghm, 
       dpi = 50, 
       width = 30, 
       height = 40
       )
