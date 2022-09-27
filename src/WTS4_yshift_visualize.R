source("src/functions_WTS4_yshift_visualize.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
# input matrix
args_input <- args[1]
# labels
args_label <- args[2]
# value type
args_type <- args[3]
# filter (label combination)
args_label_comb <- args[4]
# output 
args_output <- args[5]
    
# #### test args####
# # input matrix
# args_input <- ("output/WTS4/normalize_1/stimAfter/SBD_abs/Shift_F/SampleNumber_1.RData")
# # labels
# args_label <-("data/WTS4_Eval_behavior_ACF.xlsx")
# # value type
# args_type <- ("abs")
# # filter (label combination)
# args_label_comb <- c("ALL")
# # output 
# args_output <- ("output/WTS4/normalize_1/stimAfter/SBD_abs/Shift_F/vis_abs/ALL.png")

#### load####
# load matrix
load(args_input)
input_mat <- shift_matrix
# load label table
read.xlsx(args_label,
          sheet = "Sheet1",
          rowNames = FALSE,
          colNames =TRUE) |> 
    dplyr::rename(CellType = celltype, 
                  Classes = class) |> 
    dplyr::arrange(Classes,CellType) -> label_df

#### filter cell type####
if (args_label_comb=="No_F") {
    input_mat_F <- input_mat
} else {
    label_df %>%
        dplyr::filter(Classes %in% vis_labels[[args_label_comb]]) %>%
            .$CellType -> label_F
    which_F <- which(!is.na(match(colnames(input_mat),label_F)))
    input_mat_F <- input_mat[which_F, which_F]
}

#### sort cell type####
if (args_label_comb=="No_F") {
    label_df %>%
        dplyr::filter(CellType %in% colnames(input_mat)) %>%
        .$CellType -> hit_cell
    unhit_cell <- colnames(input_mat)[which(is.na(match(colnames(input_mat), hit_cell)))]
    No_F_order <- c(hit_cell, unhit_cell)
    input_mat_F_S <- input_mat_F[No_F_order, No_F_order]
} else {
    label_df %>%
        dplyr::filter(CellType %in% colnames(input_mat_F)) %>%
            .$CellType -> input_mat_order
    input_mat_F_S <- input_mat_F[input_mat_order, input_mat_order]
}

#### visualize ####
if (args_label_comb=="No_F") {
    args_textsize <- 25
} else {
    args_textsize <- 60
    }

ghm <- vis_ghm[[args_type]](input_mat_F_S) + 
    # theme(text = element_text(size = 40)) +
    theme(text = element_text(size = args_textsize)) +
    theme(legend.key.height = unit(1.5, "cm")) +
    theme(legend.key.width = unit(1.5, "cm"))

#### save####
ggsave(filename = args_output, 
       plot = ghm, 
       dpi = 80, 
       width = 40, 
       height = 36
       )