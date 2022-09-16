source("src/functions_WTS4_yshift_visualize.R")

#### args setting####
#### test args####
# input matrix
args_input <- ("output/WTS4/normalize_1/stimAfter/SBD_abs/Shift_F/SampleNumber_1.RData")

#### test data####
load(args_input)
shift_matrix_F <- shift_matrix_F[1:10,1:10]

#### load####
# load(args_input)
input_mat <- shift_matrix_F

#### visualize ####
input_mat |> 
    as.data.frame() |> 
    rownames_to_column("row_celltype") |> 
    pivot_longer(-row_celltype, 
                 names_to = "col_celltype", 
                 values_to = "shift_value") -> long_celltype

ghm <- ggplot_ghm(long_celltype)
#### visualize abs ####
abs(input_mat) |> 
    as.data.frame() |> 
    rownames_to_column("row_celltype") |> 
    pivot_longer(-row_celltype, 
                 names_to = "col_celltype", 
                 values_to = "shift_value") -> long_celltype_abs

ghm_abs <- ggplot_ghm(long_celltype_abs)

#### save####