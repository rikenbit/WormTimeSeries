source("src/functions_WTS4_DimReduc_cord.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
args_input_distance <- args[1]
args_DimReduc <- args[2]
args_output  <- args[3]

#### load dist object####
load(args_input_distance)
d <- merged_distance

#### Dimensionality Reduction####
set.seed(1234)
df_cord <- switch(args_DimReduc,
          "tsne" = .wts_tsne(d),
          "umap" = .wts_umap(d),
          stop("Only can use tsne,")
          )

#### save####
save(df_cord, file=args_output)