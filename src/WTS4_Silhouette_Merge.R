source("src/functions_WTS4_Silhouette_Merge.R")

#### args setting####
args <- commandArgs(trailingOnly = T)
# input merged_distance
args_input_distance <- args[1]
# input merged_cls
args_input_cls <- args[2]
# output eval_result
args_output_value <- args[3]
# output silhouette eps
args_output_plot <- args[4]
# output silhouette gg object
args_output_gg <- args[5]

#### test args####
# args_input_distance <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/Merged_distance/k_Number_4.RData")
# args_input_cls <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/Merged_cls/k_Number_4.RData")
# args_output_value <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/Eval/Silhouette/k_Number_4.RData")
# args_output_plot <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/Sil_plot/k_Number_4.eps")
# args_output_gg <- c("output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/Sil_gg/k_Number_4.RData")

#### load####
load(args_input_distance)
d <- merged_distance
load(args_input_cls)
cls <- merged_cls

#### silhouette####
sil <- silhouette(cls, d)
rownames(sil) <- attr(d, "Labels")
gg_sil <- fviz_silhouette(sil)
eval_result <- mean(gg_sil$data$sil_width)

#### save eval_result####
# save(eval_result, file=args_output_value)

#### silhouette ggplot####
gg_sil <- myfviz_silhouette(sil, gg_sil$data$cluster, label=TRUE) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.text.x = element_text(size = 8)) +
  theme(text = element_text(size = 90)) +
  theme(axis.text.x = element_blank())

#### ggsave fviz_silhouette####
ggsave(filename = args_output_plot,
       plot = gg_sil,
       dpi = 100,
       width = 30.0,
       height = 20.0,
       limitsize = FALSE,
       device="eps")

#### save gg_sil####
# patchworkで次元圧縮図に加える
# save(gg_sil, file=args_output_gg)
